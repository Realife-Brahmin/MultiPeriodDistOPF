# Is FilterDDP's (n_x+1)-column KKT solve slow because UMFPACK solves one column
# at a time? Same factors, different solve algorithm.
#
# UMFPACK's solve takes one right-hand side per call, so FilterDDP's ldiv!(F, B)
# walks the whole factor once per column (1,021 times at large10k), doing one
# multiply-add per factor entry read. A blocked solve walks the factor once per
# block of w columns and applies each entry to all w columns (a contiguous axpy),
# so it reads the factor n_cols/w times. The arithmetic is identical.
#
# Both use the SAME UMFPACK factorization:  L*U = (Rs .* K)[p, q], so
#   K X = B  <=>  L U Y = (Rs .* B)[p, :],  X[q, :] = Y.
# The blocked path pays for extracting L, U, p, q, Rs from UMFPACK (timed
# separately, since FilterDDP would pay it on every factorization) and for the
# scale/permute/transpose into and out of the block. Sequential (one thread)
# unless KKT_BLOCK_THREADS > 1, which is reported as a separate, labelled row.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/blocked_multirhs_solve_benchmark.jl \
#         <capture.jls> <system> <arm> <out.csv> [repeats]

using LinearAlgebra, SparseArrays, Serialization, Printf, Statistics, Dates

input, system, arm, out_csv = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
repeats = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 5
BLAS.set_num_threads(1)

cap = deserialize(input)
K = SparseMatrixCSC{Float64, Int64}(cap.K)
B = Matrix{Float64}(cap.rhs)
n, m = size(B)
@printf("%s %s: n=%d nnz(K)=%d rhs_cols=%d repeats=%d julia_threads=%d\n",
        system, arm, n, nnz(K), m, repeats, Threads.nthreads())

# Forward substitution with a CSC lower-triangular L on a block stored
# transposed (Bt is w x n, so Bt[:, j] is one contiguous row of the block).
function lower_solve!(Bt::Matrix{Float64}, L::SparseMatrixCSC{Float64, Int64})
    cp, rv, nz = L.colptr, L.rowval, L.nzval
    w = size(Bt, 1)
    @inbounds for j in 1:size(Bt, 2)
        k1, k2 = cp[j], cp[j + 1] - 1
        rv[k1] == j || error("L column $j: diagonal is not the first entry")
        d = nz[k1]
        if d != 1.0
            @simd for r in 1:w
                Bt[r, j] /= d
            end
        end
        for k in k1 + 1:k2
            i = rv[k]; v = nz[k]
            @simd for r in 1:w
                Bt[r, i] -= v * Bt[r, j]
            end
        end
    end
    return Bt
end

# Back substitution with a CSC upper-triangular U (diagonal is each column's last entry).
function upper_solve!(Bt::Matrix{Float64}, U::SparseMatrixCSC{Float64, Int64})
    cp, rv, nz = U.colptr, U.rowval, U.nzval
    w = size(Bt, 1)
    @inbounds for j in size(Bt, 2):-1:1
        k1, k2 = cp[j], cp[j + 1] - 1
        rv[k2] == j || error("U column $j: diagonal is not the last entry")
        d = nz[k2]
        @simd for r in 1:w
            Bt[r, j] /= d
        end
        for k in k1:k2 - 1
            i = rv[k]; v = nz[k]
            @simd for r in 1:w
                Bt[r, i] -= v * Bt[r, j]
            end
        end
    end
    return Bt
end

# Solve K X = B for all columns, w at a time, into X.
function blocked_solve!(X, B, L, U, p, q, Rs, w)
    n, m = size(B)
    Bt = Matrix{Float64}(undef, min(w, m), n)
    for c0 in 1:w:m
        cols = c0:min(c0 + w - 1, m)
        nb = length(cols)
        Bb = nb == size(Bt, 1) ? Bt : Matrix{Float64}(undef, nb, n)
        @inbounds for i in 1:n
            pi_ = p[i]; s = Rs[pi_]
            for (r, c) in enumerate(cols)
                Bb[r, i] = s * B[pi_, c]
            end
        end
        lower_solve!(Bb, L)
        upper_solve!(Bb, U)
        @inbounds for i in 1:n
            qi = q[i]
            for (r, c) in enumerate(cols)
                X[qi, c] = Bb[r, i]
            end
        end
    end
    return X
end

# The same, with column blocks spread over Julia threads (each thread its own block).
function blocked_solve_threaded!(X, B, L, U, p, q, Rs, w)
    n, m = size(B)
    blocks = collect(1:w:m)
    Threads.@threads for c0 in blocks
        cols = c0:min(c0 + w - 1, m)
        nb = length(cols)
        Bb = Matrix{Float64}(undef, nb, n)
        @inbounds for i in 1:n
            pi_ = p[i]; s = Rs[pi_]
            for (r, c) in enumerate(cols)
                Bb[r, i] = s * B[pi_, c]
            end
        end
        lower_solve!(Bb, L)
        upper_solve!(Bb, U)
        @inbounds for i in 1:n
            qi = q[i]
            for (r, c) in enumerate(cols)
                X[qi, c] = Bb[r, i]
            end
        end
    end
    return X
end

function main(K, B, n, m)
relres(X) = norm(K * X - B) / norm(B)
rows = Vector{NamedTuple}()
F = lu(K)                                        # FilterDDP's factorization, unchanged
nLU = nnz(F.L) + nnz(F.U)

# Baseline: exactly what FilterDDP does.
tb = Float64[]; local Xref
for r in 0:repeats
    GC.gc(); Xr = copy(B)
    t = @elapsed ldiv!(F, Xr)
    r == 0 ? (Xref = Xr) : push!(tb, t)
end
base = median(tb)
push!(rows, (method = "umfpack_ldiv", block = 1, threads = 1, extract_s = 0.0, solve_s = base,
             total_s = base, speedup = 1.0, relres = relres(Xref), maxreldiff = 0.0))
@printf("  %-22s solve %8.4f s  res %.2e\n", "UMFPACK ldiv! (baseline)", base, relres(Xref))

te = Float64[]
local Lf, Uf, pf, qf, Rsf
for r in 0:repeats
    GC.gc()
    t = @elapsed begin
        Lf = F.L; Uf = F.U; pf = F.p; qf = F.q; Rsf = F.Rs
    end
    r == 0 || push!(te, t)
end
ext = median(te)
@printf("  factor extraction (L, U, p, q, Rs): %.4f s   nnz(L)+nnz(U)=%d\n", ext, nLU)

widths = unique(filter(w -> w <= m, [1, 4, 8, 16, 32, 64, 128, 256, m]))
nthreads = parse(Int, get(ENV, "KKT_BLOCK_THREADS", "1"))
for (fn, label, th) in ((blocked_solve!, "blocked", 1),
                        (blocked_solve_threaded!, "blocked_threaded", Threads.nthreads()))
    th == 1 && label == "blocked_threaded" && continue
    for w in widths
        ts = Float64[]; X = zeros(n, m)
        for r in 0:repeats
            GC.gc(); fill!(X, 0.0)
            t = @elapsed fn(X, B, Lf, Uf, pf, qf, Rsf, w)
            r == 0 || push!(ts, t)
        end
        s = median(ts)
        md = maximum(abs.(X .- Xref)) / max(maximum(abs.(Xref)), eps())
        rr = relres(X)
        push!(rows, (method = label, block = w, threads = th, extract_s = ext, solve_s = s,
                     total_s = ext + s, speedup = base / (ext + s), relres = rr, maxreldiff = md))
        @printf("  %-16s w=%5d  solve %8.4f s  +extract = %8.4f s  speedup %5.2fx  res %.2e  max|dX|/max|X| %.1e\n",
                label, w, s, ext + s, base / (ext + s), rr, md)
    end
end

mkpath(dirname(out_csv))
open(out_csv, "w") do io
    println(io, "system,arm,n,nnz_K,nnz_LU,rhs_cols,method,block,threads,extract_s,solve_s,total_s,speedup_vs_umfpack_ldiv,relres,max_rel_diff_vs_umfpack")
    for r in rows
        @printf(io, "%s,%s,%d,%d,%d,%d,%s,%d,%d,%.6g,%.6g,%.6g,%.4f,%.3e,%.3e\n", system, arm, n, nnz(K),
                nLU, m, r.method, r.block, r.threads, r.extract_s, r.solve_s, r.total_s, r.speedup,
                r.relres, r.maxreldiff)
    end
end
println("wrote $out_csv  finished=$(now())")
end

main(K, B, n, m)
