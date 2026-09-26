# HSL MA57 and HSL_MA97 on the captured FilterDDP stage KKT systems, against
# FilterDDP's UMFPACK LU and the blocked multi-column solve over its factors.
#
# Same captures, timing discipline and residual checks as
# kkt_ordering_benchmark.jl: one discarded warm-up, medians of `repeats`
# repetitions of analysis/ordering, numerical factorization, a one-column solve
# and the full (n_x+1)-column solve; relative residual ||K X - B||_F / ||B||_F.
#
# The HSL libraries are licensed and are NOT in this repository. They are built
# locally from the user's HSL sources (hsl/build_hsl_windows.sh) into
# HSL_LIB_DIR: libma57.dll (MA57 3.11.3) and libhsl_ma97.dll (HSL_MA97 2.8.1
# plus the s97 shim). Both use the placeholder METIS, so METIS orderings are
# unavailable; OpenBLAS threads are set by OPENBLAS_NUM_THREADS and MA97's
# OpenMP threads by OMP_NUM_THREADS (the caller sets both).
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/hsl_kkt_benchmark.jl \
#         <capture.jls> <system> <arm> <out.csv> [repeats] [families]
# families: comma list of umfpack,blocked,ma57,ma97 (default all)

using LinearAlgebra, SparseArrays, Serialization, Printf, Statistics, Dates
using DDP4OPF

const HSL_DIR = get(ENV, "HSL_LIB_DIR", raw"C:\Users\Aryan Ritwajeet Jha\Documents\hsl\build")
const LIB57 = joinpath(HSL_DIR, "libma57.dll")
const LIB97 = joinpath(HSL_DIR, "libhsl_ma97.dll")

input, system, arm, out_csv = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
repeats = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 5
families = Set(split(length(ARGS) >= 6 ? ARGS[6] : "umfpack,blocked,ma57,ma97", ','))
BLAS.set_num_threads(1)
omp = get(ENV, "OMP_NUM_THREADS", "unset"); obl = get(ENV, "OPENBLAS_NUM_THREADS", "unset")

cap = deserialize(input)
K = SparseMatrixCSC{Float64, Int64}(cap.K)
B = Matrix{Float64}(cap.rhs)
n, m = size(B); nnzK = nnz(K)
b1 = B[:, 1:1]
Kl = tril(K)                                   # MA57/MA97 take one triangle
relres(X, R) = norm(K * X - R) / norm(R)
@printf("%s %s: n=%d nnz(K)=%d rhs=%d repeats=%d OMP=%s OPENBLAS=%s julia_threads=%d\n",
        system, arm, n, nnzK, m, repeats, omp, obl, Threads.nthreads())

rows = NamedTuple[]
function record!(solver, config, threads, to, tf, t1, tw, entries, extra, X1, XW, note)
    r = (system = system, arm = arm, n = n, nnz_K = nnzK, rhs_cols = m, solver = solver,
         config = config, threads = threads, ordering_s = median(to), factor_s = median(tf),
         solve1_s = median(t1), solvewide_s = median(tw),
         factor_plus_wide_s = median(tf) + median(tw), total_s = median(to) + median(tf) + median(tw),
         factor_entries = entries, fill_ratio = entries / nnzK, stats = extra,
         relres_1 = relres(X1, b1), relres_wide = relres(XW, B), note = note)
    push!(rows, r)
    @printf("  %-8s %-26s thr %-3s ord %8.4f fact %8.4f solve1 %8.5f wide %8.4f f+w %8.4f | fill %5.2f res %.1e/%.1e %s\n",
            solver, config, threads, r.ordering_s, r.factor_s, r.solve1_s, r.solvewide_s,
            r.factor_plus_wide_s, r.fill_ratio, r.relres_1, r.relres_wide, extra)
    flush(stdout)
end

# ---------------------------------------------------------------- UMFPACK --
const U = SparseArrays.UMFPACK
if "umfpack" in families || "blocked" in families
    to = Float64[]; tf = Float64[]; t1 = Float64[]; tw = Float64[]; tb = Float64[]
    local X1, XW, XB, F
    for r in 0:repeats
        GC.gc()
        F = U.UmfpackLU(K)
        a = @elapsed U.umfpack_symbolic!(F, nothing)
        f = @elapsed U.umfpack_numeric!(F)
        X1 = copy(b1); s1 = @elapsed ldiv!(F, X1)
        XW = copy(B);  sw = @elapsed ldiv!(F, XW)
        XB = copy(B);  sb = @elapsed DDP4OPF._blocked_umfpack_solve!(F, XB, 16)
        r == 0 && continue
        push!(to, a); push!(tf, f); push!(t1, s1); push!(tw, sw); push!(tb, sb)
    end
    nLU = nnz(F.L) + nnz(F.U)
    "umfpack" in families && record!("UMFPACK", "default (FilterDDP lu)", "1", to, tf, t1, tw, nLU, "", X1, XW, "column-by-column ldiv!")
    "blocked" in families && record!("UMFPACK", "blocked w=16", string(Threads.nthreads()), to, tf, t1, tb, nLU, "", X1, XB,
                                     "same factors, blocked solve incl. extraction")
end

# ------------------------------------------------------------------- MA57 --
function ma57_run(icntl6, icntl15, cntl1)
    cntl = zeros(Float64, 5); icntl = zeros(Cint, 20)
    ccall((:ma57id_, LIB57), Cvoid, (Ptr{Float64}, Ptr{Cint}), cntl, icntl)
    icntl[1] = -1; icntl[2] = -1; icntl[3] = -1; icntl[5] = 0
    icntl[6] = icntl6; icntl[15] = icntl15; cntl[1] = cntl1
    irn = Cint.(rowvals(Kl)); jcn = Cint[j for j in 1:n for _ in nzrange(Kl, j)]; a = copy(nonzeros(Kl))
    ne = Cint(length(a)); nn = Cint(n)
    lkeep = Cint(5n + ne + max(n, ne) + 42); keep = zeros(Cint, lkeep)
    info = zeros(Cint, 40); rinfo = zeros(Float64, 20)
    ta = @elapsed ccall((:ma57ad_, LIB57), Cvoid,
        (Ref{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}),
        nn, ne, irn, jcn, lkeep, keep, zeros(Cint, 5n), icntl, info, rinfo)
    info[1] < 0 && error("MA57AD INFO(1)=$(info[1]) INFO(2)=$(info[2])")
    ordering_used = info[36]
    lfact = Cint(ceil(1.5 * info[9])); lifact = Cint(ceil(1.5 * info[10]))
    local fact, ifact
    tf = 0.0
    for attempt in 1:4
        fact = zeros(Float64, lfact); ifact = zeros(Cint, lifact)
        tf = @elapsed ccall((:ma57bd_, LIB57), Cvoid,
            (Ref{Cint}, Ref{Cint}, Ptr{Float64}, Ptr{Float64}, Ref{Cint}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Float64}),
            nn, ne, a, fact, lfact, ifact, lifact, lkeep, keep, zeros(Cint, n), icntl, cntl, info, rinfo)
        info[1] in (-3, -4) || break
        # Delayed pivots needed more space than forecast: grow and refactorize
        # (the retry is counted in the timing, as it would be in use).
        lfact = Cint(ceil(1.5 * max(lfact, info[17]))); lifact = Cint(ceil(1.5 * max(lifact, info[18])))
    end
    info[1] < 0 && error("MA57BD INFO(1)=$(info[1]) INFO(2)=$(info[2])")
    stats = (entries = Int(info[14]), two = Int(info[22]), delayed = Int(info[23]), neg = Int(info[24]),
             rank = Int(info[25]), maxfront = Int(info[21]), ordering = Int(ordering_used))
    solve(X) = begin
        nr = Cint(size(X, 2)); w = zeros(Float64, n * size(X, 2)); inf2 = zeros(Cint, 40)
        t = @elapsed ccall((:ma57cd_, LIB57), Cvoid,
            (Ref{Cint}, Ref{Cint}, Ptr{Float64}, Ref{Cint}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ptr{Float64}, Ref{Cint},
             Ptr{Float64}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
            Cint(1), nn, fact, lfact, ifact, lifact, nr, X, nn, w, Cint(length(w)), zeros(Cint, n), icntl, inf2)
        inf2[1] < 0 && error("MA57CD INFO(1)=$(inf2[1])")
        t
    end
    X1 = copy(b1); s1 = solve(X1)
    XW = copy(B);  sw = solve(XW)
    return ta, tf, s1, sw, stats, X1, XW
end

if "ma57" in families
    for (name, o6, s15, c1, note) in (
        ("auto (ICNTL6=5)",        5, 1, 0.01, "MA57 defaults; METIS unavailable so automatic = MC47 AMD"),
        ("AMD MC47 (ICNTL6=2)",    2, 1, 0.01, ""),
        ("AMD no-dense (ICNTL6=0)",0, 1, 0.01, ""),
        ("MA27 MD (ICNTL6=3)",     3, 1, 0.01, ""),
        ("Ipopt settings",         5, 0, 1e-8, "Ipopt's ma57 defaults: pivtol 1e-8, no MC64 scaling"))
        try
            local to = Float64[]; local tf = Float64[]; local t1 = Float64[]; local tw = Float64[]
            local st, X1, XW
            for r in 0:repeats
                GC.gc()
                a, f, s1, sw, st, X1, XW = ma57_run(o6, s15, c1)
                r == 0 && continue
                push!(to, a); push!(tf, f); push!(t1, s1); push!(tw, sw)
            end
            extra = @sprintf("ordering_used=%d neg=%d delayed=%d two_by_two=%d rank=%d maxfront=%d",
                             st.ordering, st.neg, st.delayed, st.two, st.rank, st.maxfront)
            # MA57 stores L and D once; 2x makes fill comparable with LU's nnz(L)+nnz(U).
            record!("MA57", name, obl, to, tf, t1, tw, 2 * st.entries, extra, X1, XW, note)
        catch err
            @printf("  MA57     %-26s FAILED: %s\n", name, first(split(sprint(showerror, err), '\n')))
            push!(rows, (system = system, arm = arm, n = n, nnz_K = nnzK, rhs_cols = m, solver = "MA57",
                         config = name, threads = obl, ordering_s = NaN, factor_s = NaN, solve1_s = NaN,
                         solvewide_s = NaN, factor_plus_wide_s = NaN, total_s = NaN, factor_entries = 0,
                         fill_ratio = NaN, stats = "", relres_1 = NaN, relres_wide = NaN,
                         note = "FAILED: " * first(split(sprint(showerror, err), '\n'))))
        end
    end
end

# ------------------------------------------------------------------- MA97 --
if "ma97" in families
    ptr = Cint.(Kl.colptr); row = Cint.(rowvals(Kl)); val = copy(nonzeros(Kl))
    info = zeros(Float64, 12)
    for (name, ord, scal, blas3, note) in (
        ("auto, BLAS3 solve",   5, 1, 1, "METIS unavailable so automatic = AMD"),
        ("auto, BLAS2 solve",   5, 1, 0, "MA97 default solve_blas3=false"),
        ("AMD, BLAS3 solve",    1, 1, 1, ""),
        ("MD (MA27), BLAS3",    2, 1, 1, ""),
        ("auto, no scaling",    5, 0, 1, ""))
        try
            local to = Float64[]; local tf = Float64[]; local t1 = Float64[]; local tw = Float64[]
            local X1, XW
            for r in 0:repeats
                GC.gc()
                ccall((:s97_init, LIB97), Cvoid, (Cint, Cint, Cint, Cdouble), ord, scal, blas3, -1.0)
                a = @elapsed (fl = ccall((:s97_analyse, LIB97), Cint, (Cint, Ptr{Cint}, Ptr{Cint}), n, ptr, row))
                fl < 0 && error("ma97_analyse flag $fl")
                f = @elapsed (fl = ccall((:s97_factor, LIB97), Cint, (Ptr{Cint}, Ptr{Cint}, Ptr{Float64}), ptr, row, val))
                fl < 0 && error("ma97_factor flag $fl")
                ccall((:s97_info, LIB97), Cvoid, (Ptr{Float64},), info)
                X1 = copy(b1); s1 = @elapsed (fl = ccall((:s97_solve, LIB97), Cint, (Cint, Ptr{Float64}, Cint), 1, X1, n))
                fl < 0 && error("ma97_solve flag $fl")
                XW = copy(B);  sw = @elapsed (fl = ccall((:s97_solve, LIB97), Cint, (Cint, Ptr{Float64}, Cint), m, XW, n))
                fl < 0 && error("ma97_solve flag $fl")
                ccall((:s97_free, LIB97), Cvoid, ())
                r == 0 && continue
                push!(to, a); push!(tf, f); push!(t1, s1); push!(tw, sw)
            end
            extra = @sprintf("ordering_used=%d neg=%d delayed=%d two_by_two=%d rank=%d maxfront=%d",
                             Int(info[8]), Int(info[4]), Int(info[5]), Int(info[6]), Int(info[9]), Int(info[7]))
            record!("MA97", name, omp, to, tf, t1, tw, 2 * Int(info[2]), extra, X1, XW, note)
        catch err
            ccall((:s97_free, LIB97), Cvoid, ())
            @printf("  MA97     %-26s FAILED: %s\n", name, first(split(sprint(showerror, err), '\n')))
            push!(rows, (system = system, arm = arm, n = n, nnz_K = nnzK, rhs_cols = m, solver = "MA97",
                         config = name, threads = omp, ordering_s = NaN, factor_s = NaN, solve1_s = NaN,
                         solvewide_s = NaN, factor_plus_wide_s = NaN, total_s = NaN, factor_entries = 0,
                         fill_ratio = NaN, stats = "", relres_1 = NaN, relres_wide = NaN,
                         note = "FAILED: " * first(split(sprint(showerror, err), '\n'))))
        end
    end
end

cols = [:system, :arm, :n, :nnz_K, :rhs_cols, :solver, :config, :threads, :ordering_s, :factor_s,
        :solve1_s, :solvewide_s, :factor_plus_wide_s, :total_s, :factor_entries, :fill_ratio,
        :stats, :relres_1, :relres_wide, :note]
cell(v) = v isa AbstractString ? "\"" * replace(v, "\"" => "'") * "\"" :
          v isa AbstractFloat ? @sprintf("%.6g", v) : string(v)
mkpath(dirname(out_csv))
open(out_csv, "w") do io
    println(io, join(string.(cols), ","))
    for r in rows
        println(io, join((cell(getfield(r, c)) for c in cols), ","))
    end
end
println("wrote $out_csv ($(length(rows)) rows) finished=$(now())")
