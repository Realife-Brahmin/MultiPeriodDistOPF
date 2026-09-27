# Are the KKT entries significant, or are many of them tiny? (Agenda, R. Gupta,
# meeting of 2026-10-02.) On a captured stage KKT system K = [H cu'; cu 0]:
#
# 1. Relative magnitude of every entry after symmetric scaling,
#        r_ij = |K_ij| / sqrt(d_i d_j),   d_i = max_k |K_ik|,
#    so r_ij <= 1 and the diagonal scaling itself (which is harmless to a
#    factorization) does not count as "small". Distribution per block: the
#    Hessian diagonal and the Jacobian by constraint type.
# 2. Thresholding: drop every off-diagonal entry with r_ij < tau, or with
#    |K_ij| < tau (the plain reading of the question) (symmetrically;
#    the diagonal is always kept), factor and solve the thresholded matrix, and
#    report the entries dropped, the factor fill, the time, and how far the
#    (n_x+1)-column solution moves from the exact one -- a different sparsity
#    plot only matters if the solution is unchanged.
# 3. Writes the entries (i, j, log10 r_ij, log10 |K_ij|) for plot_kkt_magnitude.py.
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_magnitude_analysis.jl \
#         <capture.jls> <network_data.jls> <system> <iter> <out_dir> <pattern_dir>

using Serialization, SparseArrays, LinearAlgebra, Printf, Statistics

capfile, datafile, system, iter, outdir, patdir = ARGS[1:6]
mkpath(outdir); mkpath(patdir)
BLAS.set_num_threads(1)

c = deserialize(capfile)
K = SparseMatrixCSC{Float64, Int64}(c.K); B = Matrix{Float64}(c.rhs)
n = size(K, 1); nu, nc = c.nu, c.nc
d = deserialize(datafile)
N, L = length(d[:Nset]), length(d[:Lset]); nB = length(d[:Bset])
rowgroup = Vector{String}(undef, n)
rowgroup[1:nu] .= "H"
for (g, r) in (("Pbal", 1:N), ("Qbal", N+1:2N), ("vdrop", 2N+1:2N+L), ("SOCP", 2N+L+1:2N+2L),
               ("rootV", 2N+2L+1:2N+2L+1), ("energy", 2N+2L+2:2N+2L+1+nB))
    rowgroup[nu .+ r] .= g
end

I, J, V = findnz(K)
dmax = zeros(n)
for (i, v) in zip(I, V); dmax[i] = max(dmax[i], abs(v)); end
r = [abs(v) / sqrt(dmax[i] * dmax[j]) for (i, j, v) in zip(I, J, V)]
lr = log10.(max.(r, 1e-300))

# Block of each entry: Hessian diagonal, Hessian off-diagonal, or the Jacobian
# (either triangle) labelled by its constraint row.
block = [i == j && i <= nu ? "H.diag" : (i <= nu && j <= nu ? "H.offdiag" :
         (i > nu ? "J:" * rowgroup[i] : "J:" * rowgroup[j])) for (i, j) in zip(I, J)]

open(joinpath(patdir, "pattern_$(system)_iter$(iter).csv"), "w") do io
    println(io, "i,j,log10r,log10abs,block")
    for k in eachindex(I)
        @printf(io, "%d,%d,%.3f,%.3f,%s\n", I[k], J[k], lr[k], log10(max(abs(V[k]), 1e-300)), block[k])
    end
end

open(joinpath(outdir, "magnitude_blocks_$(system)_iter$(iter).csv"), "w") do io
    println(io, "system,iter,mu,block,entries,median_log10r,p10_log10r,p1_log10r,min_log10r," *
                "frac_below_1e-8,frac_below_1e-6,frac_below_1e-4,frac_below_1e-2")
    for b in sort(unique(block))
        x = sort(lr[block .== b])
        q(p) = x[max(1, ceil(Int, p * length(x)))]
        @printf(io, "%s,%s,%.3e,%s,%d,%.3f,%.3f,%.3f,%.3f,%.4e,%.4e,%.4e,%.4e\n",
                system, iter, c.barrier_mu, b, length(x), median(x), q(0.10), q(0.01), x[1],
                count(<(-8), x) / length(x), count(<(-6), x) / length(x),
                count(<(-4), x) / length(x), count(<(-2), x) / length(x))
    end
end

# Exact solution for reference (warm-up factorization first).
lu(K)
t0 = @elapsed F0 = lu(K)
X0 = copy(B); ts0 = @elapsed ldiv!(F0, X0)
nnz0 = nnz(F0.L) + nnz(F0.U)
rel(A, R) = norm(A - R) / norm(R)
open(joinpath(outdir, "magnitude_threshold_$(system)_iter$(iter).csv"), "w") do io
    println(io, "system,iter,mu,criterion,tau,dropped,dropped_frac,dropped_H,dropped_J,nnz_K,nnz_LU,factor_s,solve_s," *
                "dev_all,dev_primal,dev_feedforward_primal,relres_exact_K")
    @printf(io, "%s,%s,%.3e,none,0,0,0,0,0,%d,%d,%.4f,%.4f,0,0,0,%.3e\n", system, iter, c.barrier_mu,
            nnz(K), nnz0, t0, ts0, norm(K * X0 - B) / norm(B))
    # relative: r_ij < tau after symmetric scaling; absolute: |K_ij| < tau.
    for (criterion, mag) in (("relative", r), ("absolute", abs.(V))),
        tau in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2)
        keep = [(i == j) || (mag[k] >= tau) for (k, (i, j)) in enumerate(zip(I, J))]
        drop = .!keep
        dropped = count(drop)
        dH = count(drop .& startswith.(block, "H")); dJ = count(drop .& startswith.(block, "J"))
        Kt = sparse(I[keep], J[keep], V[keep], n, n)
        line = try
            lu(Kt)
            tf = @elapsed Ft = lu(Kt)
            Xt = copy(B); tsol = @elapsed ldiv!(Ft, Xt)
            @sprintf("%d,%d,%.4f,%.4f,%.3e,%.3e,%.3e,%.3e", nnz(Kt), nnz(Ft.L) + nnz(Ft.U), tf, tsol,
                     rel(Xt, X0), rel(Xt[1:nu, :], X0[1:nu, :]), rel(Xt[1:nu, 1], X0[1:nu, 1]),
                     norm(K * Xt - B) / norm(B))
        catch
            @sprintf("%d,singular,NaN,NaN,NaN,NaN,NaN,NaN", nnz(Kt))
        end
        @printf(io, "%s,%s,%.3e,%s,%.0e,%d,%.4e,%d,%d,%s\n", system, iter, c.barrier_mu, criterion, tau,
                dropped, dropped / length(I), dH, dJ, line)
        flush(io)
    end
end
println("KKT_MAGNITUDE $system iter $iter n=$n nnz=$(nnz(K)) done")
