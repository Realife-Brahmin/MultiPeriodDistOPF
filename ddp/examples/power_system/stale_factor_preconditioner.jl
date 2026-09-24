# Is a stale KKT factorisation reusable as a PRECONDITIONER rather than as the
# solve?  (Agenda point, R. Gupta 2026-09-18: "any scope of solving X = Ainv*B,
# assuming Ainv doesn't change much?")
#
# FILTERDDP_FREEZE_KKT already answered the naive version: reusing lu(K_j) as
# though it were lu(K_k) diverges, because K carries the interior-point barrier
# terms Sigma = z/s and those move by orders of magnitude.  That test throws the
# true K away.  The sound version keeps it: solve K_k x = b by iterative
# refinement preconditioned with the stale factor,
#
#     x <- x + M^{-1} (b - K_k x),     M = lu(K_j),  j < k
#
# which converges to the EXACT solution whenever rho(I - M^{-1} K_k) < 1, and
# costs one stale triangular solve per step instead of a fresh factorisation.
#
# Break-even is set by the measured factor/solve ratio.  k refinement steps beat
# a fresh factorisation when  k * solve < factor + solve,  i.e.  k < 1 + F/S.
# Measured F/S per stage: ieee123 T=3 2.71 (13.06 ms / 4.82 ms), large10k T=6
# 2.81 (12.88 s / 4.58 s).  So refinement has to converge in <= 3 steps to pay.
#
# Refinement (not GMRES) is the scheme that matters here because the backward
# pass solves against nx+1 right-hand sides at once; refinement handles the whole
# block in one application, block-GMRES would not.  Per-column GMRES is reported
# alongside as the best case any Krylov method could reach.
#
# Usage:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/stale_factor_preconditioner.jl <capture_dir> [stage] [out.csv]
# where <capture_dir> was produced by FILTERDDP_PERIODIC_CAPTURE_DIR with
# FILTERDDP_PERIODIC_CAPTURE_STRIDE=1.

using LinearAlgebra
using Printf
using Serialization
using SparseArrays

capture_dir = length(ARGS) >= 1 ? ARGS[1] : error("pass the capture directory")
stage       = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
out_csv     = length(ARGS) >= 3 ? ARGS[3] : joinpath(capture_dir, "stale_factor_preconditioner.csv")

files = sort(filter(f -> occursin(@sprintf("stage%03d.jls", stage), f), readdir(capture_dir)))
isempty(files) && error("no captures for stage $stage in $capture_dir")

snaps = map(files) do f
    p = deserialize(joinpath(capture_dir, f))
    (iter=p.iteration, K=p.K, rhs=Matrix(p.rhs_multi),
     Sigma=vcat(p.Sigma_L, p.Sigma_U), mu=p.barrier_mu)
end
sort!(snaps, by = s -> s.iter)
@printf("stage %d: %d snapshots, K is %dx%d with %d nonzeros, rhs has %d columns\n",
        stage, length(snaps), size(snaps[1].K)..., nnz(snaps[1].K), size(snaps[1].rhs, 2))

# Steps of preconditioned iterative refinement to drive the relative residual of
# the FULL multi-column system below tol.  Returns maxit+1 if it does not get
# there, and -1 if it diverges outright.
function refinement_steps(K, B, M; tol=1e-10, maxit=30)
    X = M \ B
    bn = norm(B)
    R = B - K * X
    r = norm(R) / bn
    r <= tol && return 1, r
    for k in 2:maxit
        X += M \ R
        R = B - K * X
        rnew = norm(R) / bn
        (!isfinite(rnew) || rnew > 1e6 * r) && return -1, rnew
        r = rnew
        r <= tol && return k, r
    end
    return maxit + 1, r
end

# Left-preconditioned GMRES on a single column: the best case for any Krylov
# method using this preconditioner.
function gmres_steps(K, b, M; tol=1e-10, maxit=60)
    n = length(b)
    Mb = M \ b
    bn = norm(Mb)
    bn == 0 && return 0
    V = zeros(n, maxit + 1)
    H = zeros(maxit + 1, maxit)
    V[:, 1] .= Mb ./ bn
    for j in 1:maxit
        w = M \ (K * V[:, j])
        for i in 1:j
            H[i, j] = dot(view(V, :, i), w)
            w .-= H[i, j] .* view(V, :, i)
        end
        H[j+1, j] = norm(w)
        e1 = zeros(j + 1); e1[1] = bn
        y = H[1:j+1, 1:j] \ e1
        norm(H[1:j+1, 1:j] * y - e1) / bn <= tol && return j
        H[j+1, j] < 1e-14 && return j
        V[:, j+1] .= w ./ H[j+1, j]
    end
    return maxit + 1
end

# rho(I - M^{-1} K) by power iteration -- the asymptotic refinement rate.
function reuse_spectral_radius(K, M; iters=40)
    n = size(K, 1)
    v = normalize!(ones(n) .+ 0.001 .* collect(1:n) ./ n)
    ρ = 0.0
    for _ in 1:iters
        w = v - M \ (K * v)
        ρ = norm(w)
        ρ < 1e-300 && return 0.0
        v = w ./ ρ
    end
    return ρ
end

anchors = [s.iter for s in snaps if s.iter % 10 == 0 && s.iter < snaps[end].iter]
rows = NamedTuple[]
for a in anchors
    ai = findfirst(s -> s.iter == a, snaps)
    M = lu(snaps[ai].K)
    K0, S0 = snaps[ai].K, snaps[ai].Sigma
    for s in snaps[ai:end]
        s.iter - a > 20 && break
        dK = norm(s.K - K0) / norm(K0)
        dS = norm(s.Sigma - S0) / max(norm(S0), eps())
        steps, res = refinement_steps(s.K, s.rhs, M)
        gm = gmres_steps(s.K, view(s.rhs, :, 1), M)
        ρ = reuse_spectral_radius(s.K, M)
        push!(rows, (anchor=a, iter=s.iter, lag=s.iter - a, dK_rel=dK, dSigma_rel=dS,
                     rho=ρ, refine_steps=steps, refine_res=res, gmres_steps=gm, mu=s.mu))
        @printf("anchor %3d  iter %3d  lag %2d  |dK|/|K| %.3e  |dSigma|/|Sigma| %.3e  rho %.3e  refine %3d  gmres %3d\n",
                a, s.iter, s.iter - a, dK, dS, ρ, steps, gm)
    end
end

open(out_csv, "w") do io
    println(io, "anchor,iter,lag,dK_rel,dSigma_rel,rho,refine_steps,refine_res,gmres_steps,mu")
    for r in rows
        @printf(io, "%d,%d,%d,%.6e,%.6e,%.6e,%d,%.6e,%d,%.6e\n",
                r.anchor, r.iter, r.lag, r.dK_rel, r.dSigma_rel, r.rho,
                r.refine_steps, r.refine_res, r.gmres_steps, r.mu)
    end
end
println("wrote $out_csv")
