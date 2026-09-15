# What does d2Phi_t/dP_B2 actually look like?
#
# The exact-curvature reduced-space run converges but pays nB inner solves per
# stage per backward pass, and the battery-only proxy does not converge at all.
# So the question that decides scaling is whether the Hessian of the reduced
# stage value has structure cheap enough to exploit:
#
#   * LOW RANK  -> k << nB random probes recover it (k solves, not nB).
#   * SPARSE    -> graph-colouring FD recovers it in (colours) solves.
#   * DIAGONALLY DOMINANT -> a diagonal proxy might suffice where the
#     battery-only proxy failed.
#
# Reported: SVD spectrum and effective rank at a relative-Frobenius threshold,
# off-diagonal mass, and how well a rank-k or diagonal truncation reproduces it.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#     ddp/examples/power_system/probe_phi_hessian_structure.jl [system] [T] [t]

using LinearAlgebra
using Printf
using Serialization

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
tsel   = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T$(T).jls"))
Bset = data[:Bset]; nB = length(Bset)

# evaluate at the converged reduced-space dispatch when available, else at zero
solfile = joinpath(REPO, "ddp", "results", "reduced_space",
                   "reduced_filterddp_$(system)_T$(T)_exact.jls")
have_sol = isfile(solfile)
sol = have_sol ? deserialize(solfile) : nothing
times = tsel > 0 ? [tsel] : collect(1:T)

outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
csv = open(joinpath(outdir, "phi_hessian_structure_$(system)_T$(T).csv"), "w")
println(csv, "system,horizon,time_index,at_solution,nB,fro_norm,diag_mass_frac," *
             "rank_1pct,rank_0p1pct,rank_0p01pct,diag_rel_err,rank10_rel_err," *
             "rank25_rel_err,cond_est,min_eig,max_eig")

for t in times
    pb = have_sol ? sol[:u][t][1:nB] : zeros(nB)
    r0 = inner_opf(data, t, pb)
    r0.feasible || (@printf("t=%d: base point infeasible, skipping\n", t); continue)
    g0 = r0.lambda_bal
    h = 1e-6
    H = zeros(nB, nB)
    pert = copy(pb)
    tstart = time()
    for b in 1:nB
        pert[b] = pb[b] + h
        rb = inner_opf(data, t, pert)
        H[:, b] = rb.feasible ? (rb.lambda_bal .- g0) ./ h : zeros(nB)
        pert[b] = pb[b]
    end
    H = 0.5 .* (H .+ H')
    el = time() - tstart

    fro = norm(H)
    diag_mass = norm(Diagonal(H)) / max(fro, eps())
    F = svd(H)
    s = F.S
    cum = sqrt.(reverse(cumsum(reverse(s .^ 2))))          # tail Frobenius norm
    effrank(tol) = begin
        k = findfirst(i -> cum[i] <= tol * fro, 1:length(cum))
        k === nothing ? length(s) : k - 1
    end
    r1, r01, r001 = effrank(1e-2), effrank(1e-3), effrank(1e-4)

    trunc_err(k) = k >= length(s) ? 0.0 :
        norm(F.U[:, 1:k] * Diagonal(s[1:k]) * F.Vt[1:k, :] .- H) / max(fro, eps())
    diag_err = norm(Diagonal(H) .- H) / max(fro, eps())

    ev = eigvals(Symmetric(H))
    @printf(csv, "%s,%d,%d,%s,%d,%.6e,%.4f,%d,%d,%d,%.4e,%.4e,%.4e,%.4e,%.6e,%.6e\n",
            system, T, t, have_sol, nB, fro, diag_mass, r1, r01, r001,
            diag_err, trunc_err(10), trunc_err(25),
            s[1] / max(s[end], eps()), minimum(ev), maximum(ev))
    @printf("t=%d (%.1fs, %d solves): ||H||_F=%.4g  diag mass=%.1f%%\n",
            t, el, nB, fro, 100 * diag_mass)
    @printf("   effective rank  1%%:%d   0.1%%:%d   0.01%%:%d   (of %d)\n", r1, r01, r001, nB)
    @printf("   rel err   diagonal:%.3e   rank-10:%.3e   rank-25:%.3e\n",
            diag_err, trunc_err(10), trunc_err(25))
    @printf("   eig range [%.4g, %.4g]   top-10 sv: %s\n", minimum(ev), maximum(ev),
            join((@sprintf("%.3g", x) for x in s[1:min(10, end)]), " "))
    flush(stdout)
end
close(csv)
println("wrote phi_hessian_structure_$(system)_T$(T).csv")
