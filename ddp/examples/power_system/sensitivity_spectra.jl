# Singular-value spectra of the objects in one FilterDDP stage (agenda of
# 2026-10-07): which of them, if any, is dominated by a few directions?
#
#   V      incoming value curvature V_xx (n_B x n_B), and its off-diagonal part
#   beta   the full state-sensitivity map of the controls (n_u x n_B)
#   beta_B its battery-power rows (n_B x n_B), the part the value update uses
#   Vinc   the value-curvature increment beta_B' B + omega_E' c_x,E
#
# For each: the number of singular values needed to capture the matrix to
# 10%, 1% and 0.1% in Frobenius norm, out of its full rank. Input: stage
# captures written by FILTERDDP_CAPTURE_KKT (they carry K, rhs, beta, omega
# and the incoming V_xx).
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/sensitivity_spectra.jl <capture files...>

using LinearAlgebra, Printf, Serialization, SparseArrays

function ranks(M)
    sv = svdvals(Matrix(M))
    total = sum(abs2, sv)
    total == 0 && return (0, 0, 0, length(sv), 0.0)
    tail = reverse(cumsum(reverse(abs2.(sv))))          # tail[r] = sum_{i >= r} sv_i^2
    need(tol) = (r = findfirst(i -> (i == length(sv) ? 0.0 : tail[i+1]) / total <= tol^2, 1:length(sv)); r)
    return (need(0.1), need(0.01), need(0.001), length(sv), sv[1] / max(sv[end], floatmin()))
end

report(name, M) = (r = ranks(M); @printf("  %-34s %5d x %-5d rank for 10%%/1%%/0.1%%: %4d %4d %4d of %4d   sv1/svn=%.1e\n",
                                         name, size(M, 1), size(M, 2), r[1], r[2], r[3], r[4], r[5]))

for file in ARGS
    d = deserialize(file)
    nx, nu = d.nx, d.nu
    V = d.future_Vxx
    R = d.rhs
    E = findall(i -> any(!iszero, @view R[i, 2:end]), 1:size(R, 1))
    pb = filter(<=(nu), E); en = filter(>(nu), E)
    beta = d.beta; omega = d.omega
    @printf("SPECTRA %s  stage=%d mu=%.1e nx=%d nu=%d\n", basename(file), d.stage, d.barrier_mu, nx, nu)
    offd = V - Diagonal(diag(V))
    @printf("  V_xx off-diagonal share of Frobenius norm: %.3f   (symmetric: %.1e)\n",
            norm(offd) / max(norm(V), floatmin()), norm(V - V') / max(norm(V), floatmin()))
    report("V_xx (incoming)", V)
    report("V_xx minus its diagonal", offd)
    report("beta (all controls)", beta)
    report("beta_B (battery-power rows)", beta[pb, :])
    report("beta_B minus its diagonal", beta[pb, :] - Diagonal(diag(beta[pb, :])))
    Vinc = beta[pb, :]' * (-R[pb, 2:end]) + omega[en .- nu, :]' * (-R[en, 2:end])
    report("value-curvature increment", Vinc)
    report("increment minus its diagonal", Vinc - Diagonal(diag(Vinc)))
    flush(stdout)
end
