# How far does an Ipopt-style constraint-block regularization move the stage
# KKT solution? K(δ_c) = [H cu'; cu -δ_c I] against the exact K = [H cu'; cu 0],
# both solved with UMFPACK on the captured stage-1 systems (diagonal Hessian).
# Reports the relative change of the whole (n_x+1)-column solution, of its
# primal rows (the control directions and gains FilterDDP uses) and of its
# multiplier rows, plus the spread of the Hessian diagonal that explains it.
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_dual_reg_sensitivity.jl \
#         <out.csv> <system> [system ...]

using Serialization, SparseArrays, LinearAlgebra, Printf

out = ARGS[1]
open(out, "w") do io
    println(io, "system,nu,nc,H_diag_min,H_diag_median,H_diag_max,dual_reg,dev_all,dev_primal,dev_multiplier,dev_feedforward_primal")
    for sys in ARGS[2:end]
        c = deserialize("ddp/results/kkt_ordering/captures/kkt_$(sys)_T3_diag_iter20_stage1.jls")
        K0 = SparseMatrixCSC{Float64, Int64}(c.K); B = Matrix{Float64}(c.rhs)
        n = size(K0, 1); nu = c.nu
        d = sort(diag(K0)[1:nu])
        X0 = lu(K0) \ B
        @printf("%s nu=%d nc=%d  H diag %.1e .. %.1e (median %.1e)  exact residual %.1e\n",
                sys, nu, n - nu, d[1], d[end], d[div(end, 2)], norm(K0 * X0 - B) / norm(B))
        rel(A, R) = norm(A - R) / norm(R)
        for δ in (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4)
            X = lu(K0 - δ * sparse(Diagonal([zeros(nu); ones(n - nu)]))) \ B
            r = (rel(X, X0), rel(X[1:nu, :], X0[1:nu, :]), rel(X[nu+1:end, :], X0[nu+1:end, :]),
                 rel(X[1:nu, 1], X0[1:nu, 1]))
            @printf("  δ_c=%.0e  all %.1e  primal %.1e  multiplier %.1e  feedforward primal %.1e\n", δ, r...)
            @printf(io, "%s,%d,%d,%.3e,%.3e,%.3e,%.0e,%.3e,%.3e,%.3e,%.3e\n",
                    sys, nu, n - nu, d[1], d[div(end, 2)], d[end], δ, r...)
        end
        flush(io)
    end
end
