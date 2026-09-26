# Why a constraint-block regularization -δ_c I hurts UMFPACK: it fills the
# structurally zero diagonal, and UMFPACK's automatic strategy then picks the
# symmetric strategy instead of the unsymmetric one. Factors each captured
# stage KKT (diagonal Hessian) with δ_c = 0 and 1e-8 under the automatic,
# forced-unsymmetric and forced-symmetric strategies; the automatic choice is
# identified by matching factor size. Wide-solve times are single, warm runs.
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_dual_reg_umfpack_strategy.jl <system> ...

using Serialization, SparseArrays, LinearAlgebra, Printf
const LS = SparseArrays.LibSuiteSparse
for sys in ARGS
    c = deserialize("ddp/results/kkt_ordering/captures/kkt_$(sys)_T3_diag_iter20_stage1.jls")
    K0 = SparseMatrixCSC{Float64, Int64}(c.K); n = size(K0, 1); nu = c.nu; B = Matrix{Float64}(c.rhs)
    for δ in (0.0, 1e-8), (st, name) in ((0, "auto"), (1, "unsymmetric"), (3, "symmetric"))
        K = δ == 0 ? K0 : K0 - δ * sparse(Diagonal([zeros(nu); ones(n - nu)]))
        ctl = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int64)
        ctl[LS.UMFPACK_STRATEGY + 1] = st
        F = lu(K; control=ctl)
        X = copy(B); ldiv!(F, X); X = copy(B); t = @elapsed ldiv!(F, X)
        @printf("%s δ_c=%.0e strategy=%-11s nnz(L)+nnz(U)=%8d  wide solve %.3f s  residual %.1e\n",
                sys, δ, name, nnz(F.L) + nnz(F.U), t, norm(K * X - B) / norm(B))
    end
end
