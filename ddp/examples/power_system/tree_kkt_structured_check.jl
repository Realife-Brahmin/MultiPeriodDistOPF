# Check the feeder-structured Schur complement and Woodbury solve
# (tree_kkt.jl) against MUMPS and against dense solves, on production-config
# stage captures (FILTERDDP_CAPTURE_KKT with the diagonal Hessian).
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/tree_kkt_structured_check.jl <system> <capture>

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))
include(joinpath(@__DIR__, "tree_kkt.jl"))
include(joinpath(@__DIR__, "battery_schur_hook.jl"))

system, file = ARGS[1], ARGS[2]
data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T6_periodic.jls"))
idx, nu = control_layout(data)
d = deserialize(file)
K, rhs = d.K, d.rhs
lay = tree_kkt_layout(data, idx, nu, K)
E = lay.E
RE = rhs[E, 2:end]                               # battery rows of the feedback right-hand side
# Reference: MUMPS Schur complement and the full UMFPACK solve
S_m = mumps_battery_schur(K, E)
X_full = (lu(K) \ rhs)[E, 2:end]
# Structured, with a dense check of S
fac = tree_kkt_factor(lay, K)
st = tree_kkt_structured(lay, fac, K; dense_check=true)
S_t = tree_kkt_schur_dense(st)
X_t = tree_kkt_schur_solve(st, RE)
@printf("STRUCT system=%s feeders=%d groups=%d rel_S_vs_mumps=%.3e rel_X_vs_full_solve=%.3e\n",
        system, length(lay.children[lay.root]), length(st.groups),
        norm(S_t - S_m) / norm(S_m), norm(X_t - X_full) / norm(X_full))
# Full solver against UMFPACK: one column (the feedforward) and the battery rows
sol = tree_kkt_solver(lay, K)
x1 = copy(rhs[:, 1]); ldiv!(sol, x1)
x1_ref = lu(K) \ rhs[:, 1]
@printf("SOLVER rel_full_solve_vs_umfpack=%.3e rel_battery_rows=%.3e
",
        norm(x1 - x1_ref) / norm(x1_ref),
        norm(tree_kkt_battery_rows(sol, RE) - X_full) / norm(X_full))
for _ in 1:2
    t_all = @elapsed (sol = tree_kkt_solver(lay, K))
    t_one = @elapsed ldiv!(sol, copy(rhs[:, 1]))
    t_rows = @elapsed tree_kkt_battery_rows(sol, RE)
    @printf("SOLVER_TIMING build=%.3f one_column_solve=%.3f battery_rows=%.3f
", t_all, t_one, t_rows)
end
# Timing (warm): factor, structure, Woodbury solve; against MUMPS refactor and UMFPACK factor + full solve
for _ in 1:2
    t_fac = @elapsed (fac = tree_kkt_factor(lay, K))
    t_st = @elapsed (st = tree_kkt_structured(lay, fac, K))
    t_sol = @elapsed tree_kkt_schur_solve(st, RE)
    t_m = @elapsed mumps_battery_schur(K, E)
    t_u = @elapsed (F = lu(K)); t_us = @elapsed (F \ rhs)
    @printf("TIMING tree_factor=%.3f structure=%.3f woodbury=%.3f total=%.3f | mumps_refactor=%.3f | umfpack_lu=%.3f umfpack_full_solve_1col_by_col=%.3f\n",
            t_fac, t_st, t_sol, t_fac + t_st + t_sol, t_m, t_u, t_us)
end
