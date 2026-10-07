# Check the block-tree elimination (tree_kkt.jl) against MUMPS's Schur
# complement on captured stage KKT systems, and time both informally.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/tree_kkt_check.jl <system> <capture files...>

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))      # control_layout, REPO
include(joinpath(@__DIR__, "tree_kkt.jl"))
include(joinpath(@__DIR__, "battery_schur_hook.jl"))      # mumps_battery_schur
using Statistics

system = ARGS[1]
data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T6_periodic.jls"))
idx, nu = control_layout(data)
lay = nothing; pattern = nothing
for file in ARGS[2:end]
    d = deserialize(file)
    K = d.K
    E_rhs = findall(i -> any(!iszero, @view d.rhs_multi[i, 2:end]), 1:size(K, 1))
    # The layout stores nzval positions, so it is rebuilt if the pattern
    # changes (old captures drop explicit zeros; production runs keep the
    # pattern fixed through the KKT pattern cache).
    t_lay = @elapsed if isnothing(pattern) || pattern != (K.colptr, K.rowval)
        global lay = tree_kkt_layout(data, idx, nu, K)
        global pattern = (copy(K.colptr), copy(K.rowval))
    end
    sort(lay.E) == E_rhs || error("battery rows differ from the capture's")
    tree_kkt_schur(lay, K)                                   # compile
    t_tree = @elapsed S_tree = tree_kkt_schur(lay, K)
    t_fac = @elapsed fac = tree_kkt_factor(lay, K)
    t_fc = @elapsed tree_kkt_fc_block(lay, fac)
    _SCHUR_STATE[] = nothing                                 # fresh MUMPS analysis per capture
    mumps_battery_schur(K, lay.E)
    t_mumps = @elapsed S_mumps = mumps_battery_schur(K, lay.E)  # analysis reused
    rel = norm(S_tree - S_mumps) / norm(S_mumps)
    @printf("TREE_KKT file=%s n=%d nodes=%d nE=%d nFc=%d max_own=%d max_iface=%d rel_diff_vs_mumps=%.3e tree_s=%.3f (factor %.3f, Fc block %.3f) mumps_refactor_s=%.3f layout_s=%.2f\n",
            basename(dirname(file)) * "/" * basename(file), size(K, 1), length(lay.own),
            length(lay.E), length(lay.Fc), maximum(length, lay.own), maximum(length, lay.iface),
            rel, t_tree, t_fac, t_fc, t_mumps, t_lay)
    flush(stdout)
end
