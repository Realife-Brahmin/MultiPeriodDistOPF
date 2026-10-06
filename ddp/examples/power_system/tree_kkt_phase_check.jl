# Which part of a stage's tree solve can be done before the backward sweep.
#
# The value function of stage t+1 enters the stage-t KKT matrix only on the
# battery-power diagonal, and the right-hand side only on the battery-power
# rows. Both lie in E, the block the tree solver keeps (tree_kkt.jl). This
# script checks, on a captured production stage matrix, that the network
# factorization, the Schur preparation (_tree_kkt_prepare) and the first
# network solve of the feedforward column do not read E at all: every entry of
# K_EE, and every E and energy-slack entry of the right-hand side, is replaced
# by NaN and the results must be identical, bit for bit. That is the premise of
# the one-worker-per-period accounting (FILTERDDP_PARSIM, backward_pass.jl).
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/tree_kkt_phase_check.jl <system> <capture>

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))
include(joinpath(@__DIR__, "tree_kkt.jl"))

system, file = ARGS[1], ARGS[2]
data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T6_periodic.jls"))
idx, nu = control_layout(data)
d = deserialize(file)
K, rhs = d.K, d.rhs
lay = tree_kkt_layout(data, idx, nu, K)
stat = tree_kkt_static(lay, K)
E = lay.E; inE = falses(size(K, 1)); inE[E] .= true
issubset(idx.pb, E) || error("battery powers are not all in the kept block")

# K with K_EE destroyed
Kp = copy(K)
rows = rowvals(Kp); vals = nonzeros(Kp); n_poisoned = 0
for col in axes(Kp, 2), p in nzrange(Kp, col)
    if inE[col] && inE[rows[p]]
        vals[p] = NaN; global n_poisoned += 1
    end
end

same(a, b) = isequal(a, b)
function same_blocks(A, B)
    for i in eachindex(A)
        isassigned(A, i) == isassigned(B, i) || return false
        isassigned(A, i) || continue
        isequal(A[i], B[i]) || return false
    end
    return true
end

fac = tree_kkt_factor(lay, K);  Bd, Up = _tree_kkt_prepare(lay, stat, fac, K)
facp = tree_kkt_factor(lay, Kp); Bdp, Upp = _tree_kkt_prepare(lay, stat, facp, Kp)
ok_factor = all(same(fac.F[j].A, facp.F[j].A) && same(fac.F[j].perm, facp.F[j].perm) for j in eachindex(fac.F)) &&
    same_blocks(fac.W, facp.W) && same_blocks(fac.Aio, facp.Aio) && same(fac.Mroot, facp.Mroot)
ok_prepare = same_blocks(Bd, Bdp) && same(Up, Upp)
finite = all(isfinite, Up) && all(all(isfinite, fac.F[j].A) for j in eachindex(fac.F))

# first network solve of the feedforward column, right-hand side destroyed on E and the slacks
sol = tree_kkt_solver(lay, stat, K)
b = copy(rhs[:, 1]); bp = copy(b); bp[E] .= NaN; bp[lay.iso] .= NaN
y = zeros(lay.n); yp = zeros(lay.n)
tree_kkt_network_solve!(y, sol, b); tree_kkt_network_solve!(yp, sol, bp)
net = reduce(vcat, lay.own)
ok_solve = same(y[net], yp[net]) && all(isfinite, y[net])

# and the finish step does read what was destroyed (LAPACK refuses the NaNs, or they propagate)
finish_reads_E = try
    st = _tree_kkt_finish(lay, stat, facp, Kp, Bdp, Upp, false, nothing)
    any(isnan, st.DUp) || any(F -> any(isnan, F.factors), st.Dfac)
catch err
    err isa ArgumentError || rethrow()
    true
end

@printf("PHASE_CHECK system=%s n=%d nE=%d K_EE_entries_destroyed=%d factor_identical=%s prepare_identical=%s first_network_solve_identical=%s finite=%s finish_reads_K_EE=%s\n",
        system, size(K, 1), length(E), n_poisoned, ok_factor, ok_prepare, ok_solve, finite, finish_reads_E)
(ok_factor && ok_prepare && ok_solve && finite && finish_reads_E) || error("phase check failed")
