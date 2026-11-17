# The stage constraint Jacobian cu against central differences of c: the
# constraints are at most quadratic in u, so (c(u + h d) - c(u - h d)) / 2h is
# cu * d up to rounding for any direction d. Run after a model change (here:
# the capacitor term q = qC * v of parse_capacitors!).
#   julia --project=problems/mpopf/spar_ddp .../check_jacobian.jl <system> <T>
include(joinpath(@__DIR__, "..", "..", "examples", "power_system", "ieee123c_filterddp.jl"))
using Random
system = ARGS[1]; T = parse(Int, ARGS[2])
data = deserialize(joinpath(REPO, "problems", "mpopf", "spar_ddp", "results", "network_filterddp",
                            "network_data_$(system)_T$(T)_periodic.jls"))
data[:C_B] = battery_cb(system, "system")
ENV["FILTERDDP_TYPED_EQUATIONS"] = "1"
ocp, idx, nx, nu, nc = build_model(data; gamma=gamma_terminal(system))
rng = MersenneTwister(1117); worst = 0.0
for t in 1:T, k in 1:3
    x = randn(rng, nx); u = randn(rng, nu); d = randn(rng, nu); h = 1e-3
    con = ocp.stage_constraints[t]
    fd = (con.c(x, u .+ h .* d) .- con.c(x, u .- h .* d)) ./ (2h)
    J = con.cu(x, u) * d
    global worst = max(worst, maximum(abs.(fd .- J)) / max(1.0, maximum(abs.(J))))
end
println("JACOBIAN_CHECK system=$system T=$T capacitors=$(length(get(data, :Cset, Int[]))) worst_relative_difference=$worst")
worst < 1e-8 || error("Jacobian and residual disagree")
