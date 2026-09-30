# Check that FILTERDDP_TYPED_EQUATIONS=1 changes only the speed of the stage
# residual callback, not its output: build the model with and without it and
# compare c(x, u) BITWISE at every stage, on random points and on points with
# exact zeros (which exercise signed-zero handling of the zero terms).
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/check_typed_equations.jl <system> <T>

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))
using Random, Statistics

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T$(T)_periodic.jls"))
data[:C_B] = battery_cb(system, "system")
gammaT = gamma_terminal(system)

ENV["FILTERDDP_TYPED_EQUATIONS"] = "0"
ocp_ref, idx, nx, nu, nc = build_model(data; gamma=gammaT)
ENV["FILTERDDP_TYPED_EQUATIONS"] = "1"
ocp_new, _, _, _, _ = build_model(data; gamma=gammaT)

rng = MersenneTwister(20260930)
points = Any[]
for k in 1:4
    x = randn(rng, nx); u = randn(rng, nu)
    k == 3 && (u .= 0.0; x .= 0.0)                      # all zeros
    k == 4 && (u[rand(rng, 1:nu, nu ÷ 3)] .= 0.0)      # scattered zeros
    push!(points, (x, u))
end
bits(v) = reinterpret(UInt64, Vector{Float64}(v))
mismatch = 0
for t in 1:T, (x, u) in points
    a = bits(ocp_ref.stage_constraints[t].c(x, u))
    b = bits(ocp_new.stage_constraints[t].c(x, u))
    length(a) == length(b) || error("length differs at stage $t")
    global mismatch += count(a .!= b)
end
x, u = points[1]
time_ms(f) = (f(); median([(s = time_ns(); f(); (time_ns() - s) / 1e6) for _ in 1:7]))
t_ref = time_ms(() -> ocp_ref.stage_constraints[1].c(x, u))
t_new = time_ms(() -> ocp_new.stage_constraints[1].c(x, u))
@printf("TYPED_EQUATIONS_CHECK system=%s T=%d nc=%d stages=%d points=%d bitwise_mismatches=%d untyped_ms=%.3f typed_ms=%.3f speedup=%.0fx\n",
        system, T, nc, T, length(points), mismatch, t_ref, t_new, t_ref / t_new)
mismatch == 0 || exit(1)
