# Compare one converged FilterDDP stage with a JuMP/Ipopt stage carrying the
# identical incoming quadratic value-function message. Selected state
# sensitivities are obtained by central finite differences of complete Ipopt
# stage solves and compared with FilterDDP's analytic beta/omega feedback.

using Ipopt
using JuMP
using LinearAlgebra
using Printf
using Serialization

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))

const MOI = JuMP.MOI

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
stage = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
h = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e-5

data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
    "network_data_$(system)_T$(T).jls"))
capture_path = joinpath(REPO, "ddp", "results", "network_filterddp",
    "$(system)_T$(T)_stage$(stage)_oracle_capture.jls")
capture = deserialize(capture_path)
capture.stage == stage || error("capture stage does not match requested stage")

ocp, idx, nx, nu, nc = build_model(data)
obj = ocp.stage_objectives[stage]
con = ocp.stage_constraints[stage]
dyn = ocp.dynamics
cl = ocp.control_limits
Vx = capture.future_Vx
Vxx = Symmetric((capture.future_Vxx + capture.future_Vxx') / 2)
f_reference = capture.next_state

function solve_oracle(x_parameter; start=capture.u)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", 1e-10)
    set_optimizer_attribute(model, "constr_viol_tol", 1e-10)
    set_optimizer_attribute(model, "acceptable_tol", 1e-9)
    set_optimizer_attribute(model, "max_iter", 500)

    @variable(model, u[1:nu])
    for k in 1:nu
        isfinite(cl.l[k]) && set_lower_bound(u[k], cl.l[k])
        isfinite(cl.u[k]) && set_upper_bound(u[k], cl.u[k])
        set_start_value(u[k], start[k])
    end

    delta_f = dyn.f(x_parameter, u) .- f_reference
    stage_cost = obj.l(x_parameter, u)[1]
    future_cost = sum(Vx[b] * delta_f[b] for b in 1:nx) +
        0.5 * sum(Vxx[b, j] * delta_f[b] * delta_f[j]
                  for b in 1:nx, j in 1:nx)
    @objective(model, Min, stage_cost + future_cost)

    equations = con.c(x_parameter, u)
    constraints = [@constraint(model, equations[r] == 0) for r in 1:nc]
    optimize!(model)
    status = termination_status(model)
    status in (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED) ||
        error("Ipopt stage solve failed with status $status")
    return (
        u=value.(u),
        phi=dual.(constraints),
        objective=objective_value(model),
        iterations=MOI.get(model, MOI.BarrierIterations()),
        residual=maximum(abs, con.c(x_parameter, value.(u))),
    )
end

base = solve_oracle(capture.x)
# JuMP/MOI's equality-dual convention is opposite to the multiplier sign used
# in FilterDDP's stationarity expression `Qu + cu' * phi`.
phi_sign = -1.0
base_phi = phi_sign .* base.phi

directions = Pair{String,Vector{Float64}}[]
push!(directions, "battery_1" => [b == 1 ? 1.0 : 0.0 for b in 1:nx])
middle = cld(nx, 2)
push!(directions, "battery_middle" => [b == middle ? 1.0 : 0.0 for b in 1:nx])
push!(directions, "aggregate_unit" => fill(inv(sqrt(nx)), nx))

rows = NamedTuple[]
for (name, direction) in directions
    plus = solve_oracle(capture.x .+ h .* direction; start=base.u)
    minus = solve_oracle(capture.x .- h .* direction; start=base.u)
    du = (plus.u - minus.u) / (2h)
    dphi = phi_sign .* (plus.phi - minus.phi) / (2h)
    beta_direction = capture.beta * direction
    omega_direction = capture.omega * direction
    push!(rows, (
        direction=name,
        beta_abs=norm(du - beta_direction),
        beta_rel=norm(du - beta_direction) / max(norm(beta_direction), eps()),
        beta_inf=norm(du - beta_direction, Inf),
        omega_abs=norm(dphi - omega_direction),
        omega_rel=norm(dphi - omega_direction) / max(norm(omega_direction), eps()),
        omega_inf=norm(dphi - omega_direction, Inf),
        plus_residual=plus.residual,
        minus_residual=minus.residual,
    ))
end

@printf("ORACLE system=%s T=%d stage=%d barrier_mu=%.3e\n",
    system, T, stage, capture.barrier_mu)
@printf("ORACLE base_iterations=%d base_residual=%.3e primal_inf_diff=%.3e phi_inf_diff=%.3e phi_sign=%+.0f alpha_inf=%.3e\n",
    base.iterations, base.residual, norm(base.u-capture.u, Inf),
    norm(base_phi-capture.phi, Inf), phi_sign, norm(capture.alpha, Inf))
@printf("ORACLE displacement_vs_alpha_inf=%.3e capture_phi_inf=%.3e ipopt_phi_inf=%.3e\n",
    norm((base.u-capture.u)-capture.alpha, Inf), norm(capture.phi, Inf),
    norm(base_phi, Inf))
for row in rows
    @printf("ORACLE direction=%s beta_rel=%.3e beta_inf=%.3e omega_rel=%.3e omega_inf=%.3e\n",
        row.direction, row.beta_rel, row.beta_inf, row.omega_rel, row.omega_inf)
end

output = joinpath(REPO, "ddp", "results", "network_filterddp",
    "$(system)_T$(T)_stage$(stage)_ipopt_oracle.csv")
open(output, "w") do io
    println(io, "finite_difference_h,direction,beta_abs,beta_rel,beta_inf,omega_abs,omega_rel,omega_inf,plus_residual,minus_residual")
    for row in rows
        @printf(io, "%.12e,%s,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
            h, row.direction, row.beta_abs, row.beta_rel, row.beta_inf,
            row.omega_abs, row.omega_rel, row.omega_inf,
            row.plus_residual, row.minus_residual)
    end
end
println("ORACLE wrote=$output")
