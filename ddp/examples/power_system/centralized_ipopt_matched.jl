# Centralized JuMP/Ipopt solve of EXACTLY the instance FilterDDP solves.
#
# The centralized timings in the TPEC paper (ddp/results/centralized_ipopt/,
# including the large10k T = 48 "knee" at 8262 s) were measured on the tADMM
# instance family: C_B ~ 8.8e-8 and the pre-2026-09-15 price sampling. The
# FilterDDP diagonal-Hessian results use the periodic profile and C_B = 1e-3.
# Those are different optimization problems, so neither set of timings can be
# raced against the other. This driver closes the gap: it reads the same
# network_data_<system>_T<T><tag>.jls the FilterDDP driver reads, and applies
# the same C_B override, selected by the SAME environment variables
# (REDUCED_PROFILE, REDUCED_CB) as ieee123c_filterddp.jl, so the two cannot
# drift apart.
#
# The model is ieee123c_ipopt_matrix_diagnostic.jl's -- a line-for-line
# transcription of the FilterDDP driver's formulation (cost + C_B battery term,
# BFM with the SOC relaxation as an inequality, free terminal SOC). Ipopt
# options mirror envs/tadmm/root_level/run_bf.jl, which produced the paper's
# centralized sweep: defaults plus max_iter = 5000.
#
#   REDUCED_PROFILE=periodic REDUCED_CB=1e-3 \
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/centralized_ipopt_matched.jl <system> <T> [ipopt_log]
#
# Prints one parseable CENTRAL_IPOPT line.

using Ipopt
using JuMP
using LinearAlgebra
using Printf
using Serialization

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MOI = JuMP.MOI
include(joinpath(@__DIR__, "terminal_soc_penalty.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                    "network_data_$(system)_T$(T)$(ptag).jls")
ipopt_log = length(ARGS) >= 3 ? ARGS[3] :
    joinpath(REPO, "ddp", "results", "centralized_ipopt_matched", "$(system)_T$(T)$(ptag)_ipopt.log")
mkpath(dirname(ipopt_log))

wall_start = time()
data = deserialize(datafile)
haskey(ENV, "REDUCED_CB") && (data[:C_B] = parse(Float64, ENV["REDUCED_CB"]))

Nset, Lset, Bset, Dset = data[:Nset], data[:Lset], data[:Bset], data[:Dset]
Tset = 1:T
root, nonroot = data[:substationBus], data[:Nm1set]
dt, pbase = data[:delta_t_h], data[:kVA_B]

build_start = time()
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 5)
set_optimizer_attribute(model, "max_iter", 5000)
set_optimizer_attribute(model, "print_timing_statistics", "yes")
set_optimizer_attribute(model, "output_file", ipopt_log)

@variable(model, P_Subs[Tset] >= 0)
@variable(model, Q_Subs[Tset])
@variable(model, P[Lset, Tset])
@variable(model, Q[Lset, Tset])
@variable(model, v[Nset, Tset])
@variable(model, ell[Lset, Tset] >= 0)
@variable(model, P_B[Bset, Tset])
@variable(model, B[Bset, Tset])
@variable(model, q_D[Dset, Tset])

# Soft terminal SOC, same switch and same per-system gamma as the FilterDDP
# driver (terminal_soc_penalty.jl). B[j,t] is the energy at the END of period t,
# so B[j,T] is exactly the FilterDDP driver's B^T = x_T - dt*u_T[pb].
gammaT = terminal_soc_soft() ? gamma_terminal(system) : 0.0
@objective(model, Min,
    sum(data[:LoadShapeCost][t] * pbase * dt * P_Subs[t] for t in Tset) +
    sum(data[:C_B] * pbase^2 * dt * P_B[j,t]^2 for j in Bset, t in Tset) +
    gammaT * sum((B[j,T] - data[:B0_pu][j])^2 for j in Bset))

for t in Tset
    @constraint(model, P_Subs[t] == sum(P[e,t] for e in data[:L1set]))
    @constraint(model, Q_Subs[t] == sum(Q[e,t] for e in data[:L1set]))
    for j in nonroot
        incoming = (data[:parent][j], j)
        outgoing_p = sum(P[(j,k),t] for k in data[:children][j]; init=0.0)
        outgoing_q = sum(Q[(j,k),t] for k in data[:children][j]; init=0.0)
        pL = j in data[:NLset] ? data[:p_L_pu][j,t] : 0.0
        qL = j in data[:NLset] ? data[:q_L_pu][j,t] : 0.0
        pD = j in Dset ? data[:p_D_pu][j,t] : 0.0
        pb = j in Bset ? P_B[j,t] : 0.0
        qd = j in Dset ? q_D[j,t] : 0.0
        @constraint(model, outgoing_p - P[incoming,t] +
            data[:rdict_pu][incoming] * ell[incoming,t] == pb + pD - pL)
        @constraint(model, outgoing_q - Q[incoming,t] +
            data[:xdict_pu][incoming] * ell[incoming,t] == qd - qL)
    end
    for e in Lset
        i, j = e
        r, x = data[:rdict_pu][e], data[:xdict_pu][e]
        @constraint(model, v[j,t] == v[i,t] - 2(r*P[e,t] + x*Q[e,t]) +
            (r^2+x^2)*ell[e,t])
        @constraint(model, P[e,t]^2 + Q[e,t]^2 <= v[i,t]*ell[e,t])
    end
    @constraint(model, v[root,t] == 1.05^2)
    for j in Nset
        set_lower_bound(v[j,t], data[:Vminpu][j]^2)
        set_upper_bound(v[j,t], data[:Vmaxpu][j]^2)
    end
    for j in Dset
        qmax = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
        set_lower_bound(q_D[j,t], -qmax)
        set_upper_bound(q_D[j,t], qmax)
    end
    for j in Bset
        set_lower_bound(P_B[j,t], -data[:P_B_R_pu][j])
        set_upper_bound(P_B[j,t], data[:P_B_R_pu][j])
        set_lower_bound(B[j,t], data[:soc_min][j] * data[:B_R_pu][j])
        set_upper_bound(B[j,t], data[:soc_max][j] * data[:B_R_pu][j])
        if t == 1
            @constraint(model, B[j,t] == data[:B0_pu][j] - dt*P_B[j,t])
        else
            @constraint(model, B[j,t] == B[j,t-1] - dt*P_B[j,t])
        end
    end
end
build_s = time() - build_start

optimize!(model)
wall_s = time() - wall_start

status = termination_status(model)
obj = has_values(model) ? objective_value(model) : NaN
iters = try MOI.get(model, MOI.BarrierIterations()) catch; -1 end
solve_s = try solve_time(model) catch; NaN end
@printf("CENTRAL_IPOPT system=%s T=%d profile=%s C_B=%.6e gamma=%.6e status=%s iterations=%d objective=%.12f solve_time_s=%.3f build_s=%.3f wall_s=%.3f variables=%d solver=%s\n",
        system, T, isempty(ptag) ? "default" : ptag[2:end], data[:C_B], gammaT, string(status),
        iters, obj, solve_s, build_s, wall_s, num_variables(model),
        replace(string(MOI.get(model, MOI.SolverVersion())), ' ' => '_'))

# How the batteries are used and what each objective term contributes. Printed
# after the timed solve, so it does not affect any timing.
if has_values(model)
    pb = value.(P_B); bb = value.(B)
    util = [abs(pb[j,t]) / max(data[:P_B_R_pu][j], eps()) for j in Bset, t in Tset]
    energy_term = sum(data[:LoadShapeCost][t] * pbase * dt * value(P_Subs[t]) for t in Tset)
    cb_term = sum(data[:C_B] * pbase^2 * dt * pb[j,t]^2 for j in Bset, t in Tset)
    term_term = gammaT * sum((bb[j,T] - data[:B0_pu][j])^2 for j in Bset)
    # Share of each battery's usable energy window (soc_min..soc_max of B_R)
    # that its SOC trajectory, including B0, actually spans; averaged over batteries.
    window = [begin
                  traj = vcat(data[:B0_pu][j], [bb[j,t] for t in Tset])
                  span = (data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j]
                  (maximum(traj) - minimum(traj)) / max(span, eps())
              end for j in Bset]
    @printf("CENTRAL_IPOPT_BATTERY at_power_limit=%.3f mean_utilization=%.3f energy_window_used=%.3f throughput_pu_h=%.6e energy_term=%.6e cb_term=%.6e terminal_term=%.6e\n",
            count(>=(0.99), util) / length(util), sum(util) / length(util),
            sum(window) / length(window), sum(abs.(pb)) * dt, energy_term, cb_term, term_term)
end
