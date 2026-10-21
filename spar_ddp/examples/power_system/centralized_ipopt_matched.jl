# Centralized JuMP/Ipopt solve of EXACTLY the instance FilterDDP solves.
#
# The centralized timings in the TPEC paper (spar_ddp/results/centralized_ipopt/,
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
#   julia --startup-file=no --project=envs/spar_ddp \
#         spar_ddp/examples/power_system/centralized_ipopt_matched.jl <system> <T> [ipopt_log]
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
include(joinpath(@__DIR__, "voltage_screening.jl"))
include(joinpath(@__DIR__, "loadflow_start.jl"))
include(joinpath(@__DIR__, "control_layout.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
datafile = joinpath(REPO, "spar_ddp", "results", "network_filterddp",
                    "network_data_$(system)_T$(T)$(ptag).jls")
ipopt_log = length(ARGS) >= 3 ? ARGS[3] :
    joinpath(REPO, "spar_ddp", "results", "centralized_ipopt_matched", "$(system)_T$(T)$(ptag)_ipopt.log")
mkpath(dirname(ipopt_log))

wall_start = time()
data = deserialize(datafile)
haskey(ENV, "REDUCED_CB") && (data[:C_B] = battery_cb(system, ENV["REDUCED_CB"]))  # a number, or "system"

Nset, Lset, Bset, Dset = data[:Nset], data[:Lset], data[:Bset], data[:Dset]
Tset = 1:T
root, nonroot = data[:substationBus], data[:Nm1set]
dt, pbase = data[:delta_t_h], data[:kVA_B]

build_start = time()
# CENTRAL_SOLVER=gurobi solves the same JuMP model with Gurobi instead (the
# model is a convex QCQP: quadratic cost, rotated second-order cones). Run it
# from an environment that has Gurobi (envs/tadmm). Its result lines are tagged
# CENTRAL_GUROBI, so nothing can mistake them for an Ipopt reference, and its
# log (presolve statistics included) goes to the same log path.
#   GUROBI_EXTRA_OPTIONS="Method=2;Crossover=0;Threads=1"
# Gurobi recognises the cone P^2 + Q^2 <= v ell only while v and ell are
# declared nonnegative, so with IPOPT_SCREEN it keeps those: every lower
# voltage limit (the substation's included) and ell >= 0 stay; the screened
# upper voltage limits and P_Subs >= 0 go, as for Ipopt.
central_solver = lowercase(get(ENV, "CENTRAL_SOLVER", "ipopt"))
central_solver in ("ipopt", "gurobi") || error("CENTRAL_SOLVER must be ipopt or gurobi")
result_tag = central_solver == "gurobi" ? "CENTRAL_GUROBI" : "CENTRAL_IPOPT"
if central_solver == "gurobi"
    using Gurobi
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "LogFile", ipopt_log)
    for kv in split(get(ENV, "GUROBI_EXTRA_OPTIONS", ""), ';'; keepempty=false)
        k, v = strip.(split(kv, '='; limit=2))
        val = something(tryparse(Int, v), tryparse(Float64, v), String(v))
        set_optimizer_attribute(model, String(k), val)
        println("GUROBI_OPTION ", k, "=", val)
    end
else
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 5)
    set_optimizer_attribute(model, "max_iter", 5000)
    set_optimizer_attribute(model, "print_timing_statistics", "yes")
    set_optimizer_attribute(model, "output_file", ipopt_log)
    # Extra Ipopt options, e.g. the linear solver:
    #   IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=<dir>/libma57.dll;linear_system_scaling=none"
    # Values that parse as integers or floats are passed as numbers.
    for kv in split(get(ENV, "IPOPT_EXTRA_OPTIONS", ""), ';'; keepempty=false)
        k, v = strip.(split(kv, '='; limit=2))
        val = something(tryparse(Int, v), tryparse(Float64, v), String(v))
        set_optimizer_attribute(model, String(k), val)
        println("IPOPT_OPTION ", k, "=", val)
    end
end

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
# The exact bound reductions of the FilterDDP driver (FILTERDDP_SCREEN there,
# same item names; voltage_screening.jl), so that both solvers can be given
# the identical model. Off by default.
screening_start_ns = time_ns()
ipopt_screen = Set(String.(strip.(split(get(ENV, "IPOPT_SCREEN", ""), ','; keepempty=false))))
"all" in ipopt_screen && (ipopt_screen = Set(["substation", "vupper", "vlower", "ell", "psubs"]))
if !isempty(ipopt_screen)
    isempty(setdiff(ipopt_screen, ["substation", "vupper", "vlower", "ell", "psubs"])) ||
        error("IPOPT_SCREEN: unknown item in $(collect(ipopt_screen))")
    vs = voltage_screen(data, T)
    drop_lo = 0; drop_up = 0
    keep_cone = central_solver == "gurobi"     # v >= 0 and ell >= 0 stay (see above)
    for (n, j) in enumerate(vs.buses)
        fixed = j == root && "substation" in ipopt_screen
        if !keep_cone && (fixed || (j != root && "vlower" in ipopt_screen && !any(@view vs.keep_lo[n, :])))
            for t in Tset; delete_lower_bound(v[j,t]); end
            global drop_lo += 1
        end
        if fixed || (j != root && "vupper" in ipopt_screen && !any(@view vs.keep_up[n, :]))
            for t in Tset; delete_upper_bound(v[j,t]); end
            global drop_up += 1
        end
    end
    drop_ell = "ell" in ipopt_screen && !keep_cone
    drop_ell && for e in Lset, t in Tset; delete_lower_bound(ell[e,t]); end
    drop_ps = "psubs" in ipopt_screen && minimum(vs.psubs_min) > 0
    drop_ps && for t in Tset; delete_lower_bound(P_Subs[t]); end
    @printf("BOUND_SCREEN %s: voltage limits dropped %d lower, %d upper of %d buses; ell >= 0 dropped on %d lines; P_Subs >= 0 dropped: %s\n",
            join(sort!(collect(ipopt_screen)), ","), drop_lo, drop_up, length(Nset), drop_ell ? length(Lset) : 0, drop_ps)
end
# IPOPT_LOADFLOW_START=1: the starting point the FilterDDP driver gets from
screening_s = isempty(ipopt_screen) ? 0.0 : (time_ns() - screening_start_ns) / 1e9
# FILTERDDP_LOADFLOW_START (loadflow_start.jl): batteries idle, DERs at zero
# reactive power, and the power flow of each period. Ipopt's own options are
# left at their defaults (it still pushes the point inside the bounds).
powerflow_init_s = 0.0
if get(ENV, "IPOPT_LOADFLOW_START", "0") != "0"
    powerflow_start_ns = time_ns()
    for t in Tset
        vv, PP, QQ, LL, _ = get(ENV, "IPOPT_MATCH_DDP_START", "0") == "1" ?
            loadflow_initial_values(data, t) : loadflow_sweep(data, t)
        for (n, j) in enumerate(Nset); set_start_value(v[j,t], vv[n]); end
        ps0 = 0.0; qs0 = 0.0
        for (le, e) in enumerate(Lset)
            set_start_value(P[e,t], PP[le]); set_start_value(Q[e,t], QQ[le]); set_start_value(ell[e,t], LL[le])
            if e[1] == root
                ps0 += PP[le]; qs0 += QQ[le]
            end
        end
        set_start_value(P_Subs[t], ps0); set_start_value(Q_Subs[t], qs0)
        for j in Bset; set_start_value(P_B[j,t], 0.0); set_start_value(B[j,t], data[:B0_pu][j]); end
        for j in Dset; set_start_value(q_D[j,t], 0.0); end
    end
    println("IPOPT_LOADFLOW_START applied")
    powerflow_init_s = (time_ns() - powerflow_start_ns) / 1e9
end
build_s = time() - build_start

# IPOPT_TRACE=1: record every iteration through Ipopt's intermediate callback
# (clock, objective, primal and dual infeasibility, barrier parameter) and
# print one IPOPT_ITER line per iteration after the solve. The time at which
# any looser criterion was first met can then be read from one run to the
# solver's own tolerance (ipopt_trace_posthoc.py), as near_opt_posthoc.py does
# for FilterDDP. The callback only stores six numbers; nothing is printed
# during the solve. The objective is the unscaled one of Ipopt's iteration
# table. inf_pr is Ipopt's internal primal infeasibility, which counts an
# inequality through its slack (d(x) - s); the iteration table prints the
# violation of the original constraints instead, which is zero whenever an
# inequality holds. FilterDDP's primal_inf is of the internal kind (the cone
# rows carry an explicit slack), so the two are compared like for like;
# ipopt_trace_posthoc.py reports both.
ipopt_trace = Tuple{Int, UInt64, Float64, Float64, Float64, Float64}[]
if central_solver == "ipopt" && get(ENV, "IPOPT_TRACE", "0") != "0"
    function ipopt_trace_callback(alg_mod::Cint, iter_count::Cint, obj_value::Float64, inf_pr::Float64,
                                  inf_du::Float64, mu::Float64, d_norm::Float64, regularization_size::Float64,
                                  alpha_du::Float64, alpha_pr::Float64, ls_trials::Cint)
        push!(ipopt_trace, (Int(iter_count), time_ns(), obj_value, inf_pr, inf_du, mu))
        return true
    end
    MOI.set(model, Ipopt.CallbackFunction(), ipopt_trace_callback)
    # Compile it now: its first call would otherwise cost about 50 ms inside the solve.
    ipopt_trace_callback(Cint(0), Cint(0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Cint(0)); empty!(ipopt_trace)
    println("IPOPT_TRACE on")
end

ipopt_warm_trace = copy(ipopt_trace)
ipopt_warmup_s = 0.0
if central_solver == "ipopt" && get(ENV, "IPOPT_WARMUP", "0") == "full"
    set_optimizer_attribute(model, "warm_start_init_point", "no")
    warm_begin_ns = time_ns()
    optimize!(model)
    termination_status(model) == MOI.LOCALLY_SOLVED || error("Ipopt full warm-up failed")
    global ipopt_warmup_s = (time_ns() - warm_begin_ns) / 1e9
    global ipopt_warm_trace = copy(ipopt_trace)
    empty!(ipopt_trace)
    GC.gc()
end
optimize!(model)
ipopt_end_ns = time_ns()
if !isempty(ipopt_warm_trace) && !isempty(ipopt_trace)
    delta_start = abs(ipopt_warm_trace[1][3] - ipopt_trace[1][3])
    delta_start <= 1e-10 * max(1.0, abs(ipopt_warm_trace[1][3])) || error("Ipopt starting objective mismatch")
    @printf("IPOPT_WARMUP_CHECK warmup_s=%.9f initial_objective_difference=%.9e warmup_iterations=%d timed_iterations=%d\n",
            ipopt_warmup_s, delta_start, ipopt_warm_trace[end][1], ipopt_trace[end][1])
end
wall_s = time() - wall_start

status = termination_status(model)
obj = has_values(model) ? objective_value(model) : NaN
iters = try MOI.get(model, MOI.BarrierIterations()) catch; -1 end
solve_s = try solve_time(model) catch; NaN end
# elapsed_s counts from the start of Ipopt's own solve time, so that the last
# line's elapsed_s is the reported solve time (to the return of optimize!).
for (k, tns, o, ip, id, m) in ipopt_trace
    @printf("IPOPT_ITER iteration=%d elapsed_s=%.6f objective=%.12e inf_pr=%.6e inf_du=%.6e mu=%.6e\n",
            k, solve_s - (ipopt_end_ns - tns) / 1e9, o, ip, id, m)
end
@printf("%s system=%s T=%d profile=%s C_B=%.6e gamma=%.6e status=%s iterations=%d objective=%.12f solve_time_s=%.3f build_s=%.3f wall_s=%.3f variables=%d solver=%s\n",
        result_tag, system, T, isempty(ptag) ? "default" : ptag[2:end], data[:C_B], gammaT, string(status),
        iters, obj, solve_s, build_s, wall_s, num_variables(model),
        replace(string(MOI.get(model, MOI.SolverVersion())), ' ' => '_'))
# Peak resident memory of the whole process (model build + solve), for the
# memory-knee comparison between linear solvers.
@printf("%s_MEMORY maxrss_mib=%.1f\n", result_tag, Sys.maxrss() / 2^20)
@printf("%s_PREPROCESS_TIMING screening_s=%.9f powerflow_init_s=%.9f build_excluding_preprocessing_s=%.9f matched_physical_start=%d\n",
        result_tag, screening_s, powerflow_init_s, build_s-screening_s-powerflow_init_s,
        get(ENV, "IPOPT_MATCH_DDP_START", "0") == "1")

# IPOPT_SAVE_DISPATCH=<file>: the solution in the FilterDDP control layout
# (control_layout.jl; slacks left at zero), so that opendss_check.jl checks a
# centralized dispatch exactly as it checks a FilterDDP one. Written after the
# solve, outside every timer. Any CENTRAL_SOLVER.
if haskey(ENV, "IPOPT_SAVE_DISPATCH") && has_values(model)
    lay, nu_lay = control_layout(data)
    u_save = [zeros(nu_lay) for _ in Tset]
    for t in Tset
        ut = u_save[t]
        ut[lay.ps] = value(P_Subs[t]); ut[lay.qs] = value(Q_Subs[t])
        for (k, e) in enumerate(Lset)
            ut[lay.P[k]] = value(P[e,t]); ut[lay.Q[k]] = value(Q[e,t]); ut[lay.ell[k]] = value(ell[e,t])
        end
        for (k, j) in enumerate(Nset); ut[lay.v[k]] = value(v[j,t]); end
        for (k, j) in enumerate(Bset); ut[lay.pb[k]] = value(P_B[j,t]); end
        for (k, j) in enumerate(Dset)
            qmax = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
            ut[lay.qnorm[k]] = qmax > 0 ? value(q_D[j,t]) / qmax : 0.0
        end
    end
    serialize(ENV["IPOPT_SAVE_DISPATCH"], Dict(:system => system, :T => T, :status => string(status),
        :iterations => iters, :objective => obj, :u => u_save))
    println("saved_dispatch=", ENV["IPOPT_SAVE_DISPATCH"])
end

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
    @printf("%s_BATTERY at_power_limit=%.3f mean_utilization=%.3f energy_window_used=%.3f throughput_pu_h=%.6e energy_term=%.6e cb_term=%.6e terminal_term=%.6e\n",
            result_tag, count(>=(0.99), util) / length(util), sum(util) / length(util),
            sum(window) / length(window), sum(abs.(pb)) * dt, energy_term, cb_term, term_term)
end
