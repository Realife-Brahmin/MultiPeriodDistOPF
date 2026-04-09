# ============================================================================
# run_bf.jl — Brute-Force (monolithic) SOCP solve for MPOPF
# ============================================================================
# Usage: edit config.jl, then run this script
#   julia run_bf.jl
#   OR in VS Code: Ctrl+Shift+Enter
# ============================================================================

# Fix GR backend display issues on headless/Windows environments
ENV["GKSwstype"] = "png"

# Load shared configuration
include("config.jl")

# Activate the tadmm environment
import Pkg
Pkg.activate(ENV_PATH)

using Revise
using LinearAlgebra
using JuMP
using Ipopt
if USE_GUROBI
    using Gurobi
end
using OpenDSSDirect
using Printf
using Statistics
using Parameters: @unpack
using Crayons
using Serialization

const MOI = JuMP.MOI

# Create shared Gurobi environment (suppresses repeated license messages)
const GUROBI_ENV = USE_GUROBI ? Gurobi.Env() : nothing

# Color scheme
const COLOR_SUCCESS = Crayon(foreground = :green, bold = true)
const COLOR_WARNING = Crayon(foreground = :yellow, bold = true)
const COLOR_ERROR = Crayon(foreground = :red, bold = true)
const COLOR_INFO = Crayon(foreground = :cyan, bold = true)
const COLOR_HIGHLIGHT = Crayon(foreground = :magenta, bold = true)
const COLOR_RESET = Crayon(reset = true)

# Include utilities
includet(joinpath(ENV_PATH, "parse_opendss.jl"))
includet(joinpath(ENV_PATH, "opendss_validator.jl"))
includet(joinpath(ENV_PATH, "solution_validator.jl"))
includet(joinpath(ENV_PATH, "logger.jl"))

# ============================================================================
# PARSE SYSTEM DATA
# ============================================================================

println("\n" * "="^80)
println("PARSING SYSTEM DATA")
println("="^80)

data = parse_system_from_dss(
    SYSTEM_NAME, T;
    LoadShapeLoad=LoadShapeLoad,
    LoadShapeCost=LoadShapeCost,
    LoadShapePV=LoadShapePV,
    C_B=C_B,
    delta_t_h=DELTA_T_H
)

print_powerflow_summary(data)

# ============================================================================
# BF SOLVE FUNCTION
# ============================================================================

function solve_MPOPF_with_SOCP_BruteForced(data; solver=:ipopt)

    @unpack Nset, Lset, Dset, Bset, Tset, NLset, Nm1set = data
    @unpack substationBus, L1set, Lm1set, parent, children = data
    @unpack p_L_pu, p_D_pu, q_L_pu = data
    @unpack P_B_R_pu, B_R_pu, B0_pu, S_D_R = data
    @unpack rdict_pu, xdict_pu = data
    @unpack Vminpu, Vmaxpu = data
    @unpack LoadShapeCost, C_B, delta_t_h, soc_min, soc_max = data
    @unpack kVA_B = data

    dt = delta_t_h
    P_BASE = kVA_B
    j1 = substationBus

    # ========== CREATE MODEL ==========
    model = Model()
    if solver == :gurobi
        set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
        set_optimizer_attribute(model, "NonConvex", 2)
        set_optimizer_attribute(model, "OutputFlag", 0)
        set_optimizer_attribute(model, "DualReductions", 0)
    else
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 5)
        set_optimizer_attribute(model, "max_iter", 5000)
        set_optimizer_attribute(model, "output_file", "ipopt_bf.log")
    end

    # ========== VARIABLES ==========
    @variable(model, P_Subs[t in Tset] >= 0)
    @variable(model, Q_Subs[t in Tset])
    @variable(model, P[(i, j) in Lset, t in Tset])
    @variable(model, Q[(i, j) in Lset, t in Tset])
    @variable(model, v[j in Nset, t in Tset])
    @variable(model, l[(i, j) in Lset, t in Tset] >= 0)
    @variable(model, P_B[j in Bset, t in Tset])
    @variable(model, B[j in Bset, t in Tset])
    @variable(model, q_D[j in Dset, t in Tset])

    # ========== OBJECTIVE ==========
    @expression(model, energy_cost,
        sum(LoadShapeCost[t] * P_Subs[t] * P_BASE * dt for t in Tset))
    @expression(model, battery_cost,
        sum(C_B * (P_B[j, t] * P_BASE)^2 * dt for j in Bset, t in Tset))
    @objective(model, Min, energy_cost + battery_cost)

    # ========== CONSTRAINTS ==========
    for t in Tset
        # Real power balance - substation
        @constraint(model, P_Subs[t] - sum(P[(j1, j), t] for (j1, j) in L1set) == 0,
            base_name = "RealPowerBalance_Substation_t$(t)")

        # Real power balance - non-substation
        for j in Nm1set
            i = parent[j]
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P[(j, k), t] for k in children[j])
            p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0
            p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0
            P_B_j_t = (j in Bset) ? P_B[j, t] : 0.0
            @constraint(model, sum_Pjk - P[(i, j), t] == P_B_j_t + p_D_j_t - p_L_j_t,
                base_name = "RealPowerBalance_Node$(j)_t$(t)")
        end

        # Reactive power balance - substation
        @constraint(model, Q_Subs[t] - sum(Q[(j1, j), t] for (j1, j) in L1set) == 0,
            base_name = "ReactivePowerBalance_Substation_t$(t)")

        # Reactive power balance - non-substation
        for j in Nm1set
            i = parent[j]
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
            q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
            @constraint(model, sum_Qjk - Q[(i, j), t] == q_D_j_t - q_L_j_t,
                base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
        end

        # Voltage drop (BFM-NL)
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model,
                v[j, t] == v[i, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) + (r_ij^2 + x_ij^2) * l[(i, j), t],
                base_name = "VoltageDrop_Branch_$(i)_$(j)_t$(t)")
        end

        # SOC relaxation
        for (i, j) in Lset
            @constraint(model,
                P[(i, j), t]^2 + Q[(i, j), t]^2 <= v[i, t] * l[(i, j), t],
                base_name = "SOC_Branch_$(i)_$(j)_t$(t)")
        end

        # Fixed substation voltage
        @constraint(model, v[j1, t] == 1.05^2,
            base_name = "FixedSubstationVoltage_t$(t)")

        # Voltage limits
        for j in Nset
            @constraint(model, Vminpu[j]^2 <= v[j, t] <= Vmaxpu[j]^2,
                base_name = "VoltageLimits_Node$(j)_t$(t)")
        end

        # PV reactive power limits
        for j in Dset
            p_D_val = p_D_pu[j, t]
            S_D_R_val = S_D_R[j]
            q_max_t = sqrt(S_D_R_val^2 - p_D_val^2)
            @constraint(model, -q_max_t <= q_D[j, t] <= q_max_t,
                base_name = "PVReactiveLimits_DER$(j)_t$(t)")
        end

        # Battery constraints
        for j in Bset
            if t == 1
                @constraint(model, B[j, t] == B0_pu[j] - P_B[j, t] * dt,
                    base_name = "BatterySOC_Init_$(j)_t$(t)")
            else
                @constraint(model, B[j, t] == B[j, t-1] - P_B[j, t] * dt,
                    base_name = "BatterySOC_$(j)_t$(t)")
            end
            @constraint(model, soc_min[j] * B_R_pu[j] <= B[j, t] <= soc_max[j] * B_R_pu[j],
                base_name = "BatterySOCLimits_$(j)_t$(t)")
            @constraint(model, -P_B_R_pu[j] <= P_B[j, t] <= P_B_R_pu[j],
                base_name = "BatteryPowerLimits_$(j)_t$(t)")
        end
    end

    # ========== SAVE MODEL SUMMARY & SOLVE ==========
    mkpath(SYSTEM_DIR)
    model_file = joinpath(SYSTEM_DIR, "model_summary_socp_bf.txt")
    write_model_summary(model, data, model_file)

    bf_model_stats = get_model_size_statistics(model)

    println("\n" * "="^80)
    println("STARTING BF OPTIMIZATION")
    println("="^80)
    println("[BF] Launching Ipopt solver...")

    bf_wallclock_start = time()

    if BF_TIMEOUT_SEC > 0
        set_time_limit_sec(model, BF_TIMEOUT_SEC)
        println("[BF] Timeout set to $(BF_TIMEOUT_SEC)s")
    end

    optimize!(model)
    bf_wallclock_time = time() - bf_wallclock_start

    if BF_TIMEOUT_SEC > 0 && bf_wallclock_time >= BF_TIMEOUT_SEC
        println("[BF] TIMEOUT: BF killed after $(round(bf_wallclock_time, digits=1))s")
    else
        println("[BF] Optimization completed in $(round(bf_wallclock_time, digits=1))s")
    end

    # ========== EXTRACT RESULTS ==========
    status = termination_status(model)
    obj_val = has_values(model) ? objective_value(model) : NaN
    solver_time = solve_time(model)

    bf_ipopt_iters = 0
    try
        bf_ipopt_iters = MOI.get(model, MOI.BarrierIterations())
    catch
        bf_ipopt_iters = 0
    end
    bf_time_per_iter = bf_ipopt_iters > 0 ? solver_time / bf_ipopt_iters : NaN

    result = Dict(
        :model => model,
        :status => status,
        :objective => obj_val,
        :solve_time => solver_time,
        :wallclock_time => bf_wallclock_time,
        :ipopt_iters => bf_ipopt_iters,
        :time_per_iter => bf_time_per_iter,
        :P_Subs => has_values(model) ? value.(P_Subs) : P_Subs,
        :Q_Subs => has_values(model) ? value.(Q_Subs) : Q_Subs,
        :P => has_values(model) ? value.(P) : P,
        :Q => has_values(model) ? value.(Q) : Q,
        :v => has_values(model) ? value.(v) : v,
        :l => has_values(model) ? value.(l) : l,
        :P_B => has_values(model) ? value.(P_B) : P_B,
        :B => has_values(model) ? value.(B) : B,
        :q_D => has_values(model) ? value.(q_D) : q_D,
        :model_stats => bf_model_stats,
    )

    # Validate solution
    if has_values(model)
        println("\n" * "="^80)
        println(COLOR_INFO, "VALIDATING BF SOLUTION", COLOR_RESET)
        println("="^80)
        bf_validation = validate_branch_flow_equations(result, data; tol=1e-4, verbose=false)
        print_validation_summary(bf_validation; solution_name="Brute Force SOCP")
    end

    return result
end

# ============================================================================
# SAVE RESULTS
# ============================================================================

function save_bf_results(sol, data)
    mkpath(SYSTEM_DIR)

    status = sol[:status]
    if !(status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED)
        println(COLOR_ERROR, "BF did not converge (status: $status), skipping result save", COLOR_RESET)
        return
    end

    # Validate solution
    bf_validation = validate_branch_flow_equations(sol, data; tol=1e-4, verbose=false)

    # Write results text file
    results_file = joinpath(SYSTEM_DIR, "results_socp_bf.txt")
    open(results_file, "w") do io
        println(io, "="^80)
        println(io, "SOCP BRUTE-FORCED OPTIMIZATION RESULTS")
        println(io, "="^80)
        println(io, "System: $(SYSTEM_NAME)")
        println(io, "Time horizon: T=$(T) periods")
        println(io, "Time step: $(DELTA_T_H) hours")
        println(io, "Number of buses: $(length(data[:Nset]))")
        println(io, "Number of branches: $(length(data[:Lset]))")
        println(io, "Number of batteries: $(length(data[:Bset]))")
        println(io, "Number of PV units: $(length(data[:Dset]))")
        println(io, "\n--- PROBLEM SIZE ---")
        @printf(io, "Number of variables: %d\n", sol[:model_stats][:n_variables])
        @printf(io, "Linear constraints: %d\n", sol[:model_stats][:n_linear_constraints])
        @printf(io, "Quadratic constraints (SOCP): %d\n", sol[:model_stats][:n_quadratic_constraints])
        @printf(io, "Nonlinear constraints: %d\n", sol[:model_stats][:n_nonlinear_constraints])
        @printf(io, "Total constraints: %d\n", sol[:model_stats][:total_constraints])
        println(io, "\n--- OPTIMIZATION STATUS ---")
        println(io, "Status: $(sol[:status])")
        println(io, "\n--- OBJECTIVE VALUE ---")
        @printf(io, "Total Cost: \$%.4f\n", sol[:objective])
        println(io, "\n--- COMPUTATION TIME ---")
        @printf(io, "Wall-clock time: %.4f seconds\n", sol[:wallclock_time])
        @printf(io, "Solver time: %.4f seconds\n", sol[:solve_time])
        @printf(io, "Ipopt iterations: %d\n", sol[:ipopt_iters])
        if sol[:ipopt_iters] > 0
            @printf(io, "Avg time per iteration: %.4f seconds\n", sol[:time_per_iter])
        end
        println(io, "\n--- SOLUTION FEASIBILITY ---")
        if bf_validation[:feasible]
            println(io, "Status: FEASIBLE - All constraints satisfied")
        else
            println(io, "Status: INFEASIBLE - $(bf_validation[:total_violations]) constraint violations")
            println(io, "\nViolation summary:")
            for (key, count) in bf_validation[:violations]
                if count > 0
                    max_viol = bf_validation[:max_violations][key]
                    @printf(io, "  %-30s: %6d violations (max: %.2e)\n", string(key), count, max_viol)
                end
            end
        end
        println(io, "\n--- SOLVER ---")
        solver_name = USE_GUROBI_FOR_BF ? "Gurobi" : "Ipopt"
        println(io, "Solver: $(solver_name)")
        println(io, "Formulation: SOCP (BFM-NL)")
        println(io, "="^80)
    end
    println(COLOR_SUCCESS, "Results written to $(results_file)", COLOR_RESET)

    # Save serialized solution (.jls, gitignored)
    sol_save = Dict(k => v for (k, v) in sol if k != :model)
    sol_file = joinpath(SYSTEM_DIR, "sol_socp_bf.jls")
    open(sol_file, "w") do io
        serialize(io, sol_save)
    end
    println(COLOR_SUCCESS, "Solution saved to $(sol_file)", COLOR_RESET)
end

# ============================================================================
# RUN
# ============================================================================

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH SOCP (BRUTE-FORCED)", COLOR_RESET)
println("="^80)

bf_start = time()
solver_choice = USE_GUROBI_FOR_BF ? :gurobi : :ipopt
sol_socp_bf = solve_MPOPF_with_SOCP_BruteForced(data; solver=solver_choice)

# Report
if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
    print(COLOR_SUCCESS)
    println("\nStatus: $(sol_socp_bf[:status])")
    @printf "Total Cost: \$%.2f\n" sol_socp_bf[:objective]
    @printf "Wall-clock time: %.2f seconds\n" sol_socp_bf[:wallclock_time]
    @printf "Solver time: %.2f seconds\n" sol_socp_bf[:solve_time]
    @printf "Ipopt iterations: %d\n" sol_socp_bf[:ipopt_iters]
    print(COLOR_RESET)
else
    print(COLOR_ERROR)
    println("Optimization failed: $(sol_socp_bf[:status])")
    print(COLOR_RESET)
end

save_bf_results(sol_socp_bf, data)

total_time = time() - bf_start
println("\n" * "="^80)
@printf "Total BF pipeline time: %.2f seconds\n" total_time
println("="^80)
