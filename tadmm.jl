#! Activate project environment

import Pkg; Pkg.activate(".")
using Revise
using JuMP
import OpenDSSDirect as dss
using OpenDSSDirect
using Dates
using Printf
using Random
using Statistics  # For mean, std, etc.
using Parameters: @unpack
using Gurobi
using Ipopt

# Import the ODD parser
Revise.includet("src/Parser/parseFromDSS.jl")
import .parseFromDSS as Parser


# --- USER DEFINED META VARIABLES ---
# System and simulation parameters (match copper plate example)
systemName = "ads10A_1ph"   # Change as needed
T = 24                      # Number of time steps
delta_t_h = 1.0               # Time step duration in hours

# --- COPPER PLATE VALUES (from admm_temporal_copper_plate.jl) ---
LoadShape = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2  # Normalized load shape [0, 1]
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2   # $/kWh
C_B = 1e-6 * minimum(LoadShapeCost)  # Quadratic cost coefficient for battery power: C_B * P_B^2

# If you want to keep PVShape, you can set it to zeros or as needed
# PVShape = zeros(T)

# --- Helper: update data dict with user-supplied kwargs ---
function update_data_with_kwargs!(data::Dict; kwargs...)
    for (k, v) in kwargs
        data[k] = v
    end
    return data
end

# --- PARSE SYSTEM AND APPLY USER OVERRIDES ---
data = Parser.get_system_config_from_dss(systemName, T)

user_overrides = Dict(
    :LoadShapeLoad => LoadShape,
    :LoadShapeCost => LoadShapeCost,
    :C_B => C_B,
    :delta_t_h => delta_t_h,
    # :LoadShapePV => PVShape,
)
update_data_with_kwargs!(data; user_overrides...)



# --- Print basic OpenDSS powerflow stats for quick validation ---
include("src/openDSSValidator.jl")
using .openDSSValidator
print_basic_powerflow_stats(data)


# --- SAVE DATA DICT KEYS (GROUPED BY CONTEXT WITH HEADERS) ---
const context_map = Dict(
    "Substation" => ["PSubsMax_kW", "V_Subs_pu", "substationBus"],
    "Simulation_Problem" => ["T", "Tset", "gedAppendix", "gedString", "inputForecastDescription", "objfun0", "objfun2", "objfunSense", "objfunString", "objfunPrefix", "objfunUnit", "objfunAppendix", "objfunConciseDescription", "tSOC_hard", "relax_terminal_soc_constraint", "LoadShapeLoad", "LoadShapePV", "LoadShapeCost", "C_B", "offPeakCost", "peakCost", "peakHoursFraction", "delta_t", "gedDict_ud", "systemName"],
    "Simulation_Algorithm_Params" => ["solver", "linearizedModel", "linearizedModelAppendix", "linearizedModelString", "algo_temporal_decmp", "alpha_fpi", "gamma_fpi", "numAreas", "threshold_conv_iters", "warmStart_mu", "temporal_decmp", "temporalDecmpAppendix", "temporalDecmpString", "simNatureAppendix", "simNatureString", "spatialDecAppendix", "spatialDecString"],
    "Simulation_Machine" => ["machine_ID"],
    "Simulation_Run" => ["macroItrsCompleted", "solution_time"],
    "Directories" => ["processedDataFolderPath", "rawDataFolderPath", "rootFolderPath", "srcFolderPath"],
    "System_Topology" => ["kVA_B", "kVA_B_dict", "kV_B", "kV_B_dict", "MVA_B", "MVA_B_dict", "Z_B", "Z_B_dict", "rdict", "rdict_pu", "xdict", "xdict_pu", "L1set", "LTset", "Lm1set", "LnotTset", "Lset", "N", "N1", "N1set", "NLset", "N_L", "Nc1", "Nc1set", "Nm1", "Nm1set", "Nnc1", "Nnc1set", "Nset", "children", "m", "m1", "mm1", "parent"],
    "System_Storage_Data" => ["B0", "B0_pu", "B_R", "B_R_pu", "Batt_percent", "Bref", "Bref_percent", "Bref_pu", "Bset", "P_B_R", "P_B_R_pu", "S_B_R", "S_B_R_pu", "Vmaxpu_B", "Vminpu_B", "eta_C", "eta_D", "n_B", "soc_0", "soc_max", "soc_min"],
    "System_DER_Data" => ["Dset", "n_D", "DER_percent", "p_D", "p_D_R", "p_D_R_pu", "p_D_pu", "S_D_R", "S_D_R_pu", "Vmaxpu_D", "Vminpu_D", "irrad"],
    "System_Load_Data" => ["NLset", "N_L", "p_L", "p_L_R", "p_L_R_pu", "p_L_pu", "q_L", "q_L_R", "q_L_R_pu", "q_L_pu", "Vmaxpu_L", "Vminpu_L"],
    "Voltage_Limits" => ["Vmaxpu", "Vminpu"],
    "Other_Misc" => []
)


all_keys_str = Set(string(k) for k in keys(data))
open("tadmm_data_keys.txt", "w") do io
    for (section, keyslist) in context_map
        println(io, "==== $section ====")
        for k in sort(keyslist)
            if k in all_keys_str
                println(io, k)
            end
        end
        println(io)
    end
    # Print any keys not matched
    matched_keys = reduce(vcat, values(context_map))
    unmatched = setdiff(all_keys_str, Set(matched_keys))
    if !isempty(unmatched)
        println(io, "==== Unmatched ====")
        for k in sort(collect(unmatched))
            println(io, k)
        end
    end
end

# --- SAVE ALL DATA VALUES (GROUPED BY CONTEXT WITH HEADERS) ---
open("tadmm_data_values.txt", "w") do io
    for (section, keyslist) in context_map
        println(io, "==== $section ====")
        for k in sort(keyslist)
            if k in all_keys_str
                val = get(data, Symbol(k), "<missing>")
                println(io, "[", k, "] => ", val)
            end
        end
        println(io)
    end
    # Print any keys not matched
    matched_keys = reduce(vcat, values(context_map))
    unmatched = setdiff(all_keys_str, Set(matched_keys))
    if !isempty(unmatched)
        println(io, "==== Unmatched ====")
        for k in sort(collect(unmatched))
            val = get(data, Symbol(k), "<missing>")
            println(io, "[", k, "] => ", val)
        end
    end
end

# =============================================================================
# MPOPF LinDistFlow Brute-Force Solver
# =============================================================================
# Comprehensive MPOPF solver using LinDistFlow approximation with network constraints
# Formulation matches ModelBuilder conventions but uses P_B (signed) for battery power
# and includes battery quadratic cost similar to copper plate example

"""
    solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:ipopt)

Solve Multi-Period Optimal Power Flow using LinDistFlow approximation.

# Objective
    min ∑ₜ [LoadShapeCost[t] * P_Subs[t] * Δt + C_B * P_B[b,t]² * Δt]

# Constraints
    - Nodal real power balance (substation & all nodes)
    - Nodal reactive power balance (substation & all nodes)  
    - KVL for all branches (LinDistFlow)
    - Fixed substation voltage (1.07 pu)
    - Voltage box limits for all nodes
    - PV reactive power limits based on S_D_R and p_D[t]
    - Battery SOC trajectory: B[t] = B[t-1] - P_B[t-1] * Δt
    - Battery SOC and power limits
    - No battery reactive power (for now)

# Arguments
    - data: Dictionary with all system data from parser
    - solver: :ipopt (default) or :gurobi

# Returns
    Dictionary with:
    - :model => JuMP model
    - :status => termination status
    - :objective => objective value
    - :P_Subs => substation power [t]
    - :P => branch real power [(i,j), t]
    - :Q => branch reactive power [(i,j), t]
    - :v => squared voltage [n, t]
    - :P_B => battery power [b, t] (signed: + discharging, - charging)
    - :B => battery SOC [b, t]
    - :q_D => PV reactive power [d, t]
"""
function solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:ipopt)
    
    # ========== 1. UNPACK DATA ==========
    @unpack Nset, Lset, Dset, Bset, Tset, NLset, Nm1set = data
    @unpack substationBus, L1set, Lm1set, parent, children = data
    @unpack p_L_pu, p_D_pu, q_L_pu = data
    @unpack P_B_R_pu, B_R_pu, B0_pu, S_D_R = data
    @unpack rdict_pu, xdict_pu = data
    @unpack Vminpu, Vmaxpu = data
    @unpack LoadShapeCost, C_B, delta_t_h, soc_min, soc_max = data
    @unpack kVA_B = data
    
    Δt = delta_t_h  # Time step duration in hours
    P_BASE = kVA_B  # kW base power
    j1 = substationBus  # Substation bus number
    
    # Debug: Print some basic network info
    println("\n[DEBUG] Network Info:")
    println("  Substation bus: ", j1)
    println("  Number of nodes: ", length(Nset))
    println("  Number of branches: ", length(Lset))
    println("  Number of load nodes: ", length(NLset))
    println("  Number of PV nodes: ", length(Dset))
    println("  Number of battery nodes: ", length(Bset))
    println("  Base power: ", P_BASE, " kVA")
    
    # Check total load at t=1
    total_load_t1 = sum(p_L_pu[j, 1] for j in NLset if haskey(p_L_pu, (j, 1)))
    total_pv_t1 = sum(p_D_pu[d, 1] for d in Dset if haskey(p_D_pu, (d, 1)))
    println("  Total load at t=1: ", total_load_t1 * P_BASE, " kW")
    println("  Total PV at t=1: ", total_pv_t1 * P_BASE, " kW")
    println("  Expected net substation power: ", (total_load_t1 - total_pv_t1) * P_BASE, " kW")
    
    # ========== 2. CREATE MODEL ==========
    model = Model()
    if solver == :gurobi
        set_optimizer(model, Gurobi.Optimizer)
        set_optimizer_attribute(model, "OutputFlag", 0)
    elseif solver == :ipopt
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0)
    else
        error("Unknown solver: $solver. Use :ipopt or :gurobi")
    end
    
    # ========== 3. DEFINE VARIABLES ==========
    @variable(model, P_Subs[t in Tset] >= 0)  # Substation real power (non-negative)
    @variable(model, P[(i, j) in Lset, t in Tset])  # Branch real power
    @variable(model, Q[(i, j) in Lset, t in Tset])  # Branch reactive power
    @variable(model, v[n in Nset, t in Tset] >= 0)  # Squared voltage magnitude (non-negative)
    @variable(model, P_B[b in Bset, t in Tset])  # Battery power (signed)
    @variable(model, B[b in Bset, t in Tset] >= 0)  # Battery SOC (non-negative)
    @variable(model, q_D[d in Dset, t in Tset])  # PV reactive power
    
    # ========== 4. OBJECTIVE FUNCTION ==========
    # Energy cost + Battery quadratic cost
    @expression(model, energy_cost, 
        sum(LoadShapeCost[t] * P_Subs[t] * P_BASE * Δt for t in Tset))
    @expression(model, battery_cost, 
        sum(C_B * (P_B[b, t] * P_BASE)^2 * Δt for b in Bset, t in Tset))
    @objective(model, Min, energy_cost + battery_cost)
    
    # ========== 5. CONSTRAINTS ==========
    
    for t in Tset
        # ----- 5.1 NODAL REAL POWER BALANCE -----
        # Substation node (j1)
        @constraint(model, 
            P_Subs[t] - sum(P[(j1, j), t] for (i, j) in L1set) == 0,
            base_name = "RealPowerBalance_Substation_t$(t)")
        
        # Non-substation nodes (using parent-child topology)
        for j in Nm1set
            i = parent[j]  # Parent node
            P_ij_t = P[(i, j), t]  # Power from parent
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P[(j, k), t] for k in children[j])
            
            # Net injection at node j
            p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0
            p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0
            P_B_j_t = (j in Bset) ? P_B[j, t] : 0.0
            
            # Power balance: Incoming = Outgoing + Load - PV - Battery
            # P_ij = sum_Pjk + p_L - p_D - P_B
            @constraint(model,
                P_ij_t - sum_Pjk - p_L_j_t + p_D_j_t + P_B_j_t == 0,
                base_name = "RealPowerBalance_Node$(j)_t$(t)")
        end
        
        # ----- 5.2 NODAL REACTIVE POWER BALANCE -----
        # Substation node (j1) - assume no reactive injection at substation
        @constraint(model,
            sum(Q[(j1, j), t] for (i, j) in L1set) == 0,
            base_name = "ReactivePowerBalance_Substation_t$(t)")
        
        # Non-substation nodes
        for j in Nm1set
            i = parent[j]
            Q_ij_t = Q[(i, j), t]
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            
            q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
            q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
            # No battery reactive power for now (q_B = 0)
            
            # Power balance: Incoming = Outgoing + Load - PV
            # Q_ij = sum_Qjk + q_L - q_D
            @constraint(model,
                Q_ij_t - sum_Qjk - q_L_j_t + q_D_j_t == 0,
                base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
        end
        
        # ----- 5.3 KVL CONSTRAINTS (LinDistFlow) -----
        # Substation branches
        for (i, j) in L1set
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model,
                v[i, t] - v[j, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) == 0,
                base_name = "KVL_SubstationBranch_$(i)_$(j)_t$(t)")
        end
        
        # Non-substation branches
        for (i, j) in Lm1set
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model,
                v[i, t] - v[j, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) == 0,
                base_name = "KVL_NonSubstationBranch_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.4 VOLTAGE CONSTRAINTS -----
        # Fixed substation voltage (1.07 pu, squared)
        @constraint(model, v[j1, t] == 1.07^2,
            base_name = "FixedSubstationVoltage_t$(t)")
        
        # Voltage limits (all nodes)
        for n in Nset
            @constraint(model, Vminpu[n]^2 <= v[n, t] <= Vmaxpu[n]^2,
                base_name = "VoltageLimits_Node$(n)_t$(t)")
        end
        
        # ----- 5.5 PV REACTIVE POWER LIMITS -----
        for d in Dset
            # q_D must satisfy: q_D² + p_D² ≤ S_D_R²
            # For now, use simple box: -√(S_D_R² - p_D²) ≤ q_D ≤ √(S_D_R² - p_D²)
            p_D_val = p_D_pu[d, t]
            S_D_R_val = S_D_R[d]
            q_max = sqrt(max(0, S_D_R_val^2 - p_D_val^2))
            @constraint(model, -q_max <= q_D[d, t] <= q_max,
                base_name = "PVReactiveLimits_DER$(d)_t$(t)")
        end
        
        # ----- 5.6 BATTERY CONSTRAINTS -----
        for b in Bset
            # SOC trajectory
            if t == 1
                @constraint(model, B[b, t] == B0_pu[b] - P_B[b, t] * Δt,
                    base_name = "BatterySOC_Init_$(b)_t$(t)")
            else
                @constraint(model, B[b, t] == B[b, t-1] - P_B[b, t-1] * Δt,
                    base_name = "BatterySOC_$(b)_t$(t)")
            end
            
            # SOC limits (use per-unit values)
            @constraint(model, soc_min[b] * B_R_pu[b] <= B[b, t] <= soc_max[b] * B_R_pu[b],
                base_name = "BatterySOCLimits_$(b)_t$(t)")
            
            # Power limits (use per-unit values)
            @constraint(model, -P_B_R_pu[b] <= P_B[b, t] <= P_B_R_pu[b],
                base_name = "BatteryPowerLimits_$(b)_t$(t)")
        end
    end
    
    # ========== 6. SOLVE ==========
    optimize!(model)
    
    # ========== 7. EXTRACT RESULTS ==========
    status = termination_status(model)
    obj_val = has_values(model) ? objective_value(model) : NaN
    
    # Extract solution values (return as JuMP variable refs if solved)
    result = Dict(
        :model => model,
        :status => status,
        :objective => obj_val,
        :P_Subs => has_values(model) ? value.(P_Subs) : P_Subs,
        :P => has_values(model) ? value.(P) : P,
        :Q => has_values(model) ? value.(Q) : Q,
        :v => has_values(model) ? value.(v) : v,
        :P_B => has_values(model) ? value.(P_B) : P_B,
        :B => has_values(model) ? value.(B) : B,
        :q_D => has_values(model) ? value.(q_D) : q_D,
    )
    
    return result
end

# =============================================================================
# SOLVE AND REPORT
# =============================================================================

println("\n" * "="^80)
println("SOLVING MPOPF WITH LINDISTFLOW (BRUTE-FORCED)")
println("="^80)

# Solve the problem
sol_ldf = solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:ipopt)

# Report results
println("\n--- SOLUTION STATUS ---")
println("Status: ", sol_ldf[:status])
if sol_ldf[:status] == MOI.OPTIMAL || sol_ldf[:status] == MOI.LOCALLY_SOLVED
    println("✓ Optimization successful!")
    println("\n--- OBJECTIVE VALUE ---")
    @printf "Total Cost: \$%.2f\n" sol_ldf[:objective]
    
    # Extract solution arrays
    P_Subs_vals = sol_ldf[:P_Subs]
    P_B_vals = sol_ldf[:P_B]
    B_vals = sol_ldf[:B]
    v_vals = sol_ldf[:v]
    
    # Convert to physical units for reporting
    P_BASE = data[:kVA_B]
    E_BASE = P_BASE * 1.0  # kWh per 1 hour
    
    println("\n--- POWER SUMMARY ---")
    if isa(P_Subs_vals, AbstractArray)
        P_Subs_kW = P_Subs_vals .* P_BASE
        @printf "Substation Power (kW): min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_Subs_kW) maximum(P_Subs_kW) mean(P_Subs_kW)
        # Show first few timesteps for debugging
        @printf "First 5 timesteps (kW): "
        for t in 1:min(5, length(P_Subs_kW))
            @printf "%.1f " P_Subs_kW[t]
        end
        println()
        
        # Verify power balance for t=1
        t = 1
        total_load_check = sum(data[:p_L_pu][j, t] for j in data[:NLset] if haskey(data[:p_L_pu], (j, t)))
        total_pv_check = sum(data[:p_D_pu][d, t] for d in data[:Dset] if haskey(data[:p_D_pu], (d, t)))
        battery_power_check = !isempty(data[:Bset]) ? sum(P_B_vals[b, t] for b in data[:Bset]) : 0.0
        expected_psubs = (total_load_check - total_pv_check - battery_power_check) * P_BASE
        @printf "At t=1: Expected P_Subs=%.1f kW, Got P_Subs=%.1f kW\n" expected_psubs P_Subs_kW[1]
        
        # Check branch flows from substation at t=1
        P_vals = sol_ldf[:P]
        println("[DEBUG] Type of P_vals: ", typeof(P_vals))
        if isa(P_vals, AbstractDict)
            @printf "\nBranch flows from substation at t=1 (kW):\n"
            for (i, j) in data[:L1set]
                @printf "  P[(%d,%d)] = %.2f kW\n" i j (P_vals[(i,j), 1] * P_BASE)
            end
            total_P_from_subs = sum(P_vals[(i,j), 1] for (i,j) in data[:L1set]) * P_BASE
            @printf "  Total = %.2f kW (should equal P_Subs = %.2f kW)\n" total_P_from_subs P_Subs_kW[1]
        else
            println("[DEBUG] P_vals is not a Dict, checking if it's a JuMP container")
            # Try treating it as a JuMP container
            @printf "\nBranch flows from substation at t=1 (kW):\n"
            for (i, j) in data[:L1set]
                @printf "  P[(%d,%d)] = %.2f kW\n" i j (value(P_vals[(i,j), 1]) * P_BASE)
            end
            
            # Trace power flow through node 2
            println("\nTracing power flow through node 2:")
            j = 2
            @printf "  Incoming: P[(1,2)] = %.2f kW\n" (value(P_vals[(1,2), 1]) * P_BASE)
            if !isempty(data[:children][j])
                @printf "  Outgoing:\n"
                for k in data[:children][j]
                    @printf "    P[(%d,%d)] = %.2f kW\n" j k (value(P_vals[(j,k), 1]) * P_BASE)
                end
            end
            # Check if node 2 has load, PV, or battery
            if haskey(data[:p_L_pu], (j, 1))
                @printf "  Load at node %d: %.2f kW\n" j (data[:p_L_pu][j, 1] * P_BASE)
            end
            if j in data[:Dset] && haskey(data[:p_D_pu], (j, 1))
                @printf "  PV at node %d: %.2f kW\n" j (data[:p_D_pu][j, 1] * P_BASE)
            end
            if j in data[:Bset]
                @printf "  Battery at node %d: %.2f kW\n" j (value(P_B_vals[j, 1]) * P_BASE)
            end
            
            # Trace power flow through node 6
            println("\nTracing power flow through node 6:")
            j = 6
            i = data[:parent][j]
            @printf "  Incoming: P[(%d,%d)] = %.2f kW\n" i j (value(P_vals[(i,j), 1]) * P_BASE)
            if !isempty(data[:children][j])
                @printf "  Outgoing:\n"
                for k in data[:children][j]
                    @printf "    P[(%d,%d)] = %.2f kW\n" j k (value(P_vals[(j,k), 1]) * P_BASE)
                end
            end
            if haskey(data[:p_L_pu], (j, 1))
                @printf "  Load at node %d: %.2f kW\n" j (data[:p_L_pu][j, 1] * P_BASE)
            end
            if j in data[:Dset] && haskey(data[:p_D_pu], (j, 1))
                @printf "  PV at node %d: %.2f kW\n" j (data[:p_D_pu][j, 1] * P_BASE)
            end
            if j in data[:Bset]
                @printf "  Battery at node %d: %.2f kW (discharging=positive)\n" j (value(P_B_vals[j, 1]) * P_BASE)
            end
            
            # Trace power flow through node 8
            println("\nTracing power flow through node 8:")
            j = 8
            i = data[:parent][j]
            @printf "  Incoming: P[(%d,%d)] = %.2f kW\n" i j (value(P_vals[(i,j), 1]) * P_BASE)
            if !isempty(data[:children][j])
                @printf "  Outgoing:\n"
                for k in data[:children][j]
                    @printf "    P[(%d,%d)] = %.2f kW\n" j k (value(P_vals[(j,k), 1]) * P_BASE)
                end
            end
            if haskey(data[:p_L_pu], (j, 1))
                @printf "  Load at node %d: %.2f kW\n" j (data[:p_L_pu][j, 1] * P_BASE)
            end
            if j in data[:Dset] && haskey(data[:p_D_pu], (j, 1))
                @printf "  PV at node %d: %.2f kW\n" j (data[:p_D_pu][j, 1] * P_BASE)
            end
            if j in data[:Bset]
                @printf "  Battery at node %d: %.2f kW (discharging=positive)\n" j (value(P_B_vals[j, 1]) * P_BASE)
            end
            
            # Trace power flow through node 5 (where battery is)
            println("\nTracing power flow through node 5 (Battery location):")
            j = 5
            i = data[:parent][j]
            @printf "  Incoming: P[(%d,%d)] = %.2f kW\n" i j (value(P_vals[(i,j), 1]) * P_BASE)
            if !isempty(data[:children][j])
                @printf "  Outgoing:\n"
                for k in data[:children][j]
                    @printf "    P[(%d,%d)] = %.2f kW\n" j k (value(P_vals[(j,k), 1]) * P_BASE)
                end
            end
            if haskey(data[:p_L_pu], (j, 1))
                @printf "  Load at node %d: %.2f kW\n" j (data[:p_L_pu][j, 1] * P_BASE)
            end
            if j in data[:Dset] && haskey(data[:p_D_pu], (j, 1))
                @printf "  PV at node %d: %.2f kW\n" j (data[:p_D_pu][j, 1] * P_BASE)
            end
            if j in data[:Bset]
                @printf "  Battery at node %d: %.2f kW (discharging=positive)\n" j (value(P_B_vals[j, 1]) * P_BASE)
            end
        end
    end
    
    if !isempty(data[:Bset]) && isa(P_B_vals, AbstractDict)
        for b in data[:Bset]
            P_B_kW = [P_B_vals[b, t] * P_BASE for t in data[:Tset]]
            @printf "Battery %d Power (kW): min=%.1f, max=%.1f\n" b minimum(P_B_kW) maximum(P_B_kW)
        end
    end
    
    println("\n--- BATTERY SOC SUMMARY ---")
    if !isempty(data[:Bset]) && isa(B_vals, AbstractDict)
        for b in data[:Bset]
            B_kWh = [B_vals[b, t] * E_BASE for t in data[:Tset]]
            B_R_kWh = data[:B_R][b] * E_BASE
            soc_percent = (B_kWh ./ B_R_kWh) .* 100
            @printf "Battery %d SOC (%%): min=%.1f%%, max=%.1f%%, final=%.1f%%\n" b minimum(soc_percent) maximum(soc_percent) soc_percent[end]
        end
    end
    
    println("\n--- VOLTAGE SUMMARY ---")
    if isa(v_vals, AbstractDict)
        # Extract voltage magnitudes (sqrt of squared voltages)
        v_mag = Dict()
        for n in data[:Nset], t in data[:Tset]
            v_mag[(n,t)] = sqrt(v_vals[n, t])
        end
        v_all = [v_mag[(n,t)] for n in data[:Nset], t in data[:Tset]]
        @printf "Voltage (pu): min=%.4f, max=%.4f\n" minimum(v_all) maximum(v_all)
    end
    
else
    println("⚠ Optimization failed or did not converge to optimality")
    println("Status: ", sol_ldf[:status])
end

println("\n" * "="^80)
println("MPOPF LINDISTFLOW BRUTE-FORCED SOLUTION COMPLETE")
println("="^80 * "\n")
