#=
tADMM: Temporal ADMM for Multi-Period Optimal Power Flow
Standalone implementation with minimal dependencies
=#

# %%
begin
# Activate the tadmm environment
import Pkg
env_path = @__DIR__
Pkg.activate(env_path)

using Revise
using JuMP
using Ipopt
using Gurobi
using OpenDSSDirect
using Printf
using Statistics
using Parameters: @unpack
using Plots

# Include standalone utilities
includet("parse_opendss.jl")
includet("opendss_validator.jl")

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# System and simulation parameters
# systemName = "ads10A_1ph"
systemName = "ads10_1ph"
# systemName = "ieee123_1ph"
T = 24  # Number of time steps
delta_t_h = 1.0  # Time step duration in hours

# Load shapes
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2 # $/kWh time-varying energy cost

# Solar PV profile: peaks at noon (t=12), zero at night
# t=1 is 1AM, t=12 is 12PM (noon), t=24 is 12AM (midnight)
# Solar generation from roughly 6AM (t=6) to 6PM (t=18)
LoadShapePV = [max(0.0, sin(π * (t - 6) / 12)) for t in 1:T]  # Sine curve from 6AM to 6PM
LoadShapePV = LoadShapePV ./ maximum(LoadShapePV)  # Normalize to [0, 1]

C_B = 1e-6 * minimum(LoadShapeCost)  # Battery quadratic cost coefficient

# =============================================================================
# PARSE SYSTEM DATA
# =============================================================================

println("\n" * "="^80)
println("PARSING SYSTEM DATA")
println("="^80)

data = parse_system_from_dss(
    systemName, 
    T;
    LoadShapeLoad=LoadShapeLoad,
    LoadShapeCost=LoadShapeCost,
    LoadShapePV=LoadShapePV,
    C_B=C_B,
    delta_t_h=delta_t_h
)

# Validate with OpenDSS powerflow
print_powerflow_summary(data)
end
# =============================================================================
# MPOPF SOLVER WITH LINDISTFLOW
# =============================================================================
function solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:gurobi)
    
    # ========== 1. UNPACK DATA ==========
    @unpack Nset, Lset, Dset, Bset, Tset, NLset, Nm1set = data
    @unpack substationBus, L1set, Lm1set, parent, children = data
    @unpack p_L_pu, p_D_pu, q_L_pu = data
    @unpack P_B_R_pu, B_R_pu, B0_pu, S_D_R = data
    @unpack rdict_pu, xdict_pu = data
    @unpack Vminpu, Vmaxpu = data
    @unpack LoadShapeCost, C_B, delta_t_h, soc_min, soc_max = data
    @unpack kVA_B = data
    
    Δt = delta_t_h
    P_BASE = kVA_B
    j1 = substationBus
    
    # ========== 2. CREATE MODEL ==========
    model = Model()
    if solver == :gurobi
        set_optimizer(model, Gurobi.Optimizer)
        set_optimizer_attribute(model, "NonConvex", 2)
        set_optimizer_attribute(model, "OutputFlag", 1)
    else
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 3)
        set_optimizer_attribute(model, "max_iter", 3000)
    end
    
    # ========== 3. DEFINE VARIABLES ==========
    @variable(model, P_Subs[t in Tset] >= 0)
    @variable(model, P[(i, j) in Lset, t in Tset])
    @variable(model, Q[(i, j) in Lset, t in Tset])
    @variable(model, v[j in Nset, t in Tset])
    @variable(model, P_B[j in Bset, t in Tset])
    @variable(model, B[j in Bset, t in Tset])
    @variable(model, q_D[j in Dset, t in Tset])
    
    # ========== 4. OBJECTIVE FUNCTION ==========
    @expression(model, energy_cost, 
        sum(LoadShapeCost[t] * P_Subs[t] * P_BASE * Δt for t in Tset))
    @expression(model, battery_cost, 
        sum(C_B * (P_B[j, t] * P_BASE)^2 * Δt for j in Bset, t in Tset))
    # @expression(model, battery_cost,
    #     0.0)  # Temporarily set to zero for debugging
    @objective(model, Min, energy_cost + battery_cost)
    
    # ========== 5. CONSTRAINTS ==========
    
    for t in Tset
        # ----- 5.1 NODAL REAL POWER BALANCE -----
        # Substation node
        @constraint(model, 
            P_Subs[t] - sum(P[(j1, j), t] for (j1, j) in L1set) == 0,
            base_name = "RealPowerBalance_Substation_t$(t)")
        
        # Non-substation nodes (ALL nodes in Nm1set, not just load nodes)
        for j in Nm1set
            i = parent[j]
            P_ij_t = P[(i, j), t]
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P[(j, k), t] for k in children[j])
            
            p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0
            p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0
            P_B_j_t = (j in Bset) ? P_B[j, t] : 0.0
            
            # Power balance: Incoming = Outgoing + Load - PV - Battery
            @constraint(model,
                sum_Pjk - P_ij_t == P_B_j_t + p_D_j_t - p_L_j_t,
                # P_ij_t - sum_Pjk - p_L_j_t + p_D_j_t + P_B_j_t == 0,
                base_name = "RealPowerBalance_Node$(j)_t$(t)")
        end
        
        # ----- 5.2 NODAL REACTIVE POWER BALANCE -----
        # # Substation node
        # @constraint(model,
        #     sum(Q[(j1, j), t] for (j1, j) in L1set) == 0,
        #     base_name = "ReactivePowerBalance_Substation_t$(t)")
        
        # # Non-substation nodes
        # for j in Nm1set
        #     i = parent[j]
        #     Q_ij_t = Q[(i, j), t]
        #     sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            
        #     q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
        #     q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
            
        #     @constraint(model,
        #         sum_Qjk - Q_ij_t == q_D_j_t - q_L_j_t,
        #         # Q_ij_t - sum_Qjk - q_L_j_t + q_D_j_t == 0,
        #         base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
        # end
        
        # ----- 5.3 KVL CONSTRAINTS (LinDistFlow) -----
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model,
                v[i, t] - v[j, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) == 0,
                base_name = "KVL_Branch_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.4 VOLTAGE CONSTRAINTS -----
        # Fixed substation voltage (1.05 pu, squared)
        @constraint(model, v[j1, t] == 1.00^2,
            base_name = "FixedSubstationVoltage_t$(t)")
        
        # Voltage limits (all nodes)
        for j in Nset
            @constraint(model, Vminpu[j]^2 <= v[j, t] <= Vmaxpu[j]^2,
                base_name = "VoltageLimits_Node$(j)_t$(t)")
        end
        
        # ----- 5.5 PV REACTIVE POWER LIMITS -----
        for j in Dset
            p_D_val = p_D_pu[j, t] # not DV
            S_D_R_val = S_D_R[j] # not DV
            q_max_t = sqrt(S_D_R_val^2 - p_D_val^2) # fixed as well
            @constraint(model, -q_max_t <= q_D[j, t] <= q_max_t,
                base_name = "PVReactiveLimits_DER$(j)_t$(t)")
        end
        
        # ----- 5.6 BATTERY CONSTRAINTS -----
        for j in Bset
            # SOC trajectory
            # P_B > 0: discharging (battery supplies power like generator) → SOC decreases
            # P_B < 0: charging (battery absorbs power like load) → SOC increases
            if t == 1
                @constraint(model, B[j, t] == B0_pu[j] - P_B[j, t] * Δt,
                    base_name = "BatterySOC_Init_$(j)_t$(t)")
            else
                @constraint(model, B[j, t] == B[j, t-1] - P_B[j, t] * Δt,
                    base_name = "BatterySOC_$(j)_t$(t)")
            end
            
            # SOC limits (use per-unit values)
            @constraint(model, soc_min[j] * B_R_pu[j] <= B[j, t] <= soc_max[j] * B_R_pu[j],
                base_name = "BatterySOCLimits_$(j)_t$(t)")

            # Power limits (use per-unit values)
            @constraint(model, -P_B_R_pu[j] <= P_B[j, t] <= P_B_R_pu[j],
                base_name = "BatteryPowerLimits_$(j)_t$(t)")
        end
    end
    
    # ========== 6. SOLVE ==========
    
    # Save model summary to file (with full constraint listing)
    model_file = joinpath(@__DIR__, "model_summary.txt")
    open(model_file, "w") do io
        println(io, "="^80)
        println(io, "MODEL SUMMARY")
        println(io, "="^80)
        
        # Manually print all constraints in organized order
        println(io, "\n" * "="^80)
        println(io, "ALL CONSTRAINTS (ORGANIZED BY TYPE AND TIME)")
        println(io, "="^80)
        
        # Define constraint name patterns in desired order
        constraint_groups = [
            ("SUBSTATION POWER BALANCE", r"RealPowerBalance_Substation_t\d+"),
            ("NODAL REAL POWER BALANCE", r"RealPowerBalance_Node\d+_t\d+"),
            ("NODAL REACTIVE POWER BALANCE (Substation)", r"ReactivePowerBalance_Substation_t\d+"),
            ("NODAL REACTIVE POWER BALANCE", r"ReactivePowerBalance_Node\d+_t\d+"),
            ("KVL CONSTRAINTS", r"KVL_Branch_\d+_\d+_t\d+"),
            ("FIXED SUBSTATION VOLTAGE", r"FixedSubstationVoltage_t\d+"),
            ("BATTERY SOC TRAJECTORY", r"BatterySOC.*_t\d+"),
            ("VOLTAGE LIMITS", r"VoltageLimits_Node\d+_t\d+"),
            ("PV REACTIVE POWER LIMITS", r"PVReactiveLimits_DER\d+_t\d+"),
            ("BATTERY SOC LIMITS", r"BatterySOCLimits_\d+_t\d+"),
            ("BATTERY POWER LIMITS", r"BatteryPowerLimits_\d+_t\d+"),
        ]
        
        # Collect all constraint references
        all_constraint_refs = []
        for (F, S) in list_of_constraint_types(model)
            for con in all_constraints(model, F, S)
                push!(all_constraint_refs, (name(con), con))
            end
        end
        
        # Print constraints by group
        for (group_name, pattern) in constraint_groups
            println(io, "\n" * "-"^80)
            println(io, group_name)
            println(io, "-"^80)
            
            # Filter and sort constraints matching this pattern
            matching_cons = filter(x -> occursin(pattern, x[1]), all_constraint_refs)
            sort!(matching_cons, by=x -> x[1])  # Alphabetical sort within group
            
            for (con_name, con) in matching_cons
                println(io, con_name, ": ", con)
            end
            
            println(io, "Total: $(length(matching_cons)) constraints")
        end
        
        # Print any remaining constraints not captured above
        printed_names = Set{String}()
        for (group_name, pattern) in constraint_groups
            for (con_name, con) in all_constraint_refs
                if occursin(pattern, con_name)
                    push!(printed_names, con_name)
                end
            end
        end
        
        remaining = filter(x -> !(x[1] in printed_names), all_constraint_refs)
        if !isempty(remaining)
            println(io, "\n" * "-"^80)
            println(io, "OTHER CONSTRAINTS")
            println(io, "-"^80)
            for (con_name, con) in remaining
                println(io, con_name, ": ", con)
            end
        end
        
        println(io, "\n--- VARIABLE BOUNDS CHECK ---")
        println(io, "Number of variables: ", num_variables(model))
        println(io, "Number of constraints: ", num_constraints(model; count_variable_in_set_constraints=true))
        
        # Check for potential issues
        println(io, "\n--- DATA VALIDATION ---")
        println(io, "Load nodes (NLset): ", NLset)
        println(io, "PV nodes (Dset): ", Dset)
        println(io, "Battery nodes (Bset): ", Bset)
        println(io, "All nodes (Nset): ", Nset)
        println(io, "Non-substation nodes (Nm1set): ", Nm1set)
        println(io, "Substation branches (L1set): ", L1set)
        println(io, "Other branches (Lm1set): ", Lm1set)
        
        # Check power balance feasibility
        println(io, "\n--- POWER BALANCE CHECK (t=1) ---")
        total_load_t1 = sum(p_L_pu[j, 1] for j in NLset)
        total_pv_t1 = sum(p_D_pu[j, 1] for j in Dset)
        println(io, "Total load at t=1: $(total_load_t1 * P_BASE) kW")
        println(io, "Total PV at t=1: $(total_pv_t1 * P_BASE) kW")
        println(io, "Net demand: $((total_load_t1 - total_pv_t1) * P_BASE) kW")
        
        # Check battery parameters
        if !isempty(Bset)
            println(io, "\n--- BATTERY PARAMETERS ---")
            for b in Bset
                println(io, "Battery $b:")
                println(io, "  Initial SOC: $(B0_pu[b] * P_BASE) kWh")
                println(io, "  Capacity: $(B_R_pu[b] * P_BASE) kWh")
                println(io, "  Power rating: $(P_B_R_pu[b] * P_BASE) kW")
                println(io, "  SOC min: $(soc_min[b] * 100)%")
                println(io, "  SOC max: $(soc_max[b] * 100)%")
                println(io, "  Min energy: $(soc_min[b] * B_R_pu[b] * P_BASE) kWh")
                println(io, "  Max energy: $(soc_max[b] * B_R_pu[b] * P_BASE) kWh")
            end
        end
    end
    println("Model summary saved to: $model_file")
    
    println("\n" * "="^80)
    println("STARTING OPTIMIZATION")
    println("="^80)
    
    optimize!(model)
    
    # ========== 7. EXTRACT RESULTS ==========
    status = termination_status(model)
    obj_val = has_values(model) ? objective_value(model) : NaN
    
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

sol_ldf = solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:gurobi)

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
    
    # Convert to physical units
    P_BASE = data[:kVA_B]
    E_BASE = P_BASE * 1.0
    
    println("\n--- POWER SUMMARY ---")
    if isa(P_Subs_vals, AbstractArray)
        P_Subs_kW = P_Subs_vals .* P_BASE
        @printf "Substation Power (kW): min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_Subs_kW) maximum(P_Subs_kW) mean(P_Subs_kW)
    end
    
    if !isempty(data[:Bset]) && isa(P_B_vals, AbstractDict)
        println("\n--- BATTERY POWER ---")
        for b in data[:Bset]
            P_B_kW = [P_B_vals[b, t] * P_BASE for t in data[:Tset]]
            @printf "Battery %d Power (kW): min=%.1f, max=%.1f\n" b minimum(P_B_kW) maximum(P_B_kW)
        end
    end
    
    if !isempty(data[:Bset]) && isa(B_vals, AbstractDict)
        println("\n--- BATTERY SOC ---")
        for b in data[:Bset]
            B_kWh = [B_vals[b, t] * E_BASE for t in data[:Tset]]
            B_R_kWh = data[:B_R][b]
            soc_percent = (B_kWh ./ B_R_kWh) .* 100
            @printf "Battery %d SOC (%%): min=%.1f%%, max=%.1f%%, final=%.1f%%\n" b minimum(soc_percent) maximum(soc_percent) soc_percent[end]
        end
    end
    
    println("\n--- VOLTAGE SUMMARY ---")
    if isa(v_vals, AbstractDict)
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
println("="^80)

# =============================================================================
# PLOTTING
# =============================================================================

# Include plotting utilities
include("Plotter.jl")

println("\n" * "="^80)
println("GENERATING PLOTS")
println("="^80)

# Create output directories
processedData_dir = "processedData"
system_dir = joinpath(processedData_dir, systemName)
mkpath(processedData_dir)  # Create processedData folder
mkpath(system_dir)  # Create system-specific subfolder

# Plot input curves (load, PV, cost) - save in processedData folder
input_curves_path = joinpath(processedData_dir, "input_curves.png")
plot_input_curves(data, showPlots=true, savePlots=true, filename=input_curves_path)

# Plot battery actions (only if optimization was successful and batteries exist)
# Save in system-specific subfolder
if (sol_ldf[:status] == MOI.OPTIMAL || sol_ldf[:status] == MOI.LOCALLY_SOLVED) && !isempty(data[:Bset])
    battery_actions_path = joinpath(system_dir, "battery_actions_lindistflow.png")
    plot_battery_actions(sol_ldf, data, "LinDistFlow-Gurobi", 
                        showPlots=true, savePlots=true, 
                        filename=battery_actions_path)
else
    if isempty(data[:Bset])
        println("No batteries in system, skipping battery actions plot")
    end
end

println("\n" * "="^80)
println("PLOTTING COMPLETE")
println("="^80)
