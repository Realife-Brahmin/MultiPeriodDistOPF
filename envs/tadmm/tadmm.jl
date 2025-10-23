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
using Crayons

# Include standalone utilities
includet("parse_opendss.jl")
includet("opendss_validator.jl")
includet("logger.jl")

# Define color schemes
const COLOR_SUCCESS = Crayon(foreground = :green, bold = true)
const COLOR_WARNING = Crayon(foreground = :yellow, bold = true)
const COLOR_ERROR = Crayon(foreground = :red, bold = true)
const COLOR_INFO = Crayon(foreground = :cyan, bold = true)
const COLOR_HIGHLIGHT = Crayon(foreground = :magenta, bold = true)
const COLOR_RESET = Crayon(reset = true)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# System and simulation parameters
systemName = "ads10A_1ph"
# systemName = "ads10_1ph"
# systemName = "ieee123_1ph"
T = 24  # Number of time steps
delta_t_h = 1.0  # Time step duration in hours

# Plotting settings
showPlots = false  # Set to true to display plots interactively

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
    
    # ZERO REACTIVE POWER MODEL (comment these 3 lines to enable reactive power)
    # q_L_pu = zeros(size(q_L_pu))  # All loads at unity power factor
    # Q[(i,j),t] will be fixed to 0 via constraints
    # q_D will be fixed to 0 via constraints (PV at unity power factor)
    
    # ========== 2. CREATE MODEL ==========
    model = Model()
    if solver == :gurobi
        set_optimizer(model, Gurobi.Optimizer)
        set_optimizer_attribute(model, "NonConvex", 2)
        set_optimizer_attribute(model, "OutputFlag", 1)
        set_optimizer_attribute(model, "DualReductions", 0)  # Better infeasibility diagnosis
    else
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 3)
        set_optimizer_attribute(model, "max_iter", 3000)
    end
    
    # ========== 3. DEFINE VARIABLES ==========
    @variable(model, P_Subs[t in Tset] >= 0)
    @variable(model, Q_Subs[t in Tset])  # Substation reactive power
    @variable(model, P[(i, j) in Lset, t in Tset])
    @variable(model, Q[(i, j) in Lset, t in Tset])
    @variable(model, v[j in Nset, t in Tset])
    @variable(model, P_B[j in Bset, t in Tset])
    @variable(model, B[j in Bset, t in Tset])
    # PV operates at unity power factor (q_D = 0)
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
        # Substation node
        @constraint(model,
            Q_Subs[t] - sum(Q[(j1, j), t] for (j1, j) in L1set) == 0,
            base_name = "ReactivePowerBalance_Substation_t$(t)")
        
        # Non-substation nodes
        for j in Nm1set
            i = parent[j]
            Q_ij_t = Q[(i, j), t]
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            
            q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
            # q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
            q_D_j_t = 0.0  # PV operates at unity power factor
            
            @constraint(model,
                sum_Qjk - Q_ij_t == q_D_j_t - q_L_j_t,
                # Q_ij_t - sum_Qjk - q_L_j_t + q_D_j_t == 0,
                base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
        end
        
        # ----- 5.3 KVL CONSTRAINTS (LinDistFlow) -----
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model,
                v[i, t] - v[j, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) == 0,
                base_name = "KVL_Branch_$(i)_$(j)_t$(t)")
        end
        
        # ----- 5.4 ZERO REACTIVE POWER MODEL -----
        # Fix all branch reactive power flows to zero (comment to enable reactive power)
        # for (i, j) in Lset
        #     @constraint(model, Q[(i, j), t] == 0.0,
        #         base_name = "ZeroQ_Branch_$(i)_$(j)_t$(t)")
        # end
        
        # ----- 5.5 VOLTAGE CONSTRAINTS -----
        # Fixed substation voltage (1.05 pu, squared)
        @constraint(model, v[j1, t] == 1.00^2,
            base_name = "FixedSubstationVoltage_t$(t)")
        
        # Voltage limits (all nodes)
        for j in Nset
            @constraint(model, Vminpu[j]^2 <= v[j, t] <= Vmaxpu[j]^2,
                base_name = "VoltageLimits_Node$(j)_t$(t)")
        end
        
        # ----- 5.5 PV REACTIVE POWER LIMITS -----
        # Commented out - PV operates at unity power factor (q_D = 0)
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
    
    # Create output directory for system-specific files
    processedData_dir = joinpath(@__DIR__, "processedData")
    system_folder = "$(systemName)_T$(T)"
    system_dir = joinpath(processedData_dir, system_folder)
    mkpath(system_dir)
    
    # Save model summary to file (with full constraint listing)
    model_file = joinpath(system_dir, "model_summary.txt")
    write_model_summary(model, data, model_file)
    
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
    
    # ========== 8. COMPUTE NODAL INJECTIONS ==========
    # Compute nodal injections p_j, q_j using two methods and verify they match
    if has_values(model)
        println("\n" * "="^80)
        println(COLOR_INFO, "COMPUTING AND VERIFYING NODAL INJECTIONS", COLOR_RESET)
        println("="^80)
        
        # Initialize storage for nodal injections (all nodes, all times)
        p_j_flow = Dict{Tuple{Int, Int}, Float64}()  # Flow-based: sum(P_jk) - P_ij
        p_j_gen = Dict{Tuple{Int, Int}, Float64}()   # Generation-based: P_B + p_D - p_L
        q_j_flow = Dict{Tuple{Int, Int}, Float64}()  # Flow-based: sum(Q_jk) - Q_ij
        q_j_gen = Dict{Tuple{Int, Int}, Float64}()   # Generation-based: Q_B + q_D - q_L
        
        P_vals = value.(P)
        Q_vals = value.(Q)
        P_B_vals = value.(P_B)
        q_D_vals = value.(q_D)
        
        max_p_diff = 0.0
        max_q_diff = 0.0
        
        for t in Tset
            # Substation node (j=1)
            j = 1
            # Flow-based: outgoing flows (no incoming to substation)
            p_j_flow[(j, t)] = -sum(P_vals[(j, k), t] for k in children[j])
            q_j_flow[(j, t)] = -sum(Q_vals[(j, k), t] for k in children[j])
            
            # Generation-based: substation supplies what's needed
            # For substation: P_Subs - (negative sign because it's supplying)
            p_j_gen[(j, t)] = -value(P_Subs[t])
            q_j_gen[(j, t)] = -value(Q_Subs[t])
            
            # Non-substation nodes
            for j in Nm1set
                i = parent[j]
                
                # Flow-based injection: sum(outgoing) - incoming
                # Positive injection = net power leaving node (generation > consumption)
                P_jk_sum = isempty(children[j]) ? 0.0 : sum(P_vals[(j, k), t] for k in children[j])
                Q_jk_sum = isempty(children[j]) ? 0.0 : sum(Q_vals[(j, k), t] for k in children[j])
                
                p_j_flow[(j, t)] = P_jk_sum - P_vals[(i, j), t]
                q_j_flow[(j, t)] = Q_jk_sum - Q_vals[(i, j), t]
                
                # Generation-based injection: generation - consumption
                # P_B > 0: discharging (generation), P_B < 0: charging (consumption)
                # p_D > 0: PV generation
                # p_L > 0: load consumption
                p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0
                q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0
                p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0
                q_D_j_t = (j in Dset) ? q_D_vals[j, t] : 0.0
                P_B_j_t = (j in Bset) ? P_B_vals[j, t] : 0.0
                # Q_B = 0 (batteries don't provide reactive power)
                
                p_j_gen[(j, t)] = P_B_j_t + p_D_j_t - p_L_j_t
                q_j_gen[(j, t)] = q_D_j_t - q_L_j_t
            end
        end
        
        # Verify that both methods match
        println("\nVerifying nodal injection calculations...")
        all_match = true
        tolerance = 1e-6
        
        for j in Nset, t in Tset
            p_diff = abs(p_j_flow[(j, t)] - p_j_gen[(j, t)])
            q_diff = abs(q_j_flow[(j, t)] - q_j_gen[(j, t)])
            
            max_p_diff = max(max_p_diff, p_diff)
            max_q_diff = max(max_q_diff, q_diff)
            
            if p_diff > tolerance || q_diff > tolerance
                all_match = false
                print(COLOR_WARNING)
                @printf "  ⚠ Mismatch at node %d, t=%d: Δp=%.6e, Δq=%.6e\n" j t p_diff q_diff
                print(COLOR_RESET)
            end
        end
        
        if all_match
            print(COLOR_SUCCESS)
            println("✓ All nodal injections verified! Both methods match.")
            print(COLOR_RESET)
            @printf "  Max real power difference: %.6e pu\n" max_p_diff
            @printf "  Max reactive power difference: %.6e pu\n" max_q_diff
            
            # Since they match, add to solution dictionary (use flow-based values)
            result[:p_j] = p_j_flow
            result[:q_j] = q_j_flow
            
            print(COLOR_INFO)
            println("\nNodal injections added to solution dictionary.")
            print(COLOR_RESET)
        else
            print(COLOR_ERROR)
            println("✗ WARNING: Nodal injection methods DO NOT match!")
            println("  Not adding p_j and q_j to solution dictionary.")
            print(COLOR_RESET)
            @printf "  Max real power difference: %.6e pu\n" max_p_diff
            @printf "  Max reactive power difference: %.6e pu\n" max_q_diff
        end
        
        println("="^80)
    end
    
    return result
end

# =============================================================================
# SOLVE AND REPORT
# =============================================================================

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH LINDISTFLOW (BRUTE-FORCED)", COLOR_RESET)
println("="^80)

sol_ldf_bf = solve_MPOPF_with_LinDistFlow_BruteForced(data; solver=:gurobi)

# Report results
println("\n--- SOLUTION STATUS ---")
print("Status: ")

if sol_ldf_bf[:status] == MOI.OPTIMAL || sol_ldf_bf[:status] == MOI.LOCALLY_SOLVED
    print(COLOR_SUCCESS)
    println(sol_ldf_bf[:status])
    println("✓ Optimization successful!")
    print(COLOR_RESET)
    println("\n--- OBJECTIVE VALUE ---")
    @printf "Total Cost: \$%.2f\n" sol_ldf_bf[:objective]
    
    # Extract solution arrays
    P_Subs_vals = sol_ldf_bf[:P_Subs]
    P_B_vals = sol_ldf_bf[:P_B]
    B_vals = sol_ldf_bf[:B]
    v_vals = sol_ldf_bf[:v]
    
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
    # Always run voltage diagnostics (v_vals is a JuMP container, not AbstractDict)
    v_mag = Dict()
    for n in data[:Nset], t in data[:Tset]
        v_mag[(n,t)] = sqrt(v_vals[n, t])
    end
    v_all = [v_mag[(n,t)] for n in data[:Nset], t in data[:Tset]]
    @printf "Voltage (pu): min=%.4f, max=%.4f\n" minimum(v_all) maximum(v_all)
    
    # KVL VERIFICATION: Check if KVL constraints are satisfied post-optimization
    P_vals = sol_ldf_bf[:P]
    Q_vals = sol_ldf_bf[:Q]
    
    println("\n" * "="^80)
    println(COLOR_INFO, "KVL CONSTRAINT VERIFICATION (across all time periods)", COLOR_RESET)
    println("="^80)
    local max_kvl_violation = 0.0
    local total_violations = 0
    violation_threshold = 1e-6
    
    for t in data[:Tset]
        for (i, j) in data[:Lset]
            r_ij = data[:rdict_pu][(i, j)]
            x_ij = data[:xdict_pu][(i, j)]
            # LHS: v[i,t] - v[j,t]
            lhs = v_vals[i, t] - v_vals[j, t]
            # RHS: 2*(r*P + x*Q)
            rhs = 2 * (r_ij * P_vals[(i, j), t] + x_ij * Q_vals[(i, j), t])
            # Violation: |LHS - RHS|
            violation = abs(lhs - rhs)
            
            if violation > max_kvl_violation
                max_kvl_violation = violation
            end
            
            if violation > violation_threshold
                total_violations += 1
                if total_violations <= 10  # Print first 10 violations
                    print(COLOR_WARNING)
                    @printf "  ⚠ KVL violation at t=%d, branch (%d,%d): |%.8f - %.8f| = %.2e\n" t i j lhs rhs violation
                    print(COLOR_RESET)
                end
            end
        end
    end
    
    if total_violations == 0
        print(COLOR_SUCCESS)
        println("✓ All KVL constraints satisfied (max violation: $(max_kvl_violation))")
        print(COLOR_RESET)
    else
        print(COLOR_WARNING)
        println("⚠ Total KVL violations (>$(violation_threshold)): $(total_violations)")
        println("  Maximum violation: $(max_kvl_violation)")
        print(COLOR_RESET)
    end
    println("="^80)
    
else
    print(COLOR_ERROR)
    println("⚠ Optimization failed or did not converge to optimality")
    println("Status: ", sol_ldf_bf[:status])
    print(COLOR_RESET)
end

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "MPOPF LINDISTFLOW BRUTE-FORCED SOLUTION COMPLETE", COLOR_RESET)
println("="^80)

# =============================================================================
# PLOTTING
# =============================================================================

# Include plotting utilities
include("Plotter.jl")

println("\n" * "="^80)
println(COLOR_INFO, "GENERATING PLOTS", COLOR_RESET)
println("="^80)

# Create output directories (use absolute path like model writing section)
processedData_dir = joinpath(@__DIR__, "processedData")
system_folder = "$(systemName)_T$(T)"
system_dir = joinpath(processedData_dir, system_folder)
mkpath(processedData_dir)  # Create processedData folder
mkpath(system_dir)  # Create system-specific subfolder with horizon

# Plot input curves (load, PV, cost) - save in processedData folder
input_curves_path = joinpath(processedData_dir, "input_curves.png")
plot_input_curves(data, showPlots=showPlots, savePlots=true, filename=input_curves_path)

# Plot battery actions (only if optimization was successful and batteries exist)
# Save in system-specific subfolder
if (sol_ldf_bf[:status] == MOI.OPTIMAL || sol_ldf_bf[:status] == MOI.LOCALLY_SOLVED) && !isempty(data[:Bset])
    battery_actions_path = joinpath(system_dir, "battery_actions_lindistflow.png")
    plot_battery_actions(sol_ldf_bf, data, "LinDistFlow-Gurobi", 
                        showPlots=showPlots, savePlots=true, 
                        filename=battery_actions_path)
else
    if isempty(data[:Bset])
        println("No batteries in system, skipping battery actions plot")
    end
end

# Plot substation power and cost (only if optimization was successful)
if (sol_ldf_bf[:status] == MOI.OPTIMAL || sol_ldf_bf[:status] == MOI.LOCALLY_SOLVED)
    subs_power_cost_path = joinpath(system_dir, "substation_power_cost.png")
    plot_substation_power_and_cost(sol_ldf_bf, data, "LinDistFlow-Gurobi",
                                   showPlots=showPlots, savePlots=true,
                                   filename=subs_power_cost_path)
    
    # Plot voltage profile for last bus
    voltage_one_bus_path = joinpath(system_dir, "voltage_profile_last_bus.png")
    plot_voltage_profile_one_bus(sol_ldf_bf, data, "LinDistFlow-Gurobi",
                                showPlots=showPlots, savePlots=true,
                                filename=voltage_one_bus_path)
    
    # Plot voltage profile for all buses at middle time step + create GIF
    voltage_all_buses_path = joinpath(system_dir, "voltage_profile_all_buses.png")
    voltage_gif_path = joinpath(system_dir, "voltage_animation.gif")
    plot_voltage_profile_all_buses(sol_ldf_bf, data, "LinDistFlow-Gurobi",
                                   showPlots=showPlots, savePlots=true,
                                   filename=voltage_all_buses_path,
                                   create_gif=true,
                                   gif_filename=voltage_gif_path)
    
    # Plot PV power (p_D and q_D) if PV exists
    if !isempty(data[:Dset])
        pv_power_path = joinpath(system_dir, "pv_power.png")
        plot_pv_power(sol_ldf_bf, data, "LinDistFlow-Gurobi",
                     pv_index=1,  # Plot first PV (or only PV)
                     showPlots=showPlots, savePlots=true,
                     filename=pv_power_path)
        
        # Create PV power circle GIF
        pv_circle_gif_path = joinpath(system_dir, "pv_power_circle.gif")
        plot_pv_power_circle_gif(sol_ldf_bf, data, "LinDistFlow-Gurobi",
                                pv_index=1,
                                showPlots=showPlots, savePlots=true,
                                filename=pv_circle_gif_path)
    else
        println("No PV units in system, skipping PV power plots")
    end
end

println("\n" * "="^80)
println(COLOR_SUCCESS, "PLOTTING COMPLETE", COLOR_RESET)
println("="^80)
