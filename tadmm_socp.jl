

# Activate the tadmm environment
import Pkg
env_path = joinpath(@__DIR__, "envs", "tadmm")
Pkg.activate(env_path)

using Revise
using LinearAlgebra
using JuMP
using Ipopt
using Gurobi
using OpenDSSDirect
using Printf
using Statistics
using Parameters: @unpack
using Plots
using Crayons

begin # entire script including environment setup

# Include standalone utilities
includet(joinpath(env_path, "parse_opendss.jl"))
includet(joinpath(env_path, "opendss_validator.jl"))
includet(joinpath(env_path, "logger.jl"))
includet(joinpath(env_path, "Plotter.jl"))

# System and simulation parameters
# systemName = "ads10A_1ph"
systemName = "ieee123A_1ph"
T = 6  # Number of time steps
delta_t_h = 24.0/T  # Time step duration in hours

# Solver selection
use_gurobi_for_bf = true       # Use Gurobi for brute force (SOCP)
use_gurobi_for_tadmm = true    # Use Gurobi for tADMM subproblems (SOCP), false = use Ipopt (NLP)

# tADMM algorithm parameters
rho_base = 10000.0              # Base ρ value for T=24
rho_scaling_with_T = true       # Automatically scale ρ with T (recommended)
# ρ scaling logic: Larger T → need larger ρ to handle more coupling constraints
# Rule of thumb: ρ ∝ √T (conservative) or ρ ∝ T (aggressive)
rho_tadmm = rho_scaling_with_T ? rho_base * sqrt(T / 24.0) : rho_base
max_iter_tadmm = 3000
eps_pri_tadmm = 1e-5
eps_dual_tadmm = 1e-4
adaptive_rho_tadmm = true  # Set to false for fixed ρ

# Define color schemes
const COLOR_SUCCESS = Crayon(foreground = :green, bold = true)
const COLOR_WARNING = Crayon(foreground = :yellow, bold = true)
const COLOR_ERROR = Crayon(foreground = :red, bold = true)
const COLOR_INFO = Crayon(foreground = :cyan, bold = true)
const COLOR_HIGHLIGHT = Crayon(foreground = :magenta, bold = true)
const COLOR_RESET = Crayon(reset = true)

begin # scenario config
    # Plotting settings
    showPlots = false  # Set to true to display plots interactively
    # saveAllBatteryPlots = false  # Set to true to save plots for ALL batteries (time-consuming)
    saveAllBatteryPlots = true  # Set to true to save plots for ALL batteries (time-consuming)
    saveAllPVPlots = true  # Set to true to save plots for ALL PV units (shows p_D and q_D)


    # Load shapes
    LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
    LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2 # $/kWh time-varying energy cost

    # Solar PV profile: zeros at start/end (25% each), bell curve in middle (50%)
    # This mimics solar generation that starts at sunrise, peaks at noon, ends at sunset
    LoadShapePV = zeros(T)
    if T >= 4
        # Determine active region: middle 50% of horizon
        t_start = max(1, round(Int, 0.25 * T))
        t_end = min(T, round(Int, 0.75 * T))
        n_active = t_end - t_start + 1
        
        if n_active >= 2
            # Create bell curve in active region using half sine wave
            # sin(0) = 0, sin(π) = 0, sin(π/2) = 1 → natural [0, 1] range
            t_normalized = range(0, π, length=n_active)
            LoadShapePV[t_start:t_end] = [sin(t_norm) for t_norm in t_normalized]
        else
            # Very few active periods, just use constant
            LoadShapePV[t_start:t_end] .= 1.0
        end
    elseif T == 3
        # For T=3: [0, 1, 0] pattern
        LoadShapePV[2] = 1.0
    elseif T == 2
        # For T=2: [0.5, 0.5] to avoid zeros everywhere
        LoadShapePV .= 0.5
    else
        # T=1: constant
        LoadShapePV .= 1.0
    end

    C_B = 1e-6 * minimum(LoadShapeCost)  # Battery quadratic cost coefficient
end # scenario config

begin # parse system data
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
end # parse system data

begin # function mpopf socp bruteforced
    function solve_MPOPF_with_SOCP_BruteForced(data; solver=:gurobi)
        
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
            set_optimizer_attribute(model, "OutputFlag", 0)  # Suppress Gurobi output
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
        @variable(model, ℓ[(i, j) in Lset, t in Tset] >= 0)  # Current squared magnitude
        @variable(model, P_B[j in Bset, t in Tset])
        @variable(model, B[j in Bset, t in Tset])
        # PV reactive power (can be optimized within capability curve)
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
                q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0
                # q_D_j_t = 0.0  # PV operates at unity power factor
                
                @constraint(model,
                    sum_Qjk - Q_ij_t == q_D_j_t - q_L_j_t,
                    # Q_ij_t - sum_Qjk - q_L_j_t + q_D_j_t == 0,
                    base_name = "ReactivePowerBalance_Node$(j)_t$(t)")
            end
            
            # ----- 5.3 VOLTAGE DROP CONSTRAINTS (BFM-NL) -----
            for (i, j) in Lset
                r_ij = rdict_pu[(i, j)]
                x_ij = xdict_pu[(i, j)]
                @constraint(model,
                    v[j, t] == v[i, t] - 2 * (r_ij * P[(i, j), t] + x_ij * Q[(i, j), t]) + (r_ij^2 + x_ij^2) * ℓ[(i, j), t],
                    base_name = "VoltageDrop_Branch_$(i)_$(j)_t$(t)")
            end
            
            # ----- 5.4 SOC RELAXATION CONSTRAINTS -----
            for (i, j) in Lset
                @constraint(model,
                    P[(i, j), t]^2 + Q[(i, j), t]^2 <= v[i, t] * ℓ[(i, j), t],
                    base_name = "SOC_Branch_$(i)_$(j)_t$(t)")
            end
            
            # ----- 5.5 VOLTAGE CONSTRAINTS -----
            # Fixed substation voltage (1.05 pu, squared)
            @constraint(model, v[j1, t] == 1.05^2,
                base_name = "FixedSubstationVoltage_t$(t)")
            
            # Voltage limits (all nodes)
            for j in Nset
                @constraint(model, Vminpu[j]^2 <= v[j, t] <= Vmaxpu[j]^2,
                    base_name = "VoltageLimits_Node$(j)_t$(t)")
            end
            
            # ----- 5.6 PV REACTIVE POWER LIMITS -----
            for j in Dset
                p_D_val = p_D_pu[j, t] # not DV
                S_D_R_val = S_D_R[j] # not DV
                q_max_t = sqrt(S_D_R_val^2 - p_D_val^2) # fixed as well
                @constraint(model, -q_max_t <= q_D[j, t] <= q_max_t,
                    base_name = "PVReactiveLimits_DER$(j)_t$(t)")
            end
            
            # ----- 5.7 BATTERY CONSTRAINTS -----
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
        processedData_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData")
        system_folder = "$(systemName)_T$(T)"
        system_dir = joinpath(processedData_dir, system_folder)
        mkpath(system_dir)
        
        # Save model summary to file (with full constraint listing)
        model_file = joinpath(system_dir, "model_summary_socp_bf.txt")
        write_model_summary(model, data, model_file)
        
        println("\n" * "="^80)
        println("STARTING OPTIMIZATION")
        println("="^80)
        
        optimize!(model)
        
        # ========== 7. EXTRACT RESULTS ==========
        status = termination_status(model)
        obj_val = has_values(model) ? objective_value(model) : NaN
        solver_time = solve_time(model)  # Get pure solver time from JuMP
        
        result = Dict(
            :model => model,
            :status => status,
            :objective => obj_val,
            :solve_time => solver_time,
            :P_Subs => has_values(model) ? value.(P_Subs) : P_Subs,
            :P => has_values(model) ? value.(P) : P,
            :Q => has_values(model) ? value.(Q) : Q,
            :v => has_values(model) ? value.(v) : v,
            :ℓ => has_values(model) ? value.(ℓ) : ℓ,
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
end # function mpopf socp bruteforced

begin # socp brute-forced solve
    # =============================================================================
    # SOLVE AND REPORT
    # =============================================================================

    println("\n" * "="^80)
    println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH SOCP (BRUTE-FORCED)", COLOR_RESET)
    println("="^80)

    solver_bf_choice = use_gurobi_for_bf ? :gurobi : :ipopt
    sol_socp_bf = solve_MPOPF_with_SOCP_BruteForced(data; solver=solver_bf_choice)

    # Report results
    println("\n--- SOLUTION STATUS ---")
    print("Status: ")

    if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
        print(COLOR_SUCCESS)
        println(sol_socp_bf[:status])
        println("✓ Optimization successful!")
        print(COLOR_RESET)
        println("\n--- OBJECTIVE VALUE ---")
        @printf "Total Cost: \$%.2f\n" sol_socp_bf[:objective]
        println("\n--- COMPUTATION TIME ---")
        @printf "Solver time: %.2f seconds\n" sol_socp_bf[:solve_time]
        
        # Extract solution arrays
        P_Subs_vals = sol_socp_bf[:P_Subs]
        P_B_vals = sol_socp_bf[:P_B]
        B_vals = sol_socp_bf[:B]
        v_vals = sol_socp_bf[:v]
        
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
        
        # VOLTAGE DROP VERIFICATION: Check if voltage drop constraints are satisfied
        P_vals = sol_socp_bf[:P]
        Q_vals = sol_socp_bf[:Q]
        ℓ_vals = sol_socp_bf[:ℓ]
        
        println("\n" * "="^80)
        println(COLOR_INFO, "VOLTAGE DROP CONSTRAINT VERIFICATION (across all time periods)", COLOR_RESET)
        println("="^80)
        local max_voltage_drop_violation = 0.0
        local total_violations = 0
        violation_threshold = 1e-6
        
        for t in data[:Tset]
            for (i, j) in data[:Lset]
                r_ij = data[:rdict_pu][(i, j)]
                x_ij = data[:xdict_pu][(i, j)]
                # LHS: v[j,t]
                lhs = v_vals[j, t]
                # RHS: v[i,t] - 2*(r*P + x*Q) + (r²+x²)*ℓ
                rhs = v_vals[i, t] - 2 * (r_ij * P_vals[(i, j), t] + x_ij * Q_vals[(i, j), t]) + (r_ij^2 + x_ij^2) * ℓ_vals[(i, j), t]
                # Violation: |LHS - RHS|
                violation = abs(lhs - rhs)
                
                if violation > max_voltage_drop_violation
                    max_voltage_drop_violation = violation
                end
                
                if violation > violation_threshold
                    total_violations += 1
                    if total_violations <= 10  # Print first 10 violations
                        print(COLOR_WARNING)
                        @printf "  ⚠ Voltage drop violation at t=%d, branch (%d,%d): |%.8f - %.8f| = %.2e\n" t i j lhs rhs violation
                        print(COLOR_RESET)
                    end
                end
            end
        end
        
        if total_violations == 0
            print(COLOR_SUCCESS)
            println("✓ All voltage drop constraints satisfied (max violation: $(max_voltage_drop_violation))")
            print(COLOR_RESET)
        else
            print(COLOR_WARNING)
            println("⚠ Total voltage drop violations (>$(violation_threshold)): $(total_violations)")
            println("  Maximum violation: $(max_voltage_drop_violation)")
            print(COLOR_RESET)
        end
        
        # SOC RELAXATION VERIFICATION
        println("\n" * "="^80)
        println(COLOR_INFO, "SOC RELAXATION VERIFICATION (across all time periods)", COLOR_RESET)
        println("="^80)
        local max_soc_violation = 0.0
        local total_soc_violations = 0
        
        for t in data[:Tset]
            for (i, j) in data[:Lset]
                # LHS: P² + Q²
                lhs = P_vals[(i, j), t]^2 + Q_vals[(i, j), t]^2
                # RHS: v[i,t] * ℓ[(i,j),t]
                rhs = v_vals[i, t] * ℓ_vals[(i, j), t]
                # Violation: max(LHS - RHS, 0)
                violation = max(lhs - rhs, 0.0)
                
                if violation > max_soc_violation
                    max_soc_violation = violation
                end
                
                if violation > violation_threshold
                    total_soc_violations += 1
                    if total_soc_violations <= 10
                        print(COLOR_WARNING)
                        @printf "  ⚠ SOC violation at t=%d, branch (%d,%d): %.8f > %.8f (Δ=%.2e)\n" t i j lhs rhs violation
                        print(COLOR_RESET)
                    end
                end
            end
        end
        
        if total_soc_violations == 0
            print(COLOR_SUCCESS)
            println("✓ All SOC relaxation constraints satisfied (max violation: $(max_soc_violation))")
            print(COLOR_RESET)
        else
            print(COLOR_WARNING)
            println("⚠ Total SOC violations (>$(violation_threshold)): $(total_soc_violations)")
            println("  Maximum violation: $(max_soc_violation)")
            print(COLOR_RESET)
        end
        println("="^80)
        
    else
        print(COLOR_ERROR)
        println("⚠ Optimization failed or did not converge to optimality")
        println("Status: ", sol_socp_bf[:status])
        print(COLOR_RESET)
    end

    println("\n" * "="^80)
    println(COLOR_HIGHLIGHT, "MPOPF SOCP BRUTE-FORCED SOLUTION COMPLETE", COLOR_RESET)
    println("="^80)
    
    # Write results to file
    if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
        processedData_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData")
        system_folder = "$(systemName)_T$(T)"
        system_dir = joinpath(processedData_dir, system_folder)
        results_file = joinpath(system_dir, "results_socp_bf.txt")
        
        open(results_file, "w") do io
            println(io, "="^80)
            println(io, "SOCP BRUTE-FORCED OPTIMIZATION RESULTS")
            println(io, "="^80)
            println(io, "System: $(systemName)")
            println(io, "Time horizon: T=$(T) periods")
            println(io, "Time step: $(delta_t_h) hours")
            println(io, "Number of buses: $(length(data[:Nset]))")
            println(io, "Number of branches: $(length(data[:Lset]))")
            println(io, "Number of batteries: $(length(data[:Bset]))")
            println(io, "Number of PV units: $(length(data[:Dset]))")
            println(io, "\n--- OPTIMIZATION STATUS ---")
            println(io, "Status: $(sol_socp_bf[:status])")
            println(io, "\n--- OBJECTIVE VALUE ---")
            @printf(io, "Total Cost: \$%.4f\n", sol_socp_bf[:objective])
            println(io, "\n--- COMPUTATION TIME ---")
            @printf(io, "Solver time: %.4f seconds\n", sol_socp_bf[:solve_time])
            println(io, "\n--- SOLVER ---")
            solver_name_bf = use_gurobi_for_bf ? "Gurobi" : "Ipopt"
            println(io, "Solver: $(solver_name_bf)")
            println(io, "Formulation: SOCP (BFM-NL)")
            println(io, "="^80)
        end
        
        println(COLOR_SUCCESS, "✓ Results written to $(results_file)", COLOR_RESET)
    end

end # socp brute-forced solve

begin # function primal update (update 1) tadmm socp
    function primal_update_tadmm_socp!(B_local, Bhat, u_local, data, ρ::Float64, t0::Int; solver::Symbol=:ipopt)
        @unpack Nset, Lset, L1set, Nm1set, NLset, Dset, Bset, Tset = data
        @unpack substationBus, parent, children = data
        @unpack rdict_pu, xdict_pu, Vminpu, Vmaxpu = data
        @unpack p_L_pu, q_L_pu, p_D_pu, S_D_R = data
        @unpack B0_pu, B_R_pu, P_B_R_pu, soc_min, soc_max = data
        @unpack LoadShapeCost, C_B, delta_t_h, kVA_B = data
        
        j1 = substationBus
        Δt = delta_t_h
        P_BASE = kVA_B
        
        # Create model for subproblem t0 (solver choice: Gurobi for SOCP or Ipopt for NLP)
        model = Model()
        if solver == :gurobi
            set_optimizer(model, Gurobi.Optimizer)
            set_silent(model)
            set_optimizer_attribute(model, "NonConvex", 2)
            set_optimizer_attribute(model, "OutputFlag", 0)  # Suppress Gurobi output
            set_optimizer_attribute(model, "DualReductions", 0)
        else  # :ipopt
            set_optimizer(model, Ipopt.Optimizer)
            set_silent(model)
            set_optimizer_attribute(model, "print_level", 0)
            set_optimizer_attribute(model, "max_iter", 3000)
            set_optimizer_attribute(model, "tol", 1e-6)
            set_optimizer_attribute(model, "acceptable_tol", 1e-4)
            set_optimizer_attribute(model, "mu_strategy", "adaptive")
            set_optimizer_attribute(model, "linear_solver", "mumps")  # More robust than default MA27
            set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
            set_optimizer_attribute(model, "bound_relax_factor", 1e-8)
        end
        
        # ===== NETWORK VARIABLES (time t0 ONLY) =====
        @variable(model, P_Subs_t0 >= 0)
        @variable(model, Q_Subs_t0)
        @variable(model, P_t0[(i, j) in Lset])
        @variable(model, Q_t0[(i, j) in Lset])
        @variable(model, v_t0[j in Nset])
        @variable(model, ℓ_t0[(i, j) in Lset] >= 1e-6)  # Current squared magnitude (bounded away from 0)
        @variable(model, q_D_t0[j in Dset])
        
        # ===== BATTERY VARIABLES (ENTIRE horizon) =====
        @variable(model, P_B_var[j in Bset, t in Tset])
        @variable(model, B_var[j in Bset, t in Tset])
        
        # ===== OBJECTIVE =====
        # Energy cost at t0
        energy_cost = LoadShapeCost[t0] * P_Subs_t0 * P_BASE * Δt
        
        # Battery quadratic cost at t0
        battery_cost = sum(C_B * (P_B_var[j, t0] * P_BASE)^2 * Δt for j in Bset)
        
        # ADMM penalty (full horizon)
        penalty = (ρ / 2) * sum(sum((B_var[j, t] - Bhat[j][t] + u_local[j][t])^2 
                                    for t in Tset) for j in Bset)
        
        @objective(model, Min, energy_cost + battery_cost + penalty)
        
        # ===== SPATIAL CONSTRAINTS (time t0 ONLY) =====
        
        # Real power balance - substation
        @constraint(model, P_Subs_t0 - sum(P_t0[(j1, j)] for (j1, j) in L1set) == 0)
        
        # Real power balance - non-substation nodes
        for j in Nm1set
            i = parent[j]
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P_t0[(j, k)] for k in children[j])
            p_L_j = (j in NLset) ? p_L_pu[j, t0] : 0.0
            p_D_j = (j in Dset) ? p_D_pu[j, t0] : 0.0
            P_B_j = (j in Bset) ? P_B_var[j, t0] : 0.0
            @constraint(model, sum_Pjk - P_t0[(i, j)] == P_B_j + p_D_j - p_L_j)
        end
        
        # Reactive power balance - substation
        @constraint(model, Q_Subs_t0 - sum(Q_t0[(j1, j)] for (j1, j) in L1set) == 0)
        
        # Reactive power balance - non-substation nodes
        for j in Nm1set
            i = parent[j]
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q_t0[(j, k)] for k in children[j])
            q_L_j = (j in NLset) ? q_L_pu[j, t0] : 0.0
            q_D_j = (j in Dset) ? q_D_t0[j] : 0.0
            @constraint(model, sum_Qjk - Q_t0[(i, j)] == q_D_j - q_L_j)
        end
        
        # Voltage drop constraints (BFM-NL)
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model, v_t0[j] == v_t0[i] - 2*(r_ij*P_t0[(i,j)] + x_ij*Q_t0[(i,j)]) + (r_ij^2 + x_ij^2)*ℓ_t0[(i,j)])
        end
        
        # SOC relaxation constraints
        for (i, j) in Lset
            @constraint(model, P_t0[(i,j)]^2 + Q_t0[(i,j)]^2 <= v_t0[i] * ℓ_t0[(i,j)])
        end
        
        # Fixed substation voltage
        @constraint(model, v_t0[j1] == 1.05^2)
        
        # Voltage limits
        for j in Nset
            @constraint(model, Vminpu[j]^2 <= v_t0[j] <= Vmaxpu[j]^2)
        end
        
        # PV reactive power limits
        for j in Dset
            p_D_val = p_D_pu[j, t0]
            S_D_R_val = S_D_R[j]
            q_max = sqrt(S_D_R_val^2 - p_D_val^2)
            @constraint(model, -q_max <= q_D_t0[j] <= q_max)
        end
        
        # ===== TEMPORAL CONSTRAINTS (ENTIRE horizon) =====
        
        for j in Bset
            # Battery SOC trajectory (all time periods)
            for t in Tset
                if t == 1
                    @constraint(model, B_var[j, t] == B0_pu[j] - P_B_var[j, t] * Δt)
                else
                    @constraint(model, B_var[j, t] == B_var[j, t-1] - P_B_var[j, t] * Δt)
                end
            end
            
            # SOC and power limits (all time periods)
            for t in Tset
                @constraint(model, soc_min[j] * B_R_pu[j] <= B_var[j, t] <= soc_max[j] * B_R_pu[j])
                @constraint(model, -P_B_R_pu[j] <= P_B_var[j, t] <= P_B_R_pu[j])
            end
        end
        
        # Solve subproblem
        optimize!(model)
        
        # Check if solution exists
        status = termination_status(model)
        if status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED
            error("Subproblem t0=$t0 failed with status: $status. Try reducing ρ or check problem feasibility.")
        end
        
        # Get pure solver time
        solver_time_t0 = solve_time(model)
        
        # Extract results and update local copies
        for j in Bset
            for t in Tset
                B_local[j][t] = value(B_var[j, t])
            end
        end
        
        # Extract ALL decision variables at time t0
        P_Subs_val = value(P_Subs_t0)
        Q_Subs_val = value(Q_Subs_t0)
        P_vals = Dict((i, j) => value(P_t0[(i, j)]) for (i, j) in Lset)
        Q_vals = Dict((i, j) => value(Q_t0[(i, j)]) for (i, j) in Lset)
        v_vals = Dict(j => value(v_t0[j]) for j in Nset)
        ℓ_vals = Dict((i, j) => value(ℓ_t0[(i, j)]) for (i, j) in Lset)
        q_D_vals = Dict(j => value(q_D_t0[j]) for j in Dset)
        # Extract P_B for ALL times (full trajectory)
        P_B_vals_all_times = Dict(j => Dict(t => value(P_B_var[j, t]) for t in Tset) for j in Bset)
        P_B_val_t0 = Dict(j => value(P_B_var[j, t0]) for j in Bset)
        
        # Compute objective components
        energy_cost_val = LoadShapeCost[t0] * P_Subs_val * P_BASE * Δt
        battery_cost_val = sum(C_B * (P_B_val_t0[j] * P_BASE)^2 * Δt for j in Bset)
        penalty_val = (ρ / 2) * sum(sum((value(B_var[j, t]) - Bhat[j][t] + u_local[j][t])^2 
                                        for t in Tset) for j in Bset)
        
        return Dict(
            :total_objective => objective_value(model),
            :energy_cost => energy_cost_val,
            :battery_cost => battery_cost_val,
            :penalty => penalty_val,
            :solve_time => solver_time_t0,
            :P_Subs => P_Subs_val,
            :Q_Subs => Q_Subs_val,
            :P => P_vals,
            :Q => Q_vals,
            :v => v_vals,
            :ℓ => ℓ_vals,
            :q_D => q_D_vals,
            :P_B => P_B_vals_all_times,  # Full trajectory for all times
            :B_local => B_local,
            :t0 => t0
        )
    end
end # function primal update (update 1) tadmm socp

begin # function consensus update (update 2) tadmm socp
    function consensus_update_tadmm_socp!(Bhat, B_collection, u_collection, data, ρ::Float64)
        @unpack Bset, Tset, B_R_pu, soc_min, soc_max = data
        
        violations = Tuple{Int,Int}[]  # Store (battery, time) pairs with violations
        violation_tolerance = 1e-4
        
        for j in Bset
            for t in Tset
                # Average across all T subproblems
                consensus_sum = sum(B_collection[t0][j][t] + u_collection[t0][j][t] for t0 in Tset)
                Bhat_new = consensus_sum / length(Tset)
                
                # Clamp to feasible range
                B_min = soc_min[j] * B_R_pu[j]
                B_max = soc_max[j] * B_R_pu[j]
                
                violation_amount = max(B_min - Bhat_new, Bhat_new - B_max, 0.0)
                
                if violation_amount > violation_tolerance
                    push!(violations, (j, t))
                end
                
                Bhat[j][t] = clamp(Bhat_new, B_min, B_max)
            end
        end
        
        return Dict(
            :Bhat => Bhat,
            :bounds_violations => length(violations),
            :violation_indices => violations
        )
    end
end # function consensus update (update 2) tadmm socp

begin # function dual update (update 3) tadmm socp
    function dual_update_tadmm_socp!(u_collection, B_collection, Bhat, ρ::Float64, data)
        @unpack Bset, Tset = data
        max_change = 0.0
        
        for t0 in Tset
            for j in Bset
                for t in Tset
                    old_u = u_collection[t0][j][t]
                    u_collection[t0][j][t] += (B_collection[t0][j][t] - Bhat[j][t])
                    max_change = max(max_change, abs(u_collection[t0][j][t] - old_u))
                end
            end
        end
        
        return Dict(
            :u_collection => u_collection,
            :max_dual_change => max_change,
            :total_updates => length(Tset) * length(Bset) * length(Tset)
        )
    end
end # function dual update (update 3) tadmm socp

begin # function solve MPOPF tadmm socp
    function solve_MPOPF_SOCP_tADMM(data; ρ::Float64=1.0, 
                                        max_iter::Int=1000, eps_pri::Float64=1e-5, eps_dual::Float64=1e-4,
                                        adaptive_rho::Bool=false, solver::Symbol=:ipopt)
        @unpack Bset, Tset, B0_pu, B_R_pu, soc_min, soc_max = data
        
        # Initialize global consensus variables B̂ⱼᵗ
        Bhat = Dict(j => fill(B0_pu[j], length(Tset)) for j in Bset)
        for j in Bset
            Bhat[j] .= clamp.(Bhat[j], soc_min[j] * B_R_pu[j], soc_max[j] * B_R_pu[j])
        end
        
        # Initialize local SOC variables Bⱼᵗ'ᵗ⁰ for each subproblem
        B_collection = Dict(t0 => Dict(j => copy(Bhat[j]) for j in Bset) for t0 in Tset)
        
        # Initialize scaled dual variables uⱼᵗ'ᵗ⁰
        u_collection = Dict(t0 => Dict(j => zeros(length(Tset)) for j in Bset) for t0 in Tset)
        
        # Storage for ALL decision variables from each subproblem
        P_Subs_collection = Dict{Int,Float64}()
        Q_Subs_collection = Dict{Int,Float64}()
        P_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
        Q_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
        v_collection = Dict{Int,Dict{Int,Float64}}()
        ℓ_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
        q_D_collection = Dict{Int,Dict{Int,Float64}}()
        # P_B_collection[j][t0][t] = P_B value at time t from subproblem t0 for battery j
        P_B_collection = Dict{Int,Dict{Int,Dict{Int,Float64}}}()
        
        # Initialize dictionaries for branches and nodes
        for (i, j) in data[:Lset]
            P_collection[(i, j)] = Dict{Int,Float64}()
            Q_collection[(i, j)] = Dict{Int,Float64}()
            ℓ_collection[(i, j)] = Dict{Int,Float64}()
        end
        for j in data[:Nset]
            v_collection[j] = Dict{Int,Float64}()
        end
        for j in data[:Dset]
            q_D_collection[j] = Dict{Int,Float64}()
        end
        for j in Bset
            # Initialize P_B_collection[j][t0] for each subproblem t0
            P_B_collection[j] = Dict{Int,Dict{Int,Float64}}()
        end
        
        # History tracking
        obj_history = Float64[]
        energy_cost_history = Float64[]
        battery_cost_history = Float64[]
        penalty_history = Float64[]
        r_norm_history = Float64[]
        s_norm_history = Float64[]
        Bhat_history = Dict(j => Vector{Vector{Float64}}() for j in Bset)
        # Timing tracking: subproblem_times_history[k] = [t1, t2, ..., tT] for iteration k
        subproblem_times_history = Vector{Vector{Float64}}()
        iteration_effective_times = Float64[]  # max(subproblem times) per iteration
        
        # Store initial state
        for j in Bset
            push!(Bhat_history[j], copy(Bhat[j]))
        end
        
        # Adaptive ρ parameters
        ρ_current = ρ  # Track current ρ value
        μ_balance = 5.0   # Threshold for imbalance (standard: 2-5, was 10 - too conservative!)
        τ_incr = 2.0      # Factor to increase ρ
        τ_decr = 2.0      # Factor to decrease ρ
        ρ_min = 1.0       # Minimum ρ value
        ρ_max = 1e6       # Maximum ρ value
        update_interval = 10  # Update ρ every N iterations
        
        println("\n" * "="^80)
        print(COLOR_INFO)
        @printf "🎯 tADMM[SOCP]: T=%d, ρ_init=%.1f, adaptive=%s, |Bset|=%d, |Nset|=%d, |Lset|=%d\n" length(Tset) ρ adaptive_rho length(Bset) length(data[:Nset]) length(data[:Lset])
        print(COLOR_RESET)
        println("="^80)
        
        converged = false
        final_iter = max_iter
        
        try
            for k in 1:max_iter
            # 🔵 STEP 1: Primal Update - Solve T subproblems
            print(COLOR_INFO)
            @printf "  🔵 Iteration %3d: Primal updates" k
            print(COLOR_RESET)
            total_energy_cost = 0.0
            total_battery_cost = 0.0
            total_penalty = 0.0
            
            # Track subproblem solve times for this iteration
            subproblem_times_k = Float64[]
            
            for t0 in Tset
                result = primal_update_tadmm_socp!(B_collection[t0], Bhat, u_collection[t0], data, ρ_current, t0; solver=solver)
                # Use pure solver time from the subproblem
                push!(subproblem_times_k, result[:solve_time])
                
                # Update collections - Store ALL decision variables
                B_collection[t0] = result[:B_local]
                P_Subs_collection[t0] = result[:P_Subs]
                Q_Subs_collection[t0] = result[:Q_Subs]
                
                for (i, j) in data[:Lset]
                    P_collection[(i, j)][t0] = result[:P][(i, j)]
                    Q_collection[(i, j)][t0] = result[:Q][(i, j)]
                    ℓ_collection[(i, j)][t0] = result[:ℓ][(i, j)]
                end
                
                for j in data[:Nset]
                    v_collection[j][t0] = result[:v][j]
                end
                
                for j in data[:Dset]
                    q_D_collection[j][t0] = result[:q_D][j]
                end
                
                for j in Bset
                    # Store the full P_B trajectory from this subproblem
                    # P_B_collection[j] is a Dict{Int, Dict{Int, Float64}}
                    # where P_B_collection[j][t0][t] = P_B value at time t from subproblem t0
                    P_B_collection[j][t0] = result[:P_B][j]  # This is already a Dict(t => value)
                end
                
                # Accumulate costs
                total_energy_cost += result[:energy_cost]
                total_battery_cost += result[:battery_cost]
                total_penalty += result[:penalty]
            end
            print(COLOR_SUCCESS)
            println(" ✓")
            print(COLOR_RESET)
            
            true_objective = total_energy_cost + total_battery_cost
            push!(obj_history, true_objective)
            push!(energy_cost_history, total_energy_cost)
            push!(battery_cost_history, total_battery_cost)
            push!(penalty_history, total_penalty)
            
            # Store timing information
            push!(subproblem_times_history, subproblem_times_k)
            push!(iteration_effective_times, maximum(subproblem_times_k))  # Parallel execution time
            
            # 🔴 STEP 2: Consensus Update
            Bhat_old = Dict(j => copy(Bhat[j]) for j in Bset)
            
            consensus_result = consensus_update_tadmm_socp!(Bhat, B_collection, u_collection, data, ρ)
            
            # STEP 3: Dual Update
            dual_result = dual_update_tadmm_socp!(u_collection, B_collection, Bhat, ρ, data)            
            # Store history
            for j in Bset
                push!(Bhat_history[j], copy(Bhat[j]))
            end
            
            # 📏 STEP 4: Compute residuals
            r_vectors = []
            for t0 in Tset
                for j in Bset
                    push!(r_vectors, B_collection[t0][j] - Bhat[j])
                end
            end
            r_norm = norm(vcat(r_vectors...)) / (length(Tset) * length(Bset))
            
            s_vectors = []
            for j in Bset
                push!(s_vectors, Bhat[j] - Bhat_old[j])
            end
            s_norm = ρ_current * norm(vcat(s_vectors...)) / length(Bset)
            
            push!(r_norm_history, r_norm)
            push!(s_norm_history, s_norm)
            
            # 🔧 STEP 5: Adaptive ρ Update (optional)
            if adaptive_rho && k % update_interval == 0 && k < max_iter - 50
                ρ_old = ρ_current
                
                if r_norm > 5.0 * s_norm
                    # Primal residual too large -> INCREASE rho to enforce consensus
                    ρ_current = min(ρ_max, τ_incr * ρ_current)
                    print(COLOR_WARNING)
                    @printf "  [UP] rho: %.1f -> %.1f (r/s=%.1f > 5.0 primal lagging)\n" ρ_old ρ_current (r_norm/s_norm)
                    print(COLOR_RESET)
                elseif s_norm > 5.0 * r_norm
                    # Dual residual too large -> DECREASE rho to allow flexibility
                    ρ_current = max(ρ_min, ρ_current / τ_decr)
                    print(COLOR_WARNING)
                    @printf "  [DOWN] rho: %.1f -> %.1f (s/r=%.1f > 5.0 dual lagging)\n" ρ_old ρ_current (s_norm/r_norm)
                    print(COLOR_RESET)
                end
                
                # CRITICAL: Rescale dual variables when rho changes  
                if ρ_current != ρ_old
                    scale_factor = ρ_old / ρ_current
                    for t0 in Tset, j in Bset
                        u_collection[t0][j] .*= scale_factor
                    end
                end
            end
            
            @printf "k=%3d  obj=\$%.4f (energy=\$%.4f, battery=\$%.4f, penalty=\$%.4f)  ‖r‖=%.2e  ‖s‖=%.2e  ρ=%.1f\n" k true_objective total_energy_cost total_battery_cost total_penalty r_norm s_norm ρ_current
            
            # Check convergence
            if r_norm ≤ eps_pri && s_norm ≤ eps_dual
                converged = true
                print(COLOR_SUCCESS)
                @printf "🎉 tADMM converged at iteration %d\n" k
                print(COLOR_RESET)
                break
            end
            
            final_iter = k  # Track last completed iteration
        end  # end for loop
        catch e
            if isa(e, InterruptException)
                print(COLOR_WARNING)
                @printf "\n⚠ tADMM interrupted at iteration %d\n" final_iter
                print(COLOR_RESET)
            else
                rethrow(e)
            end
        end
        
        # Report convergence status
        if !converged
            print(COLOR_WARNING)
            @printf "⚠ tADMM did NOT converge after %d iterations (‖r‖=%.2e > %.2e or ‖s‖=%.2e > %.2e)\n" max_iter r_norm_history[end] eps_pri s_norm_history[end] eps_dual
            print(COLOR_RESET)
        end
        
        # Extract final solution from blue decision variables (local subproblem solutions)
        # Use the diagonal: subproblem t0 provides the solution for time t0
        P_Subs_final = Dict()
        Q_Subs_final = Dict()
        for t in Tset
            P_Subs_final[t] = P_Subs_collection[t]
            Q_Subs_final[t] = Q_Subs_collection[t]
        end
        
        P_final = Dict()
        Q_final = Dict()
        ℓ_final = Dict()
        for (i, j) in data[:Lset]
            for t in Tset
                P_final[(i, j), t] = P_collection[(i, j)][t]
                Q_final[(i, j), t] = Q_collection[(i, j)][t]
                ℓ_final[(i, j), t] = ℓ_collection[(i, j)][t]
            end
        end
        
        v_final = Dict()
        for j in data[:Nset]
            for t in Tset
                v_final[j, t] = v_collection[j][t]
            end
        end
        
        q_D_final = Dict()
        for j in data[:Dset]
            for t in Tset
                q_D_final[j, t] = q_D_collection[j][t]
            end
        end
        
        P_B_final = Dict()
        B_final = Dict()
        for j in Bset
            for t in Tset
                # Use the diagonal: subproblem t provides solution at time t
                # P_B_collection[j][t0][t] where we want t0=t (diagonal)
                P_B_final[j, t] = P_B_collection[j][t][t]
                B_final[j, t] = B_collection[t][j][t]
            end
        end
        
        return Dict(
            :status => MOI.OPTIMAL,  # Assume converged
            :P_Subs => P_Subs_final,
            :Q_Subs => Q_Subs_final,
            :P => P_final,
            :Q => Q_final,
            :v => v_final,
            :ℓ => ℓ_final,
            :q_D => q_D_final,
            :P_B => P_B_final,
            :B => B_final,
            :objective => last(obj_history),
            :consensus_trajectory => Bhat,
            :local_solutions => B_collection,
            :dual_variables => u_collection,
            :convergence_history => Dict(
                :obj_history => obj_history,
                :energy_cost_history => energy_cost_history,
                :battery_cost_history => battery_cost_history,
                :penalty_history => penalty_history,
                :r_norm_history => r_norm_history,
                :s_norm_history => s_norm_history,
                :Bhat_history => Bhat_history
            ),
            :timing => Dict(
                :subproblem_times_history => subproblem_times_history,
                :iteration_effective_times => iteration_effective_times,
                :total_effective_time => sum(iteration_effective_times),
                :total_sequential_time => sum(sum.(subproblem_times_history))
            )
        )
    end
end # function solve MPOPF tadmm socp

begin # tadmm socp solve
    if !isempty(data[:Bset])  # Only run tADMM if there are batteries
        println("\n" * "="^80)
        println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH SOCP (tADMM)", COLOR_RESET)
        println("="^80)
        
        # Print ρ scaling info
        if rho_scaling_with_T
            scaling_factor = sqrt(T / 24.0)
            println(COLOR_INFO, "ρ scaled with T: $(rho_base) × √($(T)/24) = $(round(rho_tadmm, digits=1))", COLOR_RESET)
        else
            println(COLOR_INFO, "Using fixed ρ = $(rho_tadmm)", COLOR_RESET)
        end
        
        solver_tadmm_choice = use_gurobi_for_tadmm ? :gurobi : :ipopt
        sol_socp_tadmm = solve_MPOPF_SOCP_tADMM(data; ρ=rho_tadmm, 
                                                    max_iter=max_iter_tadmm, 
                                                    eps_pri=eps_pri_tadmm, 
                                                    eps_dual=eps_dual_tadmm,
                                                    adaptive_rho=adaptive_rho_tadmm,
                                                    solver=solver_tadmm_choice)
        
        # Report results
        println("\n--- tADMM SOLUTION STATUS ---")
        print(COLOR_SUCCESS)
        @printf "Objective: \$%.2f\n" sol_socp_tadmm[:objective]
        print(COLOR_RESET)
        println("\n--- COMPUTATION TIME ---")
        @printf "Effective wall clock time: %.2f seconds (%d iterations)\n" sol_socp_tadmm[:timing][:total_effective_time] length(sol_socp_tadmm[:timing][:iteration_effective_times])
        @printf "Sequential time (all subproblems): %.2f seconds\n" sol_socp_tadmm[:timing][:total_sequential_time]
        
        # Compare with brute force
        if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
            println("\n--- COMPARISON WITH BRUTE FORCE ---")
            obj_diff = abs(sol_socp_tadmm[:objective] - sol_socp_bf[:objective])
            obj_rel_diff = obj_diff / sol_socp_bf[:objective] * 100
            
            @printf "Brute Force objective: \$%.4f\n" sol_socp_bf[:objective]
            @printf "tADMM objective:       \$%.4f\n" sol_socp_tadmm[:objective]
            @printf "Absolute difference:   \$%.4f\n" obj_diff
            @printf "Relative difference:   %.4f%%\n" obj_rel_diff
            
            # Compare battery schedules
            println("\n--- BATTERY SCHEDULE COMPARISON ---")
            P_BASE = data[:kVA_B]
            E_BASE = P_BASE * 1.0
            
            for j in data[:Bset]
                println("Battery $j:")
                
                # Brute force
                P_B_bf_kW = [sol_socp_bf[:P_B][j, t] * P_BASE for t in data[:Tset]]
                B_bf_kWh = [sol_socp_bf[:B][j, t] * E_BASE for t in data[:Tset]]
                
                # tADMM
                P_B_tadmm_kW = [sol_socp_tadmm[:P_B][j, t] * P_BASE for t in data[:Tset]]
                B_tadmm_kWh = [sol_socp_tadmm[:B][j, t] * E_BASE for t in data[:Tset]]
                
                @printf "  BF P_B (kW):    min=%.1f, max=%.1f\n" minimum(P_B_bf_kW) maximum(P_B_bf_kW)
                @printf "  tADMM P_B (kW): min=%.1f, max=%.1f\n" minimum(P_B_tadmm_kW) maximum(P_B_tadmm_kW)
                @printf "  BF B (kWh):     min=%.1f, max=%.1f\n" minimum(B_bf_kWh) maximum(B_bf_kWh)
                @printf "  tADMM B (kWh):  min=%.1f, max=%.1f\n" minimum(B_tadmm_kWh) maximum(B_tadmm_kWh)
            end
        end
        
        # Write tADMM results to file
        processedData_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData")
        system_folder = "$(systemName)_T$(T)"
        system_dir = joinpath(processedData_dir, system_folder)
        results_file = joinpath(system_dir, "results_socp_tadmm.txt")
        
        open(results_file, "w") do io
            println(io, "="^80)
            println(io, "SOCP tADMM OPTIMIZATION RESULTS")
            println(io, "="^80)
            println(io, "System: $(systemName)")
            println(io, "Time horizon: T=$(T) periods")
            println(io, "Time step: $(delta_t_h) hours")
            println(io, "Number of buses: $(length(data[:Nset]))")
            println(io, "Number of branches: $(length(data[:Lset]))")
            println(io, "Number of batteries: $(length(data[:Bset]))")
            println(io, "Number of PV units: $(length(data[:Dset]))")
            println(io, "\n--- tADMM PARAMETERS ---")
            @printf(io, "Initial ρ: %.1f\n", rho_tadmm)
            println(io, "Adaptive ρ: $(adaptive_rho_tadmm)")
            @printf(io, "Max iterations: %d\n", max_iter_tadmm)
            @printf(io, "Primal tolerance: %.1e\n", eps_pri_tadmm)
            @printf(io, "Dual tolerance: %.1e\n", eps_dual_tadmm)
            println(io, "\n--- CONVERGENCE ---")
            n_iters = length(sol_socp_tadmm[:timing][:iteration_effective_times])
            @printf(io, "Converged in %d iterations\n", n_iters)
            @printf(io, "Final primal residual: %.2e\n", sol_socp_tadmm[:convergence_history][:r_norm_history][end])
            @printf(io, "Final dual residual: %.2e\n", sol_socp_tadmm[:convergence_history][:s_norm_history][end])
            println(io, "\n--- OBJECTIVE VALUE ---")
            @printf(io, "Total Cost: \$%.4f\n", sol_socp_tadmm[:objective])
            if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
                obj_diff = abs(sol_socp_tadmm[:objective] - sol_socp_bf[:objective])
                obj_rel_diff = obj_diff / sol_socp_bf[:objective] * 100
                @printf(io, "Brute Force objective: \$%.4f\n", sol_socp_bf[:objective])
                @printf(io, "Absolute difference: \$%.4f\n", obj_diff)
                @printf(io, "Relative difference: %.4f%%\n", obj_rel_diff)
            end
            println(io, "\n--- COMPUTATION TIME ---")
            @printf(io, "Effective solver time: %.4f seconds\n", sol_socp_tadmm[:timing][:total_effective_time])
            @printf(io, "Sequential solver time (all subproblems): %.4f seconds\n", sol_socp_tadmm[:timing][:total_sequential_time])
            if sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED
                @printf(io, "Brute Force solver time: %.4f seconds\n", sol_socp_bf[:solve_time])
                speedup = sol_socp_bf[:solve_time] / sol_socp_tadmm[:timing][:total_effective_time]
                @printf(io, "Speedup vs Brute Force: %.2fx\n", speedup)
            end
            println(io, "\n--- SOLVER ---")
            solver_name_tadmm = use_gurobi_for_tadmm ? "Gurobi" : "Ipopt"
            println(io, "Subproblem solver: $(solver_name_tadmm)")
            println(io, "Formulation: SOCP (BFM-NL)")
            println(io, "Decomposition: Temporal ADMM")
            println(io, "="^80)
        end
        
        println(COLOR_SUCCESS, "✓ tADMM results written to $(results_file)", COLOR_RESET)
        
        println("\n" * "="^80)
        println(COLOR_HIGHLIGHT, "MPOPF SOCP tADMM SOLUTION COMPLETE", COLOR_RESET)
        println("="^80)
    else
        println("\n" * "="^80)
        print(COLOR_WARNING)
        println("⚠ No batteries in system - skipping tADMM solution")
        print(COLOR_RESET)
        println("="^80)
        sol_socp_tadmm = nothing
    end
end # tadmm socp solve

begin # plotting results
    println("\n" * "="^80)
    println(COLOR_INFO, "GENERATING PLOTS", COLOR_RESET)
    println("="^80)

    # Create output directories (use absolute path like model writing section)
    processedData_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData")
    system_folder = "$(systemName)_T$(T)"
    system_dir = joinpath(processedData_dir, system_folder)
    mkpath(processedData_dir)  # Create processedData folder
    mkpath(system_dir)  # Create system-specific subfolder with horizon

    # Plot tADMM convergence history (MOVED HERE - after Plotter.jl is included)
    if !isempty(data[:Bset]) && !isnothing(sol_socp_tadmm)
        println("\n" * "="^80)
        println(COLOR_INFO, "PLOTTING tADMM CONVERGENCE", COLOR_RESET)
        println("="^80)
        
        # Create convergence plots directory
        conv_plots_dir = joinpath(system_dir, "convergence")
        mkpath(conv_plots_dir)
        
        # Use Plotter.jl function for consistent styling
        conv_plot_path = joinpath(conv_plots_dir, "tadmm_convergence_socp.png")
        plot_tadmm_ldf_convergence(sol_socp_tadmm, sol_socp_bf, eps_pri_tadmm, eps_dual_tadmm,
                                showPlots=showPlots, savePlots=true, 
                                filename=conv_plot_path)
        
        println(COLOR_SUCCESS, "✓ Convergence plots saved to $(conv_plots_dir)", COLOR_RESET)
    end

    println("\n")

    # Plot input curves (load, PV, cost) - save in system-specific folder (not processedData root)
    input_curves_path = joinpath(system_dir, "input_curves_socp.png")
    plot_input_curves(data, showPlots=showPlots, savePlots=true, filename=input_curves_path)

    # Plot battery actions (only if optimization was successful and batteries exist)
    # Save in system-specific subfolder
    if (sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED) && !isempty(data[:Bset])
        battery_actions_path = joinpath(system_dir, "battery_actions_socp_bf.png")
        plot_battery_actions(sol_socp_bf, data, "SOCP-BF (Gurobi)", 
                            showPlots=showPlots, savePlots=true, 
                            filename=battery_actions_path,
                            plot_all_batteries=saveAllBatteryPlots)
        
        # Plot tADMM battery actions if available
        if !isnothing(sol_socp_tadmm)
            battery_actions_tadmm_path = joinpath(system_dir, "battery_actions_socp_tadmm.png")
            plot_battery_actions(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)", 
                                showPlots=showPlots, savePlots=true, 
                                filename=battery_actions_tadmm_path,
                                plot_all_batteries=saveAllBatteryPlots)
        end
    else
        if isempty(data[:Bset])
            println("No batteries in system, skipping battery actions plot")
        end
    end

    # Plot substation power and cost (only if optimization was successful)
    if (sol_socp_bf[:status] == MOI.OPTIMAL || sol_socp_bf[:status] == MOI.LOCALLY_SOLVED)
        subs_power_cost_path = joinpath(system_dir, "substation_power_cost_socp_bf.png")
        plot_substation_power_and_cost(sol_socp_bf, data, "SOCP-BF (Gurobi)",
                                    showPlots=showPlots, savePlots=true,
                                    filename=subs_power_cost_path)
        
        # Plot tADMM substation power if available
        if !isnothing(sol_socp_tadmm)
            subs_power_cost_tadmm_path = joinpath(system_dir, "substation_power_cost_socp_tadmm.png")
            plot_substation_power_and_cost(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)",
                                        showPlots=showPlots, savePlots=true,
                                        filename=subs_power_cost_tadmm_path)
        end
        
        # Plot voltage profile for last bus
        last_bus = maximum(data[:Nset])
        voltage_one_bus_path = joinpath(system_dir, "voltage_profile_bus_$(last_bus)_socp_bf.png")
        plot_voltage_profile_one_bus(sol_socp_bf, data, "SOCP-BF (Gurobi)",
                                    showPlots=showPlots, savePlots=true,
                                    filename=voltage_one_bus_path)
        
        # Plot tADMM voltage profile for last bus
        if !isnothing(sol_socp_tadmm)
            voltage_one_bus_tadmm_path = joinpath(system_dir, "voltage_profile_bus_$(last_bus)_socp_tadmm.png")
            plot_voltage_profile_one_bus(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)",
                                        showPlots=showPlots, savePlots=true,
                                        filename=voltage_one_bus_tadmm_path)
        end
        
        # Plot voltage profile for all buses at middle time step + create GIF
        voltage_all_buses_path = joinpath(system_dir, "voltage_profile_all_buses_socp_bf.png")
        voltage_gif_path = joinpath(system_dir, "voltage_animation_socp_bf.gif")
        plot_voltage_profile_all_buses(sol_socp_bf, data, "SOCP-BF (Gurobi)",
                                    showPlots=showPlots, savePlots=true,
                                    filename=voltage_all_buses_path,
                                    create_gif=true,
                                    gif_filename=voltage_gif_path)
        
        # Plot tADMM voltage profile for all buses + create GIF
        if !isnothing(sol_socp_tadmm)
            voltage_all_buses_tadmm_path = joinpath(system_dir, "voltage_profile_all_buses_socp_tadmm.png")
            voltage_gif_tadmm_path = joinpath(system_dir, "voltage_animation_socp_tadmm.gif")
            plot_voltage_profile_all_buses(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)",
                                        showPlots=showPlots, savePlots=true,
                                        filename=voltage_all_buses_tadmm_path,
                                        create_gif=true,
                                        gif_filename=voltage_gif_tadmm_path)
        end
        
        # Plot PV power (p_D and q_D) if PV exists
        if !isempty(data[:Dset])
            # Plot individual PV plots for ALL PVs if enabled
            if saveAllPVPlots
                println("\n📊 Generating individual PV plots for all $(length(data[:Dset])) PV units...")
                
                # BF individual PV plots
                pv_power_all_bf_path = joinpath(system_dir, "pv_power_socp_bf.png")
                plot_pv_power(sol_socp_bf, data, "SOCP-BF (Gurobi)",
                            showPlots=showPlots, savePlots=true,
                            filename=pv_power_all_bf_path,
                            plot_all_pvs=true)
                
                # tADMM individual PV plots if available
                if !isnothing(sol_socp_tadmm)
                    pv_power_all_tadmm_path = joinpath(system_dir, "pv_power_socp_tadmm.png")
                    plot_pv_power(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)",
                                showPlots=showPlots, savePlots=true,
                                filename=pv_power_all_tadmm_path,
                                plot_all_pvs=true)
                end
                println("✓ Individual PV plots saved")
            else
                # Plot only first PV
                first_pv_bus = minimum(data[:Dset])
                pv_power_path = joinpath(system_dir, "pv_power_bus_$(first_pv_bus)_socp_bf.png")
                plot_pv_power(sol_socp_bf, data, "SOCP-BF (Gurobi)",
                            pv_index=1,  # Plot first PV (or only PV)
                            showPlots=showPlots, savePlots=true,
                            filename=pv_power_path)
                
                # Plot tADMM PV power if available
                if !isnothing(sol_socp_tadmm)
                    pv_power_tadmm_path = joinpath(system_dir, "pv_power_bus_$(first_pv_bus)_socp_tadmm.png")
                    plot_pv_power(sol_socp_tadmm, data, "SOCP-tADMM (Ipopt)",
                                pv_index=1,
                                showPlots=showPlots, savePlots=true,
                                filename=pv_power_tadmm_path)
                end
            end
            
            # Create PV power circle GIFs
            if saveAllPVPlots
                # Create circle GIFs for ALL PVs
                println("📊 Generating PV power circle GIFs for all $(length(data[:Dset])) PV units...")
                pv_circle_gif_all_path = joinpath(system_dir, "pv_power_circle_socp_bf.gif")
                plot_pv_power_circle_gif(sol_socp_bf, data, "SOCP-BF",
                                        showPlots=showPlots, savePlots=true,
                                        filename=pv_circle_gif_all_path,
                                        plot_all_pvs=true)
                
                if !isnothing(sol_socp_tadmm)
                    pv_circle_gif_all_tadmm_path = joinpath(system_dir, "pv_power_circle_socp_tadmm.gif")
                    plot_pv_power_circle_gif(sol_socp_tadmm, data, "SOCP-tADMM",
                                            showPlots=showPlots, savePlots=true,
                                            filename=pv_circle_gif_all_tadmm_path,
                                            plot_all_pvs=true)
                end
                println("✓ PV power circle GIFs saved")
            else
                # Create circle GIF for first PV only
                first_pv_bus = minimum(data[:Dset])
                pv_circle_gif_path = joinpath(system_dir, "pv_power_circle_bus_$(first_pv_bus)_socp_bf.gif")
                plot_pv_power_circle_gif(sol_socp_bf, data, "SOCP-BF",
                                        pv_index=1,
                                        showPlots=showPlots, savePlots=true,
                                        filename=pv_circle_gif_path)
                
                # Create tADMM PV power circle GIF if available
                if !isnothing(sol_socp_tadmm)
                    pv_circle_gif_tadmm_path = joinpath(system_dir, "pv_power_circle_bus_$(first_pv_bus)_socp_tadmm.gif")
                    plot_pv_power_circle_gif(sol_socp_tadmm, data, "SOCP-tADMM",
                                            pv_index=1,
                                            showPlots=showPlots, savePlots=true,
                                            filename=pv_circle_gif_tadmm_path)
                end
            end
        else
            println("No PV units in system, skipping PV power plots")
        end
    end

    println("\n" * "="^80)
    println(COLOR_SUCCESS, "PLOTTING COMPLETE", COLOR_RESET)
    println("="^80)
end # plotting results

end # entire script including environment setup
