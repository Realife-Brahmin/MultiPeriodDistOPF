"""
Solution Validator for tADMM
Validates branch flow model equations and constraint feasibility
"""

using Printf
using LinearAlgebra

"""
    validate_branch_flow_equations(sol::Dict, data::Dict; tol=1e-4, verbose=false)

Validate all branch flow model equations for a given solution.
Returns a dictionary with validation results.

# Arguments
- `sol::Dict`: Solution dictionary containing P, Q, v, ℓ, P_B, q_D, etc.
- `data::Dict`: Problem data dictionary
- `tol::Float64`: Tolerance for constraint violations (default: 1e-4)
- `verbose::Bool`: Print detailed violations (default: false)

# Returns
- `results::Dict`: Validation results with counts of violations and max violations
"""
function validate_branch_flow_equations(sol::Dict, data::Dict; tol=1e-4, verbose=false)
    @unpack Nset, Lset, L1set, Nm1set, NLset, Dset, Bset, Tset = data
    @unpack substationBus, parent, children = data
    @unpack rdict_pu, xdict_pu, Vminpu, Vmaxpu = data
    @unpack p_L_pu, q_L_pu, p_D_pu, S_D_R = data
    @unpack B_R_pu, P_B_R_pu, soc_min, soc_max, delta_t_h = data
    
    j1 = substationBus
    T = length(Tset)
    
    # Extract solution
    P_Subs = sol[:P_Subs]
    Q_Subs = sol[:Q_Subs]
    P = sol[:P]
    Q = sol[:Q]
    v = sol[:v]
    ℓ = sol[:ℓ]
    P_B = sol[:P_B]
    B = sol[:B]
    q_D = sol[:q_D]
    
    # Violation counters
    violations = Dict(
        :real_power_balance => 0,
        :reactive_power_balance => 0,
        :voltage_drop => 0,
        :soc_relaxation => 0,
        :voltage_bounds => 0,
        :battery_power_bounds => 0,
        :battery_soc_bounds => 0,
        :pv_reactive_bounds => 0,
        :substation_voltage => 0
    )
    
    max_violations = Dict(
        :real_power_balance => 0.0,
        :reactive_power_balance => 0.0,
        :voltage_drop => 0.0,
        :soc_relaxation => 0.0,
        :voltage_bounds => 0.0,
        :battery_power_bounds => 0.0,
        :battery_soc_bounds => 0.0,
        :pv_reactive_bounds => 0.0,
        :substation_voltage => 0.0
    )
    
    # Detailed violations (only stored if verbose or violation occurs)
    detailed_violations = []
    
    for t in Tset
        # 1. Real power balance - substation
        lhs = P_Subs[t] - sum(P[(j1, j), t] for (j1, j) in L1set)
        viol = abs(lhs)
        if viol > tol
            violations[:real_power_balance] += 1
            max_violations[:real_power_balance] = max(max_violations[:real_power_balance], viol)
            if verbose
                push!(detailed_violations, 
                    "Real power balance (substation) at t=$t: violation = $viol")
            end
        end
        
        # 2. Real power balance - non-substation nodes
        for j in Nm1set
            i = parent[j]
            sum_Pjk = isempty(children[j]) ? 0.0 : sum(P[(j, k), t] for k in children[j])
            p_L_j = (j in NLset) ? p_L_pu[j, t] : 0.0
            p_D_j = (j in Dset) ? p_D_pu[j, t] : 0.0
            P_B_j = (j in Bset) ? P_B[j, t] : 0.0
            
            lhs = sum_Pjk - P[(i, j), t]
            rhs = P_B_j + p_D_j - p_L_j
            viol = abs(lhs - rhs)
            
            if viol > tol
                violations[:real_power_balance] += 1
                max_violations[:real_power_balance] = max(max_violations[:real_power_balance], viol)
                if verbose
                    push!(detailed_violations,
                        "Real power balance (node $j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 3. Reactive power balance - substation
        lhs = Q_Subs[t] - sum(Q[(j1, j), t] for (j1, j) in L1set)
        viol = abs(lhs)
        if viol > tol
            violations[:reactive_power_balance] += 1
            max_violations[:reactive_power_balance] = max(max_violations[:reactive_power_balance], viol)
            if verbose
                push!(detailed_violations,
                    "Reactive power balance (substation) at t=$t: violation = $viol")
            end
        end
        
        # 4. Reactive power balance - non-substation nodes
        for j in Nm1set
            i = parent[j]
            sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])
            q_L_j = (j in NLset) ? q_L_pu[j, t] : 0.0
            q_D_j = (j in Dset) ? q_D[j, t] : 0.0
            
            lhs = sum_Qjk - Q[(i, j), t]
            rhs = q_D_j - q_L_j
            viol = abs(lhs - rhs)
            
            if viol > tol
                violations[:reactive_power_balance] += 1
                max_violations[:reactive_power_balance] = max(max_violations[:reactive_power_balance], viol)
                if verbose
                    push!(detailed_violations,
                        "Reactive power balance (node $j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 5. Voltage drop constraints (BFM-NL)
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            
            lhs = v[j, t]
            rhs = v[i, t] - 2*(r_ij*P[(i, j), t] + x_ij*Q[(i, j), t]) + (r_ij^2 + x_ij^2)*ℓ[(i, j), t]
            viol = abs(lhs - rhs)
            
            if viol > tol
                violations[:voltage_drop] += 1
                max_violations[:voltage_drop] = max(max_violations[:voltage_drop], viol)
                if verbose
                    push!(detailed_violations,
                        "Voltage drop (branch $i→$j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 6. SOC relaxation constraints
        for (i, j) in Lset
            lhs = P[(i, j), t]^2 + Q[(i, j), t]^2
            rhs = v[i, t] * ℓ[(i, j), t]
            viol = max(0.0, lhs - rhs)  # Only penalize if lhs > rhs
            
            if viol > tol
                violations[:soc_relaxation] += 1
                max_violations[:soc_relaxation] = max(max_violations[:soc_relaxation], viol)
                if verbose
                    push!(detailed_violations,
                        "SOC relaxation (branch $i→$j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 7. Voltage bounds
        for j in Nset
            v_min = Vminpu[j]^2
            v_max = Vmaxpu[j]^2
            viol_min = max(0.0, v_min - v[j, t])
            viol_max = max(0.0, v[j, t] - v_max)
            viol = max(viol_min, viol_max)
            
            if viol > tol
                violations[:voltage_bounds] += 1
                max_violations[:voltage_bounds] = max(max_violations[:voltage_bounds], viol)
                if verbose
                    push!(detailed_violations,
                        "Voltage bounds (node $j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 8. Substation voltage (fixed)
        viol = abs(v[j1, t] - 1.05^2)
        if viol > tol
            violations[:substation_voltage] += 1
            max_violations[:substation_voltage] = max(max_violations[:substation_voltage], viol)
            if verbose
                push!(detailed_violations,
                    "Substation voltage at t=$t: violation = $viol")
            end
        end
        
        # 9. Battery power bounds
        for j in Bset
            viol_min = max(0.0, -P_B_R_pu[j] - P_B[j, t])
            viol_max = max(0.0, P_B[j, t] - P_B_R_pu[j])
            viol = max(viol_min, viol_max)
            
            if viol > tol
                violations[:battery_power_bounds] += 1
                max_violations[:battery_power_bounds] = max(max_violations[:battery_power_bounds], viol)
                if verbose
                    push!(detailed_violations,
                        "Battery power bounds (battery $j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 10. Battery SOC bounds
        for j in Bset
            soc_min_val = soc_min[j] * B_R_pu[j]
            soc_max_val = soc_max[j] * B_R_pu[j]
            viol_min = max(0.0, soc_min_val - B[j, t])
            viol_max = max(0.0, B[j, t] - soc_max_val)
            viol = max(viol_min, viol_max)
            
            if viol > tol
                violations[:battery_soc_bounds] += 1
                max_violations[:battery_soc_bounds] = max(max_violations[:battery_soc_bounds], viol)
                if verbose
                    push!(detailed_violations,
                        "Battery SOC bounds (battery $j) at t=$t: violation = $viol")
                end
            end
        end
        
        # 11. PV reactive power bounds
        for j in Dset
            p_D_val = p_D_pu[j, t]
            S_D_R_val = S_D_R[j]
            q_max = sqrt(max(0.0, S_D_R_val^2 - p_D_val^2))
            viol_min = max(0.0, -q_max - q_D[j, t])
            viol_max = max(0.0, q_D[j, t] - q_max)
            viol = max(viol_min, viol_max)
            
            if viol > tol
                violations[:pv_reactive_bounds] += 1
                max_violations[:pv_reactive_bounds] = max(max_violations[:pv_reactive_bounds], viol)
                if verbose
                    push!(detailed_violations,
                        "PV reactive power bounds (PV $j) at t=$t: violation = $viol")
                end
            end
        end
    end
    
    return Dict(
        :violations => violations,
        :max_violations => max_violations,
        :detailed_violations => detailed_violations,
        :total_violations => sum(values(violations)),
        :feasible => sum(values(violations)) == 0
    )
end


"""
    print_validation_summary(validation_results::Dict; solution_name="Solution")

Print a summary of validation results.
"""
function print_validation_summary(validation_results::Dict; solution_name="Solution")
    violations = validation_results[:violations]
    max_viol = validation_results[:max_violations]
    total_viol = validation_results[:total_violations]
    feasible = validation_results[:feasible]
    
    println("\n", "="^80)
    println("SOLUTION VALIDATION: $solution_name")
    println("="^80)
    
    if feasible
        println("✓ All constraints satisfied!")
    else
        println("⚠ INFEASIBLE SOLUTION - Total constraint violations: $total_viol")
        println("\nViolation breakdown:")
        println("  Constraint Type                    | Count  | Max Violation")
        println("  " * "-"^70)
        
        # Only print constraint types with violations
        for (key, count) in violations
            if count > 0
                max_val = max_viol[key]
                println(@sprintf("  %-35s | %6d | %.2e", string(key), count, max_val))
            end
        end
    end
    
    println("="^80)
end


"""
    get_model_size_statistics(model::JuMP.Model)

Extract problem size statistics from a JuMP model.

# Returns
- `stats::Dict`: Dictionary containing number of variables, linear constraints, and nonlinear constraints
"""
function get_model_size_statistics(model::JuMP.Model)
    n_vars = num_variables(model)
    
    # Count constraint types
    n_linear = 0
    n_quadratic = 0
    n_nonlinear = 0
    
    # Get all constraint types
    constraint_types = list_of_constraint_types(model)
    
    for (F, S) in constraint_types
        n_cons = num_constraints(model, F, S)
        
        # Classify constraints
        if F <: Union{VariableRef, AffExpr}
            # Linear constraints
            n_linear += n_cons
        elseif F <: QuadExpr
            # Quadratic constraints (including SOCP)
            n_quadratic += n_cons
        else
            # Nonlinear constraints
            n_nonlinear += n_cons
        end
    end
    
    # Get number of nonlinear constraints from NLP backend if available
    if haskey(model.ext, :nlp_block)
        nlp_data = model.ext[:nlp_block]
        if !isnothing(nlp_data)
            n_nonlinear += length(nlp_data.constraints)
        end
    end
    
    return Dict(
        :n_variables => n_vars,
        :n_linear_constraints => n_linear,
        :n_quadratic_constraints => n_quadratic,
        :n_nonlinear_constraints => n_nonlinear,
        :total_constraints => n_linear + n_quadratic + n_nonlinear
    )
end


"""
    print_model_size(stats::Dict; model_name="Model")

Print model size statistics.
"""
function print_model_size(stats::Dict; model_name="Model")
    println("\n--- PROBLEM SIZE: $model_name ---")
    @printf("Number of variables: %d\n", stats[:n_variables])
    @printf("Linear constraints: %d\n", stats[:n_linear_constraints])
    @printf("Quadratic constraints (SOCP): %d\n", stats[:n_quadratic_constraints])
    @printf("Nonlinear constraints: %d\n", stats[:n_nonlinear_constraints])
    @printf("Total constraints: %d\n", stats[:total_constraints])
end
