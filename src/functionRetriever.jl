module functionRetriever

export get_battery_reactive_power,
    get_battery_real_power,
    get_battery_reactive_power_transaction_magnitude,
    get_battery_real_power_transaction_magnitude,
    get_load_reactive_power,
    get_load_real_power,
    get_loss_reactive_power,
    get_loss_real_power,
    get_model_size,
    get_pv_reactive_power,
    get_pv_real_power,
    get_scd,
    get_solution_time,
    get_substation_power_cost,
    get_substation_reactive_power,
    get_substation_real_power,
    get_substation_real_power_peak,
    get_terminal_SOC_violation,
    get_total_generation_reactive_power,
    get_total_generation_real_power,
    get_static_capacitor_reactive_power

using JuMP
using MathOptInterface
const MOI = MathOptInterface
# import JuMP: value, solve_time  # Importing JuMP's value function to extract values from the modelVals
using Parameters: @unpack  # For easier unpacking of parameters from data

#region get_loss_real_power
"""
    get_loss_real_power(modelDict; horizon::String="allT")

Calculate real power losses from modelVals.
"""
function get_loss_real_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Lset, rdict_pu, kVA_B = data
    P = modelVals[:P]
    l = modelVals[:l]

    if horizon == "1toT"
        loss_real_power_vs_t_1toT_kW = [kVA_B * sum(rdict_pu[i, j] * l[(i, j), t] for (i, j) in Lset) for t in Tset]
        return loss_real_power_vs_t_1toT_kW
    elseif horizon == "allT"
        loss_real_power_allT_kW = kVA_B * sum(rdict_pu[i, j] * l[(i, j), t] for (i, j) in Lset, t in Tset)
        return loss_real_power_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_loss_reactive_power
"""
    get_loss_reactive_power(modelDict; horizon::String="allT")

Calculate reactive power losses from modelVals.
"""
function get_loss_reactive_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Lset, xdict_pu, kVA_B = data
    Q = modelVals[:Q]
    l = modelVals[:l]

    if horizon == "1toT"
        loss_reactive_power_vs_t_1toT_kVAr = [kVA_B * sum(xdict_pu[i, j] * l[(i, j), t] for (i, j) in Lset) for t in Tset]
        return loss_reactive_power_vs_t_1toT_kVAr
    elseif horizon == "allT"
        loss_reactive_power_allT_kVAr = kVA_B * sum(xdict_pu[i, j] * l[(i, j), t] for (i, j) in Lset, t in Tset)
        return loss_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_substation_reactive_power
"""
    get_substation_reactive_power(modelDict; horizon::String="allT")

Calculate substation reactive power in kVAr from modelVals.
"""
function get_substation_reactive_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, L1set, kVA_B = data
    Q = modelVals[:Q]

    if horizon == "1toT"
        # Compute reactive power at the substation for each time step as the sum of Q values in L1set
        substation_reactive_power_vs_t_1toT_kVAr = [
            kVA_B * sum(Q[(i, j), t] for (i, j) in L1set) for t in Tset
        ]
        return substation_reactive_power_vs_t_1toT_kVAr
    elseif horizon == "allT"
        # Compute total reactive power at the substation over all time steps
        substation_reactive_power_allT_kVAr = kVA_B * sum(Q[(i, j), t] for (i, j) in L1set, t in Tset)
        return substation_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_substation_real_power
"""
    get_substation_real_power(modelDict; horizon::String="allT")

Calculate substation real power in kW from modelVals.
"""
function get_substation_real_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, kVA_B = data
    P_Subs = modelVals[:P_Subs]

    if horizon == "1toT"
        substation_real_power_vs_t_1toT_kW = [P_Subs[t] * kVA_B for t in Tset]
        return substation_real_power_vs_t_1toT_kW
    elseif horizon == "allT"
        substation_real_power_allT_kW = sum(P_Subs[t] * kVA_B for t in Tset)
        return substation_real_power_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_substation_power_cost
"""
    get_substation_power_cost(modelDict; horizon::String="allT")

Calculate substation power cost in dollars from modelVals.
"""
function get_substation_power_cost(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, LoadShapeCost, delta_t, kVA_B = data
    P_Subs = modelVals[:P_Subs]

    if horizon == "1toT"
        substation_power_cost_vs_t_1toT_dollar = [LoadShapeCost[t] * P_Subs[t] * kVA_B * delta_t for t in Tset]
        return substation_power_cost_vs_t_1toT_dollar
    elseif horizon == "allT"
        substation_power_cost_allT_dollar = sum(LoadShapeCost[t] * P_Subs[t] * kVA_B * delta_t for t in Tset)
        return substation_power_cost_allT_dollar
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_scd
"""
    get_scd(modelDict; horizon::String="allT")

Calculate Simultaneous Charging and Discharging (SCD) in kW from modelVals. 
"""
function get_scd(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Bset, kVA_B = data
    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]

    if horizon == "1toT"
        scd_vs_t_1toT_kW = [kVA_B * sum(max(P_c[j, t], P_d[j, t]) - abs(P_c[j, t] - P_d[j, t]) for j in Bset) for t in Tset]
        return scd_vs_t_1toT_kW
    elseif horizon == "allT"
        scd_allT_kW = kVA_B * sum(max(P_c[j, t], P_d[j, t]) - abs(P_c[j, t] - P_d[j, t]) for j in Bset, t in Tset)
        return scd_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_terminal_SOC_violation
"""
    get_terminal_SOC_violation(modelDict)

Calculate terminal State of Charge (SOC) violation in kWh from modelVals.
"""
function get_terminal_SOC_violation(modelDict)
    @unpack modelVals, data = modelDict
    @unpack T, Bset, Bref_pu, kVA_B = data
    B = modelVals[:B]

    soc_violation_kWh = kVA_B * sum(abs(B[j, T] - Bref_pu[j]) for j in Bset)
    return soc_violation_kWh
end
#endregion

#region get_battery_real_power
"""
    get_battery_real_power(modelDict; horizon::String="allT")

Calculate total battery real power in kW (P_d - P_c) from modelVals.
"""
function get_battery_real_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Bset, kVA_B = data
    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]

    if horizon == "1toT"
        battery_real_power_vs_t_1toT_kW = [kVA_B * sum(P_d[j, t] - P_c[j, t] for j in Bset) for t in Tset]
        return battery_real_power_vs_t_1toT_kW
    elseif horizon == "allT"
        battery_real_power_allT_kW = kVA_B * sum(P_d[j, t] - P_c[j, t] for j in Bset, t in Tset)
        return battery_real_power_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_battery_reactive_power
"""
    get_battery_reactive_power(modelDict; horizon::String="allT")

Calculate total battery reactive power in kVAr from modelVals.
"""
function get_battery_reactive_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Bset, kVA_B = data
    q_B = modelVals[:q_B]

    if horizon == "1toT"
        battery_reactive_power_vs_t_1toT_kVAr = [kVA_B * sum(q_B[j, t] for j in Bset) for t in Tset]
        return battery_reactive_power_vs_t_1toT_kVAr
    elseif horizon == "allT"
        battery_reactive_power_allT_kVAr = kVA_B * sum(q_B[j, t] for j in Bset, t in Tset)
        return battery_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

function get_model_size(modelDict)
    @unpack model = modelDict
    # Retrieve number of decision variables
    num_decvars = JuMP.num_variables(model)

    # Retrieve number of linear constraints
    num_lincons =
        JuMP.num_constraints(model, AffExpr, MOI.EqualTo{Float64}) +
        JuMP.num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) +
        JuMP.num_constraints(model, AffExpr, MOI.LessThan{Float64})

    # Retrieve number of nonlinear constraints
    num_nonlincons =
        JuMP.num_constraints(model, QuadExpr, MOI.EqualTo{Float64}) +
        JuMP.num_constraints(model, QuadExpr, MOI.GreaterThan{Float64}) +
        JuMP.num_constraints(model, QuadExpr, MOI.LessThan{Float64})

    # Pack results into a dictionary
    modelSizeDict = Dict(
        :num_decvars => num_decvars,
        :num_lincons => num_lincons,
        :num_nonlincons => num_nonlincons
    )

    return modelSizeDict
end

#region get_pv_real_power
"""
    get_pv_real_power(modelDict; horizon::String="allT")

Calculate total PV real power in kW from p_D_pu in data.
"""
function get_pv_real_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Dset, p_D_pu, kVA_B = data

    if horizon == "1toT"
        pv_real_power_vs_t_1toT_kW = [kVA_B * sum(p_D_pu[(j, t)] for j in Dset) for t in Tset]
        return pv_real_power_vs_t_1toT_kW
    elseif horizon == "allT"
        pv_real_power_allT_kW = kVA_B * sum(p_D_pu[(j, t)] for j in Dset, t in Tset)
        return pv_real_power_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_pv_reactive_power
"""
    get_pv_reactive_power(modelDict; horizon::String="allT")

Calculate total PV reactive power in kVAr (q_D) from modelVals.
"""
function get_pv_reactive_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Dset, kVA_B = data
    q_D = modelVals[:q_D]

    if horizon == "1toT"
        pv_reactive_power_vs_t_1toT_kVAr = [kVA_B * sum(q_D[j, t] for j in Dset) for t in Tset]
        return pv_reactive_power_vs_t_1toT_kVAr
    elseif horizon == "allT"
        pv_reactive_power_allT_kVAr = kVA_B * sum(q_D[j, t] for j in Dset, t in Tset)
        return pv_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_battery_real_power_transaction_magnitude
"""
    get_battery_real_power_transaction_magnitude(modelDict; horizon::String="allT")

Compute battery real power transaction magnitude |P_c - P_d|.
"""
function get_battery_real_power_transaction_magnitude(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Bset, kVA_B = data
    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]

    if horizon == "1toT"
        battery_real_power_magnitude_vs_t_1toT_kW = [kVA_B * sum(abs(P_c[j, t] - P_d[j, t]) for j in Bset) for t in Tset]
        return battery_real_power_magnitude_vs_t_1toT_kW
    elseif horizon == "allT"
        battery_real_power_magnitude_allT_kW = kVA_B * sum(abs(P_c[j, t] - P_d[j, t]) for j in Bset, t in Tset)
        return battery_real_power_magnitude_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_battery_reactive_power_transaction_magnitude
"""
    get_battery_reactive_power_transaction_magnitude(modelDict; horizon::String="allT")

Compute battery reactive power transaction magnitude |q_B|.
"""
function get_battery_reactive_power_transaction_magnitude(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict
    @unpack Tset, Bset, kVA_B = data
    q_B = modelVals[:q_B]

    if horizon == "1toT"
        battery_reactive_power_magnitude_vs_t_1toT_kVAr = [kVA_B * sum(abs(q_B[j, t]) for j in Bset) for t in Tset]
        return battery_reactive_power_magnitude_vs_t_1toT_kVAr
    elseif horizon == "allT"
        battery_reactive_power_magnitude_allT_kVAr = kVA_B * sum(abs(q_B[j, t]) for j in Bset, t in Tset)
        return battery_reactive_power_magnitude_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_static_capacitor_reactive_power
"""
    get_static_capacitor_reactive_power(modelDict; horizon::String="allT")

Compute total static capacitor reactive power generation.
"""
function get_static_capacitor_reactive_power(modelDict; horizon::String="allT")
    @unpack modelVals, data = modelDict;
    @unpack Tset, kVA_B, T = data

    if haskey(modelVals, :q_C)
        q_C = modelVals[:q_C]

        if horizon == "1toT"
            static_cap_reactive_power_vs_t_1toT_kVAr = [kVA_B * q_C[t] for t in Tset]
            return static_cap_reactive_power_vs_t_1toT_kVAr
        elseif horizon == "allT"
            static_cap_reactive_power_allT_kVAr = kVA_B * sum(q_C[t] for t in Tset)
            return static_cap_reactive_power_allT_kVAr
        else
            error("Specify either '1toT' or 'allT'")
        end
    else
        if horizon == "1toT"
            return zeros(Float64, T)  # If static capacitor does not exist, return an array of zeros
        elseif horizon == "allT"
            return 0.0  # If static capacitor does not exist, return 0
        else
            error("Specify either '1toT' or 'allT'")
        end
    end
end
#endregion

#region get_substation_real_power_peak
"""
    get_substation_real_power_peak(modelDict)

Get the peak substation real power over the horizon.
"""
function get_substation_real_power_peak(modelDict)
    @unpack modelVals, data = modelDict;
    @unpack Tset, kVA_B = data
    P_Subs = modelVals[:P_Subs]

    # Calculate the peak substation power (kW) by finding the maximum across all time steps
    peak_substation_power_kW = maximum([P_Subs[t] * kVA_B for t in Tset])

    return peak_substation_power_kW
end
#endregion

#region get_total_generation_real_power
"""
    get_total_generation_real_power(modelDict; horizon::String="allT")

Compute total real power generation.
"""
function get_total_generation_real_power(modelDict; horizon::String="allT")

    if horizon == "1toT"
        # Battery + PV + Static Capacitor (if applicable)
        battery_real_power_vs_t_1toT_kW = get_battery_real_power(modelDict, horizon="1toT")
        pv_real_power_vs_t_1toT_kW = get_pv_real_power(modelDict, horizon="1toT")
        total_gen_real_power_vs_t_1toT_kW = battery_real_power_vs_t_1toT_kW + pv_real_power_vs_t_1toT_kW
        return total_gen_real_power_vs_t_1toT_kW
    elseif horizon == "allT"
        battery_real_power_allT_kW = get_battery_real_power(modelDict, horizon="allT")
        pv_real_power_allT_kW = get_pv_real_power(modelDict, horizon="allT")
        total_gen_real_power_allT_kW = battery_real_power_allT_kW + pv_real_power_allT_kW
        return total_gen_real_power_allT_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_total_generation_reactive_power
"""
    get_total_generation_reactive_power(modelDict; horizon::String="allT")

Compute total reactive power generation.
"""
function get_total_generation_reactive_power(modelDict; horizon::String="allT")
    # Battery + PV + Static Capacitor (if applicable)

    if horizon == "1toT"
        battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(modelDict, horizon="1toT")
        pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(modelDict, horizon="1toT")
        static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="1toT")
        total_gen_reactive_power_vs_t_1toT_kVAr = battery_reactive_power_vs_t_1toT_kVAr + pv_reactive_power_vs_t_1toT_kVAr + static_cap_reactive_power_vs_t_1toT_kVAr
        return total_gen_reactive_power_vs_t_1toT_kVAr
    elseif horizon == "allT"
        battery_reactive_power_allT_kVAr = get_battery_reactive_power(modelDict, horizon="allT")
        pv_reactive_power_allT_kVAr = get_pv_reactive_power(modelDict, horizon="allT")
        static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="allT")
        total_gen_reactive_power_allT_kVAr = battery_reactive_power_allT_kVAr + pv_reactive_power_allT_kVAr + static_cap_reactive_power_allT_kVAr
        return total_gen_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_load_real_power
"""
    get_load_real_power(data; horizon::String="allT")

Get total real load over the horizon.
"""
function get_load_real_power(data; horizon::String="allT")
    @unpack Tset, NLset, p_L_pu, kVA_B = data

    if horizon == "1toT"
        load_real_power_vs_t_kW = [kVA_B * sum(p_L_pu[(j, t)] for j ∈ NLset) for t ∈ Tset]
        return load_real_power_vs_t_kW
    elseif horizon == "allT"
        load_real_power_kW_allT = kVA_B * sum(p_L_pu[(j, t)] for j ∈ NLset, t ∈ Tset)
        return load_real_power_kW_allT
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_load_reactive_power
"""
    get_load_reactive_power(data; horizon::String="allT")

Get total reactive load over the horizon.
"""
function get_load_reactive_power(data; horizon::String="allT")
    @unpack Tset, NLset, q_L_pu, kVA_B = data

    if horizon == "1toT"
        load_reactive_power_vs_t_kVAr = [kVA_B * sum(q_L_pu[(j, t)] for j ∈ NLset, t ∈ Tset)]  # Assuming p_L_pu contains reactive power
        return load_reactive_power_vs_t_kVAr
    elseif horizon == "allT"
        load_reactive_power_allT_kVAr = kVA_B * sum(q_L_pu[(j, t)] for j ∈ NLset, t in Tset)  # Assuming p_L_pu contains reactive power
        return load_reactive_power_allT_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end
#endregion

#region get_solution_time
"""
    get_solution_time(modelDict)

Compute solution time in seconds.
"""
function get_solution_time(modelDict)
    @unpack modelVals = modelDict;
    solution_time = modelVals[:solve_time]  # Retrieves solution time 

    return solution_time
end
#endregion

end # module FunctionRetriever
