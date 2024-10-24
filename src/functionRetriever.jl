module functionRetriever

# export get_real_power_loss, get_reactive_power_loss, get_substation_power, get_substation_power_cost, get_scd, get_terminal_SOC_violation

export get_real_power_loss,
    get_reactive_power_loss,
    get_substation_power,
    get_substation_power_cost,
    get_scd,
    get_terminal_SOC_violation

import JuMP: value  # Importing JuMP's value function to extract values from the model
using Parameters: @unpack  # For easier unpacking of parameters from data

# Function to get real power losses from a model
function get_real_power_loss(model, data; horizon::String="allT")
    @unpack Tset, Lset, rdict_pu, kVA_B = data
    P = model[:P]
    l = model[:l]

    if horizon == "1toT"
        loss_vs_t_kW = [kVA_B * sum(rdict_pu[i, j] * value(l[(i, j), t]) for (i, j) in Lset) for t in Tset]
        return loss_vs_t_kW
    elseif horizon == "allT"
        total_loss_kW = kVA_B * sum(rdict_pu[i, j] * value(l[(i, j), t]) for (i, j) in Lset, t in Tset)
        return total_loss_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end

# Function to get reactive power losses from a model
function get_reactive_power_loss(model, data; horizon::String="allT")
    @unpack Tset, Lset, xdict_pu, kVA_B = data
    Q = model[:Q]
    l = model[:l]

    if horizon == "1toT"
        loss_vs_t_kVAr = [kVA_B * sum(xdict_pu[i, j] * value(l[(i, j), t]) for (i, j) in Lset) for t in Tset]
        return loss_vs_t_kVAr
    elseif horizon == "allT"
        total_loss_kVAr = kVA_B * sum(xdict_pu[i, j] * value(l[(i, j), t]) for (i, j) in Lset, t in Tset)
        return total_loss_kVAr
    else
        error("Specify either '1toT' or 'allT'")
    end
end

# Function to get substation power in kW
function get_substation_power(model, data; horizon::String="allT")
    @unpack Tset, kVA_B = data
    P_Subs = model[:P_Subs]

    if horizon == "1toT"
        substation_power_vs_t_kW = [value(P_Subs[t]) * kVA_B for t in Tset]
        return substation_power_vs_t_kW
    elseif horizon == "allT"
        total_substation_power_kW = sum(value(P_Subs[t]) * kVA_B for t in Tset)
        return total_substation_power_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end

# Function to get substation power cost in dollars
function get_substation_power_cost(model, data; horizon::String="allT")
    @unpack Tset, LoadShapeCost, delta_t, kVA_B = data
    P_Subs = model[:P_Subs]

    if horizon == "1toT"
        substation_cost_vs_t_dollar = [LoadShapeCost[t] * value(P_Subs[t]) * kVA_B * delta_t for t in Tset]
        return substation_cost_vs_t_dollar
    elseif horizon == "allT"
        total_substation_cost_dollar = sum(LoadShapeCost[t] * value(P_Subs[t]) * kVA_B * delta_t for t in Tset)
        return total_substation_cost_dollar
    else
        error("Specify either '1toT' or 'allT'")
    end
end

# Function to get SCD (State of Charge Difference) in kW
function get_scd(model, data; horizon::String="allT")
    @unpack Tset, Bset, kVA_B = data
    P_c = model[:P_c]
    P_d = model[:P_d]

    if horizon == "1toT"
        scd_vs_t_kW = [kVA_B * sum(max(value(P_c[j, t]), value(P_d[j, t])) - abs(value(P_c[j, t]) - value(P_d[j, t])) for j in Bset) for t in Tset]
        return scd_vs_t_kW
    elseif horizon == "allT"
        total_scd_kW = kVA_B * sum(max(value(P_c[j, t]), value(P_d[j, t])) - abs(value(P_c[j, t]) - value(P_d[j, t])) for j in Bset, t in Tset)
        return total_scd_kW
    else
        error("Specify either '1toT' or 'allT'")
    end
end

# Function to get terminal SOC violation in kWh
function get_terminal_SOC_violation(model, data)
    @unpack T, Bset, Bref_pu, kVA_B = data
    B = model[:B]

    soc_violation_kWh = kVA_B * sum(abs(value(B[j, T]) - Bref_pu[j]) for j in Bset)
    return soc_violation_kWh
end

# Function to get any variable from the model (does not depend on horizon)
function get_variable(model, variable_name::Symbol)
    return model[variable_name]  # Fetch the variable using its symbol
end

end # module FunctionRetriever
