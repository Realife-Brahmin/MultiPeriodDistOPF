module Hyperparameters

export estimate_alpha,
    estimate_fscd,
    estimate_ftsoc,
    estimate_gamma,
    estimate_line_losses,
    estimate_substation_power_cost

using Parameters: @unpack
include("../functionRetriever.jl")
import .functionRetriever as FR

#region estimate_alpha
"""
    estimate_alpha(data)

Estimate the alpha coefficient for the SCD-term appended to the objective function.

This function calculates the alpha parameter based on the estimated 'regular' objective function value and the estimated magnitude of the maximum possible state of charge discrepancy.
"""
function estimate_alpha(data)
    @unpack func_obj_est = data;
    if func_obj_est !== nothing
        fobj_est = func_obj_est(data)
    else
        fobj_est = 0
    end
    fscd_est = estimate_fscd(data)
    alpha = fobj_est / fscd_est
    return alpha
end
#endregion

#region estimate_fscd
"""
    estimate_fscd(data)

Estimates the magnitude of the maximum possible state of charge discrepancy (fscd).

This function calculates the estimated state of charge discrepancy for all batteries over the entire time horizon.
"""
function estimate_fscd(data)
    @unpack Bset, P_B_R, eta_C, eta_D, T = data;
    η_C, η_D = eta_C, eta_D
    fscd_est = sum( (1/η_D[j] - η_C[j]) * P_B_R[j] for j ∈ Bset) * T
    return fscd_est
end
#endregion

#region estimate_ftsoc
"""
    estimate_ftsoc(data)

Estimate the terminal state of charge (ftsoc).

This function calculates the estimated terminal state of charge violation for all batteries given a particular tolerance.
"""
function estimate_ftsoc(data; tol_tSOC_violation_pu=1e-3)
    @unpack n_B = data;
    ftsoc_est = tol_tSOC_violation_pu^2 * n_B
    return ftsoc_est
end
#endregion

#region estimate_gamma
"""
    estimate_gamma(data)

Estimate the gamma parameter for the objective function.

This function calculates the gamma parameter based on the estimated terminal state of charge discrepancy.
"""
function estimate_gamma(data; relax_terminal_soc_constraint=false)
    @unpack func_obj_est = data;
    if relax_terminal_soc_constraint
        gamma = 0
        return gamma
    end
    if func_obj_est !== nothing
        fobj_est = func_obj_est(data)
    else
        fobj_est = 0
    end
    ftsoc_est = estimate_ftsoc(data)
    gamma = fobj_est / ftsoc_est
    return gamma
end
#endregion

#region estimate_substation_power_cost
"""
    estimate_substation_power_cost(data)

Estimate the substation power cost in dollars.

This function calculates the estimated substation power cost based on the load shape cost and the substation power values.
"""
function estimate_substation_power_cost(data)
    @unpack LoadShapeCost, kVA_B = data;
    C = LoadShapeCost # dollars_per_kWh
    dollars_per_kWh = C
    load_real_power_vs_t_1toT_kW = FR.get_load_real_power(data, horizon="1toT")
    line_loss_accommodation_factor = 1.01 # 1% line loss accommodation factor, just a usual value
    fcost_est = line_loss_accommodation_factor * transpose(dollars_per_kWh) * load_real_power_vs_t_1toT_kW
    
    return fcost_est
end
#endregion

#region estimate_line_losses
"""
    estimate_line_losses(data)

Estimate the line losses.

This function calculates the estimated line losses in kW for the entire network over the entire time horizon.
"""
function estimate_line_losses(data)
    @unpack kVA_B = data
    load_real_power_allT_kW = FR.get_load_real_power(data, horizon="allT")
    line_loss_accommodation_factor = 0.01 # 1% line loss accommodation factor, just a usual value
    fPLoss_est = line_loss_accommodation_factor * load_real_power_allT_kW

    return fPLoss_est
end
#endregion

end # Hyperparameters module