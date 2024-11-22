module Hyperparameters

export estimate_alpha

using Parameters: @unpack
include("../functionRetriever.jl")
import .functionRetriever as FR

function estimate_alpha(data)
    @unpack func_obj_est = data;
    if func_obj_est != nothing
        fobj_est = func_obj_est(data)
    else
        fobj_est = 0
    end
    fscd_est = estimate_fscd(data)
    alpha = fobj_est / fscd_est
    return alpha
end

function estimate_fscd(data)
    @unpack Bset, P_B_R, eta_C, eta_D, T = data;
    η_C, η_D = eta_C, eta_D
    fscd_est = sum( (1/η_D[j] - η_C[j]) * P_B_R[j] for j ∈ Bset) * T
    return fscd_est
end

function estimate_substation_power_cost(data)
    @unpack LoadShapeCost, kVA_B = data;
    C = LoadShapeCost # dollars_per_kWh
    dollars_per_kWh = C
    load_real_power_vs_t_1toT_kW = FR.get_load_real_power(data, horizon="1toT")
    line_loss_accommodation_factor = 1.1 # 10% line loss accommodation factor, just a usual value
    fcost_est = line_loss_accommodation_factor * transpose(dollars_per_kWh) * load_real_power_vs_t_1toT_kW
    
    return fcost_est
end

function estimate_line_losses(data)
    @unpack kVA_B = data
    load_real_power_allT_kW = FR.get_load_real_power(data, horizon="allT")
    line_loss_accommodation_factor = 0.1 # 10% line loss accommodation factor, just a usual value
    fPLoss_est = line_loss_accommodation_factor *  load_real_power_vs_t_1toT_kW

    return fPLoss_est
end

end # Hyperparameters module