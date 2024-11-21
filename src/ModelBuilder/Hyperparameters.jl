module Hyperparameters

export estimate_alpha

using Parameters: @unpack
include("../functionRetriever.jl")
import .functionRetriever as FR

function estimate_alpha(data)
    fcost_est = estimate_substation_power_cost(data)
    alpha = 1e-3
    return alpha
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

end # Hyperparameters module