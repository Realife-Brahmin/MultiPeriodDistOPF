module Hyperparameters

export estimate_alpha

using Parameters: @unpack
include("../functionRetriever.jl")
import .functionRetriever as FR

function estimate_alpha(data)
    @unpack LoadShapeCost, kVA_B = data;
    C = LoadShapeCost # dollars_per_kWh
    dollars_per_kWh = C
    dollars_per_pu = C * kVA_B
    load_real_power_vs_t_1toT_kW = FR.get_load_real_power(data, horizon="1toT")
    line_loss_accommodation_factor = 1.1
    fcost_est = dollars_per_pu .* load_real_power_vs_t_1toT_kW * line_loss_accommodation_factor
    println("fcost_est = $fcost_est")
    alpha = 1e-3
    return alpha
end

end