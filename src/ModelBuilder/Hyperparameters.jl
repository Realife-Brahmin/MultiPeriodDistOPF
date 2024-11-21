module Hyperparameters

export estimate_alpha

using Parameters: @unpack
include("../functionRetriever.jl")
import .functionRetriever as FR

function estimate_alpha(data)
    @unpack LoadShapeCost, kVA_B = data;
    C = LoadShapeCost # dollars_per_kWh
    dollars_per_kWh = C
    load_real_power_vs_t_1toT_kW = FR.get_load_real_power(data, horizon="1toT")
    line_loss_accommodation_factor = 1.1
    fcost_est = line_loss_accommodation_factor * transpose(dollars_per_kWh) * load_real_power_vs_t_1toT_kW
    # println("fcost_est = $fcost_est")
    alpha = 1e-3
    return alpha
end

end