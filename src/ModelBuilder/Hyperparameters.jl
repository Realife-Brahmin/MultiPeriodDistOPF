module Hyperparameters

export estimate_alpha, estimate_gamma

using Parameters: @unpack

function estimate_alpha(model, data; Tset=nothing)
    alpha = 1e-3
    return alpha
end

end