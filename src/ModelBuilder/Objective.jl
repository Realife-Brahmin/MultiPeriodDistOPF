module Objective

export define_objective_function_t_in_Tset

using JuMP
using Parameters: @unpack

include("Hyperparameters.jl")
import .Hyperparameters as HP

function define_objective_function_t_in_Tset(model, data; Tset=nothing, tSOC_hard=false)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack objfun0, objfun2 = data

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow"
        # Set the objective function to zero for powerflow
        objfun = 0
    elseif objfun0 == "subsPowerCostMin"
        # Define the base objective function (generation cost minimization)
        @unpack LoadShapeCost, delta_t, kVA_B = data
        C = LoadShapeCost
        dollars_per_pu = kVA_B
        P_Subs = model[:P_Subs]
        objfun = sum(
            dollars_per_pu * C[t] * P_Subs[t] * delta_t
            for t in Tset
        )
    elseif objfun0 == "lineLossMin"
        l = model[:l]
        @unpack Lset, rdict_pu = data;
        objfun = sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset, t in Tset
        )
    end

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd"
        @unpack Bset, eta_C, eta_D = data;
        η_C = eta_C
        η_D = eta_D
        alpha = HP.estimate_alpha(data)
        α = alpha
        P_c = model[:P_c]
        P_d = model[:P_d]
        objfun += sum(
            α * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset, t in Tset
        )
    end

    # Append the gamma term only if tSOC_hard is false
    if !tSOC_hard
        @unpack T, Bset, gamma, Bref_pu = data;
        B = model[:B]
        γ = gamma
        objfun += sum(
            γ * (B[j, T] - Bref_pu[j])^2
            for j in Bset
        )
    end

    @objective(model, Min, objfun)

    return model
end

end