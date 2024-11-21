module Objective

export define_objective_function_t_in_Tset

using JuMP
using Parameters: @unpack

function define_objective_function_t_in_Tset(model, data; Tset=nothing, tSOC_hard=false)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack objfun0, objfun2, Bset, kVA_B, rdict_pu, Lset, alpha, gamma, Bref_pu = data
    P_Subs = model[:P_Subs]
    P_c = model[:P_c]
    P_d = model[:P_d]
    l = model[:l]
    B = model[:B]

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow"
        # Set the objective function to zero for powerflow
        objfun = 0
    elseif objfun0 == "subsPowerCostMin"
        # Define the base objective function (generation cost minimization)
        @unpack LoadShapeCost, delta_t = data
        C = LoadShapeCost
        dollars_per_pu = kVA_B
        objfun = sum(
            dollars_per_pu * C[t] * P_Subs[t] * delta_t
            for t in Tset
        )
    elseif objfun0 == "lineLossMin"
        objfun = sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset, t in Tset
        )
    end

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd"
        @unpack eta_C, eta_D = data
        η_C = eta_C
        η_D = eta_D
        objfun += sum(
            alpha * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset, t in Tset
        )
    end

    # Append the gamma term only if tSOC_hard is false
    if !tSOC_hard
        @unpack T = data
        objfun += sum(
            gamma * (B[j, T] - Bref_pu[j])^2
            for j in Bset
        )
    end

    @objective(model, Min, objfun)

    return model
end

end