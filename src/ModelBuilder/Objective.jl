module Objective

export define_objective_function_t_in_Tset

using JuMP
using Parameters: @unpack, @pack!

include("Hyperparameters.jl")
import .Hyperparameters as HP

include("../helperFunctions.jl")
import .helperFunctions as HF

#region define_objective_function_t_in_Tset
"""
    define_objective_function_t_in_Tset(modelDict; Tset=nothing, tSOC_hard=false)

Define the objective function for the optimization model over a given time set.

This function sets the objective function for the optimization model stored in `modelDict`. 
It handles different objective types, including substation power cost minimization and line loss minimization, 
and optionally includes terms for state of charge (SOC) discrepancies and terminal SOC constraints.
"""
function define_objective_function_t_in_Tset(modelDict; Tset=nothing, tSOC_hard=false, relax_terminal_soc_constraint=false)
    @unpack model, data = modelDict

    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack objfun0, objfun2 = data

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow"
        # Set the objective function to zero for powerflow
        objfun = 0
        func_obj_est = nothing # no objective function (to minimize) for powerflow
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
        func_obj_est = HP.estimate_substation_power_cost 
    elseif objfun0 == "lineLossMin"
        l = model[:l]
        @unpack Lset, rdict_pu, kVA_B = data;
        objfun = kVA_B * sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset, t in Tset
        )
        func_obj_est = HP.estimate_line_losses
    end

    @pack! data = func_obj_est;

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd"
        @unpack Bset, eta_C, eta_D = data;
        η_C = eta_C
        η_D = eta_D
        alpha = HP.estimate_alpha(data)
        # println("alpha = $alpha")
        alphaAppendix = HF.trim_number_for_printing(alpha)
        @pack! data = alpha, alphaAppendix;
        α = alpha
        P_c = model[:P_c]
        P_d = model[:P_d]
        objfun += sum(
            α * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset, t in Tset
        )
    end

    # Append the gamma term only if tSOC_hard is false
    @unpack T = data;
    if !tSOC_hard && T ∈ Tset
        @unpack Bset, Bref_pu = data;
        gamma = HP.estimate_gamma(data, relax_terminal_soc_constraint=relax_terminal_soc_constraint)
        # println("gamma = $gamma")
        gammaAppendix = HF.trim_number_for_printing(gamma)
        @pack! data = gamma, gammaAppendix;
        B = model[:B]
        γ = gamma
        objfun += sum(
            1//2 * γ * (B[j, T] - Bref_pu[j])^2
            for j in Bset
        )
    end

    @objective(model, Min, objfun)

    @pack! modelDict = model, data

    return modelDict
end
#endregion

end