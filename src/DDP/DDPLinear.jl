module DDPLinear

export
    check_for_ddp_convergence,
    compute_and_store_dual_variables, 
    build_ForwardStep_1ph_NL_model_t_is_1,
    build_ForwardStep_1ph_NL_model_t_in_2toTm1,
    build_ForwardStep_1ph_NL_model_t_is_T,
    DDPModel,
    forward_pass_1ph_L,
    ForwardStep_1ph_NL_t_is_1,
    ForwardStep_1ph_NL_t_in_2toTm1,
    ForwardStep_1ph_NL_t_is_T,
    optimize_ForwardStep_1ph_NL_model_t_is_1,
    optimize_MPOPF_1ph_L_DDP


include("../computeOutputs.jl")
import .computeOutputs as CO

include("../functionRetriever.jl")
import .functionRetriever as FR

include("../ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("../ModelCopier/ModelCopier.jl")
import .ModelCopier as MC

include("../helperFunctions.jl")
using .helperFunctions

include("../SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger

include("../Exporter.jl")
import .Exporter as Exporter

using Crayons
using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function compute_and_store_dual_variables(ddpModel, model_t0; Tset=nothing)
    if Tset === nothing
        println("No Tset defined for backward pass, so using the Tset from the data")
        @unpack data = ddpModel
        Tset = data[:Tset]
    end
    @unpack k_ddp, data, mu, lambda_lo, lambda_up = ddpModel
    μ = mu
    λ_lo = lambda_lo
    λ_up = lambda_up
    @unpack T, Bset = data
    # Update mu, lambda_lo, and lambda_up values post optimization
    t_ddp = Tset[1]
    for j in Bset
        if t_ddp == 1
            constraint_name = "h_SOC_j^{t=1}_Initial_SOC_Node_j_$(j)_t1"
        elseif 2 <= t_ddp <= T
            constraint_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t_ddp)"
        else
            @error "Invalid value of t_ddp: $t_ddp"
            return
        end
        constraint_j_t0 = constraint_by_name(model_t0, constraint_name)
        # println("constraint_j_t0 = ", constraint_j_t0)
        μ[j, t_ddp, k_ddp] = -dual(constraint_j_t0)

        # Update lambda_lo and lambda_up
        lambda_lo_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t_ddp)"
        lambda_up_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t_ddp)"
        constraint_lambda_lo = constraint_by_name(model_t0, lambda_lo_name)
        constraint_lambda_up = constraint_by_name(model_t0, lambda_up_name)
        λ_lo[j, t_ddp, k_ddp] = -dual(constraint_lambda_lo)
        λ_up[j, t_ddp, k_ddp] = -dual(constraint_lambda_up)
    end

    crayon_update = Crayon(foreground=:light_blue, background=:white, bold=true)
    # println(crayon_update("Backward Pass k_ddp = $(k_ddp): μ, λ_lo, and λ_up values for t_ddp = $(t_ddp)"))
    # for j in Bset
    #     println(crayon_update("j = $j: μ_t_k = ", trim_number_for_printing(μ[j, t_ddp, k_ddp], sigdigits=2)))
    #     if t_ddp != T
    #         println(crayon_update("j = $j: μ_tp1_km1 = ", trim_number_for_printing(μ[j, t_ddp+1, k_ddp-1], sigdigits=2)))
    #     end
    #     println(crayon_update("j = $j: λ_lo_t_k = ", trim_number_for_printing(λ_lo[j, t_ddp, k_ddp], sigdigits=2)))
    #     println(crayon_update("j = $j: λ_up_t_k = ", trim_number_for_printing(λ_up[j, t_ddp, k_ddp], sigdigits=2)))
    # end

    # Print KKT balance equation
    # Todo: Insert provision of gamma term for KKT balance equation in case terminal SOC constraint is enforced
    # println(crayon_update("KKT balance equation for t_ddp = $(t_ddp):"))
    # for j in Bset
    #     if t_ddp < T
    #         balance = -λ_lo[j, t_ddp, k_ddp] + λ_up[j, t_ddp, k_ddp] + μ[j, t_ddp, k_ddp] - μ[j, t_ddp+1, k_ddp-1]
    #     else
    #         balance = -λ_lo[j, t_ddp, k_ddp] + λ_up[j, t_ddp, k_ddp] + μ[j, t_ddp, k_ddp]
    #     end
    #     println(crayon_update("j = $j: KKT balance = ", trim_number_for_printing(balance, sigdigits=2)))
    # end

    mu = μ
    lambda_lo = λ_lo
    lambda_up = λ_up
    @pack! ddpModel = mu, lambda_lo, lambda_up

    return ddpModel
end

# Define the function to compute α_fpi adaptively
function compute_alpha_fpi(alpha_fpi0, gamma_fpi, k_ddp)
    return alpha_fpi0 * gamma_fpi^(k_ddp - 2)
end

# Define the function to compute the interpolated value
function get_interpolated_value(xn, xnm1, alpha)
    diff = xn - xnm1
    return (alpha*xnm1 + xn) / (1 + alpha)
end

function build_ForwardStep_1ph_NL_model_t_is_1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, data, mu = ddpModel;

    if t_ddp != 1
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end
    verbose = true

    if k_ddp >= 1
        # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [1]
        modelDict = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0) # an unsolved model
        model_t0 = modelDict[:model]
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass)

    objfun_expr_t0_without_mu_terms = objective_function(model_t0) 
    μ = mu
    @unpack Bset, alpha_fpi, gamma_fpi = data;
    α_fpi0 = alpha_fpi
    γ_fpi = gamma_fpi
    MU = Dict()
    α_fpi = compute_alpha_fpi(α_fpi0, γ_fpi, k_ddp)
    myprintln(true, "FP$(k_ddp): α_fpi = $α_fpi")

    for j ∈ Bset
        if k_ddp == 1
            MU[j, t_ddp+1] = μ[j, t_ddp+1, k_ddp-1]
        elseif k_ddp >= 2
            MU[j, t_ddp+1] = get_interpolated_value(μ[j, t_ddp+1, k_ddp-1], μ[j, t_ddp+1, k_ddp-2], α_fpi)
            MU_used_str = trim_number_for_printing(MU[j, t_ddp+1], sigdigits=2)
            MU_not_used_str = trim_number_for_printing(μ[j, t_ddp+1, k_ddp-1], sigdigits=2)
            if j in Bset[1]
                # myprintln(true, "FP$(k_ddp): μ[$(j), $(t_ddp+1)] = $(MU_used_str) instead of $(MU_not_used_str)")
            end
        else
            @error "Invalid value of k_ddp: $k_ddp"
            return
        end
    end

    # objfun_expr_t0_k_with_mu_terms = objfun_expr_t0_without_mu_terms + 
    # sum(μ[j, t_ddp+1, k_ddp-1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
    objfun_expr_t0_k_with_mu_terms = objfun_expr_t0_without_mu_terms + 
    sum(MU[j, t_ddp+1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset)
    @objective(model_t0, Min, objfun_expr_t0_k_with_mu_terms)

    # B0 values are already set in the model so no need to fix them separately

    # Now model_t0 completely represents the (yet to be solved) Forward Step model for t = 1
    @unpack models_ddp_vs_t_vs_k = ddpModel;
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel

end

function build_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, modelVals, data, mu = ddpModel
    @unpack T = data;
    if !(2 <= t_ddp <= T-1)
        @error "t_ddp = $(t_ddp) is not in [2, T-1]"
        return
    end

    verbose = true
    # if k_ddp == 1
    if k_ddp >= 1
        # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be something like [2] or [3] or ... or [T-1]
        modelDict_t0 = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass) and previous time-step's SOC values (forward pass)

    B_model = modelVals[:B]
    @unpack Bset = data;

    crayon_light_red1 = Crayon(background=:light_red, foreground=:white, bold=true)
    # println(crayon_light_red1("Printing the previous time-step's SOC values to be used for the current time-step"))
    for j ∈ Bset
        # println(crayon_light_red1("B_model[$j, $(t_ddp - 1)] = $(B_model[j, t_ddp - 1])"))
        fix(model_t0[:B][j, t_ddp - 1], B_model[j, t_ddp - 1])
    end

    # Update the model with the future dual variables from the last iteration (backward pass)

    objfun_expr_t0_without_mu_terms = objective_function(model_t0)
    μ = mu
    @unpack Bset, alpha_fpi, gamma_fpi = data
    α_fpi0 = alpha_fpi
    γ_fpi = gamma_fpi
    MU = Dict()
    α_fpi = compute_alpha_fpi(α_fpi0, γ_fpi, k_ddp)
    for j ∈ Bset
        if k_ddp == 1
            MU[j, t_ddp+1] = μ[j, t_ddp+1, k_ddp-1]
        elseif k_ddp >= 2
            MU[j, t_ddp+1] = get_interpolated_value(μ[j, t_ddp+1, k_ddp-1], μ[j, t_ddp+1, k_ddp-2], α_fpi)
            MU_used_str = trim_number_for_printing(MU[j, t_ddp+1], sigdigits=2)
            MU_not_used_str = trim_number_for_printing(μ[j, t_ddp+1, k_ddp-1], sigdigits=2)
            if j in Bset[1]
                # myprintln(true, "FP$(k_ddp): μ[$(j), $(t_ddp+1)] = $(MU_used_str) instead of $(MU_not_used_str)")
            end
        else
            @error "Invalid value of k_ddp: $k_ddp"
            return
        end
    end

    # objfun_expr_t0_k_with_mu_terms = objfun_expr_t0_without_mu_terms + 
    # sum(μ[j, t_ddp+1, k_ddp-1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
    objfun_expr_t0_k_with_mu_terms = objfun_expr_t0_without_mu_terms +
                                     sum(MU[j, t_ddp+1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset)

    @objective(model_t0, Min, objfun_expr_t0_k_with_mu_terms)

    # Now model_t0 completely represents the (yet to be solved) Forward Step model for t = t_ddp ∈ [2, T-1]

    @unpack models_ddp_vs_t_vs_k = ddpModel
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel

end

function build_ForwardStep_1ph_NL_model_t_is_T(ddpModel;
    verbose::Bool=false)
    
    @unpack k_ddp, t_ddp, modelVals, data = ddpModel;
    @unpack T = data;
    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to T = $(T)"
        return
    end

    # if k_ddp == 1
    verbose = true
    if k_ddp >= 1
        # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [T]
        modelDict_t0 = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the previous time-step's SOC values (forward pass)
    B_model = modelVals[:B]
    @unpack Bset = data;

    crayon_light_red1 = Crayon(background=:light_red, foreground=:white, bold=true)
    # println(crayon_light_red1("Printing the previous time-step's SOC values to be used for the current time-step"))
    for j ∈ Bset
        # println(crayon_light_red1("B_model[$j, $(t_ddp - 1)] = $(B_model[j, t_ddp - 1])"))
        fix(model_t0[:B][j, t_ddp-1], B_model[j, t_ddp-1])
    end

    # No need of updating objective function, as no μ term is present in the objective function for the terminal time-step

    # Now model_t0 completely represents the (yet to be solved) Forward Step model for t = T

    @unpack models_ddp_vs_t_vs_k = ddpModel
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    return ddpModel
end

function check_for_ddp_convergence(ddpModel; 
    verbose::Bool=false,
    soc_change=false,
    print_soc=false,
    mu_change=true,
    print_mu=true,
    fval_change=true,
    print_fval=false)

    @unpack k_ddp, maxiter, models_ddp_vs_t_vs_k, data, mu = ddpModel
    @unpack Tset, Bset, T = data

    verbose = true
    

    Bset_to_print = length(Bset) > 1 ? [Bset[end]] : Bset

    # So anyway, let's see what the μ(s) look like:
    # Define crayons
    crayon_green = Crayon(foreground=:green, bold=true)
    crayon_blue = Crayon(foreground=:blue, bold=true)
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_light_blue = Crayon(foreground=:light_blue, bold=true)

    if print_mu
        # Iterate over batteries in Bset_to_print
        for (line_idx, j) in enumerate(Bset_to_print)
            # Alternate colors per line (Forward Pass label)
            crayon_fp = (line_idx % 2 == 1) ? crayon_green : crayon_blue
            print(crayon_fp("FP$(lpad(k_ddp, 2, '0')): "))

            # Iterate over time steps, alternate colors per time step
            for (t_idx, t) in enumerate(Tset)
                # Alternate colors per time step
                crayon_t = ((t_idx + line_idx) % 2 == 1) ? crayon_light_green : crayon_light_blue

                # Retrieve and format the current μ value
                mu_current = trim_number_for_printing(mu[(j, t, k_ddp)], sigdigits=2)

                # Print formatted μ value
                print(crayon_t("$mu_current "))
            end

            # Finish the current line
            println()
        end
    end

    @unpack B_R_pu = data
    if print_soc
        # Iterate over batteries in Bset_to_print
        for (line_idx, j) in enumerate(Bset_to_print)
            # Alternate colors per line (Forward Pass label)
            crayon_fp = (line_idx % 2 == 1) ? crayon_green : crayon_blue
            print(crayon_fp("FP$(lpad(k_ddp, 2, '0')): "))

            # Iterate over time steps, alternate colors per time step
            for (t_idx, t) in enumerate(Tset)
                model_current = models_ddp_vs_t_vs_k[t, k_ddp]
                # model_previous = models_ddp_vs_t_vs_k[T, k_ddp-1]

                # Create dictionaries for current and previous models
                var_dict_current = MC.create_variable_dict(model_current)
                # var_dict_previous = MC.create_variable_dict(model_previous)

                # Alternate colors per time step
                crayon_t = ((t_idx + line_idx) % 2 == 1) ? crayon_light_green : crayon_light_blue

                # Retrieve and format the current μ value
                var_name = Symbol("B[$j,$t]")
                B_j_t_current = var_dict_current[var_name]
                soc_percent_j_t_current = B_j_t_current / B_R_pu[j] * 100.0
                soc_percent_j_t_current_str = trim_number_for_printing(soc_percent_j_t_current, sigdigits=2)

                # Print formatted μ value
                print(crayon_t("$soc_percent_j_t_current_str "))
            end

            # Finish the current line
            println()
        end
    end

    # Criterion 1: Check if k_ddp has crossed the maxiter threshold
    if k_ddp >= maxiter
        println("Maximum iterations reached: $k_ddp")
        iterLimitReached = true
        shouldStop = true
        @pack! ddpModel = iterLimitReached, shouldStop
        return ddpModel
    end

    if k_ddp == 1 # No need to check for convergence at the first iteration
        return ddpModel
    end

    max_discrepancy = 0.0
    threshold_soc = 1e-3
    threshold_mu = 1e-2
    threshold_fval = 1e-2
    all_under_threshold = true

    crayon_red_neg = Crayon(foreground=:red, bold=true, negative=true)
    crayon_green = Crayon(background=:green, foreground=:white, bold=true)
    crayon_blue = Crayon(foreground=:blue, bold=true)
    crayon_blue_neg = Crayon(foreground=:blue, bold=true, negative=true)

    # Criterion 2: Check the magnitude of updates in fval

    if fval_change
        @unpack k_best, fval_best, fval_best_hit = ddpModel
        @unpack PSubsCost_dollar_vs_k = ddpModel;
        PSubsCost_current = PSubsCost_dollar_vs_k[k_ddp]
        PSubsCost_previous = PSubsCost_dollar_vs_k[k_ddp-1]
        diff_fval_best_vs_now = abs(fval_best - PSubsCost_current)
        # if diff_fval_best_vs_now <= threshold_fval
        discrepancy = abs(PSubsCost_current - PSubsCost_previous)
        if discrepancy > threshold_fval
            all_under_threshold = false
            myprintln(print_fval, "*******")
            myprintln(print_fval, crayon_red_neg("Previous value of PSubsCost_allT_dollar = $PSubsCost_previous"))
            myprintln(print_fval, crayon_red_neg("Current value of PSubsCost_allT_dollar = $PSubsCost_current"))
            myprintln(true, "*******")
        else
            myprintln(true, "FP$(k_ddp): Objective function value update is under the threshold.")
        end
    end

    # println(crayon_green("Checking convergence for B values:"))

    # Criterion 3: Check the magnitude of updates in SOC values

    if soc_change
        for t_ddp in Tset
            model_current = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
            model_previous = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1]
            # model_previous = models_ddp_vs_t_vs_k[T, k_ddp-1]

            # Create dictionaries for current and previous models
            var_dict_current = MC.create_variable_dict(model_current)
            var_dict_previous = MC.create_variable_dict(model_previous)

            for j in Bset
                var_name = Symbol("B[$j,$t_ddp]")
                value_current = var_dict_current[var_name]
                value_previous = var_dict_previous[var_name]
                discrepancy = abs(value_current - value_previous)
                max_discrepancy = max(max_discrepancy, discrepancy)

                if max_discrepancy > threshold
                    all_under_threshold = false
                    if discrepancy > threshold_soc
                        # myprintln(verbose, "Exceeding update tolerance: var_name = $var_name, discrepancy = $discrepancy")
                        if j in Bset_to_print
                            println(crayon_red_neg("Previous value of var $(var_name) = $value_previous"))
                            println(crayon_red_neg("Current value of var $(var_name) = $value_current"))
                        end
                    else
                        if j in Bset_to_print
                            # println(crayon_green("var_name = $var_name, discrepancy = $discrepancy"))
                        end
                    end
                end
            end
        end
        if all_under_threshold
            myprintln(true, "All SOC updates are under the threshold.")
        end
    end

    # Criterion 4: Check the magnitude of updates in μ values

    if mu_change

        for t in Tset
            for j in Bset[1]
                if haskey(mu, (j, t, k_ddp)) && haskey(mu, (j, t, k_ddp - 1))
                    mu_current = mu[(j, t, k_ddp)]
                    mu_previous = mu[(j, t, k_ddp - 1)]
                    discrepancy = abs(mu_current - mu_previous)
                    max_discrepancy = max(max_discrepancy, discrepancy)
                    if discrepancy > threshold_mu
                        all_under_threshold = false
                        # myprintln(verbose, "Exceeding update tolerance: mu[$j, $t, $k_ddp], discrepancy = $discrepancy")
                        if j in Bset_to_print
                        # println(crayon_blue_neg("Previous value of mu[$j, $t, $(k_ddp-1)] = $mu_previous"))
                        # println(crayon_blue_neg("Current value of mu[$j, $t, $k_ddp] = $mu_current"))
                        end
                    else
                        # if j in Bset_to_print
                        #     println(crayon_blue("mu[$j, $t, $k_ddp], discrepancy = $discrepancy"))
                        # end
                    end
                else
                    if j in Bset_to_print
                        # println(crayon_red_neg("mu[$j, $t, $k_ddp] or mu[$j, $t, $(k_ddp-1)] not found"))
                    end
                end
            end
            if all_under_threshold
                myprintln(true, "All μ updates are under the threshold.")
            end
        end
    end

    if !all_under_threshold
        return ddpModel
    end

    if all_under_threshold
        println(crayon_green("All updates are under the threshold."))
        converged = true
        shouldStop = true
        @pack! ddpModel = converged, shouldStop
    end

    return ddpModel
end
# #endregion

#region forward_pass_1ph_L
"""
    forward_pass_1ph_L(ddpModel; verbose::Bool=false)

Perform a forward pass in the DDP algorithm.

This function performs a forward pass in the Differential Dynamic Programming (DDP) algorithm by iterating over the time steps and updating the model state accordingly.

# Arguments
- `ddpModel::Dict`: A dictionary containing the current state of the DDP model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: The updated dictionary after performing the forward pass.
"""
function forward_pass_1ph_L(ddpModel; verbose::Bool=false)
    verbose = true
    @unpack k_ddp = ddpModel
    # myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
    t_ddp = 1
    @unpack data = ddpModel
    @unpack Tset, T = data

    for t_ddp ∈ Tset # Tset is assumed sorted
        @pack! ddpModel = t_ddp
        modelDict_t0_k0 = MB.build_MPOPF_1ph_L_model_t_in_Tset(ddpModel, Tset=[t_ddp], verbose=verbose)
        ddpModel = reformulate_model_as_DDP_forward_step(ddpModel, modelDict_t0_k0, Tset=[t_ddp], verbose=verbose)
        ddpModel = solve_and_store_forward_step_t_in_Tset(ddpModel, Tset=[t_ddp], verbose=verbose)
    end

    ddpModel = compstore_PSubsCost(ddpModel, verbose=verbose) # Forward pass done, so compute PSubsCost

    return ddpModel
end
#endregion

function reformulate_model_as_DDP_forward_step(ddpModel, modelDict_t0_k0, Tset=nothing, verbose=verbose)
    if isnothing(Tset) || length(Tset) != 1
        @error "Tset seems invalid: $Tset"
        return
    end
    @unpack t_ddp = Tset[1]
    @unpack k_ddp, data, mu = ddpModel
    @unpack T = data
    
    μ = mu
    @unpack Bset, alpha_fpi, gamma_fpi = data;
    α_fpi0 = alpha_fpi
    γ_fpi = gamma_fpi
    MU = Dict()
    α_fpi = compute_alpha_fpi(α_fpi0, γ_fpi, k_ddp)
    myprintln(verbose, "FP$(k_ddp): α_fpi = $α_fpi")

    MU = compute_interpolated_mu(μ, Bset, k_ddp, t_ddp, α_fpi, verbose=verbose)

    model_t0_k0 = modelDict_t0_k0[:model]
    objfun_expr_t0_k0_base = objective_function(model_t0_k0)

    # Step 1: Modify objective function required for the current time-step if applicable

    if t_ddp ∈ 1:T-1
        objfun_expr_t0_k0_appendix = sum( MU[j, t_ddp+1] * (-model_t0_k0[:B][j, t_ddp]) for j ∈ Bset )
        objfun_expr_t0_k0_ddp = objfun_expr_t0_k0_base + objfun_expr_t0_k0_appendix
    elseif t_ddp == T
        objfun_expr_t0_k0_ddp = objfun_expr_t0_k0_base
    else
        @error "Invalid value of t_ddp: $t_ddp"
        return
    end

    @objective(model_t0_k0, Min, objfun_expr_t0_k0_ddp)

    # Step 2: Fix SOC variables for the previous time-step if applicable

    if t_ddp ∈ 2:T
        # Related Todo: Disambiguate between modelVals_ddp_vs_FP and modelVals_ddp_vs_t_vs_k
        @unpack modelVals_ddp_vs_t_vs_k = ddpModel;
        modelVals_tm1_k0 = modelVals_ddp_vs_t_vs_k[t_ddp-1, k_ddp]
        B_Vals_tm1_k0 = modelVals_tm1_k0[:B]

        for j ∈ Bset
            fix(model_t0_k0[:B][j, t_ddp - 1], B_Vals_tm1_k0[j, t_ddp - 1])
        end
    end

    # Step Last: Now just save the (unsolved) model in ddpModel
    @unpack models_ddp_vs_t_vs_k = ddpModel
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0_k0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel
end

#region compute_interpolated_mu
function compute_interpolated_mu(mu, Bset, k_ddp, t_ddp, α_fpi; verbose::Bool=false)
    """
    compute_interpolated_mu(mu, Bset, k_ddp, t_ddp, α_fpi; verbose::Bool=false)

    Compute the interpolated values of μ to be used for the current forward pass's latest time-step.

    # Arguments
    - `mu::Dict`: The dictionary containing μ values.
    - `Bset::Vector`: The set of battery indices.
    - `k_ddp::Int`: The current DDP iteration.
    - `t_ddp::Int`: The current time step.
    - `α_fpi::Float64`: The computed α_fpi value.
    - `verbose::Bool`: A flag to enable verbose output (default: false).

    # Returns
    - `MU::Dict`: A dictionary containing the interpolated μ values only for usage for the latest time-step.
    """
    MU = Dict()
    for j ∈ Bset
        if k_ddp == 1
            MU[j, t_ddp+1] = mu[j, t_ddp+1, k_ddp-1]
        elseif k_ddp >= 2
            MU[j, t_ddp+1] = get_interpolated_value(mu[j, t_ddp+1, k_ddp-1], mu[j, t_ddp+1, k_ddp-2], α_fpi)
            if j in Bset[1]
                MU_used_str = trim_number_for_printing(MU[j, t_ddp+1], sigdigits=2)
                MU_not_used_str = trim_number_for_printing(mu[j, t_ddp+1, k_ddp-1], sigdigits=2)
                myprintln(verbose, "FP$(k_ddp): μ[$(j), $(t_ddp+1)] = $(MU_used_str) instead of $(MU_not_used_str)")
            end
        else
            @error "Invalid value of k_ddp: $k_ddp"
            return Dict()
        end
    end
    return MU
end
#endregion

#region compstore_PSubsCost
"""
    compstore_PSubsCost(ddpModel; verbose::Bool=false)

Compute and store the cost of power borrowed from the substation across the horizon in an dict containing values from previous DDP forward passes.

# Arguments
- `ddpModel::Dict`: A dictionary containing the current state of the DDP model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: The updated dictionary with computed and stored power substitution costs.
"""
function compstore_PSubsCost(ddpModel; verbose::Bool=false)
    # Compute output values without mutating the original modelDict
    ddpModel_k0 = deepcopy(ddpModel)
    ddpModel_k0_with_outputVals = CO.compute_output_values(ddpModel_k0, verbose=verbose, forwardPass=true)
    @unpack outputVals_vs_k = ddpModel_k0
    outputVals_vs_k[k_ddp] = ddpModel_k0_with_outputVals[:data]
    PSubsCost_dollar_vs_k = [outputVals_vs_k[k][:PSubsCost_allT_dollar] for k in 1:k_ddp]
    @pack! ddpModel = outputVals_vs_k, PSubsCost_dollar_vs_k

    return ddpModel
end
#endregion

function ForwardStep_1ph_NL_t_is_1(ddpModel;
    verbose::Bool=false)

    @unpack t_ddp = ddpModel
    if t_ddp != 1
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end

    ddpModel = build_ForwardStep_1ph_NL_model_t_is_1(ddpModel)
    ddpModel = optimize_ForwardStep_1ph_NL_model_t_is_1(ddpModel)

    return ddpModel
end

function ForwardStep_1ph_NL_t_in_2toTm1(ddpModel;
    verbose::Bool=false)

    @unpack t_ddp, data = ddpModel;
    @unpack T = data;
    if !(2 <= t_ddp <= T-1)
        @error "t_ddp = $(t_ddp) is not in [2, T-1]"
        return
    end

    ddpModel = build_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel)
    ddpModel = optimize_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel)

    return ddpModel
end

function ForwardStep_1ph_NL_t_is_T(ddpModel;
    verbose::Bool=false)
    @unpack t_ddp, data = ddpModel
    @unpack T = data;
    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to T = $(T)"
        return
    end

    ddpModel = build_ForwardStep_1ph_NL_model_t_is_T(ddpModel)
    ddpModel = optimize_ForwardStep_1ph_NL_model_t_is_T(ddpModel)

    return ddpModel
end

function optimize_ForwardStep_1ph_NL_model_t_is_1(ddpModel;
    verbose::Bool=false)

    @unpack t_ddp = ddpModel
    if t_ddp != 1
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end

    @unpack k_ddp, models_ddp_vs_t_vs_k = ddpModel;

    # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp] # unsolved model

    optimize!(model_t0)

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0 # solved
    @pack! ddpModel = models_ddp_vs_t_vs_k

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_light_red("Printing the Forward Step Battery SOC values to be used for the next time-step"))
    # for j ∈ Bset
    #     println(crayon_light_red("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    # end

    Tset = [t_ddp]
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals = ddpModel

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_light_red("Printing modelVals[:B] (should be the same)"))
    # for j ∈ Bset
    #     println(crayon_light_red("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    # end

    @unpack modelVals_ddp_vs_t_vs_k = ddpModel;
    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]

    # Now that the model_t0 is solved and updated, we can compute the dual variables associated with its soc constraints for the next iteration's forward pass
    # println("Backward Pass for Tset = $Tset")
    ddpModel = compute_and_store_dual_variables(ddpModel, model_t0, Tset=Tset)

    return ddpModel
end

function optimize_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel;
    verbose::Bool=false)

    @unpack t_ddp, data = ddpModel;
    @unpack T = data;
    if !(2 <= t_ddp <= T-1)
        @error "t_ddp = $(t_ddp) is not in [2, T-1]"
        return
    end

    @unpack k_ddp = ddpModel

    # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    @unpack models_ddp_vs_t_vs_k = ddpModel;
    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp] # unsolved model
    optimize!(model_t0)
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0 # solved
    @pack! ddpModel = models_ddp_vs_t_vs_k

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_light_red("Printing the Forward Step Battery SOC values to be used by the next time-step"))
    # for j ∈ Bset
    #     println(crayon_light_red("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    # end

    Tset = [t_ddp]
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_light_red("Printing modelVals[:B] (should be the same)"))
    # for j ∈ Bset
    #     println(crayon_light_red("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    # end

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    # println("Backward Pass for Tset = $Tset")
    ddpModel = compute_and_store_dual_variables(ddpModel, model_t0, Tset=Tset)

    return ddpModel
end

function optimize_ForwardStep_1ph_NL_model_t_is_T(ddpModel;
    verbose::Bool=false)

    @unpack t_ddp, data = ddpModel;
    @unpack T = data;
    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end

    @unpack k_ddp = ddpModel

    # myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    @unpack models_ddp_vs_t_vs_k = ddpModel;
    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp] # unsolved model
    optimize!(model_t0)
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0 # solved
    @pack! ddpModel = models_ddp_vs_t_vs_k

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        # println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    # println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Best objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_blue = Crayon(foreground=:white, background=:blue, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_blue("Printing the Terminal Battery SOC values"))
    # for j ∈ Bset
    #     println(crayon_blue("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    # end

    Tset = [t_ddp]
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel

    crayon_blue = Crayon(foreground=:white, background=:blue, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    # println(crayon_blue("Printing terminal SOC values in modelVals[:B] (should be the same)"))
    # for j ∈ Bset
    #     println(crayon_blue("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    # end

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    # println("Backward Pass for Tset = $Tset")
    ddpModel = compute_and_store_dual_variables(ddpModel, model_t0, Tset=Tset)

    return ddpModel
end

#region optimize_MPOPF_1ph_L_DDP
"""
    optimize_MPOPF_1ph_L_DDP(data; verbose::Bool=false)

Optimize the Multi-Period Optimal Power Flow (MPOPF) model for a single-phase network represented as the LinDistFlow model using DDP.

This function optimizes the MPOPF model using the Differential Dynamic Programming (DDP) algorithm. 
It initializes the DDP model, performs forward passes, checks for convergence, and retrieves the results.

# Arguments
- `data::Dict`: A dictionary containing all necessary data and parameters for the model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: A dictionary containing the optimized DDP model and its parameters.
"""
function optimize_MPOPF_1ph_L_DDP(data;
    verbose::Bool=false,
    muDict=nothing,
    maxiter::Int=7)

    ddpModel = DDPModel(data, maxiter=maxiter, muDict=muDict)

    keepForwardPassesRunning = true
    while keepForwardPassesRunning
        @unpack k_ddp = ddpModel
        # myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
        ddpModel = forward_pass_1ph_L(ddpModel,
            verbose=verbose)

        dppModel = check_for_ddp_convergence(ddpModel)

        @unpack shouldStop = ddpModel
        keepForwardPassesRunning = !shouldStop

        k_ddp += 1
        @pack! ddpModel = k_ddp

    end

    @unpack modelVals = ddpModel
    modelVals[:solve_time] = sum(modelVals[:solve_time_vs_t][t] for t ∈ data[:Tset])
    @pack! ddpModel = modelVals
    # Check solver status and retrieve results
    @unpack iterLimitReached, converged, modelVals = ddpModel

    # Define crayons for green and red text
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)
    crayon_final_green = Crayon(foreground=:white, background=:green, bold=true)
    if iterLimitReached
        println(crayon_red("Maximum iterations reached!"))
    elseif converged
        println(crayon_final_green("Optimization via DDP converged."))
    else
        @error "floc"
    end

    return ddpModel
end
#endregion

#region DDPModel
"""
    DDPModel(data; maxiter::Int=10, verbose::Bool=false)

Initialize the DDP model for a given dataset.

This function initializes the DDP (Differential Dynamic Programming) model by setting up the necessary data structures and parameters. 
It handles the initialization of dual variables, model values, and other relevant parameters for the DDP algorithm.

# Arguments
- `data::Dict`: A dictionary containing all necessary data and parameters for the model.
- `maxiter::Int`: The maximum number of iterations for the DDP algorithm (default: 10).
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: A dictionary containing the initialized DDP model and its parameters.
"""
function DDPModel(data;
    maxiter::Int=7,
    muDict=nothing,
    verbose::Bool=false)

    @unpack Tset, Bset, solver = data;
    models_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Model}()
    modelVals_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Dict}()
    mu = Dict{Tuple{Int,Int,Int},Float64}()
    lambda_lo = Dict{Tuple{Int,Int,Int},Float64}()
    lambda_up = Dict{Tuple{Int,Int,Int},Float64}()
    outputVals_vs_k = Dict{Int, Any}()
    # modelVals = Dict{Symbol,Any}()
    modelVals = MC.ModelVals(data)
    # Initialize mu[j, t_ddp, 0/-1] = 0 for all j in Bset and t_ddp in Tset
    
    for j ∈ Bset, t_ddp ∈ Tset
        if isnothing(muDict)
            mu[j, t_ddp, 0] = 0.0
        else
            mu[j, t_ddp, 0] = muDict[j, t_ddp]
        end
        lambda_lo[j, t_ddp, 0] = 0.0
        lambda_up[j, t_ddp, 0] = 0.0
    end

    ddpModel = Dict(
        :converged => false,
        :data => data,
        :iterLimitReached => false,
        :fval_best => Inf,
        :fval_best_hit => 0,
        :k_best => -1,
        :k_ddp => 1,
        :lambda_lo => lambda_lo,
        :lambda_up => lambda_up,
        :maxiter => maxiter,
        :modelVals => modelVals,
        :models_ddp_vs_t_vs_k=>models_ddp_vs_t_vs_k,
        :modelVals_ddp_vs_t_vs_k=>modelVals_ddp_vs_t_vs_k,
        :mu=>mu,
        :outputVals_vs_k => outputVals_vs_k,
        :shouldStop => false,
        :t_ddp=>0
    )
    return ddpModel
end
#endregion

end