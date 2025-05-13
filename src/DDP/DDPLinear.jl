module DDPLinear

export
    check_for_ddp_convergence,
    check_solver_status,
    compute_interpolated_mu,
    compstore_PSubsCost,
    DDPModel,
    forward_pass_1ph_L,
    get_interpolated_value,
    optimize_MPOPF_1ph_L_DDP,
    reformulate_model_as_FS,
    solve_and_store_FS,
    store_FP_k_decvar_values,
    store_FS_t_k_dual_variables

include("../computeOutputs.jl")
import .computeOutputs as CO

include("../functionRetriever.jl")
import .functionRetriever as FR

include("../ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("../ModelCopier/ModelCopier.jl")
import .ModelCopier as MC

include("../helperFunctions.jl")
import .helperFunctions as HF

include("../SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger

include("../Exporter.jl")
import .Exporter as Exporter

using Crayons
using JuMP
using Gurobi
using Ipopt
using Parameters: @unpack, @pack!

#region compute_alpha_fpi
function compute_alpha_fpi(alpha_fpi0, gamma_fpi, k_ddp)
    return alpha_fpi0 * gamma_fpi^(k_ddp - 2)
end
#endregion

#region get_interpolated_value  
function get_interpolated_value(xn, xnm1, alpha)
    diff = xn - xnm1
    return (alpha*xnm1 + xn) / (1 + alpha)
end
#endregion

#region check_for_ddp_convergence
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
                mu_current = HF.trim_number_for_printing(mu[(j, t, k_ddp)], sigdigits=2)

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
                soc_percent_j_t_current_str = HF.trim_number_for_printing(soc_percent_j_t_current, sigdigits=2)

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
    threshold_soc = 1e-2
    threshold_mu = 1e-2
    threshold_fval = 1e-1
    # threshold_conv_iters = 2
    all_under_threshold = true

    crayon_red_neg = Crayon(foreground=:red, bold=true, negative=true)
    crayon_green = Crayon(background=:green, foreground=:white, bold=true)
    crayon_blue = Crayon(foreground=:blue, bold=true)
    crayon_blue_neg = Crayon(foreground=:blue, bold=true, negative=true)

    # Criterion 2: Check the magnitude of updates in fval

    converged_fval = true
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
            HF.myprintln(print_fval, "*******")
            HF.myprintln(print_fval, crayon_red_neg("Previous value of PSubsCost_allT_dollar = $PSubsCost_previous"))
            HF.myprintln(print_fval, crayon_red_neg("Current value of PSubsCost_allT_dollar = $PSubsCost_current"))
            HF.myprintln(true, "*******")
            converged_fval = false
        else
            HF.myprintln(true, "FP$(k_ddp): Objective function value update is under the threshold.")
        end
    end

    if converged_fval
        ddpModel[:conv_fval_num] += 1
        HF.myprintln(print_fval, "FP$(k_ddp): conv_fval_num = $(ddpModel[:conv_fval_num]).")
    else
        ddpModel[:conv_fval_num] = 1
    end

    
    # println(crayon_green("Checking convergence for B values:"))

    # Criterion 3: Check the magnitude of updates in SOC values

    converged_soc = true
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
                    converged_soc = false
                    all_under_threshold = false
                    if discrepancy > threshold_soc
                        # HF.myprintln(verbose, "Exceeding update tolerance: var_name = $var_name, discrepancy = $discrepancy")
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
    end

    if converged_soc
        ddpModel[:conv_soc_num] += 1
        HF.myprintln(print_soc, "FP$(k_ddp): conv_soc_num = $(ddpModel[:conv_soc_num]).")
    else
        ddpModel[:conv_soc_num] = 1
    end

    if all_under_threshold
        HF.myprintln(true, "FP$(k_ddp): All SOC updates are under the threshold.")
    end

    # Criterion 4: Check the magnitude of updates in μ values

    converged_mu = true
    if mu_change

        for t in Tset
            for j in Bset[1]
                if haskey(mu, (j, t, k_ddp)) && haskey(mu, (j, t, k_ddp - 1))
                    mu_current = mu[(j, t, k_ddp)]
                    mu_previous = mu[(j, t, k_ddp - 1)]
                    discrepancy = abs(mu_current - mu_previous)
                    max_discrepancy = max(max_discrepancy, discrepancy)
                    if discrepancy > threshold_mu
                        converged_mu = false
                        all_under_threshold = false
                        # HF.myprintln(verbose, "Exceeding update tolerance: mu[$j, $t, $k_ddp], discrepancy = $discrepancy")
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
        end
    end
    if all_under_threshold
        HF.myprintln(true, "FP$(k_ddp): All μ updates are under the threshold.")
    end

    if converged_mu
        ddpModel[:conv_mu_num] += 1
        HF.myprintln(print_mu, "FP$(k_ddp): conv_mu_num = $(ddpModel[:conv_mu_num]).")
    else
        ddpModel[:conv_mu_num] = 1
    end

    if ddpModel[:conv_fval_num] >= ddpModel[:threshold_conv_iters] && ddpModel[:conv_soc_num] >= ddpModel[:threshold_conv_iters]
        println(crayon_green("FP$(k_ddp): SOC and fval convergence criteria repeatedly satisfied for $(ddpModel[:threshold_conv_iters]) Forward Passes."))
        converged = true
        shouldStop = true
        @pack! ddpModel = converged, shouldStop
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
#endregion

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
    # HF.myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
    t_ddp = 1
    @unpack data = ddpModel
    @unpack Tset, T = data

    for t_ddp ∈ Tset # Tset is assumed sorted
        @pack! ddpModel = t_ddp
        modelDict_t0_k0 = MB.build_MPOPF_1ph_L_model_t_in_Tset(data, Tset=[t_ddp], verbose=verbose)
        ddpModel = reformulate_model_as_FS(ddpModel, modelDict_t0_k0, Tset=[t_ddp], verbose=verbose)
        ddpModel = warm_start_FS(ddpModel, Tset=[t_ddp], verbose=verbose)
        ddpModel = solve_and_store_FS(ddpModel, Tset=[t_ddp], verbose=verbose)
    end

    ddpModel = MC.store_FP_k_decvar_values(ddpModel, verbose=verbose) # Store the decision variable values for the current forward pass now that all the FS models have been solved
    ddpModel = compstore_PSubsCost(ddpModel, verbose=verbose) # Forward pass done, so compute PSubsCost

    return ddpModel
end
#endregion

"""
    warm_start_FS(ddpModel; Tset, verbose::Bool=false)

Warm start the variables for the current forward step using the values from the last time step.

# Arguments
- `ddpModel::Dict`: The current state of the DDP model.
- `Tset::Vector`: A single-element array containing the current time-step representing the Forward Step.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: The updated DDP model with warm-started variables.
"""
function warm_start_FS(ddpModel; Tset, verbose::Bool=false)
    if isnothing(Tset) || length(Tset) != 1
        @error "Tset seems invalid: $Tset"
        return ddpModel
    end

    t_ddp = Tset[1]
    @unpack k_ddp, models_ddp_vs_t_vs_k, modelVals_ddp_vs_t_vs_k, data = ddpModel

    # Ensure we are not at the first time step
    if t_ddp == 1
        return ddpModel
    end

    # Get the current model and the variable values from the previous time step
    model_current = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
    model_previous_vals = modelVals_ddp_vs_t_vs_k[t_ddp - 1, k_ddp]

    # Warm start substation power flow variable
    P_Subs = model_current[:P_Subs]
    if haskey(model_previous_vals, :P_Subs)
        set_start_value(P_Subs[t_ddp], model_previous_vals[:P_Subs][t_ddp - 1])
    end

    # Warm start power flow variables
    @unpack Lset = data
    P = model_current[:P]
    Q = model_current[:Q]
    for (i, j) in Lset
        if haskey(model_previous_vals[:P], (i, j)) && haskey(model_previous_vals[:Q], (i, j))
            set_start_value(P[(i, j), t_ddp], model_previous_vals[:P][(i, j), t_ddp - 1])
            set_start_value(Q[(i, j), t_ddp], model_previous_vals[:Q][(i, j), t_ddp - 1])
        end
    end

    # Warm start voltage variables
    @unpack Nset = data
    v = model_current[:v]
    for j in Nset
        if haskey(model_previous_vals[:v], j)
            set_start_value(v[j, t_ddp], model_previous_vals[:v][j, t_ddp - 1])
        end
    end

    # Warm start PV inverter reactive dispatch variables
    @unpack Dset = data
    q_D = model_current[:q_D]
    for j in Dset
        if haskey(model_previous_vals[:q_D], j)
            set_start_value(q_D[j, t_ddp], model_previous_vals[:q_D][j, t_ddp - 1])
        end
    end

    # Warm start battery real and reactive dispatch variables
    @unpack Bset = data
    q_B = model_current[:q_B]
    P_c = model_current[:P_c]
    P_d = model_current[:P_d]
    for j in Bset
        if haskey(model_previous_vals[:q_B], j)
            set_start_value(q_B[j, t_ddp], model_previous_vals[:q_B][j, t_ddp - 1])
        end
        if haskey(model_previous_vals[:P_c], j)
            set_start_value(P_c[j, t_ddp], model_previous_vals[:P_c][j, t_ddp - 1])
        end
        if haskey(model_previous_vals[:P_d], j)
            set_start_value(P_d[j, t_ddp], model_previous_vals[:P_d][j, t_ddp - 1])
        end
    end

    # Warm start battery state of charge (SOC) variables
    B = model_current[:B]
    for j in Bset
        if haskey(model_previous_vals[:B], j)
            set_start_value(B[j, t_ddp], model_previous_vals[:B][j, t_ddp - 1])
        end
    end

    return ddpModel
end

#region reformulate_model_as_FS
function reformulate_model_as_FS(ddpModel, modelDict_t0_k0; 
    Tset,
    verbose=false)
    if isnothing(Tset) || length(Tset) != 1
        @error "Tset seems invalid: $Tset"
        return
    end
    t_ddp = Tset[1]
    @unpack k_ddp, data, mu = ddpModel
    @unpack T = data

    model_t0_k0 = modelDict_t0_k0[:model] # model_t0_k0 is the actual FP model for t = t_ddp, and will be directly mutated (in objective function and SOC trajectory equations) before being loaded back to models_ddp_vs_t_vs_k in ddpModel for solving

    objfun_expr_t0_k0_base = objective_function(model_t0_k0)

    # Step 1: Modify objective function required for the current time-step if applicable

    if t_ddp ∈ 1:T-1
        μ = mu
        @unpack Bset, alpha_fpi, gamma_fpi = data
        α_fpi0 = alpha_fpi
        γ_fpi = gamma_fpi
        MU = Dict()
        α_fpi = compute_alpha_fpi(α_fpi0, γ_fpi, k_ddp)
        # HF.myprintln(verbose, "FP$(k_ddp): α_fpi = $α_fpi")
        MU = compute_interpolated_mu(μ, Bset, k_ddp, t_ddp, α_fpi, verbose=verbose)

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
    @unpack Bset, algo_temporal_decmp = data;

    if algo_temporal_decmp == "tENApp"
        if t_ddp ∈ 2:T
            @unpack modelVals_ddp_vs_t_vs_k = ddpModel;
            if k_ddp >= 2
                modelVals_tm1_k0m1 = modelVals_ddp_vs_t_vs_k[t_ddp-1, k_ddp-1]
                B_Vals_tm1_k0m1 = modelVals_tm1_k0m1[:B]

                for j ∈ Bset
                    fix(model_t0_k0[:B][j, t_ddp-1], B_Vals_tm1_k0m1[j, t_ddp-1])
                end

            elseif k_ddp == 1

                @unpack B0_pu = data
                for j ∈ Bset
                    fix(model_t0_k0[:B][j, t_ddp - 1], B0_pu[j])
                end
            end
        end
    elseif algo_temporal_decmp == "DDP"
        if t_ddp ∈ 2:T
            @unpack modelVals_ddp_vs_t_vs_k = ddpModel;
            modelVals_tm1_k0 = modelVals_ddp_vs_t_vs_k[t_ddp-1, k_ddp]
            B_Vals_tm1_k0 = modelVals_tm1_k0[:B]

            for j ∈ Bset
                fix(model_t0_k0[:B][j, t_ddp - 1], B_Vals_tm1_k0[j, t_ddp - 1])
            end
        end
    end

    # Step Last: Now just save the (unsolved) model in ddpModel
    @unpack models_ddp_vs_t_vs_k = ddpModel
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0_k0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel
end
#endregion

#region solve_and_store_FS
function solve_and_store_FS(ddpModel; Tset,
    verbose=false)
    if isnothing(Tset) || length(Tset) != 1
        @error "Tset seems invalid: $Tset"
        return
    end
    t_ddp = Tset[1]

    @unpack models_ddp_vs_t_vs_k, k_ddp = ddpModel
    model_t0_k0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp] # unsolved FS model
    optimize!(model_t0_k0)
    solverStatus = check_solver_status(model_t0_k0, k_ddp=k_ddp, t_ddp=t_ddp)
    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0_k0 # solved FS model
    @pack! ddpModel = models_ddp_vs_t_vs_k

    ddpModel = MC.store_FS_t_k_decvar_values(ddpModel, Tset=Tset, optModel="Linear", verbose=verbose)    
    
    # Todo: WTF it this function?
    ddpModel = store_FS_t_k_dual_variables(ddpModel, Tset=Tset, verbose=verbose)
end
#endregion

#region check_solver_status
function check_solver_status(model; k_ddp=nothing, t_ddp=nothing)
    # Define crayons for colored output
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    # Check the solver termination status
    if termination_status(model) == LOCALLY_SOLVED
        return true
    else
        if !isnothing(k_ddp) || !isnothing(t_ddp)
            println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
        end
        println(crayon_red("Solver status: $(termination_status(model))"))
        return false
    end
end
#endregion

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
                MU_used_str = HF.trim_number_for_printing(MU[j, t_ddp+1], sigdigits=2)
                MU_not_used_str = HF.trim_number_for_printing(mu[j, t_ddp+1, k_ddp-1], sigdigits=2)
                # HF.myprintln(verbose, "FP$(k_ddp): μ[$(j), $(t_ddp+1)] = $(MU_used_str) instead of $(MU_not_used_str)")
            end
        else
            @error "Invalid value of k_ddp: $k_ddp"
            return Dict()
        end
    end
    return MU
end
#endregion

#region store_FS_t_k_dual_variables
"""
    store_FS_t_k_dual_variables(ddpModel; Tset, verbose::Bool=false)

    Store the dual variable values for the given time step and iteration in the DDP model.

    # Arguments
    - `ddpModel::Dict`: The current state of the DDP model.
    - `Tset::Vector`: 1-element array containing the current time-step representing the Forward Step.
    - `verbose::Bool`: A flag to enable verbose output (default: false).

    # Returns
    - `ddpModel::Dict`: The updated DDP model with stored dual variable values.
"""
function store_FS_t_k_dual_variables(ddpModel; Tset, verbose::Bool=false)
    if isnothing(Tset) || length(Tset) != 1
        @error "Tset seems invalid: $Tset"
        return
    end

    @unpack k_ddp, data, mu, lambda_lo, lambda_up = ddpModel
    μ = mu
    λ_lo = lambda_lo
    λ_up = lambda_up
    @unpack T, Bset = data
    t_ddp = Tset[1]

    @unpack models_ddp_vs_t_vs_k = ddpModel
    model_t_k = models_ddp_vs_t_vs_k[t_ddp, k_ddp]  # Get the solved model for this time step and iteration

    # Construct partial prefixes once per t_ddp
    if t_ddp == 1
        soc_traj_t0_k0_str_prefix = "h_SOC_j^{t=1}_Initial_SOC_Node_j_"
        soc_traj_suffix = "t1"
    elseif 2 <= t_ddp <= T
        soc_traj_t0_k0_str_prefix = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_"
        soc_traj_suffix = "t_$(t_ddp)"
    else
        @error "Invalid value of t_ddp: $t_ddp"
        return
    end

    soc_lim_lo_t0_k0_str_prefix = "g_11_j^t_MinSOC_Node_j_"
    soc_lim_up_t0_k0_str_prefix = "g_12_j^t_MaxSOC_Node_j_"

    for j in Bset
        # Build full constraint names using precomputed prefixes
        soc_traj_t0_k0_j_str = soc_traj_t0_k0_str_prefix * string(j) * "_" * soc_traj_suffix
        soc_lim_lo_t0_k0_j_str = soc_lim_lo_t0_k0_str_prefix * string(j) * "_t_$(t_ddp)"
        soc_lim_up_t0_k0_j_str = soc_lim_up_t0_k0_str_prefix * string(j) * "_t_$(t_ddp)"

        # Fetch and store duals
        soc_traj_t0_k0_j = constraint_by_name(model_t_k, soc_traj_t0_k0_j_str)
        μ[j, t_ddp, k_ddp] = -dual(soc_traj_t0_k0_j)

        soc_lim_lo_t0_k0_j = constraint_by_name(model_t_k, soc_lim_lo_t0_k0_j_str)
        soc_lim_up_t0_k0_j = constraint_by_name(model_t_k, soc_lim_up_t0_k0_j_str)
        λ_lo[j, t_ddp, k_ddp] = -dual(soc_lim_lo_t0_k0_j)
        λ_up[j, t_ddp, k_ddp] = -dual(soc_lim_up_t0_k0_j)
    end

    # Update the ddpModel with the modified dual variables
    mu = μ
    lambda_lo = λ_lo
    lambda_up = λ_up
    @pack! ddpModel = mu, lambda_lo, lambda_up

    return ddpModel
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
    @unpack outputVals_vs_k, k_ddp = ddpModel_k0
    outputVals_vs_k[k_ddp] = ddpModel_k0_with_outputVals[:data]
    PSubsCost_dollar_vs_k = [outputVals_vs_k[k][:PSubsCost_allT_dollar] for k in 1:k_ddp]
    @pack! ddpModel = outputVals_vs_k, PSubsCost_dollar_vs_k

    return ddpModel
end
#endregion

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
        # HF.myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
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
        println(crayon_final_green("Optimization via $(data[:algo_temporal_decmp]) converged."))
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
    modelVals_ddp_vs_FP = Dict{Int,Dict}()
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
        :conv_fval_num => 1,
        :conv_soc_num => 1,
        :conv_mu_num => 1,
        :threshold_conv_iters => get(data, :threshold_conv_iters, 3),

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
        :modelVals_ddp_vs_FP=>modelVals_ddp_vs_FP,
        :mu=>mu,
        :outputVals_vs_k => outputVals_vs_k,
        :shouldStop => false,
        :t_ddp=>0
    )
    return ddpModel
end
#endregion

end # module DDPLinear