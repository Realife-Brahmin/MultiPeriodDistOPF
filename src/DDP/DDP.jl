module DDP 

export
    backward_pass, 
    build_ForwardStep_1ph_NL_model_t_is_1,
    build_ForwardStep_1ph_NL_model_t_in_2toTm1,
    build_ForwardStep_1ph_NL_model_t_is_T,
    DDPModel,
    check_for_ddp_convergence,
    forward_pass,
    ForwardStep_1ph_NL_t_is_1,
    ForwardStep_1ph_NL_t_in_2toTm1,
    ForwardStep_1ph_NL_t_is_T,
    optimize_ForwardStep_1ph_NL_model_t_is_1,
    optimize_MPOPF_1ph_NL_DDP


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

function backward_pass(ddpModel, model_t0;
    Tset=nothing)
    if Tset === nothing
        println("No Tset defined for backward pass, so using the Tset from the data")
        @unpack data = ddpModel;
        Tset = data[:Tset]
    end
    @unpack k_ddp, data, mu = ddpModel;
    μ = mu;
    @unpack T, Bset = data;
    # Update mu values post optimization
    t_ddp = Tset[1]
    # for t_ddp ∈ Tset
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
        μ[j, t_ddp, k_ddp] = dual(constraint_j_t0)
    end

    crayon_update = Crayon(foreground=:light_blue, background=:white, bold=true)
    println(crayon_update("Backward Pass k_ddp = $(k_ddp): μ values for t_ddp = $(t_ddp)"))
    for j ∈ Bset
        # println("j = $j: μ = ", trim_number_for_printing(μ[j, t_ddp, k_ddp], sigdigits=2))
        println(crayon_update("j = $j: μ = ", trim_number_for_printing(μ[j, t_ddp, k_ddp], sigdigits=2)))
    end
    # end
    mu = μ
    @pack! ddpModel = mu

    return ddpModel
end

function build_ForwardStep_1ph_NL_model_t_is_1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, modelVals, models_ddp_vs_t_vs_k, data, mu = ddpModel;

    if t_ddp != 1
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end
    verbose = true

    # if k_ddp == 1
    if k_ddp >= 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [1]
        modelDict = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0) # an unsolved model
        model_t0 = modelDict[:model]
    # elseif k_ddp >= 2
    #     myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
    #     model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1] # I guess a solved model
    #     model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
    #     @unpack data = ddpModel;
    #     @unpack solver = data;
    #     model_t0 = SolverArranger.attach_solver(model_t0, solver)
    #     # @show model_t0
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0) # same as model_t0_km1 as it is a copy
    μ = mu
    @unpack Bset = data;
    # objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[j, t_ddp+1, k_ddp-1] - μ[j, t_ddp+1, k_ddp-2] ) * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum(μ[j, t_ddp+1, k_ddp-1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset)

    @objective(model_t0, Min, objfun_expr_t0_k)

    # B0 values are already set in the model so no need to fix them separately

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel

end

function build_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, modelVals, models_ddp_vs_t_vs_k, data, mu = ddpModel
    @unpack T = data;
    if !(2 <= t_ddp <= T-1)
        @error "t_ddp = $(t_ddp) is not in [2, T-1]"
        return
    end

    verbose = true
    # if k_ddp == 1
    if k_ddp >= 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be something like [2] or [3] or ... or [T-1]
        modelDict_t0 = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    # elseif k_ddp >= 2
    #     myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
    #     model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1]
    #     # model_t0 = deepcopy(model_t0_km1)
    #     model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
    #     @unpack data = ddpModel
    #     @unpack solver = data
    #     model_t0 = SolverArranger.attach_solver(model_t0, solver)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass) and previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = modelVals[:B]
    @unpack Bset = data;

    crayon_light_red1 = Crayon(background=:light_red, foreground=:white, bold=true)
    println(crayon_light_red1("Printing the previous time-step's SOC values to be used for the current time-step"))
    for j ∈ Bset
        println(crayon_light_red1("B_model[$j, $(t_ddp - 1)] = $(B_model[j, t_ddp - 1])"))
        fix(model_t0[:B][j, t_ddp - 1], B_model[j, t_ddp - 1])
    end

    # Update the model with the future dual variables from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0)
    μ = mu
    # objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[j, t_ddp+1, k_ddp-1] - μ[j, t_ddp+1, k_ddp-2] ) * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum( μ[j, t_ddp+1, k_ddp-1] * (-model_t0[:B][j, t_ddp]) for j ∈ Bset)

    @objective(model_t0, Min, objfun_expr_t0_k)

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel

end

function build_ForwardStep_1ph_NL_model_t_is_T(ddpModel;
    verbose::Bool=false)
    
    @unpack k_ddp, t_ddp, modelVals, models_ddp_vs_t_vs_k, data = ddpModel;
    @unpack T = data;
    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to T = $(T)"
        return
    end

    # if k_ddp == 1
    verbose = true
    if k_ddp >= 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [T]
        modelDict_t0 = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    # elseif k_ddp >= 2
    #     myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
    #     model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1]
    #     # model_t0 = deepcopy(model_t0_km1)
    #     model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
    #     @unpack data = ddpModel
    #     @unpack solver = data
    #     model_t0 = SolverArranger.attach_solver(model_t0, solver)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = modelVals[:B]
    @unpack Bset = data;
    # for j ∈ Bset
    #     fix(model_t0[:B][j, t_ddp - 1], B_model[j, t_ddp - 1])
    # end 

    crayon_light_red1 = Crayon(background=:light_red, foreground=:white, bold=true)
    println(crayon_light_red1("Printing the previous time-step's SOC values to be used for the current time-step"))
    for j ∈ Bset
        println(crayon_light_red1("B_model[$j, $(t_ddp - 1)] = $(B_model[j, t_ddp - 1])"))
        fix(model_t0[:B][j, t_ddp-1], B_model[j, t_ddp-1])
    end

    # No need of updating objective function, as no μ term is present in the objective function for the terminal time-step

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    return ddpModel
end

#region check_for_ddp_convergence
"""
    check_for_ddp_convergence(ddpModel; verbose::Bool=false)

Check for convergence of the DDP algorithm.

This function checks if the Differential Dynamic Programming (DDP) algorithm has converged by evaluating the iteration count and the magnitude of updates in the model decision variable values.

# Arguments
- `ddpModel::Dict`: A dictionary containing the current state of the DDP model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: The updated dictionary with convergence status and stopping criteria.
"""
function check_for_ddp_convergence(ddpModel; verbose::Bool=false)
    @unpack k_ddp, maxiter, models_ddp_vs_t_vs_k, data = ddpModel;
    @unpack Tset = data;

    verbose = true

    # Criterion 1: Check if k_ddp has crossed the maxiter threshold
    if k_ddp >= maxiter
        println("Maximum iterations reached: $k_ddp")
        iterLimitReached = true
        shouldStop = true
        @pack! ddpModel = iterLimitReached, shouldStop
        return ddpModel
    end

    # @show k_ddp
    if k_ddp == 1
        # println("No updates to check for k_ddp = 1")
        return ddpModel
    end

    # Criterion 2: Check the magnitude of updates in the model decision variable values
    
    # Compare variable values and compute discrepancies
    max_discrepancy = 0.0
    threshold = 1e-5
    all_under_threshold = true

    for t_ddp in Tset
        # @show models_ddp_vs_t_vs_k
        model_current = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
        model_previous = models_ddp_vs_t_vs_k[t_ddp, k_ddp - 1]

        # Create dictionaries for current and previous models
        var_dict_current = MC.create_variable_dict(model_current)
        var_dict_previous = MC.create_variable_dict(model_previous)

        for var_name in keys(var_dict_current)
            if haskey(var_dict_previous, var_name)
                value_current = var_dict_current[var_name]
                value_previous = var_dict_previous[var_name]
                discrepancy = abs(value_current - value_previous)
                max_discrepancy = max(max_discrepancy, discrepancy)

                crayon_red_neg = Crayon(foreground=:red, bold=true, negative=true)
                if max_discrepancy > threshold
                    all_under_threshold = false
                    myprintln(verbose, "Some updates exceed the threshold. So keep doing Forward Passes.")
                    myprintln(verbose, "Case in point: var_name = $var_name, discrepancy = $discrepancy")
                    println(crayon_red_neg("Previous value of var $(var_name) = $value_previous"))
                    println(crayon_red_neg("Current value of var $(var_name) = $value_current"))
                    # return ddpModel
                end
            end
        end
        if !all_under_threshold
            return ddpModel
        end
    end

    myprintln(verbose, "All updates are under the threshold.")
    converged = true
    shouldStop = true
    @pack! ddpModel = converged, shouldStop
    return ddpModel
    
end
#endregion

#region forward_pass
"""
    forward_pass(ddpModel; verbose::Bool=false)

Perform a forward pass in the DDP algorithm.

This function performs a forward pass in the Differential Dynamic Programming (DDP) algorithm by iterating over the time steps and updating the model state accordingly.

# Arguments
- `ddpModel::Dict`: A dictionary containing the current state of the DDP model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: The updated dictionary after performing the forward pass.

# Steps
1. **Unpack Data**: Unpacks necessary data from the DDP model dictionary.
2. **Initialize Time Step**: Initializes the time step for the forward pass.
3. **Iterate Over Time Steps**: Iterates over the sorted time steps in `Tset`.
4. **Forward Step Execution**: Executes the appropriate forward step function based on the current time step:
    - `ForwardStep_1ph_NL_t_is_1` for the first time step.
    - `ForwardStep_1ph_NL_t_in_2toTm1` for intermediate time steps.
    - `ForwardStep_1ph_NL_t_is_T` for the last time step.
5. **Export Model**: Exports the optimization model after each forward step.
6. **Return Data**: Returns the updated dictionary after performing the forward pass.
"""
function forward_pass(ddpModel;
    verbose::Bool=false)
    verbose = true
    @unpack k_ddp = ddpModel;
    myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
    t_ddp = 1
    @unpack data = ddpModel;
    @unpack Tset, T = data;
    for t_ddp ∈ Tset # Tset is assumed sorted
        @pack! ddpModel = t_ddp
        if t_ddp == 1
            ddpModel = ForwardStep_1ph_NL_t_is_1(ddpModel, verbose=verbose)
        elseif 2 <= t_ddp <= T-1
            ddpModel = ForwardStep_1ph_NL_t_in_2toTm1(ddpModel, verbose=verbose)
        elseif t_ddp == T
            ddpModel = ForwardStep_1ph_NL_t_is_T(ddpModel, verbose=verbose)
        else
            @error "Invalid value of t_ddp: $t_ddp"
            return
        end
        Exporter.export_optimization_model(ddpModel, verbose=verbose)
    end

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

    myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
    # @show get_attribute(model_t0, MOI.Silent())
    optimize!(model_t0)

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_light_red("Printing the Forward Step Battery SOC values to be used for the next time-step"))
    for j ∈ Bset
        println(crayon_light_red("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    end

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k;

    Tset = [t_ddp]
    # @unpack modelDict = ddpModel;
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_light_red("Printing modelVals[:B] (should be the same)"))
    for j ∈ Bset
        println(crayon_light_red("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    end

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    ddpModel = backward_pass(ddpModel, model_t0, Tset=Tset)

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

    @unpack k_ddp, models_ddp_vs_t_vs_k = ddpModel

    myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
    # set_silent(model_t0)
    optimize!(model_t0)

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_light_red("Printing the Forward Step Battery SOC values to be used for the next time-step"))
    for j ∈ Bset
        println(crayon_light_red("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    end

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    Tset = [t_ddp]
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel

    crayon_light_red = Crayon(foreground=:light_red, background=:white, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_light_red("Printing modelVals[:B] (should be the same)"))
    for j ∈ Bset
        println(crayon_light_red("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    end

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    ddpModel = backward_pass(ddpModel, model_t0, Tset=Tset)

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

    @unpack k_ddp, models_ddp_vs_t_vs_k = ddpModel

    myprintln(verbose, "Forward Pass k_ddp = $(k_ddp) : About to optimize Forward Step model for t = $(t_ddp)")

    model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
    # set_silent(model_t0)
    optimize!(model_t0)

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Best objective function value for t = $(t_ddp): $optimal_obj_value"))

    crayon_blue = Crayon(foreground=:white, background=:blue, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_blue("Printing the Terminal Battery SOC values"))
    for j ∈ Bset
        println(crayon_blue("B[$j, $t_ddp] =  $(value(model_t0[:B][j, t_ddp]))"))
    end

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    Tset = [t_ddp]
    ddpModel = MC.copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel

    crayon_blue = Crayon(foreground=:white, background=:blue, bold=true)
    @unpack data = ddpModel
    @unpack Bset = data
    println(crayon_blue("Printing terminal SOC values in modelVals[:B] (should be the same)"))
    for j ∈ Bset
        println(crayon_blue("modelVals[:B][$j, $t_ddp] = $(modelVals[:B][j, t_ddp])"))
    end

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    ddpModel = backward_pass(ddpModel, model_t0, Tset=Tset)

    return ddpModel
end

#region optimize_MPOPF_1ph_NL_DDP
"""
    optimize_MPOPF_1ph_NL_DDP(data; verbose::Bool=false)

Optimize the Multi-Period Optimal Power Flow (MPOPF) model for a single-phase network with nonlinear loads using DDP.

This function optimizes the MPOPF model using the Differential Dynamic Programming (DDP) algorithm. 
It initializes the DDP model, performs forward passes, checks for convergence, and retrieves the results.

# Arguments
- `data::Dict`: A dictionary containing all necessary data and parameters for the model.
- `verbose::Bool`: A flag to enable verbose output (default: false).

# Returns
- `ddpModel::Dict`: A dictionary containing the optimized DDP model and its parameters.

# Steps
1. **Initialize DDP Model**: Initializes the DDP model using the provided data.
2. **Forward Pass Loop**: Performs forward passes until convergence or iteration limit is reached.
3. **Convergence Check**: Checks for convergence after each forward pass.
4. **Update Iteration Counter**: Updates the iteration counter for the DDP algorithm.
5. **Calculate Solve Time**: Calculates the total solve time for the optimization.
6. **Check Solver Status**: Checks the solver status and prints the results.
7. **Calculate Objective Value**: Calculates and prints the optimal objective function value.
8. **Return Data**: Returns the final dictionary containing the optimized DDP model and its parameters.
"""
function optimize_MPOPF_1ph_NL_DDP(data;
    verbose::Bool=false)

    ddpModel = DDPModel(data)

    keepForwardPassesRunning = true
    while keepForwardPassesRunning
        @unpack k_ddp = ddpModel
        myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
        ddpModel = forward_pass(ddpModel,
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

    # optimal_obj_value = objective_value(model)
    @unpack Tset = data
    optimal_obj_value = sum(modelVals[:objective_value_vs_t][t] for t ∈ Tset)
    println("Optimal objective function value: ", optimal_obj_value)

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
    maxiter::Int=10,
    verbose::Bool=false)

    @unpack Tset, Bset, solver = data;
    models_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Model}()
    modelVals_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Dict}()
    mu = Dict{Tuple{Int,Int,Int},Float64}()
    # modelVals = Dict{Symbol,Any}()
    modelVals = MC.ModelVals(data)
    # Initialize mu[j, t_ddp, 0/-1] = 0 for all j in Bset and t_ddp in Tset
    for j ∈ Bset, t_ddp ∈ Tset
        mu[j, t_ddp, 0] = 0.0
        mu[j, t_ddp, -1] = 0.0
    end

    # Initialize an empty model using the configure_solver function
    # model = configure_solver(solver)

    ddpModel = Dict(
        :converged => false,
        :data => data,
        :iterLimitReached => false,
        :k_ddp => 1,
        :maxiter => maxiter,
        # :model => model,
        :modelVals => modelVals,
        :models_ddp_vs_t_vs_k=>models_ddp_vs_t_vs_k,
        :modelVals_ddp_vs_t_vs_k=>modelVals_ddp_vs_t_vs_k,
        :mu=>mu,
        :shouldStop => false,
        :t_ddp=>0
    )
    return ddpModel
end
#endregion

end