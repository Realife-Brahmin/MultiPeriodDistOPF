module DDP 

export attach_solver,
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

include("../helperFunctions.jl")
using .helperFunctions
include("../playbook_of_mpopf.jl")
using .playbook_of_mpopf

using Crayons
using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function attach_solver(model, solver_name)
    if solver_name == "Ipopt"
        optimizer = Ipopt.Optimizer
        optimizer_attributes = Dict("tol" => 1e-6, "max_iter" => 10000)
    elseif solver_name == "Gurobi"
        optimizer = Gurobi.Optimizer
        optimizer_attributes = Dict("TimeLimit" => 300)
    elseif solver_name == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
        optimizer_attributes = Dict()
    elseif solver_name == "EAGO"
        optimizer = EAGO.Optimizer
        optimizer_attributes = Dict()
    elseif solver_name == "MadNLP"
        optimizer = MadNLP.Optimizer
        optimizer_attributes = Dict()
    else
        error("Unsupported solver")
    end

    # Set the optimizer for the model
    set_optimizer(model, optimizer)

    # Set optimizer attributes
    for (attr, value) in optimizer_attributes
        set_optimizer_attribute(model, attr, value)
    end

    # Optionally, set the model to silent if needed
    set_silent(model)

    return model
end

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

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [1]
        modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0) # an unsolved model
        model_t0 = modelDict[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1] # I guess a solved model
        model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
        @unpack data = ddpModel;
        @unpack solver = data;
        model_t0 = attach_solver(model_t0, solver)
        # @show model_t0
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0) # same as model_t0_km1 as it is a copy
    μ = mu
    @unpack Bset = data;
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[j, t_ddp+1, k_ddp-1] - μ[j, t_ddp+1, k_ddp-2] ) * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
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

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be something like [2] or [3] or ... or [T-1]
        modelDict_t0 = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1]
        # model_t0 = deepcopy(model_t0_km1)
        model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
        @unpack data = ddpModel
        @unpack solver = data
        model_t0 = attach_solver(model_t0, solver)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass) and previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = modelVals[:B]
    @unpack Bset = data;
    for j ∈ Bset
        fix(model_t0[:B][j, t_ddp - 1], B_model[j, t_ddp - 1])
    end

    # Update the model with the future dual variables from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0)
    μ = mu
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[j, t_ddp+1, k_ddp-1] - μ[j, t_ddp+1, k_ddp-2] ) * (-model_t0[:B][j, t_ddp]) for j ∈ Bset )
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

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [T]
        modelDict_t0 = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict_t0[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[t_ddp, k_ddp-1]
        # model_t0 = deepcopy(model_t0_km1)
        model_t0, reference_t0_k_and_km1 = JuMP.copy_model(model_t0_km1)
        @unpack data = ddpModel
        @unpack solver = data
        model_t0 = attach_solver(model_t0, solver)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = modelVals[:B]
    @unpack Bset = data;
    for j ∈ Bset
        fix(model_t0[:B][j, t_ddp - 1], B_model[j, t_ddp - 1])
    end 

    # No need of updating objective function, as no μ term is present in the objective function for the terminal time-step

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    return ddpModel
end

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
        var_dict_current = create_variable_dict(model_current)
        var_dict_previous = create_variable_dict(model_previous)

        for var_name in keys(var_dict_current)
            if haskey(var_dict_previous, var_name)
                value_current = var_dict_current[var_name]
                value_previous = var_dict_previous[var_name]
                discrepancy = abs(value_current - value_previous)
                max_discrepancy = max(max_discrepancy, discrepancy)

                if max_discrepancy > threshold
                    all_under_threshold = false
                    # myprintln(verbose, "Some updates exceed the threshold. So keep doing Forward Passes.")
                    return ddpModel
                end
            end
        end
    end

    myprintln(verbose, "All updates are under the threshold.")
    converged = true
    shouldStop = true
    @pack! ddpModel = converged, shouldStop
    return ddpModel
    
end

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
    end

    return ddpModel
end

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

    if termination_status(model_t0) == LOCALLY_SOLVED
        myprintln(verbose, "Optimal solution found for Forward Step model for t = $(t_ddp)")
    else
        myprintln(verbose, "Optimal solution not found for Forward Step model for t = $(t_ddp)")
    end

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): ", optimal_obj_value)

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k;

    Tset = [t_ddp]
    # @unpack modelDict = ddpModel;
    ddpModel = copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel
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

    if termination_status(model_t0) == LOCALLY_SOLVED
        myprintln(verbose, "Optimal solution found for Forward Step model for t = $(t_ddp)")
    else
        myprintln(verbose, "Optimal solution not found for Forward Step model for t = $(t_ddp)")
    end

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): ", optimal_obj_value)

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    Tset = [t_ddp]
    ddpModel = copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel
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

    if termination_status(model_t0) == LOCALLY_SOLVED
        myprintln(verbose, "Optimal solution found for Forward Step model for t = $(t_ddp)")
    else
        myprintln(verbose, "Optimal solution not found for Forward Step model for t = $(t_ddp)")
    end

    # Check solver status and retrieve results
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(crayon_red("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println("Forward Pass k_ddp = $(k_ddp) : Best objective function value for t = $(t_ddp): ", optimal_obj_value)

    models_ddp_vs_t_vs_k[t_ddp, k_ddp] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    Tset = [t_ddp]
    ddpModel = copy_modelVals(ddpModel, model_t0, Tset=Tset)
    @unpack modelVals, modelVals_ddp_vs_t_vs_k = ddpModel
    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    ddpModel = backward_pass(ddpModel, model_t0, Tset=Tset)

    return ddpModel
end

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

function DDPModel(data;
    maxiter::Int=10,
    verbose::Bool=false)

    @unpack Tset, Bset, solver = data;
    models_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Model}()
    modelVals_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Dict}()
    mu = Dict{Tuple{Int,Int,Int},Float64}()
    # modelVals = Dict{Symbol,Any}()
    modelVals = ModelVals(data)
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

end