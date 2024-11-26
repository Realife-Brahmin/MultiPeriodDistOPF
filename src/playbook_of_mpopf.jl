# optimizer.jl
module Playbook_of_MPOPF

export optimize_MPOPF_1ph_NL_TemporallyBruteforced

include("./ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

using Crayons
using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function build_MPOPF_1ph_NL_model_t_in_Tset(data;
    Tset=nothing)

    @unpack solver = data

    # Define the optimization model including any specific solver settings
    model = configure_solver(solver)

    if Tset === nothing
        Tset = data[:Tset]
    end

    # Define variables
    model = MB.define_model_variables_1ph_NL_t_in_Tset(model, data, Tset=Tset)

    # ===========================
    # Constraints
    # ===========================

    # Substation node real power balance constraints
    model = MB.nodalRealPowerBalance_substation_t_in_Tset(model, data, Tset=Tset)

    # Non-substation node real power balance constraints
    model = MB.nodalRealPowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)

    # Non-substation node reactive power balance constraints
    model = MB.nodalReactivePowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)

    # KVL constraints for branches connected directly to the substation
    model = MB.KVL_substation_branches_t_in_Tset(model, data, Tset=Tset)

    # KVL constraints for branches not connected directly to the substation
    model = MB.KVL_non_substation_branches_t_in_Tset(model, data, Tset=Tset)

    # BCPF constraints for branches connected directly to the substation
    model = MB.BCPF_substation_branches_t_in_Tset(model, data, Tset=Tset)

    # BCPF constraints for branches not connected directly to the substation
    model = MB.BCPF_non_substation_branches_t_in_Tset(model, data, Tset=Tset)

    # Battery SOC trajectory equality constraints
    @unpack tSOC_hard = data;
    model = MB.battery_SOC_constraints_t_in_Tset(model, data, Tset=Tset, tSOC_hard=tSOC_hard)

    # Fixed substation voltage constraints
    model = MB.fixed_substation_voltage_constraints_t_in_Tset(model, data, Tset=Tset)

    # Voltage limits constraints
    model = MB.voltage_limits_constraints_t_in_Tset(model, data, Tset=Tset)

    # Reactive power limits for PV inverters
    model = MB.reactive_power_limits_PV_inverters_t_in_Tset(model, data, Tset=Tset)

    # Reactive power limits for battery inverters
    model = MB.reactive_power_limits_battery_inverters_t_in_Tset(model, data, Tset=Tset)

    # Charging power limits for batteries
    model = MB.charging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    # Discharging power limits for batteries
    model = MB.discharging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    # SOC limits for batteries
    model = MB.SOC_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    # Define objective function
    modelDict = MB.define_objective_function_t_in_Tset(model, data, Tset=Tset, tSOC_hard=tSOC_hard)

    @unpack model, data = modelDict;

    # Initialize variables
    model = MB.initialize_variables_1ph_NL_t_in_Tset(model, data, Tset=Tset)

    modelDict = Dict(
        :model => model,
        :data => data
    )

    return modelDict
end

function optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)

    Tset = data[:Tset] # In this case, Tset = [1, 2, ... T]
    modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset)

    @unpack model, data = modelDict
    optimize!(model)
    @pack! modelDict = model

    # Check solver status and retrieve results
    # Define crayons for green and red text
    green_crayon = Crayon(foreground=:light_green, bold=true)
    red_crayon = Crayon(foreground=:red, bold=true)

    if termination_status(model) == LOCALLY_SOLVED
        println(green_crayon("Optimal solution found."))
    else
        println(red_crayon("Optimization did not find an optimal solution."))
    end

    optimal_obj_value = objective_value(model)
    println("Optimal objective function value: ", optimal_obj_value)
    
    return modelDict

end

function optimize_MPOPF_1ph_NL_DDP(data;
    verbose::Bool=false)

    ddpModel = DDPModel(data)
    
    keepForwardPassesRunning = true
    while keepForwardPassesRunning
        @unpack k_ddp = ddpModel;
        myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
        ddpModel = ForwardPass(ddpModel,
        verbose=verbose)

        keepForwardPassesRunning = shouldStop(ddpModel)
    end

    # Check solver status and retrieve results
    @unpack iterLimitReached, converged, model = ddpModel
    
    # Define crayons for green and red text
    green_crayon = Crayon(foreground=:light_green, bold=true)
    red_crayon = Crayon(foreground=:red, bold=true)

    if iterLimitReached
        println(red_crayon("Maximum iterations reached."))
    elseif converged
        println(green_crayon("Optimization converged."))
    else
        println(red_crayon("Optimization did not converge."))
    end

    optimal_obj_value = objective_value(model)
    println("Optimal objective function value: ", optimal_obj_value)


    return ddpModel
end

function DDPModel(data;
    maxiter::Int=10,
    verbose::Bool=false)

    @unpack Tset, Bset, solver = data;
    models_ddp_vs_t_vs_k = Dict{Tuple{Int,Int},Model}()
    mu = Dict{Tuple{Int,Int,Int},Float64}()

    # Initialize mu[(j, t_ddp, 0)] = 0 for all j in Bset and t_ddp in Tset
    for j ∈ Bset, t_ddp ∈ Tset
        mu[(j, t_ddp, 0)] = 0.0
    end

    # Initialize an empty model using the configure_solver function
    model = configure_solver(solver)

    ddpModel = Dict(
        :converged => false,
        :data => data,
        :iterLimitReached => false,
        :k_ddp => 0,
        :maxiter => maxiter,
        :model => model,
        :models_ddp_vs_t_vs_k=>models_ddp_vs_t_vs_k,
        :mu=>mu,
        :t_ddp=>0
    )
    return ddpModel
end

function ForwardPass(ddpModel;
    verbose::Bool=false)
    
    @unpack k_ddp = ddpModel;
    myprintln(verbose, "Starting Forward Pass k_ddp = $(k_ddp)")
    t_ddp = 1
    @unpack data = ddpModel;
    @unpack Tset = data;
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

function shouldStop(ddpModel; verbose::Bool=false)
    @unpack k_ddp, maxiter, models_ddp_vs_t_vs_k, data = ddpModel;
    @unpack Tset = data;

    # Criterion 1: Check if k_ddp has crossed the maxiter threshold
    if k_ddp >= maxiter
        println("Maximum iterations reached: $k_ddp")
        return true
    end

    # Criterion 2: Check the magnitude of updates in the model decision variable values
    threshold = 1e-3  # Define your threshold here
    all_under_threshold = true

    for t_ddp in Tset
        model_current = models_ddp_vs_t_vs_k[(t_ddp, k_ddp)]
        model_previous = models_ddp_vs_t_vs_k[(t_ddp, k_ddp - 1)]

        max_discrepancy = 0.0

        for var in all_variables(model_current)
            if haskey(model_previous, var)
                value_current = value(var, model_current)
                value_previous = value(var, model_previous)
                discrepancy = abs(value_current - value_previous)
                max_discrepancy = max(max_discrepancy, discrepancy)

                if discrepancy > threshold
                    all_under_threshold = false
                end
            end
        end

        if verbose
            println("Max discrepancy for t_ddp = $t_ddp: $max_discrepancy")
        end
    end

    if all_under_threshold
        println("All updates are under the threshold.")
        return true
    else
        println("Some updates exceed the threshold.")
        return false
    end
end

function configure_solver(solver_name)
    if solver_name == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "max_iter", 10000)
        set_optimizer_attribute(model, "print_level", 5)
    elseif solver_name == "Gurobi"
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", 300)        # Limit time (in seconds)
    elseif solver == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
        model = Model(optimizer)
    elseif solver == "EAGO"
        model = Model(EAGO.Optimizer)
    elseif solver == "MadNLP"
        model = Model(MadNLP.Optimizer)
    else
        error("Unsupported solver")
    end

    return model
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

    model_t0 = models_ddp_vs_t_vs_k[(t_ddp, k_ddp)]
    optimize!(model_t0)

    if termination_status(model_t0) == LOCALLY_SOLVED
        myprintln(verbose, "Optimal solution found for Forward Step model for t = $(t_ddp)")
    else
        myprintln(verbose, "Optimal solution not found for Forward Step model for t = $(t_ddp)")
    end

    # Check solver status and retrieve results
    green_crayon = Crayon(foreground=:light_green, bold=true)
    red_crayon = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(green_crayon("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(red_crayon("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
    end

    optimal_obj_value = objective_value(model_t0)
    println("Forward Pass k_ddp = $(k_ddp) : Optimal objective function value for t = $(t_ddp): ", optimal_obj_value)

    models_ddp_vs_t_vs_k[(t_ddp, k_ddp)] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k;

    ddpModel = update_variables_from_ForwardStep_into_MPOPF_model(ddpModel)

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

    @unpack t_ddp = ddpModel;
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
    @unpack t_ddp = ddpModel
    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to T = $(T)"
        return
    end

    ddpModel = build_ForwardStep_1ph_NL_model_t_is_T(ddpModel)
    ddpModel = optimize_ForwardStep_1ph_NL_model_t_is_T(ddpModel)

    return ddpModel
end

function build_ForwardStep_1ph_NL_model_t_is_1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, model, models_ddp_vs_t_vs_k, data, mu = ddpModel;

    if t_ddp != 1
        @error "t_ddp = $(t_ddp) is not equal to 1"
        return
    end

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [1]
        modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[(t_ddp, k_ddp-1)]
        model_t0 = deepcopy(model_t0_km1)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0)
    μ = mu
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[(j, t_ddp+1, k_ddp-1)] - μ[(j, t_ddp+1, k_ddp-2)] ) * (-model_t0[:B][(j, t_ddp)]) for j ∈ Bset )
    @objective(model_t0, Min, objfun_expr_t0_k)

    # B0 values are already set in the model so no need to fix them separately

    models_ddp_vs_t_vs_k[(t_ddp, k_ddp)] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return ddpModel

end

function build_ForwardStep_1ph_NL_model_t_in_2toTm1(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, model, models_ddp_vs_t_vs_k, data, mu = ddpModel

    if !(2 <= t_ddp <= T-1)
        @error "t_ddp = $(t_ddp) is not in [2, T-1]"
        return
    end

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be something like [2] or [3] or ... or [T-1]
        modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[(t_ddp, k_ddp-1)]
        model_t0 = deepcopy(model_t0_km1)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the last iteration (backward pass) and previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = model[:B]
    fix(model_t0[:B][(j, t_ddp - 1)], B_model[(j, t_ddp - 1)])

    # Update the model with the future dual variables from the last iteration (backward pass)

    objfun_expr_t0_km1 = objective_function(model_t0)
    μ = mu
    objfun_expr_t0_k = objfun_expr_t0_km1 + sum( ( μ[(j, t_ddp+1, k_ddp-1)] - μ[(j, t_ddp+1, k_ddp-2)] ) * (-model_t0[:B][(j, t_ddp)]) for j ∈ Bset )
    @objective(model_t0, Min, objfun_expr_t0_k)

    models_ddp_vs_t_vs_k[(t_ddp, k_ddp)] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k
    return

end

function build_ForwardStep_1ph_NL_model_t_is_T(ddpModel;
    verbose::Bool=false)
    
    @unpack k_ddp, t_ddp, model, models_ddp_vs_t_vs_k, data = ddpModel;

    if t_ddp != T
        @error "t_ddp = $(t_ddp) is not equal to T = $(T)"
        return
    end

    if k_ddp == 1
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Building Forward Step model for t = $(t_ddp)")

        Tset_t0 = [t_ddp] # should be [T]
        modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset_t0)
        model_t0 = modelDict[:model]
    elseif k_ddp >= 2
        myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Modifying last iteration's Forward Step model for t = $(t_ddp)")
        model_t0_km1 = models_ddp_vs_t_vs_k[(t_ddp, k_ddp-1)]
        model_t0 = deepcopy(model_t0_km1)
    else
        @error "Invalid value of k_ddp: $k_ddp"
        return
    end

    # Update the model with the solutions from the previous time-step's SOC values (forward pass)

    # Previous time-step's SOC values are constant for this model's equations, which have been solved for in the previous Forward Step
    B_model = model[:B]

    for j ∈ Bset
        fix(model_t0[:B][(j, t_ddp - 1)], B_model[(j, t_ddp - 1)])
    end 

    # No need of updating objective function, as no μ term is present in the objective function for the terminal time-step

    models_ddp_vs_t_vs_k[(t_ddp, k_ddp)] = model_t0
    @pack! ddpModel = models_ddp_vs_t_vs_k

    return ddpModel
end

function update_variables_from_ForwardStep_t_is_1_into_MPOPF_model(ddpModel;
    verbose::Bool=false)

    @unpack k_ddp, t_ddp, model, models_ddp_vs_t_vs_k, data, mu = ddpModel;

    if t_ddp == 1
        @error "t_ddp = $(t_ddp) is equal to 1"
        return
    end

    myprintln(verbose, "Forward Pass k_ddp = $(k_ddp): Updating MPOPF model with Forward Step model for t = $(t_ddp)")

    model_t0 = models_ddp_vs_t_vs_k[(t_ddp, k_ddp)]
    
    # Update the MPOPF model with the solutions from the Forward Step model for t = 1

    # right now don't know how exactly to achieve that

    # Then I also need to compute the dual variables for the constraints of the Forward Step model for t = 1 and store them in μ

    # this can be done, but I'll see to it later.
    
    return ddpModel
end

end # module Playbook_of_MPOPF
