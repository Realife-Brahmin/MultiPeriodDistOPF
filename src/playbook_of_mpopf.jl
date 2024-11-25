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

function update_with_forwardStep_solutions!(model::Model,       model_ddp_t::Model;
    verbose::Bool=false)
    # Iterate over all variables in the surrogate model
    for var_ddp in all_variables(model_ddp_t)
        # Extract the variable's value
        var_value = value(var_ddp)

        if isnothing(var_value)
            @warn "Variable $(name(var_ddp)) in surrogate model has no value. Skipping."
            continue
        end

        # Get the variable name
        var_name = name(var_ddp)

        # Check if the variable exists in the overarching model
        if haskey(model, var_name)
            # Retrieve the corresponding variable from the overarching model
            var_in_main_model = model[:$(var_name)]

            # Update the value in the main model
            set_start_value(var_in_main_model, var_value)  # Set the value in the main model
        else
            @warn "Variable $(var_name) does not exist in the overarching model. Skipping."
        end
    end

    return model
end

function build_ddpMPOPF_1ph_NL_model_t_is_T(ddpModel, data;
    verbose::Bool=false)

    @unpack models_ddp_vs_t_vs_k = ddpModel # this loc assumes that models_ddp_vs_t_vs_k is at least an already defined dictionary (even if empty) in ddpModel

    @unpack k_ddp = ddpModel;
    if k_ddp == 1
        # cold start
        # copy the model locs here (maybe carve them into functions)
        # save the model as
        model_ddp_t0 = Model() # placeholder
    elseif k_ddp >= 2
        model_ddp_t0_km1 = models_ddp_vs_t_vs_k[(t0, k_ddp-1)]
        model_ddp_t0 = deepcopy(model_ddp_t0_km1)
        @unpack model = ddpModel; # because it has previous iteration's model's values saved
        B = model[:B]
        # modify hsoc equations (how to index them correctly?)
        # no need for using μ for terminal time-step, only modify objective function with f0, fscd, ftsoc
        model_ddp_t0 = Model() # placeholder
    else
        @error "Invalid value of k_ddp: $k_ddp"
    end

    models_ddp_vs_t_vs_k[t0, k_ddp] = model_ddp_t0 

    @pack! ddpModel = models_ddp_vs_t_vs_k;

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

end # module Playbook_of_MPOPF
