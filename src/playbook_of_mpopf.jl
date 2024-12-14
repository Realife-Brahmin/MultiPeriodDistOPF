# optimizer.jl
module Playbook_of_MPOPF

export optimize_MPOPF_1ph_NL_DDP,
    optimize_MPOPF_1ph_NL_TemporallyBruteforced

include("./ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("./helperFunctions.jl")
using .helperFunctions: myprintln

using Crayons
using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function ModelVals(data)
    modelVals = Dict{Symbol,Any}()

    # Extract necessary sets from data
    Tset = data[:Tset]
    Bset = data[:Bset]
    Dset = data[:Dset]
    Lset = data[:Lset]
    Nset = data[:Nset]

    # Initialize containers for each variable
    modelVals[:P_Subs] = Dict{Int,Float64}()  # t => value
    modelVals[:P] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()  # (i, j), t => value
    modelVals[:Q] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()
    modelVals[:l] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()
    modelVals[:v] = Dict{Tuple{Int,Int},Float64}() 
    modelVals[:q_D] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:q_B] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:P_c] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:P_d] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:B] = Dict{Tuple{Int,Int},Float64}()

    modelVals[:termination_status_vs_t] = Dict{Int,Any}()
    modelVals[:solve_time_vs_t] = Dict{Int,Float64}()
    modelVals[:objective_value_vs_t] = Dict{Int,Float64}()

    return modelVals
end

function copy_modelVals(modelDict, model_Tset;
    Tset=nothing) # modelDict could be ddpModel or modelDict (temporallyBruteforced)

    @unpack modelVals, data = modelDict;
    if Tset === nothing
        Tset = modelDict[:data][:Tset]
    end
    # Extract necessary sets from data
    @unpack Bset, Dset, Lset, Nset = data;

    # Retrieve variables from the model
    P_Subs_model = model_Tset[:P_Subs]
    P_model = model_Tset[:P]
    Q_model = model_Tset[:Q]
    l_model = model_Tset[:l]
    v_model = model_Tset[:v]
    q_D_model = model_Tset[:q_D]
    q_B_model = model_Tset[:q_B]
    P_c_model = model_Tset[:P_c]
    P_d_model = model_Tset[:P_d]
    B_model = model_Tset[:B]

    # Store values into modelVals using the indices from data

    # P_Subs[t]
    for t in Tset
        modelVals[:P_Subs][t] = value(P_Subs_model[t])
    end

    # P[(i,j), t], Q[(i,j), t], l[(i,j), t] for (i,j) in Lset
    for (i, j) in Lset, t in Tset
        modelVals[:P][(i, j), t] = value(P_model[(i, j), t])
        modelVals[:Q][(i, j), t] = value(Q_model[(i, j), t])
        modelVals[:l][(i, j), t] = value(l_model[(i, j), t])
    end

    # v[j, t] for j in Nset
    for j in Nset, t in Tset
        modelVals[:v][j, t] = value(v_model[j, t])
    end

    # q_D[j, t] for j in Dset
    for j in Dset, t in Tset
        modelVals[:q_D][j, t] = value(q_D_model[j, t])
    end

    # q_B[j, t], P_c[j, t], P_d[j, t], B[j, t] for j in Bset
    for j in Bset, t in Tset
        modelVals[:q_B][j, t] = value(q_B_model[j, t])
        modelVals[:P_c][j, t] = value(P_c_model[j, t])
        modelVals[:P_d][j, t] = value(P_d_model[j, t])
        modelVals[:B][j, t] = value(B_model[j, t])
    end

    @unpack T = data;
    if length(Tset) == 1
        for t ∈ Tset
            modelVals[:objective_value_vs_t][t] = objective_value(model_Tset)
            modelVals[:termination_status_vs_t][t] = termination_status(model_Tset)
            modelVals[:solve_time_vs_t][t] = solve_time(model_Tset)
        end
    elseif length(Tset) == T
        modelVals[:objective_value] = objective_value(model_Tset)
        modelVals[:termination_status] = termination_status(model_Tset)
        modelVals[:solve_time] = solve_time(model_Tset)
    else
        @error "Invalid length of Tset: $(length(Tset))"
        return
    end
    
    @pack! modelDict = modelVals;
    return modelDict
end

function build_MPOPF_1ph_NL_model_t_in_Tset(data;
    Tset=nothing)

    @unpack solver = data

    # Define the optimization model including any specific solver settings
    model = configure_solver(solver)

    if Tset === nothing
        Tset = data[:Tset]
    end

    modelDict = Dict{Symbol, Any}()
    @pack! modelDict = model, data
    # Define variables
    modelDict = MB.define_model_variables_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # ===========================
    # Constraints
    # ===========================

    # Substation node real power balance constraints
    modelDict = MB.nodalRealPowerBalance_substation_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node real power balance constraints
    modelDict = MB.nodalRealPowerBalance_non_substation_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node reactive power balance constraints
    modelDict = MB.nodalReactivePowerBalance_non_substation_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches connected directly to the substation
    modelDict = MB.KVL_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches not connected directly to the substation
    modelDict = MB.KVL_non_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches connected directly to the substation
    modelDict = MB.BCPF_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches not connected directly to the substation
    modelDict = MB.BCPF_non_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # Battery SOC trajectory equality constraints
    @unpack tSOC_hard = data;
    modelDict = MB.battery_SOC_constraints_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard)

    # Fixed substation voltage constraints
    modelDict = MB.fixed_substation_voltage_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Voltage limits constraints
    modelDict = MB.voltage_limits_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for PV inverters
    modelDict = MB.reactive_power_limits_PV_inverters_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for battery inverters
    modelDict = MB.reactive_power_limits_battery_inverters_t_in_Tset(modelDict, Tset=Tset)

    # Charging power limits for batteries
    modelDict = MB.charging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Discharging power limits for batteries
    modelDict = MB.discharging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # SOC limits for batteries
    modelDict = MB.SOC_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Define objective function
    modelDict = MB.define_objective_function_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard)

    @unpack model, data = modelDict;

    # Initialize variables
    modelDict = MB.initialize_variables_1ph_NL_t_in_Tset(modelDict, Tset=Tset)
    
    modelVals = ModelVals(data)
    @pack! modelDict = model, modelVals, data

    return modelDict
end

function optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)

    Tset = data[:Tset] # In this case, Tset = [1, 2, ... T]
    modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset)

    @unpack model, data = modelDict
    optimize!(model)
    
    @pack! modelDict = model

    # modelDict = generate_1ph_NL_model_decvar_value_dict(modelDict)
    modelDict = copy_modelVals(modelDict, model, Tset=Tset)
    # Check solver status and retrieve results
    # Define crayons for green and red text
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    red_crayon = Crayon(foreground=:red, bold=true)

    @unpack modelVals = modelDict
    termination_status = modelVals[:termination_status]

    if termination_status == LOCALLY_SOLVED
        println(crayon_light_green("Optimal solution found."))
    else
        println(red_crayon("Optimization did not find an optimal solution."))
    end

    # optimal_obj_value = objective_value(model)
    optimal_obj_value = modelVals[:objective_value]
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

        dppModel = check_for_ddp_convergence(ddpModel)
        
        @unpack shouldStop = ddpModel;
        keepForwardPassesRunning = !shouldStop

        k_ddp += 1
        @pack! ddpModel = k_ddp

    end

    @unpack modelVals = ddpModel;
    modelVals[:solve_time] = sum(modelVals[:solve_time_vs_t][t] for t ∈ data[:Tset])
    @pack! ddpModel = modelVals;
    # Check solver status and retrieve results
    @unpack iterLimitReached, converged, modelVals = ddpModel
    
    # Define crayons for green and red text
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    red_crayon = Crayon(foreground=:red, bold=true)

    if iterLimitReached
        println(red_crayon("Maximum iterations reached!"))
    elseif converged
        println(crayon_light_green("Optimization via DDP converged."))
    else
        @error "floc"
    end

    # optimal_obj_value = objective_value(model)
    @unpack Tset = data;
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

function ForwardPass(ddpModel;
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

# Function to create a dictionary of variable values
function create_variable_dict(model)
    var_dict = Dict{Symbol,Float64}()
    for v in all_variables(model)
        var_name = Symbol(name(v))
        var_dict[var_name] = value(v)
    end
    return var_dict
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

    @show k_ddp
    if k_ddp == 1
        println("No updates to check for k_ddp = 1")
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
                    myprintln(verbose, "Some updates exceed the threshold. So keep doing Forward Passes.")
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

function configure_solver(solver_name)
    if solver_name == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "max_iter", 10000)
        # set_optimizer_attribute(model, "print_level", 5)
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
    red_crayon = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(red_crayon("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
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
    # t_ddp = Tset[1]
    for t_ddp ∈ Tset
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
    end
    mu = μ
    @pack! ddpModel = mu

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
    red_crayon = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(red_crayon("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
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
    red_crayon = Crayon(foreground=:red, bold=true)

    if termination_status(model_t0) == LOCALLY_SOLVED
        println(crayon_light_green("Forward Pass k_ddp = $(k_ddp) : Optimal solution found for Forward Step model for t = $(t_ddp)"))
    else
        println(red_crayon("Forward Pass k_ddp = $(k_ddp) : Optimal solution not found for Forward Step model for t = $(t_ddp)"))
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

end # module Playbook_of_MPOPF
