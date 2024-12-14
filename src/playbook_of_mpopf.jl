# optimizer.jl
module Playbook_of_MPOPF

export optimize_MPOPF_1ph_NL_TemporallyBruteforced

include("./ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("./helperFunctions.jl")
using .helperFunctions

include("./DDP/DDP.jl")
using .DDP

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
    crayon_red = Crayon(foreground=:red, bold=true)

    @unpack modelVals = modelDict
    termination_status = modelVals[:termination_status]

    crayon_final_green = Crayon(foreground=:white,background=:green, bold=true)
    if termination_status == LOCALLY_SOLVED
        println(crayon_final_green("Optimal solution found."))
    else
        println(crayon_red("Optimization did not find an optimal solution."))
    end

    # optimal_obj_value = objective_value(model)
    optimal_obj_value = modelVals[:objective_value]
    println("Optimal objective function value: ", optimal_obj_value)
    
    return modelDict

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

end # module Playbook_of_MPOPF
