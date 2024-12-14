# optimizer.jl
module Playbook_of_MPOPF

export optimize_MPOPF_1ph_NL_TemporallyBruteforced

include("./ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("./ModelCopier/ModelCopier.jl")
import .ModelCopier as MC

include("./helperFunctions.jl")
using .helperFunctions

include("./DDP/DDP.jl")
using .DDP

include("./SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger

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
    model = SolverArranger.configure_solver(solver)

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
    
    modelVals = MC.ModelVals(data)
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
    modelDict = MC.copy_modelVals(modelDict, model, Tset=Tset)
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

end # module Playbook_of_MPOPF
