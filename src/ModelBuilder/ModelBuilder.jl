module ModelBuilder

export BCPF_non_substation_branches_t_in_Tset,
    BCPF_substation_branches_t_in_Tset,
    KVL_non_substation_branches_t_in_Tset,
    KVL_substation_branches_t_in_Tset,
    SOC_limits_batteries_t_in_Tset,
    battery_SOC_constraints_t_in_Tset,
    charging_power_limits_batteries_t_in_Tset,
    define_model_variables_1ph_NL_t_in_Tset,
    define_objective_function_t_in_Tset,
    discharging_power_limits_batteries_t_in_Tset,
    estimate_alpha,
    fixed_substation_voltage_constraints_t_in_Tset,
    initialize_variables_1ph_NL_t_in_Tset,
    nodalReactivePowerBalance_non_substation_t_in_Tset,
    nodalRealPowerBalance_non_substation_t_in_Tset,
    nodalRealPowerBalance_substation_t_in_Tset,
    reactive_power_limits_PV_inverters_t_in_Tset,
    reactive_power_limits_battery_inverters_t_in_Tset,
    voltage_limits_constraints_t_in_Tset
    
include("Variables.jl")
using .Variables

include("Constraints.jl")
using .Constraints

include("Hyperparameters.jl")
using .Hyperparameters

include("Objective.jl")
using .Objective

include("Initialization.jl")
using .Initialization

include("../ModelCopier/ModelCopier.jl")
import .ModelCopier as MC

include("../SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger

using Parameters

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
    modelDict = define_model_variables_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # ===========================
    # Constraints
    # ===========================

    # Substation node real power balance constraints
    modelDict = nodalRealPowerBalance_substation_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node real power balance constraints
    modelDict = nodalRealPowerBalance_non_substation_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node reactive power balance constraints
    modelDict = nodalReactivePowerBalance_non_substation_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches connected directly to the substation
    modelDict = KVL_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches not connected directly to the substation
    modelDict = KVL_non_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches connected directly to the substation
    modelDict = BCPF_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches not connected directly to the substation
    modelDict = BCPF_non_substation_branches_t_in_Tset(modelDict, Tset=Tset)

    # Battery SOC trajectory equality constraints
    @unpack tSOC_hard = data;
    modelDict = battery_SOC_constraints_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard)

    # Fixed substation voltage constraints
    modelDict = fixed_substation_voltage_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Voltage limits constraints
    modelDict = voltage_limits_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for PV inverters
    modelDict = reactive_power_limits_PV_inverters_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for battery inverters
    modelDict = reactive_power_limits_battery_inverters_t_in_Tset(modelDict, Tset=Tset)

    # Charging power limits for batteries
    modelDict = charging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Discharging power limits for batteries
    modelDict = discharging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # SOC limits for batteries
    modelDict = SOC_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Define objective function
    modelDict = define_objective_function_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard)

    @unpack model, data = modelDict;

    # Initialize variables
    modelDict = initialize_variables_1ph_NL_t_in_Tset(modelDict, Tset=Tset)
    
    modelVals = MC.ModelVals(data)
    @pack! modelDict = model, modelVals, data

    return modelDict
end

end