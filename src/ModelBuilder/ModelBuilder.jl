module ModelBuilder

export 
    battery_SOC_constraints_t_in_Tset,
    build_MPOPF_1ph_NL_model_t_in_Tset,
    build_MPOPF_1ph_L_model_t_in_Tset,
    BCPF_non_substation_branches_1ph_NL_t_in_Tset,
    BCPF_substation_branches_1ph_NL_t_in_Tset,
    charging_power_limits_batteries_t_in_Tset,
    discharging_power_limits_batteries_t_in_Tset,
    fixed_substation_voltage_constraints_t_in_Tset,
    KVL_non_substation_branches_1ph_NL_t_in_Tset,
    KVL_non_substation_branches_1ph_L_t_in_Tset,
    KVL_substation_branches_1ph_NL_t_in_Tset,
    KVL_substation_branches_1ph_L_t_in_Tset,
    nodalReactivePowerBalance_non_substation_1ph_NL_t_in_Tset,
    nodalReactivePowerBalance_non_substation_1ph_L_t_in_Tset,
    nodalRealPowerBalance_non_substation_1ph_NL_t_in_Tset,
    nodalRealPowerBalance_non_substation_1ph_L_t_in_Tset,
    nodalRealPowerBalance_substation_t_in_Tset,
    reactive_power_limits_PV_inverters_t_in_Tset,
    reactive_power_limits_battery_inverters_1ph_NL_t_in_Tset,
    reactive_power_limits_battery_inverters_1ph_L_t_in_Tset,
    SOC_limits_batteries_t_in_Tset,
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

#region build_MPOPF_1ph_NL_model_t_in_Tset
"""
    build_MPOPF_1ph_NL_model_t_in_Tset(data; Tset=nothing)

Build the Multi-Period Optimal Power Flow (MPOPF) model for a single-phase network with nonlinear loads over a given time set.

This function constructs the MPOPF model by defining the optimization variables, constraints, and objective function. 
It handles various aspects of the model, including solver configuration, variable initialization, and constraint setup.

# Arguments
- `data::Dict`: A dictionary containing all necessary data and parameters for the model.
- `Tset::Union{Nothing, Vector{Int}}`: An optional vector specifying the time steps to consider. If not provided, it defaults to the full time set in `data`.

# Returns
- `modelDict::Dict`: A dictionary containing the constructed model, data, and model values.

# Steps
1. **Solver Configuration**: Configures the solver based on the provided data.
2. **Variable Definition**: Defines the optimization variables for the model.
3. **Constraint Setup**: Sets up various constraints, including:
    - Substation node real power balance
    - Non-substation node real power balance
    - Non-substation node reactive power balance
    - KVL constraints for substation and non-substation branches
    - BCPF constraints for substation and non-substation branches
    - Battery SOC trajectory equality constraints
    - Fixed substation voltage constraints
    - Voltage limits constraints
    - Reactive power limits for PV and battery inverters
    - Charging and discharging power limits for batteries
    - SOC limits for batteries
4. **Objective Function Definition**: Defines the objective function for the model, considering substation power cost minimization, line loss minimization, and optionally SOC discrepancies.
5. **Variable Initialization**: Initializes the variables for the model.

# Nuances
- The function handles different sets of variables and constraints based on the provided data.
- It includes optional terminal SOC constraints if `tSOC_hard` is true.
- The function ensures that the model is properly packed and returned with all necessary components.
"""
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
    modelDict = nodalRealPowerBalance_non_substation_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node reactive power balance constraints
    modelDict = nodalReactivePowerBalance_non_substation_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches connected directly to the substation
    modelDict = KVL_substation_branches_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches not connected directly to the substation
    modelDict = KVL_non_substation_branches_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches connected directly to the substation
    modelDict = BCPF_substation_branches_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # BCPF constraints for branches not connected directly to the substation
    modelDict = BCPF_non_substation_branches_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

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
    modelDict = reactive_power_limits_battery_inverters_1ph_NL_t_in_Tset(modelDict, Tset=Tset)

    # Charging power limits for batteries
    modelDict = charging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Discharging power limits for batteries
    modelDict = discharging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # SOC limits for batteries
    modelDict = SOC_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Define objective function
    @unpack tSOC_hard, relax_terminal_soc_constraint = data;
    modelDict = define_objective_function_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard, relax_terminal_soc_constraint=relax_terminal_soc_constraint)

    @unpack model, data = modelDict;

    # Initialize variables
    modelDict = initialize_variables_1ph_NL_t_in_Tset(modelDict, Tset=Tset)
    
    modelVals = MC.ModelVals(data)
    @pack! modelDict = model, modelVals, data

    return modelDict
end
#endregion

#region build_MPOPF_1ph_L_model_t_in_Tset
"""
    build_MPOPF_1ph_L_model_t_in_Tset(data; Tset=nothing)

Build the Multi-Period Optimal Power Flow (MPOPF) model for a single-phase network with linearized loads over a given time set.

This function constructs the MPOPF model by defining the optimization variables, constraints, and objective function. 
It handles various aspects of the model, including solver configuration, variable initialization, and constraint setup.

# Arguments
- `data::Dict`: A dictionary containing all necessary data and parameters for the model.
- `Tset::Union{Nothing, Vector{Int}}`: An optional vector specifying the time steps to consider. If not provided, it defaults to the full time set in `data`.

# Returns
- `modelDict::Dict`: A dictionary containing the constructed model, data, and model values.

# Steps
1. **Solver Configuration**: Configures the solver based on the provided data.
2. **Variable Definition**: Defines the optimization variables for the model.
3. **Constraint Setup**: Sets up various constraints, including:
    - Substation node real power balance
    - Non-substation node real power balance
    - Non-substation node reactive power balance
    - KVL constraints for substation and non-substation branches
    - Battery SOC trajectory equality constraints
    - Fixed substation voltage constraints
    - Voltage limits constraints
    - Reactive power limits for PV and battery inverters
    - Charging and discharging power limits for batteries
    - SOC limits for batteries
4. **Objective Function Definition**: Defines the objective function for the model, considering substation power cost minimization, line loss minimization, and optionally SOC discrepancies.
5. **Variable Initialization**: Initializes the variables for the model.

# Nuances
- The function handles different sets of variables and constraints based on the provided data.
- It includes optional terminal SOC constraints if `tSOC_hard` is true.
- The function ensures that the model is properly packed and returned with all necessary components.
"""
function build_MPOPF_1ph_L_model_t_in_Tset(data;
    Tset=nothing)

    @unpack solver = data

    # Define the optimization model including any specific solver settings
    model = SolverArranger.configure_solver(solver)

    if Tset === nothing
        Tset = data[:Tset]
    end

    modelDict = Dict{Symbol,Any}()
    @pack! modelDict = model, data
    # Define variables
    modelDict = define_model_variables_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # ===========================
    # Constraints
    # ===========================

    # Substation node real power balance constraints
    modelDict = nodalRealPowerBalance_substation_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node real power balance constraints
    modelDict = nodalRealPowerBalance_non_substation_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # Non-substation node reactive power balance constraints
    modelDict = nodalReactivePowerBalance_non_substation_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches connected directly to the substation
    modelDict = KVL_substation_branches_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # KVL constraints for branches not connected directly to the substation
    modelDict = KVL_non_substation_branches_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # Battery SOC trajectory equality constraints
    @unpack tSOC_hard = data
    modelDict = battery_SOC_constraints_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard)

    # Fixed substation voltage constraints
    modelDict = fixed_substation_voltage_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Voltage limits constraints
    modelDict = voltage_limits_constraints_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for PV inverters
    modelDict = reactive_power_limits_PV_inverters_t_in_Tset(modelDict, Tset=Tset)

    # Reactive power limits for battery inverters
    modelDict = reactive_power_limits_battery_inverters_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    # Charging power limits for batteries
    modelDict = charging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Discharging power limits for batteries
    modelDict = discharging_power_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # SOC limits for batteries
    modelDict = SOC_limits_batteries_t_in_Tset(modelDict, Tset=Tset)

    # Define objective function
    @unpack tSOC_hard, relax_terminal_soc_constraint = data
    modelDict = define_objective_function_t_in_Tset(modelDict, Tset=Tset, tSOC_hard=tSOC_hard, relax_terminal_soc_constraint=relax_terminal_soc_constraint)

    @unpack model, data = modelDict

    # Initialize variables
    modelDict = initialize_variables_1ph_L_t_in_Tset(modelDict, Tset=Tset)

    modelVals = MC.ModelVals(data)
    @pack! modelDict = model, modelVals, data

    return modelDict
end
#endregion

end # module ModelBuilder