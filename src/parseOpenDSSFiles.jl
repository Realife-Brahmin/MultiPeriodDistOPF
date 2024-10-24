# parseOpenDSSFiles.jl

module parseOpenDSSFiles

export myprintln, parse_all_data, parse_system_simulation_data, parse_branch_data, parse_load_data, parse_pv_data, parse_battery_data, evaluate_voltage_limits, generateBinaryLoadShape

include("parseSystemSimulationData.jl")
include("parseBranchData.jl")
include("parseLoadData.jl")
include("parsePVData.jl")
include("parseBatteryData.jl")
include("evaluateVoltageLimits.jl")
# include("helperFunctions.jl")

include("./helperFunctions.jl")
# using .helperFunctions

# ... include other parsing scripts as needed

using .parseSystemSimulationData: parse_system_simulation_data
using .parseBranchData: parse_branch_data
using .parseLoadData: parse_load_data
using .parsePVData: parse_pv_data
using .parseBatteryData: parse_battery_data
using .evaluateVoltageLimits: evaluate_voltage_limits
using .helperFunctions: generateBinaryLoadShape, myprintln
# ... using other parsing modules as needed

function parse_all_data(systemName::String, T::Int;
    numAreas=1,
    alpha=1e-3,
    objfun0="genCostMin",
    objfun2="scd")

    # Parse system simulation data
    sysSimData = parse_system_simulation_data(systemName, T,
    numAreas=numAreas, alpha=alpha, objfun0=objfun0, objfun2=objfun2)
    # Parse branch data
    branch_data = parse_branch_data(systemName)
    # Parse load data
    load_data = parse_load_data(systemName, T)
    N_L = load_data[:N_L]
    # Parse PV data
    pv_data = parse_pv_data(systemName, T, N_L=N_L)
    # Parse Battery data
    battery_data = parse_battery_data(systemName, N_L=N_L)
    # Evaluate Voltage Limits for every bus based on various components (load, pv, battery) attached to it
    component_data = evaluate_voltage_limits(load_data, pv_data, battery_data)
    # Parse substation real power cost data
    cost_data = generateBinaryLoadShape(T)
    # Merge dictionaries
    data = merge(sysSimData, branch_data, load_data, pv_data, battery_data, component_data, cost_data)
    
    return data
end

end # module
