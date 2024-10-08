# parseOpenDSSFiles.jl

module parseOpenDSSFiles

export parse_all_data

include("parseSystemSimulationData.jl")
include("parseBranchData.jl")
include("parseLoadData.jl")
include("parsePVData.jl")
include("parseBatteryData.jl")
include("helperFunctions.jl")
# ... include other parsing scripts as needed

using .parseSystemSimulationData: parse_system_simulation_data
using .parseBranchData: parse_branch_data
using .parseLoadData: parse_load_data
using .parsePVData: parse_pv_data
using .parseBatteryData: parse_battery_data
using .helperFunctions: generateBinaryLoadShape
# ... using other parsing modules as needed

function parse_all_data(systemName::String, T::Int)

    # Parse system simulation data
    sysSimData = parse_system_simulation_data(systemName)
    # Parse branch data
    branch_data = parse_branch_data(systemName)
    # Parse load data
    load_data = parse_load_data(systemName, T)
    # Parse PV data
    pv_data = parse_pv_data(systemName, T)
    # Parse Battery data
    battery_data = parse_battery_data(systemName)
    # Parse substation real power cost data
    cost_data = generateBinaryLoadShape(T)
    # Merge dictionaries
    data = merge(sysSimData, branch_data, load_data, pv_data, battery_data, cost_data)

    Tset = Set(1:T)
    data[:Tset] = Tset
    
    return data
end

end # module
