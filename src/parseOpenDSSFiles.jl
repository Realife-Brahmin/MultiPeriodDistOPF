# parseOpenDSSFiles.jl

module parseOpenDSSFiles

export parse_all_data

include("parseSystemSimulationData.jl")
include("parseBranchData.jl")
include("parseLoadData.jl")
# ... include other parsing scripts as needed

using .parseSystemSimulationData: parse_system_simulation_data
using .parseBranchData: parse_branch_data
using .parseLoadData: parse_load_data
# ... using other parsing modules as needed

function parse_all_data(systemName::String, T::Int)

    # Parse system simulation data
    sysSimData = parse_system_simulation_data(systemName)
    # Parse branch data
    branch_data = parse_branch_data(systemName)
    # Parse load data
    load_data = parse_load_data(systemName, T)
    # Merge dictionaries
    data = merge(sysSimData, branch_data, load_data)

    return data
end

end # module