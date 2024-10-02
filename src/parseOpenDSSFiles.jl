# parseOpenDSSFiles.jl

# Include individual parsing scripts
include("parseBranchData.jl")
include("parseLoadData.jl")  # Updated to include bus data parsing
include("parsePVData.jl")
include("parseBatteryData.jl")
include("parseSystemSimulationData.jl")

# Import functions from the parsing scripts
using .parseBranchData: parse_branch_data
using .parseLoadData: parse_load_data  # Updated
using .parsePVData: parse_pv_data
using .parseBatteryData: parse_battery_data
using .parseSystemSimulationData: parse_system_simulation_data

function parseOpenDSSFiles()
    # Paths to the .dss files (update with actual paths)
    branch_data_file = "BranchData.dss"
    load_data_file = "Loads.dss"
    pv_data_file = "PVSystem.dss"
    battery_data_file = "Storage.dss"
    system_sim_data_file = "SysSim.dss"

    substationBus, V_Subs, V_base, Δt = parse_system_simulation_data(systemName)

    # Parse branch data
    Lset, L1set, Lm1set, r, x, Parent, Children = parse_branch_data(branch_data_file)

    # Parse load data
    Nset, p_L_R, q_L_R, V_minpu, V_maxpu, p_L, q_L = parse_load_data(systemName, T)

    # Parse PV data
    Dset, p_D = parse_pv_data(pv_data_file, T)

    # Parse battery data
    Bset, battery_params = parse_battery_data(battery_data_file)

    return nothing
end
