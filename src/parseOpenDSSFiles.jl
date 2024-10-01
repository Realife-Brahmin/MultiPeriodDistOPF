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

    # Parse system simulation data first to get T
    T, Tset, C, η_C, η_D, V_base, V_Subs, v_min, v_max = parse_system_simulation_data(system_sim_data_file)

    # Parse branch data
    Lset, L1set, Lm1set, r, x, Parent, Children = parse_branch_data(branch_data_file)

    # Parse load data (includes bus data)
    N, Nset, Nm1set, p_L, q_L = parse_load_data(load_data_file, T)

    # Parse PV data
    Dset, p_D = parse_pv_data(pv_data_file, T)

    # Parse battery data
    Bset, battery_params = parse_battery_data(battery_data_file)

    return N, Nset, Nm1set, Lset, L1set, Lm1set, r, x, Parent, Children,
    T, Tset, C, η_C, η_D, V_base, V_Subs, v_min, v_max, p_L, q_L, Dset, p_D, Bset, battery_params
end
