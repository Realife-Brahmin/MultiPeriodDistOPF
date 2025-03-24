module parseOpenDSSFiles

using Parameters

export 
    evaluate_voltage_limits, 
    parse_all_data,
    parse_battery_data,
    parse_branch_data, 
    parse_load_data, 
    parse_pv_data,  
    post_process_data,
    parse_system_simulation_data

include("parseSystemSimulationData.jl")
using .parseSystemSimulationData

include("parseBranchData.jl")
using .parseBranchData

include("parseLoadData.jl")
using .parseLoadData

include("parsePVData.jl")
using .parsePVData

include("parseBatteryData.jl")
using .parseBatteryData

include("evaluateVoltageLimits.jl")
using .evaluateVoltageLimits

include("../helperFunctions.jl")
import .helperFunctions as HF

#region parse_all_data
"""
    parse_all_data(systemName::String, T::Int; numAreas=1, alpha=1e-3, gamma=1e-3, objfun0="genCostMin", objfun2="scd", temporal_decmp=false, PSubsMax_kW=Inf, inputForecastDescription="nonspecific", solver="Ipopt", tSOC_hard=true)

Parse all relevant data from OpenDSS files for a given system.

This function reads and parses various OpenDSS files for the specified system, extracting relevant data and initializing various parameters. 
It handles the extraction of branch, load, PV, and battery data, evaluates voltage limits, and merges all data into a single dictionary.

# Arguments
- `systemName::String`: The name of the system for which to parse data.
- `T::Int`: The number of time steps for the simulation.
- `numAreas::Int`: The number of areas in the system (default: 1).
- `alpha::Float64`: The alpha parameter for the objective function (default: 1e-3).
- `gamma::Float64`: The gamma parameter for the objective function (default: 1e-3).
- `objfun0::String`: The primary objective function type (default: "genCostMin").
- `objfun2::String`: The secondary objective function type (default: "scd").
- `temporal_decmp::Bool`: A flag to enable temporal decomposition (default: false).
- `PSubsMax_kW::Float64`: The maximum substation power in kW (default: Inf).
- `inputForecastDescription::String`: A description of the input forecast (default: "nonspecific").
- `solver::String`: The solver to use for optimization (default: "Ipopt").
- `tSOC_hard::Bool`: A flag to enable hard terminal SOC constraints (default: true).

# Returns
- `data::Dict`: A dictionary containing all parsed and processed data.

# Steps
1. **System Simulation Data**: Parses system simulation data.
2. **Branch Data**: Parses branch data from the BranchData.dss file.
3. **Load Data**: Parses load data from the Loads.dss file.
4. **PV Data**: Parses PV data from the PVData.dss file.
5. **Battery Data**: Parses battery data from the Storage.dss file.
6. **Voltage Limits Evaluation**: Evaluates voltage limits for each bus based on load, PV, and battery data.
7. **Cost Data**: Generates substation real power cost data.
8. **Data Merging**: Merges all parsed data into a single dictionary.
9. **Folder Paths**: Sets up folder paths for raw and processed data.
10. **Post-Processing**: Performs post-processing on the merged data.
11. **Return Data**: Returns the final dictionary containing all parsed and processed data.
"""
function parse_all_data(systemName::String, T::Int;
    numAreas=1,
    objfun0="subsPowerCostMin",
    objfun2="scd",
    temporal_decmp=false,
    PSubsMax_kW=Inf,
    inputForecastDescription="nonspecific",
    solver="Ipopt",
    tSOC_hard=false, 
    relax_terminal_soc_constraint=false,
    linearizedModel=false,
    gedDict_ud=nothing,
    alpha_fpi=0.43,
    warmStart_mu=nothing)

    # Parse system simulation data
    sysSimData = parse_system_simulation_data(systemName, T,
        numAreas=numAreas, objfun0=objfun0, objfun2=objfun2, temporal_decmp=temporal_decmp,
        PSubsMax_kW=PSubsMax_kW,
        inputForecastDescription=inputForecastDescription, solver=solver,
        tSOC_hard=tSOC_hard,
        relax_terminal_soc_constraint=relax_terminal_soc_constraint,
        linearizedModel=linearizedModel,
        gedDict_ud=gedDict_ud,
        alpha_fpi=alpha_fpi,
        warmStart_mu=warmStart_mu)

    @unpack kVA_B, kV_B = sysSimData
    # Parse branch data
    branch_data = parse_branch_data(systemName, kVA_B=kVA_B, kV_B=kV_B)
    @unpack kVA_B_dict, MVA_B_dict, kV_B_dict = branch_data;
    baseValuesDict = Dict(:kVA_B_dict=>kVA_B_dict, :kV_B_dict=>kV_B_dict, :MVA_B_dict=>MVA_B_dict)
    # Parse load data
    load_data = parse_load_data(systemName, T, baseValuesDict=baseValuesDict)
    @unpack N_L = load_data
    # Parse PV data
    pv_data = parse_pv_data(systemName, T, N_L=N_L, gedDict_ud=gedDict_ud, baseValuesDict=baseValuesDict)
    # Parse Battery data
    battery_data = parse_battery_data(systemName, N_L=N_L, gedDict_ud=gedDict_ud, baseValuesDict=baseValuesDict)
    # Evaluate Voltage Limits for every bus based on various components (load, pv, battery) attached to it
    component_data = evaluate_voltage_limits(load_data, pv_data, battery_data)
    # Parse substation real power cost data
    cost_data = HF.generateBinaryLoadShape(T)
    # Merge dictionaries
    data = merge(sysSimData, branch_data, load_data, pv_data, battery_data, component_data, cost_data)
    
    rootFolderPath = dirname(dirname(@__DIR__)) # This really assumes that this specific file lies in root/src/Parser/
    rawDataFolderPath = joinpath(rootFolderPath, "rawData")
    processedDataFolderPath = joinpath(rootFolderPath, "processedData")
    srcFolderPath = joinpath(rootFolderPath, "src")

    @pack! data = processedDataFolderPath, rawDataFolderPath, rootFolderPath, srcFolderPath
    data = post_process_data(data)

    return data
end
#endregion

end # module parseOpenDSSFiles
