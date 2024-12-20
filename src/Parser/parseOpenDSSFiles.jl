module parseOpenDSSFiles

using Parameters: @pack!

export 
    evaluate_voltage_limits, 
    myprintln,
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
using .helperFunctions

function parse_all_data(systemName::String, T::Int;
    numAreas=1,
    alpha=1e-3,
    gamma=1e-3,
    objfun0="genCostMin",
    objfun2="scd",
    temporal_decmp=false,
    PSubsMax_kW=Inf,
    inputForecastDescription="nonspecific",
    solver="Ipopt",
    tSOC_hard=true)

    # Parse system simulation data
    sysSimData = parse_system_simulation_data(systemName, T,
        numAreas=numAreas, alpha=alpha, gamma=gamma, objfun0=objfun0, objfun2=objfun2, temporal_decmp=temporal_decmp,
        PSubsMax_kW=PSubsMax_kW,
        inputForecastDescription=inputForecastDescription, solver=solver,
        tSOC_hard=tSOC_hard)
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
    
    rootFolderPath = dirname(dirname(@__DIR__)) # This really assumes that this specific file lies in root/src/Parser/
    rawDataFolderPath = joinpath(rootFolderPath, "rawData")
    processedDataFolderPath = joinpath(rootFolderPath, "processedData")
    srcFolderPath = joinpath(rootFolderPath, "src")

    @pack! data = processedDataFolderPath, rawDataFolderPath, rootFolderPath, srcFolderPath
    data = post_process_data(data)

    return data
end

end # module parseOpenDSSFiles
