import os
from src.Parser_Py.parseSystemSimulationData import parse_system_simulation_data, post_process_data
from src.Parser_Py.parseBranchData import parse_branch_data
from src.Parser_Py.parseLoadData import parse_load_data
from src.Parser_Py.parsePVData import parse_pv_data
from src.Parser_Py.parseBatteryData import parse_battery_data
from src.evaluateVoltageLimits import evaluate_voltage_limits
from src.helperFunctions import generateBinaryLoadShape, myprintln

def parse_all_data(
    systemName: str,
    T: int,
    numAreas=1,
    alpha=1e-3,
    gamma=1e-3,
    objfun0="genCostMin",
    objfun2="scd",
    temporal_decmp=False,
    PSubsMax_kW=float("inf"),
    inputForecastDescription="nonspecific",
    solver="Ipopt",
    tSOC_hard=True
    ):
    # Parse system simulation data
    sysSimData = parse_system_simulation_data(
        systemName,
        T,
        numAreas=numAreas,
        alpha=alpha,
        gamma=gamma,
        objfun0=objfun0,
        objfun2=objfun2,
        temporal_decmp=temporal_decmp,
        PSubsMax_kW=PSubsMax_kW,
        inputForecastDescription=inputForecastDescription,
        solver=solver,
        tSOC_hard=tSOC_hard,
    )

    # Parse branch data
    branch_data = parse_branch_data(systemName)

    # Parse load data
    load_data = parse_load_data(systemName, T)
    N_L = load_data["N_L"]

    # Parse PV data
    pv_data = parse_pv_data(systemName, T, N_L=N_L)

    # Parse battery data
    battery_data = parse_battery_data(systemName, N_L=N_L)

    # Evaluate voltage limits for every bus based on various components (load, PV, battery) attached to it
    component_data = evaluate_voltage_limits(load_data, pv_data, battery_data)

    # Parse substation real power cost data
    cost_data = generateBinaryLoadShape(T)

    # Merge dictionaries
    data = {
        **sysSimData,
        **branch_data,
        **load_data,
        **pv_data,
        **battery_data,
        **component_data,
        **cost_data,
    }

    rootFolderPath = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    rawDataFolderPath = os.path.join(rootFolderPath, "rawData")
    processedDataFolderPath = os.path.join(rootFolderPath, "processedData")
    srcFolderPath = os.path.join(rootFolderPath, "src")

    data.update(
        {
            "processedDataFolderPath": processedDataFolderPath,
            "rawDataFolderPath": rawDataFolderPath,
            "rootFolderPath": rootFolderPath,
            "srcFolderPath": srcFolderPath,
        }
    )

    data = post_process_data(data)

    return data
