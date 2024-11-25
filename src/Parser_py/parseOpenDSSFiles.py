# 

## New parser for parsing opendss files

import os
from .parseSystemSimulationData import parse_system_simulation_data, post_process_data
from .parseBranchData import parse_branch_data
from .parseLoadData import parse_load_data
from .parsePVData import parse_pv_data
from .parseBatteryData import parse_battery_data
from evaluateVoltageLimits import evaluate_voltage_limits
from helperFunctions import generate_binary_load_shape

def parse_all_data(system_name, T,
                   num_areas=1,
                   alpha=1e-3,
                   objfun0="genCostMin",
                   objfun2="scd",
                   temporal_decmp=False,
                   PSubsMax_kW=float('inf'),
                   input_forecast_description="nonspecific",
                   solver="Ipopt"):

    # Parse system simulation data
    sys_sim_data = parse_system_simulation_data(
        system_name, T,
        num_areas=num_areas,
        alpha=alpha,
        objfun0=objfun0,
        objfun2=objfun2,
        temporal_decmp=temporal_decmp,
        PSubsMax_kW=PSubsMax_kW,
        input_forecast_description=input_forecast_description,
        solver=solver
    )

    # Parse branch data
    branch_data = parse_branch_data(system_name)

    # Parse load data
    load_data = parse_load_data(system_name, T)
    N_L = load_data["N_L"]

    # Parse PV data
    pv_data = parse_pv_data(system_name, T, N_L=N_L)

    # Parse battery data
    battery_data = parse_battery_data(system_name, N_L=N_L)

    # Evaluate Voltage Limits for every bus based on various components (load, pv, battery) attached to it
    component_data = evaluate_voltage_limits(load_data, pv_data, battery_data)

    # Parse substation real power cost data
    cost_data = generate_binary_load_shape(T)

    # Merge dictionaries
    data = {**sys_sim_data, **branch_data, **load_data, **pv_data, **battery_data, **component_data, **cost_data}

    root_folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    raw_data_folder_path = os.path.join(root_folder_path, "rawData")
    processed_data_folder_path = os.path.join(root_folder_path, "processedData")
    src_folder_path = os.path.join(root_folder_path, "src")

    # Add paths to the data dictionary
    data.update({
        "processed_data_folder_path": processed_data_folder_path,
        "raw_data_folder_path": raw_data_folder_path,
        "root_folder_path": root_folder_path,
        "src_folder_path": src_folder_path
    })

    # Post-process data
    data = post_process_data(data)

    return data
