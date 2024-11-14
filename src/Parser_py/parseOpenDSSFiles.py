# parse_open_dss_files/__init__.py

from src.Parser_py.parseSystemSimulationData import parse_system_simulation_data, post_process_data
from src.Parser_py.parseBranchData import parse_branch_data
from src.Parser_py.parseLoadData import parse_load_data
from src.Parser_py.parsePVData import parse_pv_data
from src.Parser_py.parseBatteryData import parse_battery_data
from src.evaluateVoltageLimits import evaluate_voltage_limits
from src.helperFunctions import generateLoadShape, myprintln, generateBinaryLoadShape

__all__ = [
    "myprintln",
    "parse_all_data",
    "parse_system_simulation_data",
    "parse_branch_data",
    "parse_load_data",
    "parse_pv_data",
    "parse_battery_data",
    "post_process_data",
    "evaluate_voltage_limits",
    "generateLoadShape"
]

def parse_all_data(system_name, T):
    # Parse system simulation data
    sys_sim_data = parse_system_simulation_data(system_name, T)
    
    # Parse branch data
    branch_data = parse_branch_data(system_name)
    
    # Parse load data
    load_data = parse_load_data(system_name, T)
    N_L = load_data["N_L"]
    
    # Parse PV data
    pv_data = parse_pv_data(system_name, T, N_L=N_L)
    
    # Parse battery data
    battery_data = parse_battery_data(system_name, N_L=N_L)
    
    # Evaluate voltage limits for each bus based on various components
    component_data = evaluate_voltage_limits(load_data, pv_data, battery_data)
    
    # Parse substation real power cost data
    cost_data = generateBinaryLoadShape(T)
    
    # Merge dictionaries
    data = {**sys_sim_data, **branch_data, **load_data, **pv_data, **battery_data, **component_data, **cost_data}
    
    # Post-process data
    data = post_process_data(data)
    
    return data
