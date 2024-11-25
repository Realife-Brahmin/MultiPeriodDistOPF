# import os
# import re
# from math import inf
# from collections import defaultdict

# def parse_system_simulation_data(system_name, T,kVA_B=1000):
#     """
#     Parses system simulation data from a `.dss` file and sets up simulation parameters.
#     """
    
#     # Default parameters
#     substation_bus = 1         # Default substation bus number
#     V_Subs = 1.0               # Default per-unit voltage at substation
#     kV_B = 1.0                 # Default base voltage in kV (line-to-ground)
#     delta_t = min(1.0, 24.0 / T)  # Default timestep

#     # Get the working directory of this script
#     wd = os.path.dirname(os.path.abspath(__file__))
#     # Construct the file path for SysSim.dss
#     filename = os.path.join(wd, "..", "..", "rawData", system_name, "SysSim.dss")

#     # Open and read the SysSim.dss file
#     with open(filename, "r") as file:
#         for line in file:
#             # Remove comments and strip whitespace
#             line = line.split("!")[0].strip()
#             if not line:
#                 continue

#             # Parse the 'Edit "Vsource.source"' line
#             if line.startswith('Edit "Vsource.source"'):
#                 tokens = line.split()
#                 for token in tokens[2:]:
#                     if "=" in token:
#                         key, value = token.split("=")
#                         key, value = key.strip(), value.strip()
#                         if key == "bus1":
#                             substation_bus = int(value)
#                         elif key == "pu":
#                             V_Subs = float(value)
#                         elif key == "basekv":
#                             kV_B = float(value)

#             # Parse the 'Set stepsize' line
#             elif line.startswith("Set stepsize"):
#                 match = re.search(r"Set stepsize\s*=\s*(.*)", line)
#                 if match:
#                     stepsize_str = match.group(1).strip()
#                     if "h" in stepsize_str:
#                         delta_t = float(stepsize_str.replace("h", ""))
#                     elif "m" in stepsize_str:
#                         delta_t = float(stepsize_str.replace("m", "")) / 60.0
#                     elif "s" in stepsize_str:
#                         delta_t = float(stepsize_str.replace("s", "")) / 3600.0
#                     else:
#                         delta_t = float(stepsize_str)

#     MVA_B = kVA_B / 1000
#     Z_B = (kV_B) ** 2 / MVA_B

#     Tset = sorted(list(range(1, T + 1)))

#     sys_sim_data = {
#         "system_name": system_name,
#         "substation_bus": substation_bus,
#         "V_Subs": V_Subs,
#         "kV_B": kV_B,
#         "kVA_B": kVA_B,
#         "Z_B": Z_B,
#         "delta_t": delta_t,
#         "T": T,
#         "Tset": Tset
#     }

#     return sys_sim_data

# def post_process_data(data):
#     """
#     Post-process the system simulation data to add additional descriptive fields.
#     """
#     DER_percent = data.get("DER_percent", 0)
#     Batt_percent = data.get("Batt_percent", 0)
#     ged_string = f"{DER_percent}% PVs and {Batt_percent}% Batteries"
#     ged_appendix = f"pv_{DER_percent}_batt_{Batt_percent}"

#     # Add to data dictionary
#     data["ged_string"] = ged_string
#     data["ged_appendix"] = ged_appendix

#     return data


## New parser for parsesystemsimulationdata

import os
import socket

def parse_system_simulation_data(system_name, T,
                                  num_areas=1,
                                  alpha=1e-3,
                                  kVA_B=1000,
                                  objfun0="genCostMin",
                                  objfun2="scd",
                                  temporal_decmp=False,
                                  PSubsMax_kW=float('inf'),
                                  input_forecast_description="nonspecific",
                                  solver="Ipopt"):

    # Initialize parameters with default values
    substation_bus = 1  # Default substation bus number
    V_Subs = 1.0  # Default per-unit voltage at substation
    kV_B = 1.0  # Default base voltage in kV (line-to-ground)
    delta_t = min(1.0, 24.0 / T)  # Default step size

    # Get working directory
    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path for SysSim.dss
    filename = os.path.join(wd, "..", "..", "rawData", system_name, "SysSim.dss")

    # Open and read the SysSim.dss file
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()

            # Skip empty lines
            if not line:
                continue

            # Parse the 'Edit "Vsource.source"' line
            if line.startswith("Edit \"Vsource.source\""):
                tokens = line.split()
                for token in tokens[2:]:
                    if "=" in token:
                        key, value = map(str.strip, token.split("="))
                        if key == "bus1":
                            substation_bus = int(value)
                        elif key == "pu":
                            V_Subs = float(value)
                        elif key == "basekv":
                            kV_B = float(value)

            # Parse the 'Set stepsize' line
            elif line.startswith("Set stepsize"):
                stepsize_str = line.split("=")[-1].strip()
                if "h" in stepsize_str:
                    delta_t = float(stepsize_str.replace("h", ""))
                elif "m" in stepsize_str:
                    delta_t = float(stepsize_str.replace("m", "")) / 60.0
                elif "s" in stepsize_str:
                    delta_t = float(stepsize_str.replace("s", "")) / 3600.0
                else:
                    delta_t = float(stepsize_str)

    MVA_B = kVA_B / 1000
    Z_B = (kV_B ** 2) / MVA_B

    Tset = sorted(range(1, T + 1))

    if objfun0 == "subsPowerCostMin":
        objfun_string = "Cost of Substation Power"
        objfun_sense = "Min"
        objfun_prefix = "subsPowerCost_min"
        objfun_unit = "$"
    elif objfun0 == "lineLossMin":
        objfun_string = "Line Losses"
        objfun_sense = "Min"
        objfun_prefix = "lineLoss_min"
        objfun_unit = "kW"
    elif objfun0 == "subsPowerMin":
        objfun_string = "Substation Power"
        objfun_sense = "Min"
        objfun_prefix = "subsPower_min"
        objfun_unit = "kW"
    else:
        objfun_string = "unknown objective"
        objfun_sense = "Min"
        objfun_prefix = "unknown_obj"
        objfun_unit = None

    objfun_appendix = "with_scd" if objfun2 == "scd" else ""
    objfun_concise_description = f"{objfun_prefix}_{objfun_appendix}"

    if temporal_decmp:
        temporal_decmp_string = "Temporally Decomposed via DDP"
        temporal_decmp_appendix = "tmrpl_dcmpsd"
    else:
        temporal_decmp_string = "Temporally Brute-forced"
        temporal_decmp_appendix = "tmprl_bruteforced"

    if num_areas > 1:
        spatial_decomp_string = f"Spatially Decomposed into {num_areas} areas"
        spatial_decomp_appendix = f"spat_dcmpsd_{num_areas}_areas"
    else:
        spatial_decomp_string = "Spatially Centralized"
        spatial_decomp_appendix = "spat_centr_system"

    machine_id = socket.gethostname()
    sim_nature_string = f"{temporal_decmp_string}, {spatial_decomp_string}"
    sim_nature_appendix = f"{temporal_decmp_appendix}_{spatial_decomp_appendix}"

    sys_sim_data = {
        "alpha": alpha,
        "input_forecast_description": input_forecast_description,
        "machine_id": machine_id,
        "macro_iterations_completed": 0,
        "system_name": system_name,
        "num_areas": num_areas,
        "solution_time": -1,
        "substationBus": substation_bus,
        "V_Subs": V_Subs,
        "kV_B": kV_B,
        "kVA_B": kVA_B,
        "Z_B": Z_B,
        "delta_t": delta_t,
        "objfun0": objfun0,
        "objfun2": objfun2,
        "objfun_string": objfun_string,
        "objfun_sense": objfun_sense,
        "objfun_prefix": objfun_prefix,
        "objfun_appendix": objfun_appendix,
        "objfun_concise_description": objfun_concise_description,
        "objfun_unit": objfun_unit,
        "PSubsMax_kW": PSubsMax_kW,
        "sim_nature_appendix": sim_nature_appendix,
        "sim_nature_string": sim_nature_string,
        "solver": solver,
        "spatial_decomp_appendix": spatial_decomp_appendix,
        "spatial_decomp_string": spatial_decomp_string,
        "T": T,
        "temporal_decmp": temporal_decmp,
        "temporal_decmp_string": temporal_decmp_string,
        "temporal_decmp_appendix": temporal_decmp_appendix,
        "Tset": Tset
    }

    return sys_sim_data

def post_process_data(data):
    der_percent = data.get("DER_percent", 0)
    batt_percent = data.get("Batt_percent", 0)
    ged_string = f"{der_percent}% PVs and {batt_percent}% Batteries"
    ged_appendix = f"pv_{der_percent}_batt_{batt_percent}"
    data.update({
        "ged_string": ged_string,
        "ged_appendix": ged_appendix
    })
    return data
