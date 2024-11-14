import os
import re
from math import inf
from collections import defaultdict

def parse_system_simulation_data(system_name, T,kVA_B=1000):
    """
    Parses system simulation data from a `.dss` file and sets up simulation parameters.
    """
    
    # Default parameters
    substation_bus = 1         # Default substation bus number
    V_Subs = 1.0               # Default per-unit voltage at substation
    kV_B = 1.0                 # Default base voltage in kV (line-to-ground)
    delta_t = min(1.0, 24.0 / T)  # Default timestep

    # Get the working directory of this script
    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path for SysSim.dss
    filename = os.path.join(wd, "..", "..", "rawData", system_name, "SysSim.dss")

    # Open and read the SysSim.dss file
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()
            if not line:
                continue

            # Parse the 'Edit "Vsource.source"' line
            if line.startswith('Edit "Vsource.source"'):
                tokens = line.split()
                for token in tokens[2:]:
                    if "=" in token:
                        key, value = token.split("=")
                        key, value = key.strip(), value.strip()
                        if key == "bus1":
                            substation_bus = int(value)
                        elif key == "pu":
                            V_Subs = float(value)
                        elif key == "basekv":
                            kV_B = float(value)

            # Parse the 'Set stepsize' line
            elif line.startswith("Set stepsize"):
                match = re.search(r"Set stepsize\s*=\s*(.*)", line)
                if match:
                    stepsize_str = match.group(1).strip()
                    if "h" in stepsize_str:
                        delta_t = float(stepsize_str.replace("h", ""))
                    elif "m" in stepsize_str:
                        delta_t = float(stepsize_str.replace("m", "")) / 60.0
                    elif "s" in stepsize_str:
                        delta_t = float(stepsize_str.replace("s", "")) / 3600.0
                    else:
                        delta_t = float(stepsize_str)

    MVA_B = kVA_B / 1000
    Z_B = (kV_B) ** 2 / MVA_B

    Tset = sorted(list(range(1, T + 1)))

    sys_sim_data = {
        "system_name": system_name,
        "substation_bus": substation_bus,
        "V_Subs": V_Subs,
        "kV_B": kV_B,
        "kVA_B": kVA_B,
        "Z_B": Z_B,
        "delta_t": delta_t,
        "T": T,
        "Tset": Tset
    }

    return sys_sim_data

def post_process_data(data):
    """
    Post-process the system simulation data to add additional descriptive fields.
    """
    DER_percent = data.get("DER_percent", 0)
    Batt_percent = data.get("Batt_percent", 0)
    ged_string = f"{DER_percent}% PVs and {Batt_percent}% Batteries"
    ged_appendix = f"pv_{DER_percent}_batt_{Batt_percent}"

    # Add to data dictionary
    data["ged_string"] = ged_string
    data["ged_appendix"] = ged_appendix

    return data
