import os
import re
from math import ceil
from collections import defaultdict

# Assuming `generateLoadShape` is defined in a helper module
from src.helperFunctions import generateLoadShape

def parse_pv_data(system_name, T, N_L=None, kVA_B=1000, LoadShape=None, filename_load_shape=None):
    """
    Parses PV data from a `.dss` file for the specified system.
    """
    
    # Get the working directory of this script
    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path for the PVsystem.dss file
    filename = os.path.join(wd, "..", "..", "rawData", system_name, "PVsystem.dss")

    # Initialize data structures for PV systems
    Dset = set()
    p_D_R, p_D_R_pu = {}, {}
    S_D_R, S_D_R_pu = {}, {}
    irrad, p_D, p_D_pu = {}, defaultdict(list), defaultdict(list)
    Vminpu_D, Vmaxpu_D = {}, {}

    # Use provided LoadShape or generate it if not provided
    LoadShapePV = LoadShape if LoadShape is not None else generateLoadShape(T, filename_load_shape=filename_load_shape)

    # Regular expression to match key=value pairs
    kv_pattern = re.compile(r"(\w+)\s*=\s*([\S]+)")

    # Open and read the PVsystem.dss file
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()

            # Skip empty lines
            if not line:
                continue

            # Parse lines starting with "New PVsystem."
            if line.startswith("New PVsystem."):
                # Extract parameters from the line
                pv_info = {}
                for match in kv_pattern.finditer(line):
                    key, value = match.groups()
                    pv_info[key.strip()] = value.strip()

                # Extract bus number (e.g., Bus1=6.1)
                if "Bus1" in pv_info:
                    bus = int(pv_info["Bus1"].split(".")[0])
                    Dset.add(bus)
                else:
                    raise ValueError("Bus1 not specified for a PV in PVsystem.dss")

                # Extract rated active power (Pmpp)
                p_D_R[bus] = float(pv_info.get("Pmpp", 0.0))
                p_D_R_pu[bus] = p_D_R[bus] / kVA_B  # Convert to per-unit based on system base kVA

                # Extract rated apparent power (kVA)
                S_D_R[bus] = float(pv_info.get("kVA", 0.0))
                S_D_R_pu[bus] = S_D_R[bus] / kVA_B  # Convert to per-unit based on system base kVA

                # Extract irradiance (irrad)
                irrad[bus] = float(pv_info.get("irrad", 1.0))  # Default irradiance if not specified

                # Force bus voltage limits (default values since OpenDSS does not include these for PVs)
                Vminpu_D[bus] = float(pv_info.get("Vminpu", 0.95))
                Vmaxpu_D[bus] = float(pv_info.get("Vmaxpu", 1.05))

                # Calculate active power profile over time using the user-provided or generated LoadShapePV
                p_D[bus] = [p_D_R[bus] * LoadShapePV[t] for t in range(T)]
                p_D_pu[bus] = [p_D_R_pu[bus] * LoadShapePV[t] for t in range(T)]

    # Calculate DER_percent
    n_D = len(Dset)
    if N_L is None:
        DER_percent = 100
    else:
        DER_percent = max(1, int(ceil(n_D / N_L * 100)))

    Dset = sorted(Dset)

    # Return the extracted data as a dictionary
    pv_data = {
        "n_D": n_D,
        "DER_percent": DER_percent,
        "Dset": Dset,
        "p_D_R": p_D_R,
        "p_D_R_pu": p_D_R_pu,  # Per-unit rated active power
        "S_D_R": S_D_R,
        "S_D_R_pu": S_D_R_pu,  # Per-unit rated apparent power
        "irrad": irrad,
        "p_D": p_D,
        "p_D_pu": p_D_pu,  # Per-unit active power profile over time
        "Vminpu_D": Vminpu_D,
        "Vmaxpu_D": Vmaxpu_D,
        "LoadShapePV": LoadShapePV  # Store the LoadShapePV used (either user-supplied or generated)
    }

    return pv_data
