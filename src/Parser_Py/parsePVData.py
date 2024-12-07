import os
from collections import defaultdict
from src.helperFunctions import generateLoadShape

def parse_pv_data(systemName: str, T: int, N_L=None, kVA_B=1000, LoadShape=None, filenameLoadShape=None):
    # Get wd: the path of <this> file
    wd = os.path.dirname(__file__)
    # Construct the file path for PVsystem.dss using wd
    filename_pv = os.path.join(wd, "..", "..", "rawData", systemName, "PVsystem.dss")

    # Initialize data structures for PV systems
    Dset = set()                     # Set of nodes with PVs
    p_D_R = {}            # Rated PV active power (Pmpp, kW)
    p_D_R_pu = {}         # Rated PV active power [pu]
    S_D_R = {}            # Rated PV apparent power (kVA)
    S_D_R_pu = {}         # Rated PV apparent power [pu]
    irrad = {}            # Irradiance values for each PV
    p_D = {}   # PV active power profile over time
    p_D_pu = {}  # PV active power profile over time [pu]
    Vminpu_D = {}         # PV bus voltage lower limit (artificially created data entry)
    Vmaxpu_D = {}         # PV bus voltage upper limit (artificially created data entry)

    LoadShapePV = LoadShape

    # If the user doesn't provide a LoadShape, generate it using the helper function
    if LoadShape is None:
        LoadShapePV = generateLoadShape(T, filenameLoadShape=filenameLoadShape)

    # Open and read the PVsystem.dss file
    with open(filename_pv, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()

            # Skip empty lines
            if not line:
                continue

            # Parse lines starting with "New PVsystem."
            if line.startswith("New PVsystem."):
                # Extract parameters from the line
                tokens = line.split()
                pv_info = {}

                for token in tokens[1:]:
                    if "=" in token:
                        key, value = map(str.strip, token.split("="))
                        pv_info[key] = value

                # Extract bus number (e.g., Bus1=6.1)
                if "Bus1" in pv_info:
                    bus_str = pv_info["Bus1"]
                    # Extract the integer part before the decimal point
                    bus_parts = bus_str.split(".")
                    j = int(bus_parts[0])  # Node number
                    Dset.add(j)
                else:
                    raise ValueError("Bus not specified for a PV in PVsystem.dss")

                # Extract rated active power (Pmpp, kW)
                if "Pmpp" in pv_info:
                    p_D_R[j] = float(pv_info["Pmpp"])
                    p_D_R_pu[j] = p_D_R[j] / kVA_B
                else:
                    p_D_R[j] = 0.0  # Default to zero if not specified
                    p_D_R_pu[j] = p_D_R[j] / kVA_B

                # Extract rated apparent power (kVA)
                if "kVA" in pv_info:
                    S_D_R[j] = float(pv_info["kVA"])
                    S_D_R_pu[j] = S_D_R[j] / kVA_B
                else:
                    S_D_R[j] = 0.0  # Default to zero if not specified
                    S_D_R_pu[j] = S_D_R[j] / kVA_B

                # Extract irradiance values
                if "irradiance" in pv_info:
                    irrad[j] = float(pv_info["irradiance"])
                else:
                    irrad[j] = 1.0  # Default to 1.0 if not specified

                # Extract minimum per-unit voltage (Vminpu)
                if "Vminpu" in pv_info:
                    Vminpu_D[j] = float(pv_info["Vminpu"])
                else:
                    Vminpu_D[j] = 0.95  # Default value if not specified

                # Extract maximum per-unit voltage (Vmaxpu)
                if "Vmaxpu" in pv_info:
                    Vmaxpu_D[j] = float(pv_info["Vmaxpu"])
                else:
                    Vmaxpu_D[j] = 1.05  # Default value if not specified

                # Initialize PV profiles using the provided or generated LoadShapePV
                for t in range(1, T + 1):
                    p_D[(j, t)] = p_D_R[j] * LoadShapePV[t - 1]
                    p_D_pu[(j, t)] = p_D_R_pu[j] * LoadShapePV[t - 1]

    n_D = len(Dset)
    DER_percent = 0

    # If the user doesn't provide a N_L, set DER_percent to 100, else compute an actual percentage
    if N_L is None:
        DER_percent = 100
    else:  # 1% if it is less than that (but nonzero)
        DER_percent = int((n_D / N_L) * 100)

    Dset = sorted(Dset)

    pvData = {
        "n_D": n_D,
        "DER_percent": DER_percent,
        "Dset": Dset,
        "p_D_R": p_D_R,
        "p_D_R_pu": p_D_R_pu,
        "S_D_R": S_D_R,
        "S_D_R_pu": S_D_R_pu,
        "irrad": irrad,
        "p_D": p_D,
        "p_D_pu": p_D_pu,
        "Vminpu_D": Vminpu_D,
        "Vmaxpu_D": Vmaxpu_D,
        "LoadShapePV": LoadShapePV  # Store the LoadShapePV used
    }

    # Return the extracted data as a dictionary
    return pvData
