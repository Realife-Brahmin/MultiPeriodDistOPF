import os
from collections import defaultdict
from src.helperFunctions import generateLoadShape

def parse_load_data(systemName: str, T: int, kVA_B=1000, LoadShapeLoad=None, filenameLoadShape="LoadShapeDefault.dss"):
    # Get wd: the path of <this> file
    wd = os.path.dirname(__file__)
    # Construct the file path for Loads.dss using wd
    filename_load = os.path.join(wd, "..", "..", "rawData", systemName, "Loads.dss")

    # Initialize data structures
    NLset = set()                 # Set of nodes with loads
    p_L_R = {}      # Rated active power for each load (kW)
    p_L_R_pu = {}      # Rated active power for each load [pu]
    q_L_R = {}      # Rated reactive power for each load (kVAR)
    q_L_R_pu = {}      # Rated reactive power for each load [pu]
    Vminpu_L = {}    # Minimum per-unit voltage for each node
    Vmaxpu_L = {}    # Maximum per-unit voltage for each node

    # Placeholder for load shapes (to be updated later)
    p_L = {}  # Active power profile over time
    p_L_pu = {}  # Active power profile over time [pu]
    q_L = {}  # Reactive power profile over time
    q_L_pu = {}  # Reactive power profile over time [pu]

    # If user does not provide a LoadShapeLoad, generate one using the helper function
    if LoadShapeLoad is None:
        LoadShapeLoad = generateLoadShape(T, filenameLoadShape=filenameLoadShape)

    # Open and read the Loads.dss file
    with open(filename_load, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()

            # Skip empty lines
            if not line:
                continue

            # Parse lines starting with "New load."
            if line.startswith("New load."):
                # Extract parameters from the line
                tokens = line.split()
                load_info = {}

                for token in tokens[1:]:
                    if "=" in token:
                        key, value = map(str.strip, token.split("="))
                        load_info[key] = value

                # Extract bus number (e.g., Bus=2.1)
                if "Bus" in load_info:
                    bus_str = load_info["Bus"]
                    # Extract the integer part before the decimal point
                    bus_parts = bus_str.split(".")
                    j = int(bus_parts[0])  # Node number
                    NLset.add(j)
                else:
                    raise ValueError("Bus not specified for a load in Loads.dss")

                # Extract rated active power (kw)
                if "kw" in load_info:
                    p_L_R[j] = float(load_info["kw"])
                    p_L_R_pu[j] = p_L_R[j] / kVA_B
                else:
                    p_L_R[j] = 0.0  # Default to zero if not specified
                    p_L_R_pu[j] = p_L_R[j] / kVA_B

                # Extract rated reactive power (kvar)
                if "kvar" in load_info:
                    q_L_R[j] = float(load_info["kvar"])
                    q_L_R_pu[j] = q_L_R[j] / kVA_B
                else:
                    q_L_R[j] = 0.0  # Default to zero if not specified
                    q_L_R_pu[j] = q_L_R[j] / kVA_B

                # Extract minimum per-unit voltage (Vminpu)
                if "Vminpu" in load_info:
                    Vminpu_L[j] = float(load_info["Vminpu"])
                else:
                    Vminpu_L[j] = 0.95  # Default value if not specified

                # Extract maximum per-unit voltage (Vmaxpu)
                if "Vmaxpu" in load_info:
                    Vmaxpu_L[j] = float(load_info["Vmaxpu"])
                else:
                    Vmaxpu_L[j] = 1.05  # Default value if not specified

                # Initialize load profiles using the provided or generated LoadShapeLoad
                for t in range(1, T + 1):
                    p_L[(j, t)] = p_L_R[j] * LoadShapeLoad[t - 1]
                    p_L_pu[(j, t)] = p_L_R_pu[j] * LoadShapeLoad[t - 1]
                    q_L[(j, t)] = q_L_R[j] * LoadShapeLoad[t - 1]
                    q_L_pu[(j, t)] = q_L_R_pu[j] * LoadShapeLoad[t - 1]

    N_L = len(NLset)

    NLset = sorted(NLset)

    loadData = {
        "N_L": N_L,
        "NLset": NLset,
        "p_L_R": p_L_R,
        "p_L_R_pu": p_L_R_pu,
        "q_L_R": q_L_R,
        "q_L_R_pu": q_L_R_pu,
        "Vminpu_L": Vminpu_L,
        "Vmaxpu_L": Vmaxpu_L,
        "p_L": p_L,
        "p_L_pu": p_L_pu,
        "q_L": q_L,
        "q_L_pu": q_L_pu,
        "LoadShapeLoad": LoadShapeLoad  # Store the LoadShapeLoad used
    }

    # Return the extracted data as a dictionary
    return loadData
