import os
from collections import defaultdict

# Importing the helper function for load shape generation
from src.helperFunctions import generateLoadShape

def parse_load_data(system_name, T, kVA_B=1000, LoadShapeLoad=None, filename_load_shape="LoadShapeDefault.dss"):
    # Set working directory to the path of the current file
    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path for Loads.dss
    filename_load = os.path.join(wd, "..", "..", "rawData", system_name, "Loads.dss")
    
    # Initialize data structures
    NLset = set()  # Set of nodes with loads
    p_L_R = {}  # Rated active power for each load (kW)
    p_L_R_pu = {}  # Rated active power for each load [pu]
    q_L_R = {}  # Rated reactive power for each load (kVAR)
    q_L_R_pu = {}  # Rated reactive power for each load [pu]
    Vminpu_L = {}  # Minimum per-unit voltage for each node
    Vmaxpu_L = {}  # Maximum per-unit voltage for each node
    
    # Placeholder for load shapes (to be updated later)
    p_L = defaultdict(list)  # Active power profile over time
    p_L_pu = defaultdict(list)  # Active power profile over time [pu]
    q_L = defaultdict(list)  # Reactive power profile over time
    q_L_pu = defaultdict(list)  # Reactive power profile over time [pu]

    # If LoadShapeLoad is not provided, generate one using the helper function
    if LoadShapeLoad is None:
        LoadShapeLoad = generateLoadShape(T, filename_load_shape=filename_load_shape)
    
    # Open and read the Loads.dss file
    with open(filename_load, 'r') as file:
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
                        key, value = token.split("=")
                        load_info[key.strip()] = value.strip()
                
                # Extract bus number (e.g., Bus=2.1)
                if "Bus" in load_info:
                    bus_str = load_info["Bus"]
                    bus_parts = bus_str.split(".")
                    j = int(bus_parts[0])  # Node number
                    NLset.add(j)
                else:
                    raise ValueError("Bus not specified for a load in Loads.dss")

                # Extract rated active power (kw)
                p_L_R[j] = float(load_info.get("kw", 0.0))
                p_L_R_pu[j] = p_L_R[j] / kVA_B

                # Extract rated reactive power (kvar)
                q_L_R[j] = float(load_info.get("kvar", 0.0))
                q_L_R_pu[j] = q_L_R[j] / kVA_B

                # Extract minimum per-unit voltage (Vminpu)
                Vminpu_L[j] = float(load_info.get("Vminpu", 0.95))

                # Extract maximum per-unit voltage (Vmaxpu)
                Vmaxpu_L[j] = float(load_info.get("Vmaxpu", 1.05))

                # Initialize load profiles using the provided or generated LoadShapeLoad
                p_L[j] = [p_L_R[j] * LoadShapeLoad[t] for t in range(T)]
                p_L_pu[j] = [p_L_R_pu[j] * LoadShapeLoad[t] for t in range(T)]
                q_L[j] = [q_L_R[j] * LoadShapeLoad[t] for t in range(T)]
                q_L_pu[j] = [q_L_R_pu[j] * LoadShapeLoad[t] for t in range(T)]

    # Sorting NLset
    NLset = sorted(NLset)
    
    # Collect data in a dictionary
    load_data = {
        "N_L": len(NLset),
        "NLset": NLset,
        "p_L_R": p_L_R,
        "p_L_R_pu": p_L_R_pu,
        "q_L_R": q_L_R,
        "q_L_R_pu": q_L_R_pu,
        "Vminpu_L": Vminpu_L,
        "Vmaxpu_L": Vmaxpu_L,
        "p_L": dict(p_L),
        "p_L_pu": dict(p_L_pu),
        "q_L": dict(q_L),
        "q_L_pu": dict(q_L_pu),
        "LoadShapeLoad": LoadShapeLoad  # Store the LoadShapeLoad used
    }
    
    # Return the extracted data as a dictionary
    return load_data
