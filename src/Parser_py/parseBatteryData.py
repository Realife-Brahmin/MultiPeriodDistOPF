import os
import re
from collections import defaultdict
from math import ceil

# Assuming `myprintln` is defined in a helper module
from src.helperFunctions import myprintln

def parse_battery_data(system_name, kVA_B=1000, N_L=None, verbose=False):
    """
    Parses battery data from a `.dss` file for the specified system.
    """
    
    # Get the working directory of this script
    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path for the Storage.dss file
    filename = os.path.join(wd, "..", "..", "rawData", system_name, "Storage.dss")

    myprintln(verbose, f"Reading Storage file from: {filename}")

    # Initialize output data structures
    Bset = set()
    B0, B0_pu, B_R, B_R_pu = {}, {}, {}, {}
    eta_C, eta_D = {}, {}
    P_B_R, P_B_R_pu, S_B_R, S_B_R_pu = {}, {}, {}, {}
    Vminpu_B, Vmaxpu_B = {}, {}
    soc_min, soc_max, soc_0 = {}, {}, {}

    # Regular expression to match key-value pairs
    kv_pattern = re.compile(r"(\w+)\s*=\s*([\S]+)")

    # Read and parse the Storage.dss file
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()
            myprintln(verbose, f"Processing line: {line}")

            # Skip empty lines
            if not line:
                myprintln(verbose, "Skipping empty line.")
                continue

            # Parse lines starting with "New Storage."
            if line.startswith("New Storage."):
                myprintln(verbose, "Found a New Storage entry.")

                # Parse storage information
                storage_info = {}
                for match in kv_pattern.finditer(line):
                    key, value = match.groups()
                    storage_info[key.strip()] = value.strip()
                
                myprintln(verbose, f"Parsed storage info: {storage_info}")

                # Extract bus number
                if "Bus1" in storage_info:
                    bus = int(storage_info["Bus1"].split(".")[0])
                    Bset.add(bus)
                    myprintln(verbose, f"Battery located at bus {bus}")
                else:
                    raise ValueError("Bus1 not specified for a storage in Storage.dss")

                # Extract rated power and energy values
                P_B_R[bus] = float(storage_info.get("kWrated", 0.0))
                P_B_R_pu[bus] = P_B_R[bus] / kVA_B
                S_B_R[bus] = float(storage_info.get("kVA", 0.0))
                S_B_R_pu[bus] = S_B_R[bus] / kVA_B
                B_R[bus] = float(storage_info.get("kWhrated", 0.0))
                B_R_pu[bus] = B_R[bus] / kVA_B
                myprintln(verbose, f"P_B_R: {P_B_R[bus]}, S_B_R: {S_B_R[bus]}, B_R: {B_R[bus]}")

                # Extract efficiencies
                eta_C[bus] = float(storage_info.get("EffCharge", 95)) / 100
                eta_D[bus] = float(storage_info.get("EffDischarge", 95)) / 100
                myprintln(verbose, f"eta_C: {eta_C[bus]}, eta_D: {eta_D[bus]}")

                # Extract minimum and maximum SOC
                soc_max[bus] = 0.95
                soc_min[bus] = float(storage_info.get("reserve", 30)) / 100
                myprintln(verbose, f"soc_max: {soc_max[bus]}, soc_min: {soc_min[bus]}")

                # Extract initial SOC percentage (soc_0)
                soc_0[bus] = float(storage_info.get("stored", 62.5)) / 100

                # Compute initial SOC (B0)
                B0[bus] = soc_0[bus] * B_R[bus]
                B0_pu[bus] = B0[bus] / kVA_B
                myprintln(verbose, f"Initial SOC B0: {B0[bus]}")

                # Extract voltage limits
                Vminpu_B[bus] = float(storage_info.get("Vminpu", 0.95))
                Vmaxpu_B[bus] = float(storage_info.get("Vmaxpu", 1.05))
                myprintln(verbose, f"Vminpu: {Vminpu_B[bus]}, Vmaxpu: {Vmaxpu_B[bus]}")

    # Compute the cardinality of Bset
    n_B = len(Bset)

    # Calculate Batt_percent
    if N_L is None:
        Batt_percent = 100
    else:
        Batt_percent = max(1, int(ceil(n_B / N_L * 100)))

    # Ensuring batteries return to initial SOC at the end of the optimization
    Bref = B0.copy()
    Bref_pu = B0_pu.copy()
    Bref_percent = {j: Bref_pu[j] / B_R_pu[j] for j in Bset}

    Bset = sorted(Bset)

    # Create a dictionary with all outputs
    storage_data = {
        "Bset": Bset,
        "n_B": n_B,
        "Batt_percent": Batt_percent,
        "B0": B0,
        "B0_pu": B0_pu,
        "Bref": Bref,
        "Bref_percent": Bref_percent,
        "Bref_pu": Bref_pu,
        "B_R": B_R,
        "B_R_pu": B_R_pu,
        "eta_C": eta_C,
        "eta_D": eta_D,
        "P_B_R": P_B_R,
        "P_B_R_pu": P_B_R_pu,
        "S_B_R": S_B_R,
        "S_B_R_pu": S_B_R_pu,
        "Vminpu_B": Vminpu_B,
        "Vmaxpu_B": Vmaxpu_B,
        "soc_min": soc_min,
        "soc_max": soc_max,
        "soc_0": soc_0
    }

    myprintln(verbose, f"Final parsed storage data: {storage_data}")

    return storage_data
