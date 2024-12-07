import os
import re
from collections import defaultdict

def myprintln(verbose, message):
    if verbose:
        print(message)

def parse_battery_data(systemName,
    kVA_B=1000,
    N_L=None,
    verbose=False):

    # Get the working directory of this script
    wd = os.path.dirname(os.path.abspath(__file__))

    # Construct the full file path for the Storage.dss file
    filename = os.path.join(wd, "..", "..", "rawData", systemName, "Storage.dss")

    myprintln(verbose, f"Reading Storage file from: {filename}")

    # Initialize output data structures
    Bset = set()
    B0 = {}
    B0_pu = {}
    B_R = {}
    B_R_pu = {}
    eta_C = {}
    eta_D = {}
    P_B_R = {}
    P_B_R_pu = {}
    S_B_R = {}
    S_B_R_pu = {}
    Vminpu_B = {}
    Vmaxpu_B = {}
    soc_min = {}
    soc_max = {}
    soc_0 = {}

    # Regular expression to match key-value pairs
    kv_pattern = re.compile(r"(\w+)\s*=\s*([\S]+)")

    # Read and parse the Storage.dss file
    with open(filename, "r") as file:
        for line in file:
            line = line.split("!")[0].strip()
            myprintln(verbose, f"Processing line: {line}")

            if not line:
                myprintln(verbose, "Skipping empty line.")
                continue

            if line.startswith("New Storage."):
                myprintln(verbose, "Found a New Storage entry.")

                storage_info = {}
                for match in kv_pattern.finditer(line):
                    key, value = match.groups()
                    storage_info[key.strip()] = value.strip()

                myprintln(verbose, f"Parsed storage info: {storage_info}")

                if "Bus1" in storage_info:
                    bus_str = storage_info["Bus1"]
                    bus_parts = bus_str.split(".")
                    bus = int(bus_parts[0])
                    Bset.add(bus)
                    myprintln(verbose, f"Battery located at bus {bus}.")
                else:
                    raise ValueError("Bus1 not specified for a storage in Storage.dss")

                P_B_R[bus] = float(storage_info.get("kWrated", 0.0))
                P_B_R_pu[bus] = P_B_R[bus] / kVA_B
                S_B_R[bus] = float(storage_info.get("kVA", 0.0))
                S_B_R_pu[bus] = S_B_R[bus] / kVA_B
                B_R[bus] = float(storage_info.get("kWhrated", 0.0))
                B_R_pu[bus] = B_R[bus] / kVA_B

                myprintln(verbose, f"P_B_R: {P_B_R[bus]}, S_B_R: {S_B_R[bus]}, B_R: {B_R[bus]}")

                eta_C[bus] = float(storage_info.get("EffCharge", 95)) / 100
                eta_D[bus] = float(storage_info.get("EffDischarge", 95)) / 100
                myprintln(verbose, f"eta_C: {eta_C[bus]}, eta_D: {eta_D[bus]}")

                soc_max[bus] = 0.95
                soc_min[bus] = float(storage_info.get("reserve", 30)) / 100
                myprintln(verbose, f"soc_max: {soc_max[bus]}, soc_min: {soc_min[bus]}")

                soc_0[bus] = float(storage_info.get("stored", 62.5)) / 100

                B0[bus] = soc_0[bus] * B_R[bus]
                B0_pu[bus] = B0[bus] / kVA_B
                myprintln(verbose, f"Initial SOC B0: {B0[bus]}")

                Vminpu_B[bus] = float(storage_info.get("Vminpu", 0.95))
                Vmaxpu_B[bus] = float(storage_info.get("Vmaxpu", 1.05))
                myprintln(verbose, f"Vminpu: {Vminpu_B[bus]}, Vmaxpu: {Vmaxpu_B[bus]}")

    n_B = len(Bset)
    if N_L is None:
        Batt_percent = 100
    else:
        Batt_percent = int(-(-n_B // N_L * 100))

    Bref = B0
    Bref_pu = B0_pu
    Bref_percent = {j: Bref_pu[j] / B_R_pu[j] for j in Bset}

    Bset = sorted(Bset)

    storageData = {
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
        "soc_0": soc_0,
    }

    myprintln(verbose, f"Final parsed storage data: {storageData}")

    return storageData
