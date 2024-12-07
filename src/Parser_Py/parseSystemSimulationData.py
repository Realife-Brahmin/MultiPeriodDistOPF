import os
from socket import gethostname

def parse_system_simulation_data(
    systemName: str,
    T: int,
    numAreas=1,
    alpha=1e-3,
    gamma=1e-3,
    kVA_B=1000,
    objfun0="genCostMin",
    objfun2="scd",
    temporal_decmp=False,
    PSubsMax_kW=float("inf"),
    inputForecastDescription="nonspecific",
    solver="Ipopt",
    tSOC_hard=True,
):
    # Initialize parameters with default values
    substationBus = 1       # Default substation bus number
    V_Subs = 1.0            # Default per-unit voltage at substation
    kV_B = 1.0              # Default base voltage in kV (line-to-ground)
    delta_t = min(1.0, 24.0 / T)  # Default time step size in hours

    wd = os.path.dirname(__file__)
    # Construct the file path for SysSim.dss
    filename = os.path.join(wd, "..", "..", "rawData", systemName, "SysSim.dss")

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
                for token in tokens[1:]:
                    if "=" in token:
                        key, value = map(str.strip, token.split("="))
                        if key == "bus1":
                            substationBus = int(value)
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

    Tset = list(range(1, T + 1))

    if objfun0 == "subsPowerCostMin":
        objfunString = "Cost of Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPowerCost_min"
        objfunUnit = "$"
    elif objfun0 == "lineLossMin":
        objfunString = "Line Losses"
        objfunSense = "Min"
        objfunPrefix = "lineLoss_min"
        objfunUnit = "kW"
    elif objfun0 == "subsPowerMin":
        objfunString = "Substation Power"
        objfunSense = "Min"
        objfunPrefix = "subsPower_min"
        objfunUnit = "kW"
    else:
        objfunString = "unknown objective"
        objfunSense = "Min"
        objfunPrefix = "unknown_obj"
        objfunUnit = ""

    objfunAppendix = "with_scd" if objfun2 == "scd" else ""
    objfunConciseDescription = f"{objfunPrefix}_{objfunAppendix}"

    if temporal_decmp:
        temporalDecmpString = "Temporally Decomposed via DDP"
        temporalDecmpAppendix = "tmprl_dcmpsd"
    else:
        temporalDecmpString = "Temporally Brute-forced"
        temporalDecmpAppendix = "tmprl_bruteforced"

    if numAreas > 1:
        spatialDecString = f"Spatially Decomposed into {numAreas} areas"
        spatialDecAppendix = f"spat_dcmpsd_{numAreas}_areas"
    else:
        spatialDecString = "Spatially Centralized"
        spatialDecAppendix = "spat_centr_system"

    machine_ID = gethostname()

    simNatureString = f"{temporalDecmpString}, {spatialDecString}"
    simNatureAppendix = f"{temporalDecmpAppendix}_{spatialDecAppendix}"

    sysSimData = {
        "alpha": alpha,
        "gamma": gamma,
        "inputForecastDescription": inputForecastDescription,
        "machine_ID": machine_ID,
        "macroItrsCompleted": 0,
        "systemName": systemName,
        "numAreas": numAreas,
        "solution_time": -1,
        "substationBus": substationBus,
        "V_Subs": V_Subs,
        "kV_B": kV_B,
        "kVA_B": kVA_B,
        "Z_B": Z_B,
        "delta_t": delta_t,
        "objfun0": objfun0,
        "objfun2": objfun2,
        "objfunString": objfunString,
        "objfunSense": objfunSense,
        "objfunPrefix": objfunPrefix,
        "objfunAppendix": objfunAppendix,
        "objfunConciseDescription": objfunConciseDescription,
        "objfunUnit": objfunUnit,
        "PSubsMax_kW": PSubsMax_kW,
        "simNatureAppendix": simNatureAppendix,
        "simNatureString": simNatureString,
        "solver": solver,
        "spatialDecAppendix": spatialDecAppendix,
        "spatialDecString": spatialDecString,
        "tSOC_hard": tSOC_hard,
        "T": T,
        "temporal_decmp": temporal_decmp,
        "temporalDecmpString": temporalDecmpString,
        "temporalDecmpAppendix": temporalDecmpAppendix,
        "Tset": Tset,
    }

    return sysSimData

def post_process_data(data):
    DER_percent = data.get("DER_percent", 0)
    Batt_percent = data.get("Batt_percent", 0)

    gedString = f"{DER_percent}% PVs and {Batt_percent}% Batteries"
    gedAppendix = f"pv_{DER_percent}_batt_{Batt_percent}"

    data.update({
        "gedString": gedString,
        "gedAppendix": gedAppendix,
    })

    return data
