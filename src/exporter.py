import os
import csv
from pathlib import Path

def myprintln(verbose, message):
    if verbose:
        print(message)

def export_optimization_model(modelDict, verbose=False):
    data = modelDict["data"]
    model = modelDict["model"]

    T = data["T"]
    systemName = data["systemName"]
    numAreas = data["numAreas"]
    gedAppendix = data["gedAppendix"]
    machine_ID = data["machine_ID"]
    objfunAppendix = data["objfunAppendix"]
    simNatureAppendix = data["simNatureAppendix"]

    base_dir = os.path.join(
        "processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}"
    )

    if not os.path.isdir(base_dir):
        print(f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_ID}_optimizationModel_{gedAppendix}_for_{objfunAppendix}_via_{simNatureAppendix}.txt",
    )

    if os.path.isfile(filename):
        os.remove(filename)

    with open(filename, "w") as file:
        file.write(str(model))

    myprintln(verbose, f"Model successfully written to {filename}")

def export_decision_variables(modelDict, filename="decision_variables.csv", verbose=False):
    data = modelDict["data"]
    modelVals = modelDict["modelVals"]

    T = data["T"]
    systemName = data["systemName"]
    numAreas = data["numAreas"]
    gedAppendix = data["gedAppendix"]
    machine_ID = data["machine_ID"]
    objfunAppendix = data["objfunAppendix"]
    simNatureAppendix = data["simNatureAppendix"]
    solver = data["solver"]

    base_dir = os.path.join(
        "processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}"
    )

    if not os.path.isdir(base_dir):
        print(f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    ext = ".csv"
    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_ID}_{solver}_decisionVariables_{gedAppendix}_for_{objfunAppendix}_via_{simNatureAppendix}{ext}",
    )

    myprintln(verbose, f"Saving to filename: {filename}")

    Tset = sorted(data["Tset"])
    LoadShapeLoad = data["LoadShapeLoad"]
    LoadShapePV = data["LoadShapePV"]
    LoadShapeCost = [x * 100 for x in data["LoadShapeCost"]]

    Lset = sorted(data["Lset"])
    Nset = sorted(data["Nset"])
    Bset = sorted(data["Bset"])
    Dset = sorted(data["Dset"])

    P_Subs = modelVals["P_Subs"]
    P = modelVals["P"]
    Q = modelVals["Q"]
    l = modelVals["l"]
    v = modelVals["v"]
    q_D = modelVals["q_D"]
    q_B = modelVals["q_B"]
    P_c = modelVals["P_c"]
    P_d = modelVals["P_d"]
    B = modelVals["B"]

    data_matrix = []

    data_matrix.append(["t"] + Tset)
    data_matrix.append(["lambda"] + LoadShapeLoad)
    data_matrix.append(["Irrad"] + LoadShapePV)
    data_matrix.append(["cents/kWh"] + LoadShapeCost)
    data_matrix.append(["P_Subs"] + [P_Subs[t] for t in Tset])

    for (i, j) in Lset:
        data_matrix.append([f"P_ij_{i}_{j}"] + [P[(i, j), t] for t in Tset])

    for (i, j) in Lset:
        data_matrix.append([f"Q_ij_{i}_{j}"] + [Q[(i, j), t] for t in Tset])

    for j in Nset:
        data_matrix.append([f"v_j_{j}"] + [v[j, t] for t in Tset])

    for j in Dset:
        data_matrix.append([f"q_D_j_{j}"] + [q_D[j, t] for t in Tset])

    for j in Bset:
        data_matrix.append([f"q_B_j_{j}"] + [q_B[j, t] for t in Tset])

    for j in Bset:
        data_matrix.append([f"P_c_j_{j}"] + [P_c[j, t] for t in Tset])

    for j in Bset:
        data_matrix.append([f"P_d_j_{j}"] + [P_d[j, t] for t in Tset])

    for j in Bset:
        data_matrix.append([f"B_j_{j}"] + [B[j, t] for t in Tset])

    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(data_matrix)

    myprintln(verbose, f"Decision variables exported to {filename}")

def export_simulation_key_results_txt(modelDict, filename="simulation_results.txt", verbose=False):
    """Export simulation key results to a text file."""
    data = modelDict["data"]

    T = data["T"]
    systemName = data["systemName"]
    numAreas = data["numAreas"]
    gedAppendix = data["gedAppendix"]
    machine_ID = data["machine_ID"]
    objfunConciseDescription = data["objfunConciseDescription"]
    simNatureAppendix = data["simNatureAppendix"]
    solver = data["solver"]

    base_dir = os.path.join("processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}")

    if not os.path.isdir(base_dir):
        myprintln(verbose, f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_ID}_{solver}_results_{gedAppendix}_for_{objfunConciseDescription}_via_{simNatureAppendix}.txt",
    )

    macroItrsCompleted = data.get("macroItrsCompleted", 0)
    solution_time = data.get("solution_time", -1)

    with open(filename, "w") as file:
        item_counter = 1

        # Header Section
        file.write("---------------------------------------------\n")
        file.write(f"{item_counter}. Machine ID: {data['machine_ID']}\n")
        item_counter += 1
        file.write(f"{item_counter}. Solver Used: {data['solver']}\n")
        item_counter += 1
        file.write(f"{item_counter}. System Name: {data['systemName']}\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Duration: {data['T']}\n")
        item_counter += 1
        file.write(f"{item_counter}. Nature of Optimization Simulation: {data['simNatureString']}\n")
        item_counter += 1
        file.write(f"{item_counter}. Objective: {data['objfunString']}\n")
        item_counter += 1
        file.write(f"{item_counter}. GED Configuration: {data['gedAppendix']}\n")
        item_counter += 1
        file.write(f"{item_counter}. Maximum Substation Power Allowed: {data['PSubsMax_kW']} kW\n")
        item_counter += 1
        file.write("---------------------------------------------\n")

        # Horizon Results
        file.write(f"Full {data['T']} Hour Horizon\n")
        file.write(f"{item_counter}. Horizon Total Cost of Substation Power: $ {round(data['PSubsCost_allT_dollar'], 2)}\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Line Loss: {round(data['PLoss_allT_kW'], 2)} kW\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Substation Power: {round(data['PSubs_allT_kW'], 2)} kW + {round(data['QSubs_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Load: {round(data['load_real_power_allT_kW'], 2)} kW + {round(data['load_reactive_power_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Generation: {round(data['total_gen_real_power_allT_kW'], 2)} kW + {round(data['total_gen_reactive_power_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Static Capacitor Reactive Power Generation: {round(data['static_cap_reactive_power_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total PV Generation: {round(data['pv_real_power_allT_kW'], 2)} kW + {round(data['pv_reactive_power_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Battery Generation: {round(data['battery_real_power_allT_kW'], 2)} kW + {round(data['battery_reactive_power_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total Battery Transaction Magnitude: {round(data['battery_real_power_transaction_magnitude_allT_kW'], 2)} kW + {round(data['battery_reactive_power_transaction_magnitude_allT_kVAr'], 2)} kVAr\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon Total SCD Observed: {round(data['scd_allT_kW'], 2)} kW\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon-end Battery Energy Deviation from Reference: {round(data['terminal_soc_violation_kWh'], 2)} kWh\n")
        item_counter += 1
        file.write(f"{item_counter}. Horizon-Total All Time Substation Power Peak: {round(data['substation_real_power_peak_allT_kW'], 2)} kW\n")
        item_counter += 1

        # Additional Metadata
        file.write("---------------------------------------------\n")
        file.write(f"{item_counter}. Number of Macro-Iterations: {macroItrsCompleted + 1}\n")
        item_counter += 1
        # file.write(f"{item_counter}. Simulation Time: {round(solution_time, 2)} seconds\n")
        # item_counter += 1
        file.write("---------------------------------------------\n")

    myprintln(verbose, f"Simulation key results exported to {filename}")


import os
import csv

def export_validation_decision_variables(modelDict, verbose=False):
    """Export validation decision variables to a CSV file."""
    valdVals = modelDict["valdVals"]
    data = modelDict["data"]

    # Unpack necessary variables
    T = data["T"]
    systemName = data["systemName"]
    numAreas = data["numAreas"]
    gedAppendix = data["gedAppendix"]
    machine_ID = data["machine_ID"]
    objfunConciseDescription = data["objfunConciseDescription"]
    processedDataFolderPath = data["processedDataFolderPath"]
    simNatureAppendix = data["simNatureAppendix"]
    solver = data["solver"]

    # Define base directory path
    base_dir = os.path.join(processedDataFolderPath, systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}")

    # Create the directory if it doesn't exist
    if not os.path.isdir(base_dir):
        if verbose:
            print(f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    # Additional parameters for the filename
    alphaAppendix = data["alphaAppendix"]
    gammaAppendix = data["gammaAppendix"]

    # Define the filename with the appropriate structure
    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_ID}_{solver}_validationDecisionVariables_{gedAppendix}_for_{objfunConciseDescription}_via_{simNatureAppendix}_alpha_{alphaAppendix}_gamma_{gammaAppendix}.csv"
    )

    # Write the `valdVals` dictionary to a CSV file
    try:
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            for key, value in valdVals.items():
                writer.writerow([key, value])
        if verbose:
            print(f"Validation decision variables written to {filename}")
    except Exception as e:
        if verbose:
            print(f"Error occurred while writing to {filename}: {e}")


def export_validation_key_results(
    modelDict,
    filename="validation_results.txt",
    verbose=False,
    printEveryTimeStepPowerflow=True
):
    valdVals = modelDict["valdVals"]
    data = modelDict["data"]

    # Define the path and filename based on the specified structure
    T = data["T"]
    systemName = data["systemName"]
    numAreas = data["numAreas"]
    gedAppendix = data["gedAppendix"]
    machine_ID = data["machine_ID"]
    simNatureAppendix = data["simNatureAppendix"]

    base_dir = os.path.join(
        "processedData", systemName, gedAppendix, f"Horizon_{T}", f"numAreas_{numAreas}"
    )

    if not os.path.isdir(base_dir):
        if verbose:
            print(f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    solver = data["solver"]
    alphaAppendix = data["alphaAppendix"]
    gammaAppendix = data["gammaAppendix"]
    objfunConciseDescription = data["objfunConciseDescription"]

    # Adjust filename based on `printEveryTimeStepPowerflow`
    if printEveryTimeStepPowerflow:
        filename = os.path.join(
            base_dir,
            f"Horizon_{T}_{machine_ID}_{solver}_valdResults_{gedAppendix}_for_{objfunConciseDescription}_via_{simNatureAppendix}_alpha_{alphaAppendix}_gamma_{gammaAppendix}_full.txt"
        )
    else:
        filename = os.path.join(
            base_dir,
            f"Horizon_{T}_{machine_ID}_{solver}_valdResults_{gedAppendix}_for_{objfunConciseDescription}_via_{simNatureAppendix}_alpha_{alphaAppendix}_gamma_{gammaAppendix}.txt"
        )

    # Open the file and write each section
    try:
        with open(filename, "w") as f:
            # Initialize output item counter
            item_counter = 1

            # Header Section
            f.write("---------------------------------------------\n")
            f.write(f"{item_counter}. Machine ID: {data['machine_ID']}\n")
            item_counter += 1
            f.write(f"{item_counter}. Solver Used: {data['solver']}\n")
            item_counter += 1
            f.write(f"{item_counter}. System Name: {data['systemName']}\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Duration: {data['T']}\n")
            item_counter += 1
            f.write(f"{item_counter}. Nature of Validation Simulation: {data['simNatureString']}\n")
            item_counter += 1
            f.write("---------------------------------------------\n")

            # Horizon Results
            f.write(f"Full {data['T']} Hour Horizon Validation Results\n")
            f.write(f"{item_counter}. Horizon Total Substation Power Cost: ${round(valdVals['vald_PSubsCost_allT_dollar'], 2)}\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Line Loss: {round(valdVals['vald_PLoss_allT_kW'], 2)} kW\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Substation Power: {round(valdVals['vald_PSubs_allT_kW'], 2)} kW + {round(valdVals['vald_QSubs_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Load: {round(valdVals['vald_load_real_power_allT_kW'], 2)} kW + {round(valdVals['vald_load_reactive_power_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Generation: {round(valdVals['vald_total_gen_real_power_allT_kW'], 2)} kW + {round(valdVals['vald_total_gen_reactive_power_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Static Capacitor Reactive Power Generation: {round(valdVals['vald_static_cap_reactive_power_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total PV Generation: {round(valdVals['vald_pv_real_power_allT_kW'], 2)} kW + {round(valdVals['vald_pv_reactive_power_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Battery Generation: {round(valdVals['vald_battery_real_power_allT_kW'], 2)} kW + {round(valdVals['vald_battery_reactive_power_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total Battery Transaction Magnitude: {round(valdVals['vald_battery_real_power_transaction_magnitude_allT_kW'], 2)} kW + {round(valdVals['vald_battery_reactive_power_transaction_magnitude_allT_kVAr'], 2)} kVAr\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon Total SCD Observed: N/A\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon-end Battery Energy Deviation from Reference: {round(valdVals['vald_terminal_soc_violation_kWh'], 2)} kWh\n")
            item_counter += 1
            f.write(f"{item_counter}. Horizon-Total All Time Substation Power Peak: {round(valdVals['vald_substation_real_power_peak_allT_kW'], 2)} kW\n")
            item_counter += 1

            # Discrepancies
            f.write("---------------------------------------------\n")
            f.write("Discrepancies (Maximum All Time):\n")
            f.write(f"{item_counter}. Maximum All Time Voltage Discrepancy: {round(valdVals['disc_voltage_all_time_pu'], 6)} pu\n")
            item_counter += 1
            f.write(f"{item_counter}. Maximum All Time Line Loss Discrepancy: {round(valdVals['disc_line_loss_all_time_kW'], 6)} kW\n")
            item_counter += 1
            f.write(f"{item_counter}. Maximum All Time Substation Borrowed Real Power Discrepancy: {round(valdVals['disc_PSubs_all_time_kW'], 6)} kW\n")
            item_counter += 1
            f.write(f"{item_counter}. Maximum All Time Substation Borrowed Reactive Power Discrepancy: {round(valdVals['disc_QSubs_all_time_kVAr'], 6)} kVAr\n")
            item_counter += 1

            # Optional: Print per-timestep power flow results if enabled
            if printEveryTimeStepPowerflow:
                f.write("\nPer-Timestep Power Flow Results:\n")
                for t in data["Tset"]:
                    f.write("\n" + "*" * 30 + "\n")
                    f.write(f"   Time Step: {t}\n")
                    f.write("*" * 30 + "\n")
                    f.write(f"   Power Loss              : {valdVals['vald_PLoss_vs_t_1toT_kW'][t]} kW\n")
                    f.write(f"   Substation Power        : {valdVals['vald_PSubs_vs_t_1toT_kW'][t]} kW\n")
                    f.write(f"   Reactive Power          : {valdVals['vald_QSubs_vs_t_1toT_kVAr'][t]} kVAr\n")
                    f.write(f"   Total Load Power        : {valdVals['vald_load_real_power_vs_t_1toT_kW'][t]} kW\n")
                    f.write(f"   Total Load Reactive Power: {valdVals['vald_load_reactive_power_vs_t_1toT_kVAr'][t]} kVAr\n")
                    f.write(f"   Total PV Power          : {valdVals['vald_pv_real_power_vs_t_1toT_kW'][t]} kW\n")
                    f.write(f"   Total PV Reactive Power : {valdVals['vald_pv_reactive_power_vs_t_1toT_kVAr'][t]} kVAr\n")
                    f.write(f"   Total Battery Power     : {valdVals['vald_battery_real_power_vs_t_1toT_kW'][t]} kW\n")
                    f.write(f"   Total Battery Reactive Power: {valdVals['vald_battery_reactive_power_vs_t_1toT_kVAr'][t]} kVAr\n")
                    f.write("*" * 30 + "\n")

        if verbose:
            print(f"Validation key results exported to {filename}")
    except Exception as e:
        if verbose:
            print(f"Error writing to file {filename}: {e}")

