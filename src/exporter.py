import os
import csv
from pyomo.environ import value


def myprintln(verbose, message):
    if verbose:
        print(message)

def create_directory_if_not_exists(directory, verbose=False):
    if not os.path.isdir(directory):
        myprintln(verbose, f"Creating directory: {directory}")
        os.makedirs(directory, exist_ok=True)

def export_optimization_model(model, data, verbose=False):
    T = data['T']
    system_name = data['system_name']
    num_areas = data['num_areas']
    ged_appendix = data['ged_appendix']
    machine_id = data['machine_id']
    objfun_appendix = data['objfun_appendix']
    sim_nature_appendix = data['sim_nature_appendix']

    base_dir = os.path.join(
        "processedData", system_name, ged_appendix, f"Horizon_{T}", f"numAreas_{num_areas}"
    )
    create_directory_if_not_exists(base_dir, verbose)

    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_id}_optimizationModel_{ged_appendix}_for_{objfun_appendix}_via_{sim_nature_appendix}.txt"
    )

    if os.path.isfile(filename):
        os.remove(filename)

    with open(filename, "w") as f:
        f.write(str(model))

    myprintln(verbose, f"Model successfully written to {filename}")

def export_decision_variables(model, data, filename="decision_variables.csv", verbose=False):
    T = data['T']
    system_name = data['system_name']
    num_areas = data['num_areas']
    ged_appendix = data['ged_appendix']
    machine_id = data['machine_id']
    objfun_appendix = data['objfun_appendix']
    sim_nature_appendix = data['sim_nature_appendix']
    solver = data['solver']

    base_dir = os.path.join(
        "processedData", system_name, ged_appendix, f"Horizon_{T}", f"numAreas_{num_areas}"
    )
    create_directory_if_not_exists(base_dir, verbose)

    ext = ".csv"
    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_id}_{solver}_decisionVariables_{ged_appendix}_for_{objfun_appendix}_via_{sim_nature_appendix}{ext}"
    )

    myprintln(verbose, f"Current working directory: {os.getcwd()}")
    myprintln(verbose, f"Saving to filename: {filename}")

    # Unpack data
    Tset = sorted(data['Tset'])
    LoadShapeLoad = data['LoadShapeLoad']
    LoadShapePV = data['LoadShapePV']
    LoadShapeCost = [cost * 100 for cost in data['LoadShapeCost']]  # Convert from $/kWh to cents/kWh

    Lset = sorted(data['Lset'])
    Nset = sorted(data['Nset'])
    Bset = sorted(data['Bset'])
    Dset = sorted(data['Dset'])

    # Model components
    P_Subs = model.P_Subs
    P = model.P
    Q = model.Q
    l = model.l
    v = model.v
    q_D = model.q_D
    q_B = model.q_B
    P_c = model.P_c
    P_d = model.P_d
    B = model.B

    # Prepare data for CSV
    data_matrix = []

    data_matrix.append(["t"] + Tset)
    data_matrix.append(["lambda"] + LoadShapeLoad)
    data_matrix.append(["Irrad"] + LoadShapePV)
    data_matrix.append(["cents/kWh"] + LoadShapeCost)

    data_matrix.append(["P_Subs"] + [P_Subs[t].value for t in Tset])

    for i, j in Lset:
        data_matrix.append([f"P_ij_{i}_{j}"] + [P[i, j, t].value for t in Tset])

    for i, j in Lset:
        data_matrix.append([f"Q_ij_{i}_{j}"] + [Q[i, j, t].value for t in Tset])

    for i, j in Lset:
        data_matrix.append([f"l_ij_{i}_{j}"] + [l[i, j, t].value for t in Tset])

    for j in Nset:
        data_matrix.append([f"v_j_{j}"] + [v[j, t].value for t in Tset])

    for j in Dset:
        data_matrix.append([f"q_D_j_{j}"] + [q_D[j, t].value for t in Tset])

    for j in Bset:
        data_matrix.append([f"q_B_j_{j}"] + [q_B[j, t].value for t in Tset])

    for j in Bset:
        data_matrix.append([f"P_c_j_{j}"] + [P_c[j, t].value for t in Tset])

    for j in Bset:
        data_matrix.append([f"P_d_j_{j}"] + [P_d[j, t].value for t in Tset])

    for j in Bset:
        data_matrix.append([f"B_j_{j}"] + [B[j, t].value for t in Tset])

    try:
        myprintln(verbose, f"Opening file: {filename}")
        with open(filename, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(data_matrix)
        myprintln(verbose, "Data written successfully.")
    except Exception as e:
        myprintln(verbose, f"Error occurred: {e}")

    if os.path.isfile(filename):
        myprintln(verbose, f"File exists: {filename}")
    else:
        myprintln(verbose, f"File does not exist: {filename}")

    myprintln(verbose, f"Decision variables exported to {filename}")

def export_simulation_key_results_txt(model, data, filename="simulation_results.txt", verbose=False):
    T = data['T']
    system_name = data['system_name']
    num_areas = data['num_areas']
    ged_appendix = data['ged_appendix']
    machine_id = data['machine_id']
    objfun_appendix = data['objfun_appendix']
    sim_nature_appendix = data['sim_nature_appendix']
    solver = data['solver']

    base_dir = os.path.join(
        "processedData", system_name, ged_appendix, f"Horizon_{T}", f"numAreas_{num_areas}"
    )
    create_directory_if_not_exists(base_dir, verbose)

    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_id}_{solver}_results_{ged_appendix}_for_{objfun_appendix}_via_{sim_nature_appendix}.txt"
    )

    macro_iters_completed = data.get('macroItrsCompleted', 0)
    solution_time = data.get('solution_time', -1)

    with open(filename, "w") as f:
        item_counter = 1

        # Header Section
        f.write("---------------------------------------------\n")
        f.write(f"{item_counter}. Machine ID: {machine_id}\n")
        item_counter += 1
        f.write(f"{item_counter}. Solver Used: {solver}\n")
        item_counter += 1
        f.write(f"{item_counter}. System Name: {system_name}\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Duration: {T}\n")
        item_counter += 1
        f.write(f"{item_counter}. Nature of Optimization Simulation: {data['sim_nature_string']}\n")
        item_counter += 1
        f.write(f"{item_counter}. Objective: {data['objfun_string']}\n")
        item_counter += 1
        f.write(f"{item_counter}. GED Configuration: {ged_appendix}\n")
        item_counter += 1
        f.write(f"{item_counter}. Maximum Substation Power Allowed: {data['PSubsMax_kW']} kW\n")
        item_counter += 1
        f.write("---------------------------------------------\n")

        # Horizon Results
        f.write(f"Full {T} Hour Horizon\n")
        f.write(f"{item_counter}. Horizon Total Cost of Substation Power: $ {data['PSubsCost_allT_dollar']:.2f}\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Line Loss: {data['PLoss_allT_kW']:.2f} kW\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Substation Power: {data['PSubs_allT_kW']:.2f} kW + {data['QSubs_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Load: {data['load_real_power_allT_kW']:.2f} kW + {data['load_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Generation: {value(data['total_gen_real_power_allT_kW']):.2f} kW + {value(data['total_gen_reactive_power_allT_kVAr']):.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Static Capacitor Reactive Power Generation: {data['static_cap_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total PV Generation: {value(data['pv_real_power_allT_kW']):.2f} kW + {value(data['pv_reactive_power_allT_kVAr']):.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Battery Generation: {value(data['battery_real_power_allT_kW']):.2f} kW + {value(data['battery_reactive_power_allT_kVAr']):.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Battery Transaction Magnitude: {value(data['battery_real_power_transaction_magnitude_allT_kW']):.2f} kW + {value(data['battery_reactive_power_transaction_magnitude_allT_kVAr']):.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total SCD Observed: {data['scd_allT_kW']:.2f} kW\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon-end Battery Energy Deviation from Reference: {data['terminal_soc_violation_kWh']:.2f} kWh\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon-Total All Time Substation Power Peak: {data['substation_real_power_peak_allT_kW']:.2f} kW\n")
        item_counter += 1

        # Additional Simulation Metadata
        f.write("---------------------------------------------\n")
        f.write(f"{item_counter}. Number of Macro-Iterations: {macro_iters_completed+1}\n")
        item_counter += 1
        f.write(f"{item_counter}. Simulation Time: {solution_time:.2f} s\n")
        item_counter += 1
        f.write(f"{item_counter}. Time to solve with sequential (non-parallel) computation: {solution_time:.2f} s\n")
        item_counter += 1
        f.write(f"{item_counter}. Time to solve if OPF computation parallelized: {solution_time:.2f} s\n")
        f.write("---------------------------------------------\n")

    myprintln(verbose, f"Simulation key results exported to {filename}")

def export_key_validation_results(vald, data, filename="validation_results.txt", verbose=False, print_every_time_step_powerflow=True):
    T = data['T']
    system_name = data['system_name']
    num_areas = data['num_areas']
    ged_appendix = data['ged_appendix']
    machine_id = data['machine_id']
    sim_nature_appendix = data['sim_nature_appendix']
    solver = data['solver']

    base_dir = os.path.join(
        "processedData", system_name, ged_appendix, f"Horizon_{T}", f"numAreas_{num_areas}"
    )
    create_directory_if_not_exists(base_dir, verbose)

    if print_every_time_step_powerflow:
        filename = os.path.join(
            base_dir,
            f"Horizon_{T}_{machine_id}_{solver}_valdResults_{ged_appendix}_via_{sim_nature_appendix}_full.txt"
        )
    else:
        filename = os.path.join(
            base_dir,
            f"Horizon_{T}_{machine_id}_{solver}_valdResults_{ged_appendix}_via_{sim_nature_appendix}.txt"
        )

    with open(filename, "w") as f:
        item_counter = 1

        # Header Section
        f.write("---------------------------------------------\n")
        f.write(f"{item_counter}. Machine ID: {machine_id}\n")
        item_counter += 1
        f.write(f"{item_counter}. Solver Used: {solver}\n")
        item_counter += 1
        f.write(f"{item_counter}. System Name: {system_name}\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Duration: {T}\n")
        item_counter += 1
        f.write(f"{item_counter}. Nature of Validation Simulation: {data['simNatureString']}\n")
        item_counter += 1
        f.write("---------------------------------------------\n")

        # Horizon Results
        f.write(f"Full {T} Hour Horizon Validation Results\n")
        f.write(f"{item_counter}. Horizon Total Substation Power Cost: ${vald['vald_PSubsCost_allT_dollar']:.2f}\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Line Loss: {vald['vald_PLoss_allT_kW']:.2f} kW\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Substation Power: {vald['vald_PSubs_allT_kW']:.2f} kW + {vald['vald_QSubs_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Load: {vald['vald_load_real_power_allT_kW']:.2f} kW + {vald['vald_load_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Generation: {vald['vald_total_gen_real_power_allT_kW']:.2f} kW + {vald['vald_total_gen_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Static Capacitor Reactive Power Generation: {vald['vald_static_cap_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total PV Generation: {vald['vald_pv_real_power_allT_kW']:.2f} kW + {vald['vald_pv_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Battery Generation: {vald['vald_battery_real_power_allT_kW']:.2f} kW + {vald['vald_battery_reactive_power_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total Battery Transaction Magnitude: {vald['vald_battery_real_power_transaction_magnitude_allT_kW']:.2f} kW + {vald['vald_battery_reactive_power_transaction_magnitude_allT_kVAr']:.2f} kVAr\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon Total SCD Observed: N/A\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon-end Battery Energy Deviation from Reference: {vald['vald_terminal_soc_violation_kWh']:.2f} kWh\n")
        item_counter += 1
        f.write(f"{item_counter}. Horizon-Total All Time Substation Power Peak: {vald['vald_substation_real_power_peak_allT_kW']:.2f} kW\n")
        item_counter += 1

        # Discrepancies
        f.write("---------------------------------------------\n")
        f.write("Discrepancies (Maximum All Time):\n")
        f.write(f"{item_counter}. Maximum All Time Voltage Discrepancy: {vald['disc_voltage_all_time_pu']:.6f} pu\n")
        item_counter += 1
        f.write(f"{item_counter}. Maximum All Time Line Loss Discrepancy: {vald['disc_line_loss_all_time_kW']:.6f} kW\n")
        item_counter += 1
        f.write(f"{item_counter}. Maximum All Time Substation Borrowed Real Power Discrepancy: {vald['disc_PSubs_all_time_kW']:.6f} kW\n")
        item_counter += 1
        f.write(f"{item_counter}. Maximum All Time Substation Borrowed Reactive Power Discrepancy: {vald['disc_QSubs_all_time_kVAr']:.6f} kVAr\n")
        item_counter += 1

        # Additional Metadata
        f.write("---------------------------------------------\n")
        f.write(f"{item_counter}. Solution Time: Small\n")
        item_counter += 1

        # Optional: Print per-timestep power flow results if enabled
        if print_every_time_step_powerflow:
            f.write("\nPer-Timestep Power Flow Results:\n")
            for t in data['Tset']:
                f.write("\n" + "*" * 30 + "\n")
                f.write(f"   Time Step: {t}\n")
                f.write("*" * 30 + "\n")
                f.write(f"   Power Loss              : {vald['vald_PLoss_vs_t_1toT_kW'][t]:.2f} kW\n")
                f.write(f"   Substation Power        : {vald['vald_PSubs_vs_t_1toT_kW'][t]:.2f} kW\n")
                f.write(f"   Reactive Power          : {vald['vald_QSubs_vs_t_1toT_kVAr'][t]:.2f} kVAr\n")
                f.write(f"   Total Load Power        : {vald['vald_load_real_power_vs_t_1toT_kW'][t]:.2f} kW\n")
                f.write(f"   Total Load Reactive Power: {vald['vald_load_reactive_power_vs_t_1toT_kVAr'][t]:.2f} kVAr\n")
                f.write(f"   Total PV Power          : {vald['vald_pv_real_power_vs_t_1toT_kW'][t]:.2f} kW\n")
                f.write(f"   Total PV Reactive Power : {vald['vald_pv_reactive_power_vs_t_1toT_kVAr'][t]:.2f} kVAr\n")
                f.write(f"   Total Battery Power     : {vald['vald_battery_real_power_vs_t_1toT_kW'][t]:.2f} kW\n")
                f.write(f"   Total Battery Reactive Power: {vald['vald_battery_reactive_power_vs_t_1toT_kVAr'][t]:.2f} kVAr\n")
                f.write("*" * 30 + "\n")

    myprintln(verbose, f"Simulation key results exported to {filename}")

