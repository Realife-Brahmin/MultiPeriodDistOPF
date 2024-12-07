
# openDSSValidator.py

import os
import numpy as np
import opendssdirect as dss
from src.helperFunctions import myprintln  # Assuming it exists and behaves like Julia's myprintln

def compute_highest_allTime_voltage_discrepancy(modelDict, valdVals):
    modelVals = modelDict['modelVals']
    data = modelDict['data']

    disc_voltage_all_time_pu = 0.0
    v = modelVals['v']
    Tset = data['Tset']

    # Tset assumed 1-based in Julia. If Tset is 1-based, we adjust indexing carefully.
    # We'll assume Tset is 1-based for consistency. vald_voltages_vs_t_1toT_pu is zero-based in Python arrays.
    for t in Tset:
        # t-1 indexing since Python arrays are zero-based
        vald_voltages_dict = valdVals['vald_voltages_vs_t_1toT_pu'][t-1]
        for bus_index, vald_voltage in vald_voltages_dict.items():
            model_voltage = np.sqrt(v[(bus_index, t)])
            discrepancy = abs(vald_voltage - model_voltage)
            if discrepancy > disc_voltage_all_time_pu:
                disc_voltage_all_time_pu = discrepancy

    return disc_voltage_all_time_pu

def get_battery_powers_opendss_powerflow_for_timestep_t(verbose=False):
    vald_battery_real_power_t_kW = 0.0
    vald_battery_reactive_power_t_kVAr = 0.0
    vald_battery_real_power_transaction_magnitude_t_kW = 0.0
    vald_battery_reactive_power_transaction_magnitude_t_kVAr = 0.0

    battery_names = dss.Storages.AllNames()
    for battery_name in battery_names:
        dss.Circuit.SetActiveElement(f"Storage.{battery_name}")
        battery_powers = dss.CktElement.Powers()
        # Powers format: [P1,Q1,P2,Q2,...], first two entries for first terminal
        real_power = -battery_powers[0]  # Negate to match injection convention
        reactive_power = -battery_powers[1]

        vald_battery_real_power_t_kW += real_power
        vald_battery_reactive_power_t_kVAr += reactive_power
        vald_battery_real_power_transaction_magnitude_t_kW += abs(real_power)
        vald_battery_reactive_power_transaction_magnitude_t_kVAr += abs(reactive_power)

        myprintln(verbose, f"Battery {battery_name} | Real Power: {real_power} kW, Reactive Power: {reactive_power} kVAr")

    batteryPowersDict_t = {
        'vald_battery_real_power_t_kW': vald_battery_real_power_t_kW,
        'vald_battery_reactive_power_t_kVAr': vald_battery_reactive_power_t_kVAr,
        'vald_battery_real_power_transaction_magnitude_t_kW': vald_battery_real_power_transaction_magnitude_t_kW,
        'vald_battery_reactive_power_transaction_magnitude_t_kVAr': vald_battery_reactive_power_transaction_magnitude_t_kVAr
    }

    return batteryPowersDict_t

def get_load_powers_opendss_powerflow_for_timestep_t():
    total_load_t_kW = 0.0
    total_load_t_kVAr = 0.0

    load_names = dss.Loads.AllNames()
    for load_name in load_names:
        dss.Circuit.SetActiveElement(f"Load.{load_name}")
        load_powers = dss.CktElement.Powers()
        total_load_t_kW += load_powers[0]
        total_load_t_kVAr += load_powers[1]

    loadPowersDict_t = {
        'total_load_t_kW': total_load_t_kW,
        'total_load_t_kVAr': total_load_t_kVAr
    }

    return loadPowersDict_t

def get_pv_powers_opendss_powerflow_for_timestep_t():
    total_pv_t_kW = 0.0
    total_pv_t_kVAr = 0.0

    pv_names = dss.PVsystems.AllNames()
    for pv_name in pv_names:
        dss.Circuit.SetActiveElement(f"PVSystem.{pv_name}")
        pv_powers = dss.CktElement.Powers()
        total_pv_t_kW -= pv_powers[0]
        total_pv_t_kVAr -= pv_powers[1]

    pvPowersDict_t = {
        'total_pv_t_kW': total_pv_t_kW,
        'total_pv_t_kVAr': total_pv_t_kVAr
    }

    return pvPowersDict_t

def get_source_bus():
    dss.Vsources.First()
    vsource_name = dss.Vsources.Name()
    dss.Circuit.SetActiveElement(f"Vsource.{vsource_name}")
    source_bus = dss.CktElement.BusNames()[0]
    return source_bus

def get_substation_lines(substation_bus:str):
    substation_lines = []
    line_id = dss.Lines.First()
    while line_id > 0:
        line_name = dss.Lines.Name()
        dss.Circuit.SetActiveElement(f"Line.{line_name}")
        busnames = dss.CktElement.BusNames()
        bus1, bus2 = busnames[0], busnames[1]
        if bus1 == substation_bus or bus2 == substation_bus:
            substation_lines.append(line_name)
        line_id = dss.Lines.Next()
    return substation_lines

def get_substation_powers_opendss_powerflow_for_timestep_t(data, useVSourcePower=True):
    substationBus = data['substationBus']
    P_substation_total_t_kW = 0.0
    Q_substation_total_t_kVAr = 0.0

    if useVSourcePower:
        dss.Circuit.SetActiveElement("Vsource.source")
        vsource_powers = dss.CktElement.Powers()
        # negate the entire array
        vsource_powers = [-x for x in vsource_powers]
        P_substation_total_t_kW = vsource_powers[0]
        Q_substation_total_t_kVAr = vsource_powers[1]
    else:
        s_lines = get_substation_lines(substationBus)
        for line in s_lines:
            dss.Circuit.SetActiveElement(f"Line.{line}")
            line_powers = dss.CktElement.Powers()
            P_line = line_powers[0]
            Q_line = line_powers[1]
            P_substation_total_t_kW += P_line
            Q_substation_total_t_kVAr += Q_line

    return {
        'P_substation_total_t_kW': P_substation_total_t_kW,
        'Q_substation_total_t_kVAr': Q_substation_total_t_kVAr
    }

def get_terminal_soc_values_opendss_powerflow(data):
    Bref_percent = data['Bref_percent']
    B_R = data['B_R']

    terminalSOCDict = {
        'vald_Bj_T_pu': {},
        'soc_violation_j_kWh': {},
        'vald_terminal_soc_violation_kWh': 0.0
    }

    storage_id = dss.Storages.First()
    while storage_id > 0:
        storage_name = dss.Storages.Name()
        storage_number = int(storage_name.replace("battery",""))
        vald_Bj_T_pu = dss.Storages.puSOC()
        terminalSOCDict['vald_Bj_T_pu'][storage_number] = vald_Bj_T_pu
        soc_violation = abs(vald_Bj_T_pu - Bref_percent[storage_number])*B_R[storage_number]
        terminalSOCDict['soc_violation_j_kWh'][storage_number] = soc_violation
        terminalSOCDict['vald_terminal_soc_violation_kWh'] += soc_violation

        storage_id = dss.Storages.Next()

    return terminalSOCDict

def get_voltages_opendss_powerflow_for_timestep_t():
    vald_voltage_dict_t_pu = {}
    bus_names = dss.Circuit.AllBusNames()
    bus_voltages = dss.Circuit.AllBusMagPu()
    for i, bus_name in enumerate(bus_names):
        bus_number = int(bus_name)
        vald_voltage_dict_t_pu[bus_number] = bus_voltages[i]
    return vald_voltage_dict_t_pu

def set_custom_load_shape_(LoadShapeArray, verbose=False):
    loadshape_command = "New Loadshape.LoadShapeLoad npts={} interval=1 mult=[{}]".format(len(LoadShapeArray), " ".join(map(str, LoadShapeArray)))
    dss.Text.Command(loadshape_command)
    myprintln(verbose, "Defined LoadShapeLoad with provided LoadShapeArray")

    load_id = dss.Loads.First()
    while load_id > 0:
        dss.Loads.Daily("LoadShapeLoad")
        load_id = dss.Loads.Next()
    myprintln(verbose, "Applied LoadShapeLoad to all loads")

def set_battery_controls_opendss_powerflow_for_timestep_t(modelDict, t, verbose=False):
    modelVals = modelDict['modelVals']
    data = modelDict['data']

    P_c = modelVals['P_c']
    P_d = modelVals['P_d']
    q_B = modelVals['q_B']
    kVA_B = data['kVA_B']

    storage_id = dss.Storages.First()
    while storage_id > 0:
        storage_name = dss.Storages.Name()
        storage_number = int(storage_name.replace("battery",""))
        charge_power_kW = P_c[(storage_number, t)] * kVA_B
        discharge_power_kW = P_d[(storage_number, t)] * kVA_B
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = q_B[(storage_number, t)] * kVA_B

        command_str = f"Edit Storage.Battery{storage_number} kW={net_power_kW} kvar={reactive_power_kVAr}"
        dss.Text.Command(command_str)
        if verbose:
            print(f"Time Step {t}: Setting battery {storage_number} with command: {command_str}")

        storage_id = dss.Storages.Next()

def set_pv_controls_opendss_powerflow_for_timestep_t(modelDict, t, verbose=False):
    modelVals = modelDict['modelVals']
    data = modelDict['data']
    kVA_B = data['kVA_B']
    p_D_pu = data['p_D_pu']
    q_D = modelVals['q_D']

    pv_id = dss.PVsystems.First()
    while pv_id > 0:
        pv_name = dss.PVsystems.Name()
        pv_number = int(pv_name.replace("pv",""))
        p_D_t_kW = p_D_pu[(pv_number, t)]*kVA_B
        q_D_t_kVAr = q_D[(pv_number, t)]*kVA_B

        # Setting P and Q for PV
        dss.PVsystems.Pmpp(p_D_t_kW)
        dss.PVsystems.kvar(q_D_t_kVAr)

        if verbose:
            print(f"Setting PV for bus {pv_number} at t={t}: p_D_t_kW={p_D_t_kW}, q_D_t_kVAr={q_D_t_kVAr}")

        pv_id = dss.PVsystems.Next()

def validate_opf_against_opendss(modelDict, verbose=False):
    modelVals = modelDict['modelVals']
    data = modelDict['data']

    T = data['T']
    LoadShapeLoad = data['LoadShapeLoad']

    # Initialize valdVals
    valdVals = {
        'vald_PLoss_vs_t_1toT_kW': np.zeros(T),
        'vald_PSubs_vs_t_1toT_kW': np.zeros(T),
        'vald_QLoss_vs_t_1toT_kVAr': np.zeros(T),
        'vald_QSubs_vs_t_1toT_kVAr': np.zeros(T),
        'vald_load_real_power_vs_t_1toT_kW': np.zeros(T),
        'vald_load_reactive_power_vs_t_1toT_kVAr': np.zeros(T),
        'vald_pv_real_power_vs_t_1toT_kW': np.zeros(T),
        'vald_pv_reactive_power_vs_t_1toT_kVAr': np.zeros(T),
        'vald_battery_real_power_vs_t_1toT_kW': np.zeros(T),
        'vald_battery_reactive_power_vs_t_1toT_kVAr': np.zeros(T),
        'vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW': np.zeros(T),
        'vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr': np.zeros(T),
        'vald_static_cap_reactive_power_vs_t_1toT_kVAr': np.zeros(T),
        'vald_total_gen_reactive_power_vs_t_1toT_kVAr': np.zeros(T),
        'vald_total_gen_real_power_vs_t_1toT_kW': np.zeros(T),
        'vald_voltages_vs_t_1toT_pu': [None]*T,
        'vald_PSubsCost_vs_t_1toT_dollar': np.zeros(T),
        'vald_battery_reactive_power_allT_kVAr': 0.0,
        'vald_battery_reactive_power_transaction_magnitude_allT_kVAr': 0.0,
        'vald_battery_real_power_allT_kW': 0.0,
        'vald_battery_real_power_transaction_magnitude_allT_kW': 0.0,
        'vald_pv_reactive_power_allT_kVAr': 0.0,
        'vald_pv_real_power_allT_kW': 0.0,
        'vald_PLoss_allT_kW': 0.0,
        'vald_PSubs_allT_kW': 0.0,
        'vald_QLoss_allT_kVAr': 0.0,
        'vald_QSubs_allT_kVAr': 0.0,
        'vald_static_cap_reactive_power_allT_kVAr': 0.0,
        'vald_total_gen_real_power_allT_kW': 0.0,
        'vald_total_gen_reactive_power_allT_kVAr': 0.0,
        'vald_load_real_power_allT_kW': 0.0,
        'vald_load_reactive_power_allT_kVAr': 0.0,
        'vald_PSubsCost_allT_dollar': 0.0,
        'vald_solution_time': 0.0,
        'vald_substation_real_power_peak_allT_kW': 0.0,
        'vald_terminal_soc_violation_kWh': 0.0
    }

    systemName = data['systemName']
    rawDataFolderPath = data['rawDataFolderPath']
    dss_dir = os.path.join(rawDataFolderPath, systemName)
    dss_file = os.path.join(dss_dir, "Master.dss")
    myprintln(verbose, f"Master.dss file path: {dss_file}")

    dss.Text.Command("Clear")
    dss.Text.Command(f"Redirect \"{dss_file}\"")

    set_custom_load_shape_(LoadShapeLoad, verbose=verbose)

    for t in range(1, T+1):
        set_pv_controls_opendss_powerflow_for_timestep_t(modelDict, t, verbose=verbose)
        set_battery_controls_opendss_powerflow_for_timestep_t(modelDict, t, verbose=verbose)
        dss.Solution.Solve()

        totalLosses = dss.Circuit.Losses() # returns [Watts, var]
        totalLosses_t_kW = totalLosses[0]/1000.0
        totalLosses_t_kVAr = -totalLosses[1]/1000.0
        valdVals['vald_PLoss_vs_t_1toT_kW'][t-1] = totalLosses_t_kW
        valdVals['vald_PLoss_allT_kW'] += totalLosses_t_kW
        valdVals['vald_QLoss_vs_t_1toT_kVAr'][t-1] = totalLosses_t_kVAr
        valdVals['vald_QLoss_allT_kVAr'] += totalLosses_t_kVAr

        substationPowersDict_t = get_substation_powers_opendss_powerflow_for_timestep_t(data, useVSourcePower=True)
        P_substation_total_t_kW = substationPowersDict_t['P_substation_total_t_kW']
        Q_substation_total_t_kVAr = substationPowersDict_t['Q_substation_total_t_kVAr']

        valdVals['vald_PSubs_vs_t_1toT_kW'][t-1] = P_substation_total_t_kW
        valdVals['vald_PSubs_allT_kW'] += P_substation_total_t_kW
        valdVals['vald_QSubs_vs_t_1toT_kVAr'][t-1] = Q_substation_total_t_kVAr
        valdVals['vald_QSubs_allT_kVAr'] += Q_substation_total_t_kVAr
        if P_substation_total_t_kW > valdVals['vald_substation_real_power_peak_allT_kW']:
            valdVals['vald_substation_real_power_peak_allT_kW'] = P_substation_total_t_kW

        # Cost calculation
        LoadShapeCost = data['LoadShapeCost']
        delta_t = data['delta_t']
        valdVals['vald_PSubsCost_vs_t_1toT_dollar'][t-1] = LoadShapeCost[t-1]*P_substation_total_t_kW*delta_t
        valdVals['vald_PSubsCost_allT_dollar'] += valdVals['vald_PSubsCost_vs_t_1toT_dollar'][t-1]

        loadPowersDict_t = get_load_powers_opendss_powerflow_for_timestep_t()
        total_load_t_kW = loadPowersDict_t['total_load_t_kW']
        total_load_t_kVAr = loadPowersDict_t['total_load_t_kVAr']
        valdVals['vald_load_real_power_vs_t_1toT_kW'][t-1] = total_load_t_kW
        valdVals['vald_load_reactive_power_vs_t_1toT_kVAr'][t-1] = total_load_t_kVAr
        valdVals['vald_load_real_power_allT_kW'] += total_load_t_kW
        valdVals['vald_load_reactive_power_allT_kVAr'] += total_load_t_kVAr

        pvPowersDict_t = get_pv_powers_opendss_powerflow_for_timestep_t()
        total_pv_t_kW = pvPowersDict_t['total_pv_t_kW']
        total_pv_t_kVAr = pvPowersDict_t['total_pv_t_kVAr']
        valdVals['vald_pv_real_power_vs_t_1toT_kW'][t-1] = total_pv_t_kW
        valdVals['vald_pv_real_power_allT_kW'] += total_pv_t_kW
        valdVals['vald_pv_reactive_power_vs_t_1toT_kVAr'][t-1] = total_pv_t_kVAr
        valdVals['vald_pv_reactive_power_allT_kVAr'] += total_pv_t_kVAr

        batteryPowersDict_t = get_battery_powers_opendss_powerflow_for_timestep_t(verbose=verbose)
        vald_battery_real_power_t_kW = batteryPowersDict_t['vald_battery_real_power_t_kW']
        vald_battery_reactive_power_t_kVAr = batteryPowersDict_t['vald_battery_reactive_power_t_kVAr']
        vald_battery_real_power_transaction_magnitude_t_kW = batteryPowersDict_t['vald_battery_real_power_transaction_magnitude_t_kW']
        vald_battery_reactive_power_transaction_magnitude_t_kVAr = batteryPowersDict_t['vald_battery_reactive_power_transaction_magnitude_t_kVAr']

        valdVals['vald_battery_reactive_power_vs_t_1toT_kVAr'][t-1] = vald_battery_reactive_power_t_kVAr
        valdVals['vald_battery_reactive_power_allT_kVAr'] += vald_battery_reactive_power_t_kVAr
        valdVals['vald_battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr'][t-1] = vald_battery_reactive_power_transaction_magnitude_t_kVAr
        valdVals['vald_battery_reactive_power_transaction_magnitude_allT_kVAr'] += vald_battery_reactive_power_transaction_magnitude_t_kVAr

        valdVals['vald_battery_real_power_vs_t_1toT_kW'][t-1] = vald_battery_real_power_t_kW
        valdVals['vald_battery_real_power_allT_kW'] += vald_battery_real_power_t_kW
        valdVals['vald_battery_real_power_transaction_magnitude_vs_t_1toT_kW'][t-1] = vald_battery_real_power_transaction_magnitude_t_kW
        valdVals['vald_battery_real_power_transaction_magnitude_allT_kW'] += vald_battery_real_power_transaction_magnitude_t_kW

        valdVals['vald_total_gen_reactive_power_vs_t_1toT_kVAr'][t-1] = valdVals['vald_pv_reactive_power_vs_t_1toT_kVAr'][t-1] + valdVals['vald_battery_reactive_power_vs_t_1toT_kVAr'][t-1]
        valdVals['vald_total_gen_reactive_power_allT_kVAr'] += valdVals['vald_total_gen_reactive_power_vs_t_1toT_kVAr'][t-1]

        valdVals['vald_total_gen_real_power_vs_t_1toT_kW'][t-1] = valdVals['vald_pv_real_power_vs_t_1toT_kW'][t-1] + valdVals['vald_battery_real_power_vs_t_1toT_kW'][t-1]
        valdVals['vald_total_gen_real_power_allT_kW'] += valdVals['vald_total_gen_real_power_vs_t_1toT_kW'][t-1]

        valdVals['vald_voltages_vs_t_1toT_pu'][t-1] = get_voltages_opendss_powerflow_for_timestep_t()

    terminalSOCDict = get_terminal_soc_values_opendss_powerflow(data)
    valdVals['vald_terminal_soc_violation_kWh'] = terminalSOCDict['vald_terminal_soc_violation_kWh']

    # Discrepancy calculations
    disc_voltage_all_time_pu = compute_highest_allTime_voltage_discrepancy(modelDict, valdVals)

    PLoss_vs_t_model = valdVals['vald_PLoss_vs_t_1toT_kW']
    PLoss_vs_t_data = data['PLoss_vs_t_1toT_kW']
    line_loss_discrepancies = np.abs(PLoss_vs_t_model - PLoss_vs_t_data)
    disc_line_loss_all_time_kW = np.max(line_loss_discrepancies)

    PSubs_vs_t_model = valdVals['vald_PSubs_vs_t_1toT_kW']
    PSubs_vs_t_data = data['PSubs_vs_t_1toT_kW']
    disc_PSubs_all_time_kW = np.max(np.abs(PSubs_vs_t_model - PSubs_vs_t_data))

    QSubs_vs_t_model = valdVals['vald_QSubs_vs_t_1toT_kVAr']
    QSubs_vs_t_data = data['QSubs_vs_t_1toT_kVAr']
    disc_QSubs_all_time_kVAr = np.max(np.abs(QSubs_vs_t_model - QSubs_vs_t_data))

    valdVals['disc_voltage_all_time_pu'] = disc_voltage_all_time_pu
    valdVals['disc_line_loss_all_time_kW'] = disc_line_loss_all_time_kW
    valdVals['disc_PSubs_all_time_kW'] = disc_PSubs_all_time_kW
    valdVals['disc_QSubs_all_time_kVAr'] = disc_QSubs_all_time_kVAr

    modelDict['valdVals'] = valdVals
    return modelDict

def export_validation_decision_variables(modelDict, verbose=False):
    # Placeholder for function to export validation decision variables
    # Implement file writing logic as needed
    myprintln(verbose, "Exporting validation decision variables...")




