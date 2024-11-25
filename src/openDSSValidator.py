__all__ = [
    "export_validation_decision_variables",
    "get_load_powers_opendss_powerflow_for_timestep_t",
    "get_source_bus",
    "get_substation_lines",
    "get_substation_powers_opendss_powerflow_for_timestep_t",
    "get_voltages_opendss_powerflow_for_timestep_t",
    "set_custom_load_shape",
    "set_battery_controls_opendss_powerflow_for_timestep_t",
    "set_pv_controls_opendss_powerflow_for_timestep_t",
]

import os
import pandas as pd
from opendssdirect import Circuit, Solution, Text, Vsources, CktElement, Lines, Loads, Storages, PVsystems

# Utility functions
def export_validation_decision_variables(vald, data, verbose=False):
    T = data['T']
    system_name = data['systemName']
    num_areas = data['numAreas']
    ged_appendix = data['gedAppendix']
    machine_id = data['machine_ID']
    objfun_description = data['objfunConciseDescription']
    processed_data_path = data['processedDataFolderPath']
    sim_nature_appendix = data['simNatureAppendix']
    solver = data['solver']

    base_dir = os.path.join(
        processed_data_path, system_name, ged_appendix, f"Horizon_{T}", f"numAreas_{num_areas}"
    )
    if not os.path.isdir(base_dir):
        if verbose:
            print(f"Creating directory: {base_dir}")
        os.makedirs(base_dir)

    filename = os.path.join(
        base_dir,
        f"Horizon_{T}_{machine_id}_{solver}_validationDecisionVariables_{ged_appendix}_for_{objfun_description}_via_{sim_nature_appendix}.txt",
    )

    pd.DataFrame(vald).to_csv(filename, index=False)

    if verbose:
        print(f"Validation decision variables written to {filename}")

def get_load_powers_opendss_powerflow_for_timestep_t():
    total_load_kW = 0.0
    total_load_kVAr = 0.0

    load_names = Loads.AllNames()
    for load_name in load_names:
        Circuit.SetActiveElement(f"Load.{load_name}")
        load_powers = CktElement.Powers()
        total_load_kW += load_powers[0].real
        total_load_kVAr += load_powers[0].imag

    return {
        "total_load_t_kW": total_load_kW,
        "total_load_t_kVAr": total_load_kVAr,
    }

def get_source_bus():
    Vsources.First()
    vsource_name = Vsources.Name()
    Circuit.SetActiveElement(f"Vsource.{vsource_name}")
    source_bus = CktElement.BusNames()[0]
    return source_bus

def get_substation_lines(substation_bus):
    substation_lines = []

    line_id = Lines.First()
    while line_id > 0:
        line_name = Lines.Name()
        Circuit.SetActiveElement(f"Line.{line_name}")
        bus1, bus2 = CktElement.BusNames()[:2]

        if bus1 == substation_bus or bus2 == substation_bus:
            substation_lines.append(line_name)

        line_id = Lines.Next()

    return substation_lines

def get_substation_powers_opendss_powerflow_for_timestep_t(data, use_vsource_power=True):
    substation_bus = data['substationBus']
    P_substation_total_kW = 0.0
    Q_substation_total_kVAr = 0.0

    if use_vsource_power:
        Circuit.SetActiveElement("Vsource.source")
        vsource_powers = CktElement.Powers()
        P_substation_total_kW = -vsource_powers[0].real
        Q_substation_total_kVAr = -vsource_powers[0].imag
    else:
        substation_lines = get_substation_lines(substation_bus)
        for line in substation_lines:
            Circuit.SetActiveElement(f"Line.{line}")
            line_powers = CktElement.Powers()
            P_substation_total_kW += sum(line_powers[0::2])
            Q_substation_total_kVAr += sum(line_powers[1::2])

    return {
        "P_substation_total_t_kW": P_substation_total_kW,
        "Q_substation_total_t_kVAr": Q_substation_total_kVAr,
    }

def get_voltages_opendss_powerflow_for_timestep_t():
    bus_names = Circuit.AllBusNames()
    bus_voltages = Circuit.AllBusMagPu()

    vald_voltage_dict = {
        int(bus): voltage for bus, voltage in zip(bus_names, bus_voltages)
    }
    return vald_voltage_dict

def set_custom_load_shape(load_shape_array, verbose=False):
    loadshape_command = (
        f"New Loadshape.LoadShapeLoad npts={len(load_shape_array)} interval=1 mult=["
        + " ".join(map(str, load_shape_array))
        + "]"
    )
    Text.Command(loadshape_command)
    if verbose:
        print("Defined LoadShapeLoad with provided LoadShapeArray")

    load_id = Loads.First()
    while load_id > 0:
        Loads.Daily("LoadShapeLoad")
        load_id = Loads.Next()
    if verbose:
        print("Applied LoadShapeLoad to all loads")

def set_battery_controls_opendss_powerflow_for_timestep_t(model, data, t, verbose=False):
    P_c = model['P_c']
    P_d = model['P_d']
    q_B = model['q_B']
    kVA_B = data['kVA_B']

    storage_id = Storages.First()
    while storage_id > 0:
        storage_name = Storages.Name()
        storage_number = int(storage_name.split("battery")[1])

        charge_power_kW = P_c[storage_number][t] * kVA_B
        discharge_power_kW = P_d[storage_number][t] * kVA_B
        net_power_kW = discharge_power_kW - charge_power_kW
        reactive_power_kVAr = q_B[storage_number][t] * kVA_B

        command_str = f"Edit Storage.Battery{storage_number} kW={net_power_kW} kvar={reactive_power_kVAr}"
        Text.Command(command_str)

        if verbose:
            print(f"Time Step {t}: Setting battery {storage_number} with command: {command_str}")

        storage_id = Storages.Next()

def set_pv_controls_opendss_powerflow_for_timestep_t(model, data, t, verbose=False):
    kVA_B = data['kVA_B']
    p_D_pu = data['p_D_pu']
    q_D = model['q_D']

    pv_id = PVsystems.First()
    while pv_id > 0:
        pv_name = PVsystems.Name()
        pv_number = int(pv_name.split("pv")[1])

        p_D_t_kW = p_D_pu[pv_number][t] * kVA_B
        q_D_t_kVAr = q_D[pv_number][t] * kVA_B

        PVsystems.Pmpp(p_D_t_kW)
        PVsystems.kvar(q_D_t_kVAr)

        if verbose:
            print(f"Setting PV for bus {pv_number} at t = {t}: p_D_t_kW = {p_D_t_kW}, q_D_t_kVAr = {q_D_t_kVAr}")

        pv_id = PVsystems.Next()
