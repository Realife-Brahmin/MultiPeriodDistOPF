# functionRetriever.py

import numpy as np


def get_loss_real_power(model, data, horizon="allT"):
    Tset = data['Tset']
    Lset = data['Lset']
    rdict_pu = data['rdict_pu']
    kVA_B = data['kVA_B']
    l = model.l

    if horizon == "1toT":
        return [
            kVA_B * sum(rdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(
            rdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset for t in Tset
        )
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_loss_reactive_power(model, data, horizon="allT"):
    Tset = data['Tset']
    Lset = data['Lset']
    xdict_pu = data['xdict_pu']
    kVA_B = data['kVA_B']
    l = model.l

    if horizon == "1toT":
        return [
            kVA_B * sum(xdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(
            xdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset for t in Tset
        )
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_substation_reactive_power(model, data, horizon="allT"):
    Tset = data['Tset']
    L1set = data['L1set']
    kVA_B = data['kVA_B']
    Q = model.Q

    if horizon == "1toT":
        return [
            kVA_B * sum(Q[(i, j), t].value for (i, j) in L1set)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(Q[(i, j), t].value for (i, j) in L1set for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_substation_real_power(model, data, horizon="allT"):
    Tset = data['Tset']
    kVA_B = data['kVA_B']
    P_Subs = model.P_Subs

    if horizon == "1toT":
        return [P_Subs[t].value * kVA_B for t in Tset]
    elif horizon == "allT":
        return sum(P_Subs[t].value * kVA_B for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_substation_power_cost(model, data, horizon="allT"):
    Tset = data['Tset']
    LoadShapeCost = data['LoadShapeCost']
    delta_t = data['delta_t']
    kVA_B = data['kVA_B']
    P_Subs = model.P_Subs

    if horizon == "1toT":
        return [
            LoadShapeCost[t-1] * P_Subs[t].value * kVA_B * delta_t for t in Tset
        ]
    elif horizon == "allT":
        return sum(
            LoadShapeCost[t-1] * P_Subs[t].value * kVA_B * delta_t for t in Tset
        )
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_scd(model, data, horizon="allT"):
    Tset = data['Tset']
    Bset = data['Bset']
    kVA_B = data['kVA_B']
    P_c = model.P_c
    P_d = model.P_d

    if horizon == "1toT":
        return [
            kVA_B * sum(
                max(P_c[j, t].value, P_d[j, t].value)
                - abs(P_c[j, t].value - P_d[j, t].value)
                for j in Bset
            )
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(
            max(P_c[j, t].value, P_d[j, t].value)
            - abs(P_c[j, t].value - P_d[j, t].value)
            for j in Bset for t in Tset
        )
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_terminal_SOC_violation(model, data):
    T = data['T']
    Bset = data['Bset']
    Bref_pu = data['Bref_pu']
    kVA_B = data['kVA_B']
    B = model.B

    return kVA_B * sum(abs(B[j, T].value - Bref_pu[j]) for j in Bset)

def get_variable(model, variable_name):
    return model[variable_name]

def get_battery_real_power(model, data, horizon="allT"):
    Tset = data['Tset']
    Bset = data['Bset']
    kVA_B = data['kVA_B']
    P_c = model.P_c
    P_d = model.P_d

    if horizon == "1toT":
        return [
            kVA_B * sum(P_d[j, t].value - P_c[j, t].value for j in Bset)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(
            P_d[j, t].value - P_c[j, t].value for j in Bset for t in Tset
        )
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_battery_reactive_power(model, data, horizon="allT"):
    Tset = data['Tset']
    Bset = data['Bset']
    kVA_B = data['kVA_B']
    q_B = model.q_B

    if horizon == "1toT":
        return [
            kVA_B * sum(q_B[j, t].value for j in Bset)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(q_B[j, t].value for j in Bset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

def get_pv_real_power(model, data, horizon="allT"):
    Tset = data['Tset']
    Dset = data['Dset']
    p_D_pu = data['p_D_pu']
    kVA_B = data['kVA_B']

    if horizon == "1toT":
        return [
            kVA_B * sum(p_D_pu[j][t-1] for j in Dset)
            for t in Tset
        ]
    elif horizon == "allT":
        return kVA_B * sum(p_D_pu[j][t-1] for j in Dset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to get total PV reactive power in kVAr (q_D from the model)
def get_pv_reactive_power(model, data, horizon="allT"):
    Tset, Dset, kVA_B = data['Tset'], data['Dset'], data['kVA_B']
    q_D = model.q_D

    if horizon == "1toT":
        pv_reactive_power_vs_t_1toT_kVAr = [kVA_B * sum(q_D[j, t] for j in Dset) for t in Tset]
        return pv_reactive_power_vs_t_1toT_kVAr
    elif horizon == "allT":
        pv_reactive_power_allT_kVAr = kVA_B * sum(q_D[j, t] for j in Dset for t in Tset)
        return pv_reactive_power_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to compute battery real power transaction magnitude |P_c - P_d|
def get_battery_real_power_transaction_magnitude(model, data, horizon="allT"):
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    P_c = model.P_c
    P_d = model.P_d

    if horizon == "1toT":
        battery_real_power_magnitude_vs_t_1toT_kW = [
            kVA_B * sum(abs(P_c[j, t] - P_d[j, t]) for j in Bset) for t in Tset
        ]
        return battery_real_power_magnitude_vs_t_1toT_kW
    elif horizon == "allT":
        battery_real_power_magnitude_allT_kW = kVA_B * sum(abs(P_c[j, t] - P_d[j, t]) for j in Bset for t in Tset)
        return battery_real_power_magnitude_allT_kW
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to compute battery reactive power transaction magnitude |q_B|
def get_battery_reactive_power_transaction_magnitude(model, data, horizon="allT"):
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    q_B = model.q_B

    if horizon == "1toT":
        battery_reactive_power_magnitude_vs_t_1toT_kVAr = [
            kVA_B * sum(abs(q_B[j, t]) for j in Bset) for t in Tset
        ]
        return battery_reactive_power_magnitude_vs_t_1toT_kVAr
    elif horizon == "allT":
        battery_reactive_power_magnitude_allT_kVAr = kVA_B * sum(abs(q_B[j, t]) for j in Bset for t in Tset)
        return battery_reactive_power_magnitude_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to compute total static capacitor reactive power generation (if exists)
def get_static_capacitor_reactive_power(model, data, horizon="allT"):
    Tset, kVA_B, T = data['Tset'], data['kVA_B'], data['T']

    if 'q_C' in model:
        q_C = model.q_C

        if horizon == "1toT":
            static_cap_reactive_power_vs_t_1toT_kVAr = [kVA_B * q_C[t] for t in Tset]
            return static_cap_reactive_power_vs_t_1toT_kVAr
        elif horizon == "allT":
            static_cap_reactive_power_allT_kVAr = kVA_B * sum(q_C[t] for t in Tset)
            return static_cap_reactive_power_allT_kVAr
        else:
            raise ValueError("Specify either '1toT' or 'allT'")
    else:
        if horizon == "1toT":
            return np.zeros(T)  # If static capacitor does not exist, return an array of zeros
        elif horizon == "allT":
            return 0.0  # If static capacitor does not exist, return 0
        else:
            raise ValueError("Specify either '1toT' or 'allT'")

# Function to get the peak substation real power over the horizon
def get_substation_real_power_peak(model, data):
    Tset, kVA_B = data['Tset'], data['kVA_B']
    P_Subs = model.P_Subs

    # Calculate the peak substation power (kW) by finding the maximum across all time steps
    peak_substation_power_kW = max([P_Subs[t].value * kVA_B for t in Tset])

    return peak_substation_power_kW

# Function to compute total real power generation
def get_total_generation_real_power(model, data, horizon="allT"):
    if horizon == "1toT":
        battery_real_power_vs_t_1toT_kW = get_battery_real_power_transaction_magnitude(model, data, horizon="1toT")
        pv_real_power_vs_t_1toT_kW = get_pv_reactive_power(model, data, horizon="1toT")
        total_gen_real_power_vs_t_1toT_kW = [
            b + p for b, p in zip(battery_real_power_vs_t_1toT_kW, pv_real_power_vs_t_1toT_kW)
        ]
        return total_gen_real_power_vs_t_1toT_kW
    elif horizon == "allT":
        battery_real_power_allT_kW = get_battery_real_power_transaction_magnitude(model, data, horizon="allT")
        pv_real_power_allT_kW = get_pv_reactive_power(model, data, horizon="allT")
        total_gen_real_power_allT_kW = battery_real_power_allT_kW + pv_real_power_allT_kW
        return total_gen_real_power_allT_kW
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to compute total reactive power generation
def get_total_generation_reactive_power(model, data, horizon="allT"):
    if horizon == "1toT":
        battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power_transaction_magnitude(model, data, horizon="1toT")
        pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(model, data, horizon="1toT")
        static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(model, data, horizon="1toT")
        total_gen_reactive_power_vs_t_1toT_kVAr = [
            b + p + s for b, p, s in zip(
                battery_reactive_power_vs_t_1toT_kVAr,
                pv_reactive_power_vs_t_1toT_kVAr,
                static_cap_reactive_power_vs_t_1toT_kVAr
            )
        ]
        return total_gen_reactive_power_vs_t_1toT_kVAr
    elif horizon == "allT":
        battery_reactive_power_allT_kVAr = get_battery_reactive_power_transaction_magnitude(model, data, horizon="allT")
        pv_reactive_power_allT_kVAr = get_pv_reactive_power(model, data, horizon="allT")
        static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(model, data, horizon="allT")
        total_gen_reactive_power_allT_kVAr = (
            battery_reactive_power_allT_kVAr +
            pv_reactive_power_allT_kVAr +
            static_cap_reactive_power_allT_kVAr
        )
        return total_gen_reactive_power_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to get total real load over the horizon
def get_load_real_power(model, data, horizon="allT"):
    Tset, p_L_pu, kVA_B = data['Tset'], data['p_L_pu'], data['kVA_B']

    if horizon == "1toT":
        load_real_power_vs_t_kW = [
            kVA_B * sum(p_L_pu[j][t-1] for j in p_L_pu.keys()) for t in Tset
        ]
        return load_real_power_vs_t_kW
    elif horizon == "allT":
        load_real_power_kW_allT = kVA_B * sum(p_L_pu[j][t-1] for j in p_L_pu.keys() for t in Tset)
        return load_real_power_kW_allT
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to get total reactive load over the horizon
def get_load_reactive_power(model, data, horizon="allT"):
    Tset, q_L_pu, kVA_B = data['Tset'], data['q_L_pu'], data['kVA_B']

    if horizon == "1toT":
        load_reactive_power_vs_t_kVAr = [
            kVA_B * sum(q_L_pu[j][t-1] for j in q_L_pu.keys()) for t in Tset
        ]
        return load_reactive_power_vs_t_kVAr
    elif horizon == "allT":
        load_reactive_power_allT_kVAr = kVA_B * sum(q_L_pu[j][t-1] for j in q_L_pu.keys() for t in Tset)
        return load_reactive_power_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")

# Function to compute solution time (in seconds)
def get_solution_time(model):
    # Assumes `solve_time` is part of the model object
    return model.solve_time
