__all__ = [
    "get_battery_reactive_power",
    "get_battery_real_power",
    "get_battery_reactive_power_transaction_magnitude",
    "get_battery_real_power_transaction_magnitude",
    "get_load_reactive_power",
    "get_load_real_power",
    "get_pv_reactive_power",
    "get_pv_real_power",
    "get_loss_reactive_power",
    "get_loss_real_power",
    "get_scd",
    "get_solution_time",
    "get_substation_power_cost",
    "get_substation_reactive_power",
    "get_substation_real_power",
    "get_substation_real_power_peak",
    "get_terminal_SOC_violation",
    "get_total_generation_reactive_power",
    "get_total_generation_real_power",
    "get_static_capacitor_reactive_power",
    "get_load_real_power",
    "get_load_reactive_power",
]

def get_loss_real_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Lset, rdict_pu, kVA_B = data['Tset'], data['Lset'], data['rdict_pu'], data['kVA_B']
    l = model.l

    if horizon == "1toT":
        return [kVA_B * sum(rdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(rdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_loss_reactive_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Lset, xdict_pu, kVA_B = data['Tset'], data['Lset'], data['xdict_pu'], data['kVA_B']
    l = model.l

    if horizon == "1toT":
        return [kVA_B * sum(xdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(xdict_pu[i, j] * l[(i, j), t].value for (i, j) in Lset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_substation_reactive_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, L1set, kVA_B = data['Tset'], data['L1set'], data['kVA_B']
    Q = model.Q

    if horizon == "1toT":
        return [kVA_B * sum(Q[(i, j), t].value for (i, j) in L1set) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(Q[(i, j), t].value for (i, j) in L1set for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_substation_real_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, kVA_B = data['Tset'], data['kVA_B']
    P_Subs = model.P_Subs

    if horizon == "1toT":
        return [P_Subs[t].value * kVA_B for t in Tset]
    elif horizon == "allT":
        return sum(P_Subs[t].value * kVA_B for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_substation_power_cost(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, LoadShapeCost, delta_t, kVA_B = data['Tset'], data['LoadShapeCost'], data['delta_t'], data['kVA_B']
    P_Subs = model.P_Subs

    if horizon == "1toT":
        return [LoadShapeCost[t-1] * P_Subs[t].value * kVA_B * delta_t for t in Tset]
    elif horizon == "allT":
        return sum(LoadShapeCost[t-1] * P_Subs[t].value * kVA_B * delta_t for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_scd(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    P_c = model.P_c
    P_d = model.P_d

    if horizon == "1toT":
        return [kVA_B * sum(max(P_c[j, t].value, P_d[j, t].value) - abs(P_c[j, t].value - P_d[j, t].value) for j in Bset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(max(P_c[j, t].value, P_d[j, t].value) - abs(P_c[j, t].value - P_d[j, t].value) for j in Bset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_terminal_SOC_violation(modelDict):
    model = modelDict['model']
    data = modelDict['data']
    T, Bset, Bref_pu, kVA_B = data['T'], data['Bset'], data['Bref_pu'], data['kVA_B']
    B = model.B

    return kVA_B * sum(abs(B[j, T].value - Bref_pu[j]) for j in Bset)


def get_battery_real_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    P_c, P_d = model.P_c, model.P_d

    if horizon == "1toT":
        return [kVA_B * sum(P_d[j, t].value - P_c[j, t].value for j in Bset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(P_d[j, t].value - P_c[j, t].value for j in Bset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_battery_reactive_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    q_B = model.q_B

    if horizon == "1toT":
        return [kVA_B * sum(q_B[j, t].value for j in Bset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(q_B[j, t].value for j in Bset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_pv_real_power(modelDict, horizon="allT"):
    data = modelDict['data']
    Tset, Dset, p_D_pu, kVA_B = data['Tset'], data['Dset'], data['p_D_pu'], data['kVA_B']

    if horizon == "1toT":
        return [kVA_B * sum(p_D_pu[(j, t)] for j in Dset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(p_D_pu[(j, t)] for j in Dset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_pv_reactive_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Dset, kVA_B = data['Tset'], data['Dset'], data['kVA_B']
    q_D = model.q_D

    if horizon == "1toT":
        return [kVA_B * sum(q_D[j, t].value for j in Dset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(q_D[j, t].value for j in Dset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")
    

import numpy as np

def get_battery_real_power_transaction_magnitude(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    P_c = model.P_c
    P_d = model.P_d

    if horizon == "1toT":
        battery_real_power_magnitude_vs_t_1toT_kW = [
            kVA_B * sum(abs(P_c[j, t].value - P_d[j, t].value) for j in Bset) for t in Tset
        ]
        return battery_real_power_magnitude_vs_t_1toT_kW
    elif horizon == "allT":
        battery_real_power_magnitude_allT_kW = kVA_B * sum(
            abs(P_c[j, t].value - P_d[j, t].value) for j in Bset for t in Tset
        )
        return battery_real_power_magnitude_allT_kW
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_battery_reactive_power_transaction_magnitude(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, Bset, kVA_B = data['Tset'], data['Bset'], data['kVA_B']
    q_B = model.q_B

    if horizon == "1toT":
        battery_reactive_power_magnitude_vs_t_1toT_kVAr = [
            kVA_B * sum(abs(q_B[j, t].value) for j in Bset) for t in Tset
        ]
        return battery_reactive_power_magnitude_vs_t_1toT_kVAr
    elif horizon == "allT":
        battery_reactive_power_magnitude_allT_kVAr = kVA_B * sum(
            abs(q_B[j, t].value) for j in Bset for t in Tset
        )
        return battery_reactive_power_magnitude_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_static_capacitor_reactive_power(modelDict, horizon="allT"):
    model = modelDict['model']
    data = modelDict['data']
    Tset, kVA_B, T = data['Tset'], data['kVA_B'], data['T']

    if 'q_C' in model:
        q_C = model.q_C

        if horizon == "1toT":
            static_cap_reactive_power_vs_t_1toT_kVAr = [
                kVA_B * q_C[t].value for t in Tset
            ]
            return static_cap_reactive_power_vs_t_1toT_kVAr
        elif horizon == "allT":
            static_cap_reactive_power_allT_kVAr = kVA_B * sum(q_C[t].value for t in Tset)
            return static_cap_reactive_power_allT_kVAr
        else:
            raise ValueError("Specify either '1toT' or 'allT'")
    else:
        if horizon == "1toT":
            return [0.0] * T  # If static capacitor does not exist, return an array of zeros
        elif horizon == "allT":
            return 0.0  # If static capacitor does not exist, return 0
        else:
            raise ValueError("Specify either '1toT' or 'allT'")


def get_substation_real_power_peak(modelDict):
    model = modelDict['model']
    data = modelDict['data']
    Tset, kVA_B = data['Tset'], data['kVA_B']
    P_Subs = model.P_Subs

    # Calculate the peak substation power (kW) by finding the maximum across all time steps
    peak_substation_power_kW = max([P_Subs[t].value * kVA_B for t in Tset])
    return peak_substation_power_kW


def get_total_generation_real_power(modelDict, horizon="allT"):
    if horizon == "1toT":
        battery_real_power_vs_t_1toT_kW = get_battery_real_power(modelDict, horizon="1toT")
        pv_real_power_vs_t_1toT_kW = get_pv_real_power(modelDict, horizon="1toT")
        total_gen_real_power_vs_t_1toT_kW = [
            b + p for b, p in zip(battery_real_power_vs_t_1toT_kW, pv_real_power_vs_t_1toT_kW)
        ]
        return total_gen_real_power_vs_t_1toT_kW
    elif horizon == "allT":
        battery_real_power_allT_kW = get_battery_real_power(modelDict, horizon="allT")
        pv_real_power_allT_kW = get_pv_real_power(modelDict, horizon="allT")
        total_gen_real_power_allT_kW = battery_real_power_allT_kW + pv_real_power_allT_kW
        return total_gen_real_power_allT_kW
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_total_generation_reactive_power(modelDict, horizon="allT"):
    if horizon == "1toT":
        battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(
            modelDict, horizon="1toT"
        )
        pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(
            modelDict, horizon="1toT"
        )
        static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(
            modelDict, horizon="1toT"
        )
        total_gen_reactive_power_vs_t_1toT_kVAr = [
            b + p + s
            for b, p, s in zip(
                battery_reactive_power_vs_t_1toT_kVAr,
                pv_reactive_power_vs_t_1toT_kVAr,
                static_cap_reactive_power_vs_t_1toT_kVAr,
            )
        ]
        return total_gen_reactive_power_vs_t_1toT_kVAr
    elif horizon == "allT":
        battery_reactive_power_allT_kVAr = get_battery_reactive_power(
            modelDict, horizon="allT"
        )
        pv_reactive_power_allT_kVAr = get_pv_reactive_power(
            modelDict, horizon="allT"
        )
        static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(
            modelDict, horizon="allT"
        )
        total_gen_reactive_power_allT_kVAr = (
            battery_reactive_power_allT_kVAr
            + pv_reactive_power_allT_kVAr
            + static_cap_reactive_power_allT_kVAr
        )
        return total_gen_reactive_power_allT_kVAr
    else:
        raise ValueError("Specify either '1toT' or 'allT'")



def get_load_real_power(data, horizon="allT"):
    Tset, NLset, p_L_pu, kVA_B = data['Tset'], data['NLset'], data['p_L_pu'], data['kVA_B']

    if horizon == "1toT":
        return [kVA_B * sum(p_L_pu[(j, t)] for j in NLset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(p_L_pu[(j, t)] for j in NLset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_load_reactive_power(data, horizon="allT"):
    Tset, NLset, q_L_pu, kVA_B = data['Tset'], data['NLset'], data['q_L_pu'], data['kVA_B']

    if horizon == "1toT":
        return [kVA_B * sum(q_L_pu[(j, t)] for j in NLset) for t in Tset]
    elif horizon == "allT":
        return kVA_B * sum(q_L_pu[(j, t)] for j in NLset for t in Tset)
    else:
        raise ValueError("Specify either '1toT' or 'allT'")


def get_solution_time(modelDict):
    modelVals = modelDict['modelVals']
    return modelVals['solve_time']
