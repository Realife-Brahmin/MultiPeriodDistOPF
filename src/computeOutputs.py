from functools import reduce
import numpy as np

from src.functionRetriever import *

def compute_output_values(modelDict, verbose=False):
    model = modelDict['model']
    data = modelDict['data']

    # Initialize arrays of size T
    T = data['T']
    fval_vs_t_1toT = np.zeros(T)

    # Scalars to store cumulative values
    fval_allT = 0.0

    # Retrieve required values using helper functions
    PLoss_vs_t_1toT_kW = get_loss_real_power(modelDict, horizon="1toT")
    PLoss_allT_kW = get_loss_real_power(modelDict, horizon="allT")

    QSubs_vs_t_1toT_kVAr = get_substation_reactive_power(modelDict, horizon="1toT")
    QSubs_allT_kVAr = get_substation_reactive_power(modelDict, horizon="allT")

    PSubs_vs_t_1toT_kW = get_substation_real_power(modelDict, horizon="1toT")
    PSubs_allT_kW = get_substation_real_power(modelDict, horizon="allT")

    PSubsCost_vs_t_1toT_dollar = get_substation_power_cost(modelDict, horizon="1toT")
    PSubsCost_allT_dollar = get_substation_power_cost(modelDict, horizon="allT")

    scd_vs_t_1toT_kW = get_scd(modelDict, horizon="1toT")
    scd_allT_kW = get_scd(modelDict, horizon="allT")

    terminal_soc_violation_kWh = get_terminal_SOC_violation(modelDict)

    QLoss_vs_t_1toT_kVAr = get_loss_reactive_power(modelDict, horizon="1toT")
    QLoss_allT_kVAr = get_loss_reactive_power(modelDict, horizon="allT")

    battery_real_power_vs_t_1toT_kW = get_battery_real_power(modelDict, horizon="1toT")
    battery_real_power_allT_kW = get_battery_real_power(modelDict, horizon="allT")

    battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(modelDict, horizon="1toT")
    battery_reactive_power_allT_kVAr = get_battery_reactive_power(modelDict, horizon="allT")

    pv_real_power_vs_t_1toT_kW = get_pv_real_power(modelDict, horizon="1toT")
    pv_real_power_allT_kW = get_pv_real_power(modelDict, horizon="allT")

    pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(modelDict, horizon="1toT")
    pv_reactive_power_allT_kVAr = get_pv_reactive_power(modelDict, horizon="allT")

    load_real_power_vs_t_1toT_kW = get_load_real_power(data, horizon="1toT")
    load_real_power_allT_kW = get_load_real_power(data, horizon="allT")

    load_reactive_power_vs_t_1toT_kVAr = get_load_reactive_power(data, horizon="1toT")
    load_reactive_power_allT_kVAr = get_load_reactive_power(data, horizon="allT")

    static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="1toT")
    static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(modelDict, horizon="allT")

    substation_real_power_peak_allT_kW = get_substation_real_power_peak(modelDict)

    total_gen_real_power_vs_t_1toT_kW = get_total_generation_real_power(modelDict, horizon="1toT")
    total_gen_real_power_allT_kW = get_total_generation_real_power(modelDict, horizon="allT")

    total_gen_reactive_power_vs_t_1toT_kVAr = get_total_generation_reactive_power(modelDict, horizon="1toT")
    total_gen_reactive_power_allT_kVAr = get_total_generation_reactive_power(modelDict, horizon="allT")

    battery_real_power_transaction_magnitude_vs_t_1toT_kW = get_battery_real_power_transaction_magnitude(modelDict, horizon="1toT")
    battery_real_power_transaction_magnitude_allT_kW = get_battery_real_power_transaction_magnitude(modelDict, horizon="allT")

    battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr = get_battery_reactive_power_transaction_magnitude(modelDict, horizon="1toT")
    battery_reactive_power_transaction_magnitude_allT_kVAr = get_battery_reactive_power_transaction_magnitude(modelDict, horizon="allT")

    solution_time = get_solution_time(modelDict)

    # Loop over time steps to compute all required values
    Tset = data['Tset']
    for t in Tset:
        kVA_B = data['kVA_B']

        PLoss_vs_t_1toT = [p / kVA_B for p in PLoss_vs_t_1toT_kW]
        PSubsCost_vs_t_1toT = [c / kVA_B for c in PSubsCost_vs_t_1toT_dollar]
        scd_vs_t_1toT = [s / kVA_B for s in scd_vs_t_1toT_kW]

        objfun0 = data['objfun0']
        objfun2 = data['objfun2']

        if objfun0 == "subsPowerCostMin":
            fval_vs_t_1toT[t - 1] = PSubsCost_vs_t_1toT[t - 1]
        elif objfun0 == "lineLossMin":
            fval_vs_t_1toT[t - 1] = PLoss_vs_t_1toT[t - 1]
        elif objfun0 == "subsPowerMin":
            fval_vs_t_1toT[t - 1] = PSubs_vs_t_1toT_kW[t - 1]
        else:
            raise ValueError("Unknown objfun0 value")

        if objfun2 == "scd":
            fval_vs_t_1toT[t - 1] += scd_vs_t_1toT[t - 1]
        elif objfun2 == "none":
            pass
        else:
            raise ValueError("Unknown objfun2 value")

        fval_allT += fval_vs_t_1toT[t - 1]

    # Pack data back into the dictionary
    data.update({
        "battery_reactive_power_allT_kVAr": battery_reactive_power_allT_kVAr,
        "battery_reactive_power_transaction_magnitude_allT_kVAr": battery_reactive_power_transaction_magnitude_allT_kVAr,
        "battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr": battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr,
        "battery_reactive_power_vs_t_1toT_kVAr": battery_reactive_power_vs_t_1toT_kVAr,
        "battery_real_power_allT_kW": battery_real_power_allT_kW,
        "battery_real_power_transaction_magnitude_allT_kW": battery_real_power_transaction_magnitude_allT_kW,
        "battery_real_power_transaction_magnitude_vs_t_1toT_kW": battery_real_power_transaction_magnitude_vs_t_1toT_kW,
        "battery_real_power_vs_t_1toT_kW": battery_real_power_vs_t_1toT_kW,
        "fval_allT": fval_allT,
        "fval_vs_t_1toT": fval_vs_t_1toT,
        "load_reactive_power_vs_t_1toT_kVAr": load_reactive_power_vs_t_1toT_kVAr,
        "load_reactive_power_allT_kVAr": load_reactive_power_allT_kVAr,
        "load_real_power_vs_t_1toT_kW": load_real_power_vs_t_1toT_kW,
        "load_real_power_allT_kW": load_real_power_allT_kW,
        "PSubsCost_allT_dollar": PSubsCost_allT_dollar,
        "PSubsCost_vs_t_1toT_dollar": PSubsCost_vs_t_1toT_dollar,
        "PSubs_allT_kW": PSubs_allT_kW,
        "PSubs_vs_t_1toT_kW": PSubs_vs_t_1toT_kW,
        "PLoss_allT_kW": PLoss_allT_kW,
        "PLoss_vs_t_1toT_kW" : PLoss_vs_t_1toT_kW,
        "pv_real_power_allT_kW" : pv_real_power_allT_kW,
        "pv_real_power_vs_t_1toT_kW" : pv_real_power_vs_t_1toT_kW,
        "pv_reactive_power_allT_kVAr" : pv_reactive_power_allT_kVAr,
        "pv_reactive_power_vs_t_1toT_kVAr" : pv_reactive_power_vs_t_1toT_kVAr,
        "QLoss_allT_kVAr" : QLoss_allT_kVAr,
        "QLoss_vs_t_1toT_kVAr" : QLoss_vs_t_1toT_kVAr,
        "QSubs_allT_kVAr" : QSubs_allT_kVAr,
        "QSubs_vs_t_1toT_kVAr" : QSubs_vs_t_1toT_kVAr,
        "scd_allT_kW" : scd_allT_kW,
        "scd_vs_t_1toT_kW" : scd_vs_t_1toT_kW,
        "solution_time" : solution_time,
        "static_cap_reactive_power_allT_kVAr" : static_cap_reactive_power_allT_kVAr,
        "static_cap_reactive_power_vs_t_1toT_kVAr" : static_cap_reactive_power_vs_t_1toT_kVAr,
        "substation_real_power_peak_allT_kW" : substation_real_power_peak_allT_kW,
        "terminal_soc_violation_kWh" : terminal_soc_violation_kWh,
        "total_gen_reactive_power_allT_kVAr" : total_gen_reactive_power_allT_kVAr,
        "total_gen_reactive_power_vs_t_1toT_kVAr" : total_gen_reactive_power_vs_t_1toT_kVAr,
        "total_gen_real_power_allT_kW" : total_gen_real_power_allT_kW,
        "total_gen_real_power_vs_t_1toT_kW" : total_gen_real_power_vs_t_1toT_kW
        })
    
    # Update modelDict with the updated data dictionary
    modelDict['data'] = data

    # Return the updated modelDict
    return modelDict
