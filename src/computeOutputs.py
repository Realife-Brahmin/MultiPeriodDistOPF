import numpy as np
from functionRetriever import (
    get_loss_real_power, get_substation_reactive_power, get_substation_real_power,
    get_substation_power_cost, get_scd, get_terminal_SOC_violation,
    get_loss_reactive_power, get_battery_real_power, get_battery_reactive_power,
    get_pv_real_power, get_pv_reactive_power, get_load_real_power,
    get_load_reactive_power, get_static_capacitor_reactive_power,
    get_substation_real_power_peak, get_total_generation_real_power,
    get_total_generation_reactive_power,
    get_battery_real_power_transaction_magnitude,
    get_battery_reactive_power_transaction_magnitude, get_solution_time
)

def compute_output_values(model, data, verbose=False):
    # Initialize arrays of size T
    T = data['T']
    fval_vs_t_1toT = np.zeros(T)

    # Scalars to store cumulative values
    fval_allT = 0.0

    PLoss_vs_t_1toT_kW = get_loss_real_power(model, data, horizon="1toT")
    PLoss_allT_kW = get_loss_real_power(model, data, horizon="allT")

    QSubs_vs_t_1toT_kVAr = get_substation_reactive_power(model, data, horizon="1toT")
    QSubs_allT_kVAr = get_substation_reactive_power(model, data, horizon="allT")

    PSubs_vs_t_1toT_kW = get_substation_real_power(model, data, horizon="1toT")
    PSubs_allT_kW = get_substation_real_power(model, data, horizon="allT")

    PSubsCost_vs_t_1toT_dollar = get_substation_power_cost(model, data, horizon="1toT")
    PSubsCost_allT_dollar = get_substation_power_cost(model, data, horizon="allT")

    scd_vs_t_1toT_kW = get_scd(model, data, horizon="1toT")
    scd_allT_kW = get_scd(model, data, horizon="allT")

    terminal_soc_violation_kWh = get_terminal_SOC_violation(model, data)

    QLoss_vs_t_1toT_kVAr = get_loss_reactive_power(model, data, horizon="1toT")
    QLoss_allT_kVAr = get_loss_reactive_power(model, data, horizon="allT")

    battery_real_power_vs_t_1toT_kW = get_battery_real_power(model, data, horizon="1toT")
    battery_real_power_allT_kW = get_battery_real_power(model, data, horizon="allT")

    battery_reactive_power_vs_t_1toT_kVAr = get_battery_reactive_power(model, data, horizon="1toT")
    battery_reactive_power_allT_kVAr = get_battery_reactive_power(model, data, horizon="allT")

    pv_real_power_vs_t_1toT_kW = get_pv_real_power(model, data, horizon="1toT")
    pv_real_power_allT_kW = get_pv_real_power(model, data, horizon="allT")

    pv_reactive_power_vs_t_1toT_kVAr = get_pv_reactive_power(model, data, horizon="1toT")
    pv_reactive_power_allT_kVAr = get_pv_reactive_power(model, data, horizon="allT")

    load_real_power_vs_t_1toT_kW = get_load_real_power(model, data, horizon="1toT")
    load_real_power_allT_kW = get_load_real_power(model, data, horizon="allT")

    load_reactive_power_vs_t_1toT_kVAr = get_load_reactive_power(model, data, horizon="1toT")
    load_reactive_power_allT_kVAr = get_load_reactive_power(model, data, horizon="allT")

    static_cap_reactive_power_vs_t_1toT_kVAr = get_static_capacitor_reactive_power(model, data, horizon="1toT")
    static_cap_reactive_power_allT_kVAr = get_static_capacitor_reactive_power(model, data, horizon="allT")

    substation_real_power_peak_allT_kW = get_substation_real_power_peak(model, data)

    total_gen_real_power_vs_t_1toT_kW = get_total_generation_real_power(model, data, horizon="1toT")
    total_gen_real_power_allT_kW = get_total_generation_real_power(model, data, horizon="allT")

    total_gen_reactive_power_vs_t_1toT_kVAr = get_total_generation_reactive_power(model, data, horizon="1toT")
    total_gen_reactive_power_allT_kVAr = get_total_generation_reactive_power(model, data, horizon="allT")

    battery_real_power_transaction_magnitude_vs_t_1toT_kW = get_battery_real_power_transaction_magnitude(model, data, horizon="1toT")
    battery_real_power_transaction_magnitude_allT_kW = get_battery_real_power_transaction_magnitude(model, data, horizon="allT")

    battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr = get_battery_reactive_power_transaction_magnitude(model, data, horizon="1toT")
    battery_reactive_power_transaction_magnitude_allT_kVAr = get_battery_reactive_power_transaction_magnitude(model, data, horizon="allT")

    # solution_time = get_solution_time(model)

    # Loop over time steps to compute all required values
    Tset = data['Tset']
    for time in Tset:
        t=time-1
        kVA_B = data['kVA_B']

        # PLoss_vs_t_1toT = PLoss_vs_t_1toT_kW / kVA_B
        # PSubsCost_vs_t_1toT = PSubsCost_vs_t_1toT_dollar / kVA_B
        # scd_vs_t_1toT = scd_vs_t_1toT_kW / kVA_B

        PLoss_vs_t_1toT = [loss / kVA_B for loss in PLoss_vs_t_1toT_kW]
        PSubsCost_vs_t_1toT = [cost / kVA_B for cost in PSubsCost_vs_t_1toT_dollar]
        scd_vs_t_1toT = [scd / kVA_B for scd in scd_vs_t_1toT_kW]

        objfun0 = data['objfun0']
        objfun2 = data['objfun2']

        if objfun0 == "subsPowerCostMin":
            fval_vs_t_1toT[t] = PSubsCost_vs_t_1toT[t]
        elif objfun0 == "lineLossMin":
            fval_vs_t_1toT[t] = PLoss_vs_t_1toT[t]
        elif objfun0 == "subsPowerMin":
            fval_vs_t_1toT[t] = PSubs_vs_t_1toT_kW[t]
        else:
            raise ValueError("Invalid objfun0")

        if objfun2 == "scd":
            fval_vs_t_1toT[t] += scd_vs_t_1toT[t]
        elif objfun2 == "none":
            pass
        else:
            raise ValueError("Invalid objfun2")

        fval_allT += fval_vs_t_1toT[t]

    # Update data dictionary with computed values
    data.update({
        'battery_reactive_power_allT_kVAr': battery_reactive_power_allT_kVAr,
        'battery_reactive_power_transaction_magnitude_allT_kVAr': battery_reactive_power_transaction_magnitude_allT_kVAr,
        'battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr': battery_reactive_power_transaction_magnitude_vs_t_1toT_kVAr,
        'battery_reactive_power_vs_t_1toT_kVAr': battery_reactive_power_vs_t_1toT_kVAr,
        'battery_real_power_allT_kW': battery_real_power_allT_kW,
        'battery_real_power_transaction_magnitude_allT_kW': battery_real_power_transaction_magnitude_allT_kW,
        'battery_real_power_transaction_magnitude_vs_t_1toT_kW': battery_real_power_transaction_magnitude_vs_t_1toT_kW,
        'battery_real_power_vs_t_1toT_kW': battery_real_power_vs_t_1toT_kW,
        'fval_allT': fval_allT,
        'fval_vs_t_1toT': fval_vs_t_1toT,
        'load_reactive_power_vs_t_1toT_kVAr': load_reactive_power_vs_t_1toT_kVAr,
        'load_reactive_power_allT_kVAr': load_reactive_power_allT_kVAr,
        'load_real_power_vs_t_1toT_kW': load_real_power_vs_t_1toT_kW,
        'load_real_power_allT_kW': load_real_power_allT_kW,
        'PSubsCost_allT_dollar': PSubsCost_allT_dollar,
        'PSubsCost_vs_t_1toT_dollar': PSubsCost_vs_t_1toT_dollar,
        'PSubs_allT_kW': PSubs_allT_kW,
        'PSubs_vs_t_1toT_kW': PSubs_vs_t_1toT_kW,
        'PLoss_allT_kW': PLoss_allT_kW,
        'PLoss_vs_t_1toT_kW': PLoss_vs_t_1toT_kW,
        'pv_real_power_allT_kW': pv_real_power_allT_kW,
        'pv_real_power_vs_t_1toT_kW': pv_real_power_vs_t_1toT_kW,
        'pv_reactive_power_allT_kVAr': pv_reactive_power_allT_kVAr,
        'pv_reactive_power_vs_t_1toT_kVAr': pv_reactive_power_vs_t_1toT_kVAr,
        'QLoss_allT_kVAr': QLoss_allT_kVAr,
        'QLoss_vs_t_1toT_kVAr': QLoss_vs_t_1toT_kVAr,
        'QSubs_allT_kVAr': QSubs_allT_kVAr,
        'QSubs_vs_t_1toT_kVAr': QSubs_vs_t_1toT_kVAr,
        'scd_allT_kW': scd_allT_kW,
        'scd_vs_t_1toT_kW': scd_vs_t_1toT_kW,
        # 'solution_time': solution_time,
        'static_cap_reactive_power_allT_kVAr': static_cap_reactive_power_allT_kVAr,
        'static_cap_reactive_power_vs_t_1toT_kVAr': static_cap_reactive_power_vs_t_1toT_kVAr,
        'substation_real_power_peak_allT_kW': substation_real_power_peak_allT_kW,
        'terminal_soc_violation_kWh': terminal_soc_violation_kWh,
        'total_gen_reactive_power_allT_kVAr': total_gen_reactive_power_allT_kVAr,
        'total_gen_reactive_power_vs_t_1toT_kVAr': total_gen_reactive_power_vs_t_1toT_kVAr,
        'total_gen_real_power_allT_kW': total_gen_real_power_allT_kW,
        'total_gen_real_power_vs_t_1toT_kW': total_gen_real_power_vs_t_1toT_kW
    })

    return data  # Return the updated data dictionary

