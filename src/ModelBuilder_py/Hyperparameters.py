from src.functionRetriever import * # Assuming this is the equivalent of `FR.get_load_real_power`
import numpy as np

def estimate_alpha(data):
    func_obj_est = data.get('func_obj_est', None)
    if func_obj_est is not None:
        fobj_est = func_obj_est(data)
    else:
        fobj_est = 0
    fscd_est = estimate_fscd(data)
    alpha = fobj_est / fscd_est
    return alpha

def estimate_fscd(data):
    Bset = data['Bset']
    P_B_R = data['P_B_R']
    eta_C = data['eta_C']
    eta_D = data['eta_D']
    T = data['T']

    # (1/η_D[j] - η_C[j]) * P_B_R[j], summed over j in Bset, then multiplied by T
    fscd_terms = [(1/eta_D[j] - eta_C[j]) * P_B_R[j] for j in Bset]
    fscd_est = sum(fscd_terms) * T
    return fscd_est

def estimate_ftsoc(data, tol_tSOC_violation_pu=1e-3):
    n_B = data['n_B']
    ftsoc_est = (tol_tSOC_violation_pu**2) * n_B
    return ftsoc_est

def estimate_gamma(data):
    func_obj_est = data.get('func_obj_est', None)
    if func_obj_est is not None:
        fobj_est = func_obj_est(data)
    else:
        fobj_est = 0
    ftsoc_est = estimate_ftsoc(data)
    gamma = fobj_est / ftsoc_est
    return gamma

def estimate_substation_power_cost(data):
    LoadShapeCost = data['LoadShapeCost']    # dollars_per_kWh as array
    kVA_B = data['kVA_B']
    # Assuming LoadShapeCost is a np.array of shape (T,)
    load_real_power_vs_t_1toT_kW = get_load_real_power(data, horizon="1toT")

    line_loss_accommodation_factor = 1.01
    # In Julia: transpose(dollars_per_kWh) * load_real_power_vs_t_1toT_kW
    # Assuming both are 1D arrays, this is just a dot product:
    fcost_est = line_loss_accommodation_factor * np.dot(LoadShapeCost,load_real_power_vs_t_1toT_kW)
    return fcost_est

def estimate_line_losses(data):
    kVA_B = data['kVA_B']
    load_real_power_allT_kW = get_load_real_power(data, horizon="allT")

    line_loss_accommodation_factor = 0.01
    fPLoss_est = line_loss_accommodation_factor * load_real_power_allT_kW
    return fPLoss_est
