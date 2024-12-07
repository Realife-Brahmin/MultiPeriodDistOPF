# from pyomo.environ import *
# from math import sqrt

# def nodalRealPowerBalance_substation_t_in_Tset(model, t):
#     substationBus = model.data['substationBus']
#     L1set = model.data['L1set']
#     return model.P_Subs[t] == sum(model.P[(substationBus, j), t] for (substationBus, j) in L1set)

# model.nodalRealPowerBalance_substation_t_in_Tset = Constraint(model.T, rule=nodalRealPowerBalance_substation_t_in_Tset)

# def nodalRealPowerBalance_non_substation_t_in_Tset(model, t, j):
#     children = model.data['children']
#     parent = model.data['parent']
#     rdict_pu = model.data['rdict_pu']

#     P = model.P
#     l = model.l

#     sum_Pjk = sum(P[(j, k), t] for k in children.get(j, [])) if j in children else 0.0
#     i = parent[j]
#     P_ij_t = P[(i, j), t]
#     r_ij = rdict_pu[(i, j)]
#     l_ij_t = l[(i, j), t]

#     line_loss = r_ij * l_ij_t
#     p_L_j_t = model.data['p_L_pu'].get((j, t), 0.0)
#     p_D_j_t = model.data['p_D_pu'].get((j, t), 0.0)

#     P_c = model.P_c
#     P_d = model.P_d

#     P_c_j_t = P_c[j, t] if j in model.data['Bset'] else 0.0
#     P_d_j_t = P_d[j, t] if j in model.data['Bset'] else 0.0

#     return sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0

# model.nodalRealPowerBalance_non_substation_t_in_Tset = Constraint(model.T, model.Nm1set, rule=nodalRealPowerBalance_non_substation_t_in_Tset)

# def nodalReactivePowerBalance_non_substation_t_in_Tset(model, t, j):
#     children = model.data['children']
#     parent = model.data['parent']
#     xdict_pu = model.data['xdict_pu']

#     Q = model.Q
#     l = model.l

#     sum_Qjk = sum(Q[(j, k), t] for k in children.get(j, [])) if j in children else 0.0
#     i = parent[j]
#     Q_ij_t = Q[(i, j), t]
#     x_ij = xdict_pu[(i, j)]
#     l_ij_t = l[(i, j), t]

#     line_reactive_loss = x_ij * l_ij_t
#     q_L_j_t = model.data['q_L_pu'].get((j, t), 0.0)
#     q_D_j_t = model.data['q_D'].get((j, t), 0.0)
#     q_B_j_t = model.q_B[j, t] if j in model.data['Bset'] else 0.0

#     return sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0

# model.nodalReactivePowerBalance_non_substation_t_in_Tset = Constraint(model.T, model.Nm1set, rule=nodalReactivePowerBalance_non_substation_t_in_Tset)

# def KVL_substation_branches_t_in_Tset(model, t, i, j):
#     rdict_pu = model.data['rdict_pu']
#     xdict_pu = model.data['xdict_pu']

#     r_ij = rdict_pu[(i, j)]
#     x_ij = xdict_pu[(i, j)]

#     P_ij_t = model.P[(i, j), t]
#     Q_ij_t = model.Q[(i, j), t]
#     l_ij_t = model.l[(i, j), t]
#     v_i_t = model.v[i, t]
#     v_j_t = model.v[j, t]

#     return v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij**2 + x_ij**2) * l_ij_t == 0

# model.KVL_substation_branches_t_in_Tset = Constraint(model.T, model.L1set, rule=KVL_substation_branches_t_in_Tset)

# def KVL_non_substation_branches_t_in_Tset(model, t, i, j):
#     rdict_pu = model.data['rdict_pu']
#     xdict_pu = model.data['xdict_pu']

#     r_ij = rdict_pu[(i, j)]
#     x_ij = xdict_pu[(i, j)]

#     P_ij_t = model.P[(i, j), t]
#     Q_ij_t = model.Q[(i, j), t]
#     l_ij_t = model.l[(i, j), t]
#     v_i_t = model.v[i, t]
#     v_j_t = model.v[j, t]

#     return v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij**2 + x_ij**2) * l_ij_t == 0

# model.KVL_non_substation_branches_t_in_Tset = Constraint(model.T, model.Lm1set, rule=KVL_non_substation_branches_t_in_Tset)

# def BCPF_substation_branches_t_in_Tset(model, t, i, j):
#     P_ij_t = model.P[(i, j), t]
#     Q_ij_t = model.Q[(i, j), t]
#     v_i_t = model.v[i, t]
#     l_ij_t = model.l[(i, j), t]

#     return P_ij_t**2 + Q_ij_t**2 - v_i_t * l_ij_t == 0

# model.BCPF_substation_branches_t_in_Tset = Constraint(model.T, model.L1set, rule=BCPF_substation_branches_t_in_Tset)

# def BCPF_non_substation_branches_t_in_Tset(model, t, i, j):
#     P_ij_t = model.P[(i, j), t]
#     Q_ij_t = model.Q[(i, j), t]
#     v_i_t = model.v[i, t]
#     l_ij_t = model.l[(i, j), t]

#     return P_ij_t**2 + Q_ij_t**2 - v_i_t * l_ij_t == 0

# model.BCPF_non_substation_branches_t_in_Tset = Constraint(model.T, model.Lm1set, rule=BCPF_non_substation_branches_t_in_Tset)

# def battery_SOC_constraints_t_in_Tset(model, t, j):
#     delta_t = model.data['delta_t']
#     eta_C = model.data['eta_C']
#     eta_D = model.data['eta_D']

#     P_c = model.P_c
#     P_d = model.P_d
#     B = model.B

#     if t == 1:
#         B0_pu = model.data['B0_pu']
#         return B[j, t] == (B0_pu[j] + delta_t * eta_C[j] * P_c[j, t] - delta_t * (1 / eta_D[j]) * P_d[j, t])
#     else:
#         return B[j, t] == (B[j, t-1] + delta_t * eta_C[j] * P_c[j, t] - delta_t * (1 / eta_D[j]) * P_d[j, t])

# model.battery_SOC_constraints_t_in_Tset = Constraint(model.T, model.Bset, rule=battery_SOC_constraints_t_in_Tset)

# def SOC_limits_batteries_t_in_Tset(model, t, j):
#     B = model.B
#     B_R_pu = model.data['B_R_pu']
#     soc_min = model.data['soc_min']
#     soc_max = model.data['soc_max']

#     return (soc_min[j] * B_R_pu[j] <= B[j, t]) & (B[j, t] <= soc_max[j] * B_R_pu[j])

# model.SOC_limits_batteries_t_in_Tset = Constraint(model.T, model.Bset, rule=SOC_limits_batteries_t_in_Tset)

# from pyomo.environ import *
# from math import sqrt

# def battery_SOC_constraints_t_in_Tset(model, t, j):
#     delta_t = model.data['delta_t']
#     eta_C = model.data['eta_C']
#     eta_D = model.data['eta_D']
#     B = model.B
#     P_c = model.P_c
#     P_d = model.P_d

#     # Constraint h_SOC_j^{t=1}: Initial SOC constraint
#     if t == 1:
#         B0_pu = model.data['B0_pu']
#         return B[j, t] == B0_pu[j] + delta_t * eta_C[j] * P_c[j, t] - delta_t * (1 / eta_D[j]) * P_d[j, t]
    
#     # Constraint h_SOC_j^{t=2 to T}: SOC trajectory for middle and final time periods
#     else:
#         return B[j, t] == B[j, t - 1] + delta_t * eta_C[j] * P_c[j, t] - delta_t * (1 / eta_D[j]) * P_d[j, t]

# model.battery_SOC_constraints_t_in_Tset = Constraint(model.T, model.Bset, rule=battery_SOC_constraints_t_in_Tset)

# def battery_SOC_terminal_constraint(model, j):
#     T = max(model.T)
#     tSOC_hard = model.data.get('tSOC_hard', False)
#     Bref_pu = model.data['Bref_pu']

#     # Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = Bref_j)
#     if tSOC_hard:
#         return model.B[j, T] == Bref_pu[j]
#     else:
#         return Constraint.Skip

# model.battery_SOC_terminal_constraint = Constraint(model.Bset, rule=battery_SOC_terminal_constraint)

# def fixed_substation_voltage_constraints_t_in_Tset(model, t):
#     substationBus = model.data['substationBus']
#     V_Subs = model.data['V_Subs']
#     return model.v[substationBus, t] == V_Subs**2

# model.fixed_substation_voltage_constraints_t_in_Tset = Constraint(model.T, rule=fixed_substation_voltage_constraints_t_in_Tset)

# def voltage_limits_constraints_t_in_Tset(model, t, j):
#     Vminpu_Comp = model.data['Vminpu_Comp']
#     Vmaxpu_Comp = model.data['Vmaxpu_Comp']
#     substationBus = model.data['substationBus']

#     # Voltage limit constraints
#     if j == substationBus:
#         return model.v[substationBus, t] == model.data['V_Subs']**2
#     else:
#         V_min_j_sq = Vminpu_Comp[j]**2
#         V_max_j_sq = Vmaxpu_Comp[j]**2
#         return (V_min_j_sq <= model.v[j, t]) & (model.v[j, t] <= V_max_j_sq)

# model.voltage_limits_constraints_t_in_Tset = Constraint(model.T, model.Compset, rule=voltage_limits_constraints_t_in_Tset)

# def reactive_power_limits_PV_inverters_t_in_Tset(model, t, j):
#     p_D_R_pu = model.data['p_D_R_pu']
#     Dset = model.data['Dset']
#     q_D = model.q_D

#     if j in Dset:
#         p_D_j_t = model.data['p_D_pu'][(j, t)]
#         q_D_Max_j_t = sqrt((1.2 * p_D_R_pu[j])**2 - p_D_j_t**2)
#         return (-q_D_Max_j_t <= q_D[j, t]) & (q_D[j, t] <= q_D_Max_j_t)
#     else:
#         return Constraint.Skip

# model.reactive_power_limits_PV_inverters_t_in_Tset = Constraint(model.T, model.Dset, rule=reactive_power_limits_PV_inverters_t_in_Tset)

# def reactive_power_limits_battery_inverters_t_in_Tset(model, t, j):
#     Bset = model.data['Bset']
#     P_B_R_pu = model.data['P_B_R_pu']
#     q_B = model.q_B

#     if j in Bset:
#         q_B_Max_j = sqrt((1.2 * P_B_R_pu[j])**2 - P_B_R_pu[j]**2)
#         return (-q_B_Max_j <= q_B[j, t]) & (q_B[j, t] <= q_B_Max_j)
#     else:
#         return Constraint.Skip

# model.reactive_power_limits_battery_inverters_t_in_Tset = Constraint(model.T, model.Bset, rule=reactive_power_limits_battery_inverters_t_in_Tset)

# def charging_power_limits_batteries_t_in_Tset(model, t, j):
#     Bset = model.data['Bset']
#     P_B_R_pu = model.data['P_B_R_pu']
#     P_c = model.P_c

#     if j in Bset:
#         return (0 <= P_c[j, t]) & (P_c[j, t] <= P_B_R_pu[j])
#     else:
#         return Constraint.Skip

# model.charging_power_limits_batteries_t_in_Tset = Constraint(model.T, model.Bset, rule=charging_power_limits_batteries_t_in_Tset)

# def discharging_power_limits_batteries_t_in_Tset(model, t, j):
#     Bset = model.data['Bset']
#     P_B_R_pu = model.data['P_B_R_pu']
#     P_d = model.P_d

#     if j in Bset:
#         return (0 <= P_d[j, t]) & (P_d[j, t] <= P_B_R_pu[j])
#     else:
#         return Constraint.Skip

# model.discharging_power_limits_batteries_t_in_Tset = Constraint(model.T, model.Bset, rule=discharging_power_limits_batteries_t_in_Tset)

# def SOC_limits_batteries_t_in_Tset(model, t, j):
#     B = model.B
#     B_R_pu = model.data['B_R_pu']
#     soc_min = model.data['soc_min']
#     soc_max = model.data['soc_max']

#     return (soc_min[j] * B_R_pu[j] <= B[j, t]) & (B[j, t] <= soc_max[j] * B_R_pu[j])

# model.SOC_limits_batteries_t_in_Tset = Constraint(model.T, model.Bset, rule=SOC_limits_batteries_t_in_Tset)

# Constraints.py

from pyomo.environ import Constraint, sqrt

def nodalRealPowerBalance_substation_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    # Assign relevant data to model attributes if needed
    model.substationBus = data['substationBus']
    model.L1set = data['L1set']  # list of (j1, j) tuples
    j1 = model.substationBus

    def substation_real_power_balance_rule(m, t):
        return m.P_Subs[t] - sum(m.P[(j1, j), t] for (j1, j) in m.L1set) == 0

    model.SubstationRealPowerBalance = Constraint(Tset, rule=substation_real_power_balance_rule)
    return model


def nodalRealPowerBalance_non_substation_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    # Assume Nm1set, children, parent, etc. are available on model or through data
    Nm1set = data['Nm1set']
    children = data['children']
    parent = data['parent']
    rdict_pu = data['rdict_pu']
    NLset = data['NLset']
    p_L_pu = data['p_L_pu']
    Dset = data['Dset']
    p_D_pu = data['p_D_pu']
    Bset = data['Bset']

    # We create a double-indexed constraint over (t,j) in Tset x Nm1set
    # For clarity, we'll define an indexing set in the code. Otherwise, ensure sets are on model.
    model.Nm1set = Nm1set

    def non_substation_power_balance_rule(m, t, j):
        sum_Pjk = sum(m.P[(j, k), t] for k in children[j]) if len(children[j])>0 else 0.0
        i = parent[j]
        P_ij_t = m.P[(i, j), t]
        r_ij = rdict_pu[(i, j)]
        l_ij_t = m.l[(i, j), t]
        line_loss = r_ij * l_ij_t
        p_L_j_t = p_L_pu[(j, t)] if j in NLset else 0.0
        p_D_j_t = p_D_pu[(j, t)] if j in Dset else 0.0
        P_d_j_t = m.P_d[j, t] if j in Bset else 0.0
        P_c_j_t = m.P_c[j, t] if j in Bset else 0.0

        return sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0

    model.NodeRealPowerBalanceNonSubs = Constraint(Tset, Nm1set, rule=non_substation_power_balance_rule)
    return model


def nodalReactivePowerBalance_non_substation_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    Nm1set = data['Nm1set']
    children = data['children']
    parent = data['parent']
    xdict_pu = data['xdict_pu']
    NLset = data['NLset']
    q_L_pu = data['q_L_pu']
    Dset = data['Dset']
    Bset = data['Bset']

    def reactive_balance_rule(m, t, j):
        sum_Qjk = sum(m.Q[(j, k), t] for k in children[j]) if len(children[j])>0 else 0.0
        i = parent[j]
        Q_ij_t = m.Q[(i, j), t]
        x_ij = xdict_pu[(i, j)]
        l_ij_t = m.l[(i, j), t]
        line_reactive_loss = x_ij * l_ij_t
        q_L_j_t = q_L_pu[(j, t)] if j in NLset else 0.0
        q_D_j_t = m.q_D[j, t] if j in Dset else 0.0
        q_B_j_t = m.q_B[j, t] if j in Bset else 0.0

        return sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0

    model.NodeReactivePowerBalanceNonSubs = Constraint(Tset, Nm1set, rule=reactive_balance_rule)
    return model


def KVL_substation_branches_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']
    L1set = data['L1set']
    rdict_pu = data['rdict_pu']
    xdict_pu = data['xdict_pu']

    # We'll define an indexed constraint over Tset x L1set
    # L1set is a list of (i,j) tuples
    def kvl_substation_rule(m, t, i, j):
        # i, j = i_j
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P_ij_t = m.P[(i, j), t]
        Q_ij_t = m.Q[(i, j), t]
        l_ij_t = m.l[(i, j), t]
        v_i_t = m.v[i, t]
        v_j_t = m.v[j, t]
        return v_i_t - v_j_t - 2*(r_ij*P_ij_t + x_ij*Q_ij_t) + (r_ij**2 + x_ij**2)*l_ij_t == 0

    model.KVLSubstation = Constraint(Tset, L1set, rule=kvl_substation_rule)
    return model


def KVL_non_substation_branches_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']
    Lm1set = data['Lm1set']
    rdict_pu = data['rdict_pu']
    xdict_pu = data['xdict_pu']

    def kvl_non_sub_rule(m, t, i,j):
        # i, j = i_j
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P_ij_t = m.P[(i, j), t]
        Q_ij_t = m.Q[(i, j), t]
        l_ij_t = m.l[(i, j), t]
        v_i_t = m.v[i, t]
        v_j_t = m.v[j, t]
        return v_i_t - v_j_t - 2*(r_ij*P_ij_t + x_ij*Q_ij_t) + (r_ij**2 + x_ij**2)*l_ij_t == 0

    model.KVLNonSubstation = Constraint(Tset, Lm1set, rule=kvl_non_sub_rule)
    return model


def BCPF_substation_branches_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']
    L1set = data['L1set']

    def bcpf_sub_rule(m, t, i, j):
        # i, j = i_j
        P_ij_t = m.P[(i, j), t]
        Q_ij_t = m.Q[(i, j), t]
        v_i_t = m.v[i, t]
        l_ij_t = m.l[(i, j), t]
        return (P_ij_t)**2 + (Q_ij_t)**2 - v_i_t*l_ij_t == 0

    model.BCPFSubstation = Constraint(Tset, L1set, rule=bcpf_sub_rule)
    return model


def BCPF_non_substation_branches_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']
    Lm1set = data['Lm1set']

    def bcpf_non_sub_rule(m, t, i, j):
        # i, j = i_j
        P_ij_t = m.P[(i, j), t]
        Q_ij_t = m.Q[(i, j), t]
        v_i_t = m.v[i, t]
        l_ij_t = m.l[(i, j), t]
        return (P_ij_t)**2 + (Q_ij_t)**2 - v_i_t*l_ij_t == 0

    model.BCPFNonSubstation = Constraint(Tset, Lm1set, rule=bcpf_non_sub_rule)
    return model


def battery_SOC_constraints_t_in_Tset(model, data, Tset=None, tSOC_hard=False):
    if Tset is None:
        Tset = data['Tset']

    Bset = data['Bset']
    delta_t = data['delta_t']
    eta_C = data['eta_C']
    eta_D = data['eta_D']
    B0_pu = data['B0_pu']
    Bref_pu = data['Bref_pu']
    T = data['T']

    # Initial SOC constraint only if t=1 is in Tset
    if 1 in Tset:
        def initial_soc_rule(m, j):
            return m.B[j, 1] - (B0_pu[j] + delta_t*eta_C[j]*m.P_c[j, 1] - delta_t*(1/eta_D[j])*m.P_d[j, 1]) == 0
        model.InitialSOC = Constraint(Bset, rule=initial_soc_rule)

    # SOC trajectory for t>1
    def soc_trajectory_rule(m, j, t):
        if t > 1:
            return m.B[j, t] - (m.B[j, t-1] + delta_t*eta_C[j]*m.P_c[j, t] - delta_t*(1/eta_D[j])*m.P_d[j, t]) == 0
        return Constraint.Skip

    model.SOCTrajectory = Constraint(Bset, Tset, rule=soc_trajectory_rule)

    # Final SOC constraint if tSOC_hard and max(Tset)=T
    if tSOC_hard and max(Tset) == T:
        def final_soc_rule(m, j):
            return m.B[j, T] == Bref_pu[j]
        model.FinalSOC = Constraint(Bset, rule=final_soc_rule)

    return model


def fixed_substation_voltage_constraints_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']
    substationBus = data['substationBus']
    V_Subs = data['V_Subs']

    def fixed_substation_voltage_rule(m, t):
        return m.v[substationBus, t] == (V_Subs)**2

    model.FixedSubstationVoltage = Constraint(Tset, rule=fixed_substation_voltage_rule)
    return model


def voltage_limits_constraints_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    Compset = data['Compset']
    Vminpu_Comp = data['Vminpu_Comp']
    Vmaxpu_Comp = data['Vmaxpu_Comp']
    substationBus = data['substationBus']
    V_Subs = data['V_Subs']

    # If you need special handling for the substationBus, you can do so by a separate constraint or within the rule.

    def voltage_limits_rule(m, t, j):
        # If j == substationBus, fix voltage
        if j == substationBus:
            return m.v[substationBus, t] == (V_Subs)**2

        V_min_j_sq = (Vminpu_Comp[j])**2
        V_max_j_sq = (Vmaxpu_Comp[j])**2

        # We must return a single expression. Since we have two inequalities,
        # we can either create two separate constraints or combine them differently.
        # Let's create two separate constraints for clarity.

        return Constraint.Skip  # This main rule won't return both constraints together.

    # We'll define two separate constraints for upper and lower bounds:
    def lower_voltage_bound_rule(m, t, j):
        if j == substationBus:
            return Constraint.Skip
        V_min_j_sq = (Vminpu_Comp[j])**2
        return V_min_j_sq - m.v[j, t] <= 0

    def upper_voltage_bound_rule(m, t, j):
        if j == substationBus:
            return Constraint.Skip
        V_max_j_sq = (Vmaxpu_Comp[j])**2
        return m.v[j, t] - V_max_j_sq <= 0

    # Substation fixed voltage constraint (if needed)
    def fixed_voltage_substation_rule(m, t):
        return m.v[substationBus, t] == (V_Subs)**2

    # Create separate constraints:
    model.LowerVoltageBound = Constraint(Tset, Compset, rule=lower_voltage_bound_rule)
    model.UpperVoltageBound = Constraint(Tset, Compset, rule=upper_voltage_bound_rule)
    model.FixedSubstationVoltage2 = Constraint(Tset, rule=fixed_voltage_substation_rule)

    return model


def reactive_power_limits_PV_inverters_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    p_D_R_pu = data['p_D_R_pu']
    Dset = data['Dset']
    p_D_pu = data['p_D_pu']

    def pv_reactive_limits_rule_lower(m, t, j):
        p_D_R_j = p_D_R_pu[j]
        p_D_j_t = p_D_pu[(j, t)]
        q_D_Max_j_t = sqrt((1.2*p_D_R_j)**2 - (p_D_j_t)**2)
        return -q_D_Max_j_t - m.q_D[j, t] <= 0

    def pv_reactive_limits_rule_upper(m, t, j):
        p_D_R_j = p_D_R_pu[j]
        p_D_j_t = p_D_pu[(j, t)]
        q_D_Max_j_t = sqrt((1.2*p_D_R_j)**2 - (p_D_j_t)**2)
        return m.q_D[j, t] - q_D_Max_j_t <= 0

    model.PVReactiveLower = Constraint(Tset, Dset, rule=pv_reactive_limits_rule_lower)
    model.PVReactiveUpper = Constraint(Tset, Dset, rule=pv_reactive_limits_rule_upper)
    return model


def reactive_power_limits_battery_inverters_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    Bset = data['Bset']
    P_B_R_pu = data['P_B_R_pu']

    # Precompute q_B_Max once:
    q_B_Max = {j: sqrt((1.2*P_B_R_pu[j])**2 - (1.0*P_B_R_pu[j])**2) for j in Bset}

    def battery_reactive_lower_rule(m, t, j):
        return -q_B_Max[j] - m.q_B[j, t] <= 0

    def battery_reactive_upper_rule(m, t, j):
        return m.q_B[j, t] - q_B_Max[j] <= 0

    model.BatteryReactiveLower = Constraint(Tset, Bset, rule=battery_reactive_lower_rule)
    model.BatteryReactiveUpper = Constraint(Tset, Bset, rule=battery_reactive_upper_rule)
    return model


def charging_power_limits_batteries_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    P_B_R_pu = data['P_B_R_pu']
    Bset = data['Bset']

    def charging_nonneg_rule(m, t, j):
        return -m.P_c[j, t] <= 0

    def charging_max_rule(m, t, j):
        return m.P_c[j, t] - P_B_R_pu[j] <= 0

    model.ChargingNonNeg = Constraint(Tset, Bset, rule=charging_nonneg_rule)
    model.ChargingMax = Constraint(Tset, Bset, rule=charging_max_rule)
    return model


def discharging_power_limits_batteries_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    P_B_R_pu = data['P_B_R_pu']
    Bset = data['Bset']

    def discharging_nonneg_rule(m, t, j):
        return -m.P_d[j, t] <= 0

    def discharging_max_rule(m, t, j):
        return m.P_d[j, t] - P_B_R_pu[j] <= 0

    model.DischargingNonNeg = Constraint(Tset, Bset, rule=discharging_nonneg_rule)
    model.DischargingMax = Constraint(Tset, Bset, rule=discharging_max_rule)
    return model


def SOC_limits_batteries_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    Bset = data['Bset']
    B_R_pu = data['B_R_pu']
    soc_min = data['soc_min']
    soc_max = data['soc_max']

    def soc_min_rule(m, t, j):
        return soc_min[j]*B_R_pu[j] - m.B[j, t] <= 0

    def soc_max_rule(m, t, j):
        return m.B[j, t] - soc_max[j]*B_R_pu[j] <= 0

    model.SOCMin = Constraint(Tset, Bset, rule=soc_min_rule)
    model.SOCMax = Constraint(Tset, Bset, rule=soc_max_rule)
    return model
