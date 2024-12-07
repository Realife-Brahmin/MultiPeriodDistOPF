from pyomo.environ import Param, Var  # Assuming Pyomo is used as an analog to JuMP

def initialize_variables_1ph_NL_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data.get('Tset')

    # Initialize substation power flow variable
    P_Subs = model.P_Subs
    for t in Tset:
        P_Subs[t].value = 0.0

    # Initialize power flow variables
    Lset = data.get('Lset')
    P = model.P
    Q = model.Q
    l = model.l
    for (i, j) in Lset:
        for t in Tset:
            P[(i, j), t].value = 0.0
            Q[(i, j), t].value = 0.0
            l[(i, j), t].value = 0.0

    # Initialize voltage variables
    Compset = data.get('Compset')
    V_Subs = data.get('V_Subs')
    v = model.v
    for j in Compset:
        for t in Tset:
            v[j, t].value = V_Subs ** 2

    # Initialize PV inverter reactive dispatch variables
    Dset = data.get('Dset')
    q_D = model.q_D
    for j in Dset:
        for t in Tset:
            q_D[j, t].value = 0.0

    # Initialize battery real and reactive dispatch variables
    Bset = data.get('Bset')
    q_B = model.q_B
    P_c = model.P_c
    P_d = model.P_d
    for j in Bset:
        for t in Tset:
            q_B[j, t].value = 0.0
            P_c[j, t].value = 0.0
            P_d[j, t].value = 0.0

    # Initialize battery state of charge (SOC) variables
    B0_pu = data.get('B0_pu')
    B = model.B
    for j in Bset:
        for t in Tset:
            B[j, t].value = B0_pu[j]

    return model
