from pyomo.environ import ConcreteModel, Var, Constraint, NonNegativeReals, Reals, sqrt, value, Objective, minimize
from pyomo.opt import TerminationCondition, SolverFactory



def optimize_MPOPF_1ph_NL(data):
    """
    Python implementation of the Julia function `optimize_MPOPF_1ph_NL`.
    """

    # Unpack the solver name
    solver_name = data['solver']

    # Define the function to configure the solver
    def configure_solver(solver_name):
        if solver_name == "Ipopt":
            model = ConcreteModel()
            solver_options = {
                "tol": 1e-6,
                "max_iter": 10000,
                "print_level": 5
            }
        elif solver_name == "Gurobi":
            model = ConcreteModel()
            solver_options = {
                "TimeLimit": 300  # Set time limit in seconds
            }
        elif solver_name == "Juniper":
            raise NotImplementedError("Juniper solver configuration is not implemented.")
        elif solver_name == "EAGO":
            raise NotImplementedError("EAGO solver configuration is not implemented.")
        elif solver_name == "MadNLP":
            raise NotImplementedError("MadNLP solver configuration is not implemented.")
        else:
            raise ValueError(f"Unsupported solver: {solver_name}")
        return model, solver_options

    # ===========================
    # Variables
    # ===========================

    # Configure the solver and create the model
    model, solver_options = configure_solver(solver_name)

    # Unpack data
    Tset = data['Tset']
    Nset = data['Nset']
    Lset = data['Lset']
    Dset = data['Dset']
    Bset = data['Bset']
    PSubsMax_kW = data['PSubsMax_kW']
    kVA_B = data['kVA_B']
    substationBus = data['substationBus']
    NLset = data['NLset']
    Nm1set = data['Nm1set']
    children = data['children']
    parent = data['parent']
    rdict_pu = data['rdict_pu']
    xdict_pu = data['xdict_pu']
    p_L_pu = data['p_L_pu']
    p_D_pu = data['p_D_pu']

    # Derived data
    PSubsMax_pu = PSubsMax_kW / kVA_B

    # Define variables
    model.P_Subs = Var(Tset, domain=NonNegativeReals, bounds=(0, PSubsMax_pu))
    model.P = Var(Lset, Tset, domain=Reals, bounds=(-PSubsMax_pu, PSubsMax_pu))
    model.Q = Var(Lset, Tset, domain=Reals)
    model.l = Var(Lset, Tset, domain=NonNegativeReals)
    model.v = Var(Nset, Tset, domain=NonNegativeReals)
    model.q_D = Var(Dset, Tset, domain=Reals)
    model.q_B = Var(Bset, Tset, domain=Reals)
    model.P_c = Var(Bset, Tset, domain=NonNegativeReals)
    model.P_d = Var(Bset, Tset, domain=NonNegativeReals)
    model.B = Var(Bset, Tset, domain=NonNegativeReals)

    # Substation node
    j1 = substationBus

    # Real power balance at substation node
    def substation_power_balance_rule(model, t):
        return model.P_Subs[t] - sum(model.P[j1, j, t] for (j1, j) in Lset if j1 == substationBus) == 0
    model.substation_power_balance = Constraint(Tset, rule=substation_power_balance_rule)

    # Real power balance at non-substation nodes
    def node_power_balance_rule(model, t, j):
        sum_Pjk = sum(model.P[j, k, t] for k in children.get(j, [])) if j in children else 0
        i = parent[j]
        P_ij_t = model.P[i, j, t]
        r_ij = rdict_pu[(i, j)]
        line_loss = r_ij * model.l[i, j, t]
        p_L_j_t = p_L_pu[j][t-1] if j in NLset else 0.0
        p_D_j_t = p_D_pu[j][t-1] if j in Dset else 0.0
        P_d_j_t = model.P_d[j, t] if j in Bset else 0.0
        P_c_j_t = model.P_c[j, t] if j in Bset else 0.0
        return sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0
    
    model.active_power_balance = Constraint(
        [(t, j) for t in Tset for j in Nm1set], rule=node_power_balance_rule
    )

    # Reactive Power Balance Constraints
    Tset = data['Tset']
    Nm1set = data['Nm1set']
    NLset = data['NLset']
    Dset = data['Dset']
    Bset = data['Bset']
    children = data['children']
    parent = data['parent']
    xdict_pu = data['xdict_pu']
    q_L_pu = data['q_L_pu']
    Q = model.Q
    l = model.l
    q_D = model.q_D
    q_B = model.q_B

    def reactive_power_balance_rule(model, t, j):
        sum_Qjk = sum(Q[j, k, t] for k in children.get(j, [])) if j in children else 0.0
        i = parent[j]
        Q_ij_t = Q[i, j, t]
        x_ij = xdict_pu[(i, j)]
        l_ij_t = l[i, j, t]
        line_reactive_loss = x_ij * l_ij_t
        q_L_j_t = q_L_pu[j][t-1] if j in NLset else 0.0
        q_D_j_t = q_D[j, t] if j in Dset else 0.0
        q_B_j_t = q_B[j, t] if j in Bset else 0.0

        return (
            sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t
            == 0
        )

    model.reactive_power_balance = Constraint(
        [(t, j) for t in Tset for j in Nm1set], rule=reactive_power_balance_rule
    )

    Tset = data['Tset']
    L1set = data['L1set']
    Lm1set = data['Lm1set']
    rdict_pu = data['rdict_pu']
    xdict_pu = data['xdict_pu']
    P = model.P
    Q = model.Q
    l = model.l
    v = model.v

    def kvl_substation_rule(model, t, i, j):
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P_ij_t = P[i, j, t]
        Q_ij_t = Q[i, j, t]
        l_ij_t = l[i, j, t]
        v_i_t = v[i, t]
        v_j_t = v[j, t]

        return (
            v_i_t
            - v_j_t
            - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t)
            + (r_ij**2 + x_ij**2) * l_ij_t
            == 0
        )

    model.kvl_substation = Constraint(
        [(t, i, j) for t in Tset for (i, j) in L1set], rule=kvl_substation_rule
    )

    def kvl_non_substation_rule(model, t, i, j):
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P_ij_t = P[i, j, t]
        Q_ij_t = Q[i, j, t]
        l_ij_t = l[i, j, t]
        v_i_t = v[i, t]
        v_j_t = v[j, t]

        return (
            v_i_t
            - v_j_t
            - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t)
            + (r_ij**2 + x_ij**2) * l_ij_t
            == 0
        )

    model.kvl_non_substation = Constraint(
        [(t, i, j) for t in Tset for (i, j) in Lm1set], rule=kvl_non_substation_rule
    )

    # Branch Complex Power Flow Equations (BCPF)
    Tset = data['Tset']
    L1set = data['L1set']
    Lm1set = data['Lm1set']
    P = model.P
    Q = model.Q
    v = model.v
    l = model.l

    def bcpf_substation_rule(model, t, i, j):
        P_ij_t = P[i, j, t]
        Q_ij_t = Q[i, j, t]
        v_i_t = v[i, t]
        l_ij_t = l[i, j, t]

        return P_ij_t**2 + Q_ij_t**2 - v_i_t * l_ij_t == 0

    model.bcpf_substation = Constraint(
        [(t, i, j) for t in Tset for (i, j) in L1set], rule=bcpf_substation_rule
    )

    def bcpf_non_substation_rule(model, t, i, j):
        P_ij_t = P[i, j, t]
        Q_ij_t = Q[i, j, t]
        v_i_t = v[i, t]
        l_ij_t = l[i, j, t]

        return P_ij_t**2 + Q_ij_t**2 - v_i_t * l_ij_t == 0

    model.bcpf_non_substation = Constraint(
        [(t, i, j) for t in Tset for (i, j) in Lm1set], rule=bcpf_non_substation_rule
    )

    # Battery SOC Trajectory Equality Constraints

    Bset = data['Bset']
    delta_t = data['delta_t']
    eta_C = data['eta_C']
    eta_D = data['eta_D']
    B0_pu = data['B0_pu']
    Bref_pu = data['Bref_pu']
    T = data['T']
    Tset = data['Tset']
    substationBus = data['substationBus']
    V_Subs = data['V_Subs']
    Compset = data['Compset']
    Vminpu_Comp = data['Vminpu_Comp']
    Vmaxpu_Comp = data['Vmaxpu_Comp']
    p_D_R_pu = data['p_D_R_pu']
    P_B_R_pu = data['P_B_R_pu']
    soc_min = data['soc_min']
    soc_max = data['soc_max']
    B_R_pu = data['B_R_pu']
    q_D = model.q_D
    q_B = model.q_B
    P_c = model.P_c
    P_d = model.P_d
    B = model.B
    v = model.v

    Δt = delta_t
    η_C = eta_C
    η_D = eta_D

    # Initial SOC constraint
    def initial_soc_rule(model, j):
        return (
            B[j, 1] - (B0_pu[j] + Δt * η_C[j] * P_c[j, 1] - Δt * (1 / η_D[j]) * P_d[j, 1])
            == 0
        )
    model.initial_soc = Constraint(Bset, rule=initial_soc_rule)

    # SOC trajectory for middle and final time periods
    def soc_trajectory_rule(model, j, t):
        return (
            B[j, t]
            - (B[j, t - 1] + Δt * η_C[j] * P_c[j, t] - Δt * (1 / η_D[j]) * P_d[j, t])
            == 0
        )
    model.soc_trajectory = Constraint(
        [(j, t) for j in Bset for t in range(2, T + 1)],
        rule=soc_trajectory_rule,
    )

    # Final SOC constraint
    def final_soc_rule(model, j):
        return B[j, T] == Bref_pu[j]
    model.final_soc = Constraint(Bset, rule=final_soc_rule)

    # Fixed substation node voltage
    def fixed_substation_voltage_rule(model, t):
        return v[substationBus, t] == V_Subs**2
    model.fixed_substation_voltage = Constraint(Tset, rule=fixed_substation_voltage_rule)

    # Voltage limits
    def voltage_lower_bound_rule(model, t, j):
        V_min_j_sq = Vminpu_Comp[j] ** 2
        return V_min_j_sq - v[j, t] <= 0
    model.voltage_lower_bound = Constraint(
        [(t, j) for t in Tset for j in Compset],
        rule=voltage_lower_bound_rule,
    )

    def voltage_upper_bound_rule(model, t, j):
        V_max_j_sq = Vmaxpu_Comp[j] ** 2
        return v[j, t] - V_max_j_sq <= 0
    model.voltage_upper_bound = Constraint(
        [(t, j) for t in Tset for j in Compset],
        rule=voltage_upper_bound_rule,
    )

    # Reactive Power Limits for PV Inverters
    def pv_reactive_power_lower_rule(model, t, j):
        p_D_R_j = p_D_R_pu[j]
        p_D_j_t = data['p_D_pu'][j][t-1]
        q_D_Max_j_t = sqrt((1.2 * p_D_R_j)**2 - (p_D_j_t)**2)
        return -q_D_Max_j_t - q_D[j, t] <= 0
    model.pv_reactive_power_lower = Constraint(
        [(t, j) for t in Tset for j in data['Dset']],
        rule=pv_reactive_power_lower_rule,
    )

    def pv_reactive_power_upper_rule(model, t, j):
        p_D_R_j = p_D_R_pu[j]
        p_D_j_t = data['p_D_pu'][j][t-1]
        q_D_Max_j_t = sqrt((1.2 * p_D_R_j)**2 - (p_D_j_t)**2)
        return q_D[j, t] - q_D_Max_j_t <= 0
    model.pv_reactive_power_upper = Constraint(
        [(t, j) for t in Tset for j in data['Dset']],
        rule=pv_reactive_power_upper_rule,
    )

    # Reactive Power Limits for Battery Inverters
    q_B_Max = {
        j: sqrt((1.2 * P_B_R_pu[j])**2 - (1.0 * P_B_R_pu[j])**2) for j in Bset
    }

    def battery_reactive_power_lower_rule(model, t, j):
        return -q_B_Max[j] - q_B[j, t] <= 0
    model.battery_reactive_power_lower = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=battery_reactive_power_lower_rule,
    )

    def battery_reactive_power_upper_rule(model, t, j):
        return q_B[j, t] - q_B_Max[j] <= 0
    model.battery_reactive_power_upper = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=battery_reactive_power_upper_rule,
    )

    # Charging Power Limits
    def charging_power_nonnegativity_rule(model, t, j):
        return -P_c[j, t] <= 0
    model.charging_power_nonnegativity = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=charging_power_nonnegativity_rule,
    )

    def charging_power_upper_rule(model, t, j):
        return P_c[j, t] - P_B_R_pu[j] <= 0
    model.charging_power_upper = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=charging_power_upper_rule,
    )

    # Discharging Power Limits
    def discharging_power_nonnegativity_rule(model, t, j):
        return -P_d[j, t] <= 0
    model.discharging_power_nonnegativity = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=discharging_power_nonnegativity_rule,
    )

    def discharging_power_upper_rule(model, t, j):
        return P_d[j, t] - P_B_R_pu[j] <= 0
    model.discharging_power_upper = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=discharging_power_upper_rule,
    )

    # SOC Limits for Batteries
    def soc_min_rule(model, t, j):
        return soc_min[j] * B_R_pu[j] - B[j, t] <= 0
    model.soc_min = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=soc_min_rule,
    )

    def soc_max_rule(model, t, j):
        return B[j, t] - soc_max[j] * B_R_pu[j] <= 0
    model.soc_max = Constraint(
        [(t, j) for t in Tset for j in Bset],
        rule=soc_max_rule,
    )

    # Unpack data
    objfun0 = data['objfun0']
    objfun2 = data['objfun2']
    Tset = data['Tset']
    Bset = data['Bset']
    eta_C = data['eta_C']
    eta_D = data['eta_D']
    LoadShapeCost = data['LoadShapeCost']
    kVA_B = data['kVA_B']
    delta_t = data['delta_t']
    alpha = data.get('alpha', 0)
    rdict_pu = data['rdict_pu']
    Lset = data['Lset']

    # Define constants
    C = LoadShapeCost
    η_C = eta_C
    η_D = eta_D

    # Initialize the objective function
    if objfun0 == "powerflow":
        # Set the objective function to zero for powerflow
        objfun = 0
    elif objfun0 == "subsPowerCostMin":
        # Define the base objective function (generation cost minimization)
        dollars_per_pu = kVA_B
        objfun = sum(
            dollars_per_pu * C[t-1] * model.P_Subs[t] * delta_t
            for t in Tset
        )
    elif objfun0 == "lineLossMin":
        objfun = sum(
            rdict_pu[(i, j)] * model.l[i, j, t]
            for (i, j) in Lset for t in Tset
        )
    else:
        raise ValueError("Unsupported objfun0 value")

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd":
        objfun += sum(
            alpha * ((1 - η_C[j]) * model.P_c[j, t] + (1 / η_D[j] - 1) * model.P_d[j, t])
            for j in Bset for t in Tset
        )

    # Add the objective function to the model
    model.obj = Objective(expr=objfun, sense=minimize)



    """
    Initialize variables with starting values.
    """
    Tset = data['Tset']
    Lset = data['Lset']
    Compset = data['Compset']
    Dset = data['Dset']
    Bset = data['Bset']
    B0_pu = data['B0_pu']
    V_Subs = data['V_Subs']

    # Initialize substation power flow variable
    for t in Tset:
        model.P_Subs[t].value = 0.0

    # Initialize power flow variables
    for (i, j) in Lset:
        for t in Tset:
            model.P[i, j, t].value = 0.0
            model.Q[i, j, t].value = 0.0
            model.l[i, j, t].value = 0.0

    # Initialize voltage variables
    for j in Compset:
        for t in Tset:
            model.v[j, t].value = V_Subs**2

    # Initialize PV inverter reactive dispatch variables
    for j in Dset:
        for t in Tset:
            model.q_D[j, t].value = 0.0

    # Initialize battery real and reactive dispatch variables
    for j in Bset:
        for t in Tset:
            model.q_B[j, t].value = 0.0
            model.P_c[j, t].value = 0.0
            model.P_d[j, t].value = 0.0

    # Initialize battery state of charge (SOC) variables
    for j in Bset:
        for t in Tset:
            model.B[j, t].value = B0_pu[j]

    """
    Solve the optimization model and retrieve results.
    """
    solver = SolverFactory(solver_name.lower())

    # Set solver options (example for IPOPT)
    if solver_name == "ipopt":
        solver.options['tol'] = 1e-6
        solver.options['max_iter'] = 10000
        solver.options['print_level'] = 5

    # Solve the model
    results = solver.solve(model, tee=True)

    # Check solver status
    if results.solver.termination_condition == TerminationCondition.optimal:
        print("Optimal solution found.")
        # Explicitly load results into the model
        model.solutions.load_from(results)
    else:
        print("Optimization did not find an optimal solution.")

    # Retrieve optimal objective value
    optimal_obj_value = value(model.obj)
    print("Optimal objective function value:", optimal_obj_value)

    return model






