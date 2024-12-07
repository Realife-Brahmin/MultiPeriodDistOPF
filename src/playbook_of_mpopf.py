# Playbook_of_MPOPF.py

from pyomo.environ import *
from pyomo.opt import TerminationCondition
from src.ModelBuilder_py.ModelBuilder import *  # Assuming model_builder.py is the equivalent of ModelBuilder.jl
from src.ModelBuilder_py.Hyperparameters import estimate_alpha, estimate_gamma  # Placeholder for additional imports
from src.ModelBuilder_py.Initialization import initialize_variables_1ph_NL_t_in_Tset  # Placeholder for initialization
import crayons  # For colored terminal output (e.g., green, red text)


def generate_1ph_NL_model_decvar_value_dict(modelDict,results):
    model = modelDict['model']
    data = modelDict['data']
    modelVals = {}

    # Extract necessary sets from data
    Tset = data['Tset']
    Bset = data['Bset']
    Dset = data['Dset']
    Lset = data['Lset']
    Nset = data['Nset']

    # Initialize containers for each variable
    modelVals['P_Subs'] = {}  # t => value
    modelVals['P'] = {}  # (i, j), t => value
    modelVals['Q'] = {}
    modelVals['l'] = {}
    modelVals['v'] = {}  # (j, t) => value
    modelVals['q_D'] = {}
    modelVals['q_B'] = {}
    modelVals['P_c'] = {}
    modelVals['P_d'] = {}
    modelVals['B'] = {}

    # Store variable values
    for t in Tset:
        modelVals['P_Subs'][t] = value(model.P_Subs[t])

    for (i, j) in Lset:
        for t in Tset:
            modelVals['P'][(i, j), t] = value(model.P[(i, j), t])
            modelVals['Q'][(i, j), t] = value(model.Q[(i, j), t])
            modelVals['l'][(i, j), t] = value(model.l[(i, j), t])

    for j in Nset:
        for t in Tset:
            modelVals['v'][(j, t)] = value(model.v[(j, t)])

    for j in Dset:
        for t in Tset:
            modelVals['q_D'][(j, t)] = value(model.q_D[(j, t)])

    for j in Bset:
        for t in Tset:
            modelVals['q_B'][(j, t)] = value(model.q_B[(j, t)])
            modelVals['P_c'][(j, t)] = value(model.P_c[(j, t)])
            modelVals['P_d'][(j, t)] = value(model.P_d[(j, t)])
            modelVals['B'][(j, t)] = value(model.B[(j, t)])

    modelVals['objective_value'] = value(model.obj)
    modelVals['termination_status'] = results.solver.status
    modelVals['termination_condition'] = results.solver.termination_condition

    if 'solve_time' not in modelVals:
        modelVals['solve_time'] = results.solver.wallclock_time
    else:
        modelVals['solve_time'] += results.solver.wallclock_time

    modelDict['modelVals'] = modelVals
    return modelDict


def build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=None):
    solver = data['solver']

    # Define the optimization model with specific solver settings
    model = configure_solver(solver)

    if Tset is None:
        Tset = data['Tset']

    # Define variables
    model = define_model_variables_1ph_NL_t_in_Tset(model, data, Tset=Tset)

    # Constraints
    model = nodalRealPowerBalance_substation_t_in_Tset(model, data, Tset=Tset)
    model = nodalRealPowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)
    model = nodalReactivePowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)
    model = KVL_substation_branches_t_in_Tset(model, data, Tset=Tset)
    model = KVL_non_substation_branches_t_in_Tset(model, data, Tset=Tset)
    model = BCPF_substation_branches_t_in_Tset(model, data, Tset=Tset)
    model = BCPF_non_substation_branches_t_in_Tset(model, data, Tset=Tset)
    model = battery_SOC_constraints_t_in_Tset(model, data, Tset=Tset, tSOC_hard=data['tSOC_hard'])
    model = fixed_substation_voltage_constraints_t_in_Tset(model, data, Tset=Tset)
    model = voltage_limits_constraints_t_in_Tset(model, data, Tset=Tset)
    model = reactive_power_limits_PV_inverters_t_in_Tset(model, data, Tset=Tset)
    model = reactive_power_limits_battery_inverters_t_in_Tset(model, data, Tset=Tset)
    model = charging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)
    model = discharging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)
    model = SOC_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    # Objective function
    modelDict = define_objective_function_t_in_Tset(model, data, Tset=Tset, tSOC_hard=data['tSOC_hard'])

    # Initialize variables
    modelDict['model'] = initialize_variables_1ph_NL_t_in_Tset(modelDict['model'], data, Tset=Tset)

    return modelDict


def optimize_MPOPF_1ph_NL_TemporallyBruteforced(data):
    Tset = data['Tset']  # In this case, Tset = [1, 2, ..., T]
    modelDict = build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset)

    model = modelDict['model']
    solver = SolverFactory(data['solver'])
    results = solver.solve(model)

    modelDict = generate_1ph_NL_model_decvar_value_dict(modelDict,results)

    # Check solver status and retrieve results
    # green_crayon = crayons(foreground="light_green", bold=True)
    # red_crayon = crayons(foreground="red", bold=True)

    green_text = crayons.green
    red_text = crayons.red

    termination_status = modelDict['modelVals']['termination_status']
    termination_condition = modelDict['modelVals']['termination_condition']

    if termination_status == 'ok' and termination_condition == 'optimal':
        print(green_text("Optimal solution found.", bold=True))
    else:
        print(red_text("Optimization did not find an optimal solution.", bold=True))

    optimal_obj_value = modelDict['modelVals']['objective_value']
    print("Optimal objective function value: ", optimal_obj_value)

    return modelDict

def configure_solver(solver_name):
    model = None  # Initialize model as None for now

    if solver_name == "Ipopt":
        # For Ipopt, we use Pyomo with Ipopt as the solver
        model = ConcreteModel()
        solver = SolverFactory('ipopt')
        solver.options['tol'] = 1e-6  # Set tolerance
        solver.options['max_iter'] = 10000  # Set max iterations
        solver.options['print_level'] = 5  # Set print level
    elif solver_name == "Gurobi":
        # For Gurobi, we use Pyomo with Gurobi as the solver
        model = ConcreteModel()
        solver = SolverFactory('gurobi')
        solver.options['TimeLimit'] = 300  # Limit time in seconds
    elif solver_name == "Juniper":
        # Juniper-specific configuration
        # Using Ipopt as a sub-solver for Juniper optimization
        model = ConcreteModel()
        solver = SolverFactory('juniper')
        solver.options['nl_solver'] = SolverFactory('ipopt')
        solver.options['nl_solver.options']['print_level'] = 0  # Suppress output for Ipopt
    elif solver_name == "EAGO":
        # EAGO solver (hypothetical, no direct Python equivalent)
        model = ConcreteModel()
        solver = SolverFactory('eago')  # Assuming EAGO has a Python interface
    elif solver_name == "MadNLP":
        # MadNLP solver (hypothetical, no direct Python equivalent)
        model = ConcreteModel()
        solver = SolverFactory('madnlp')  # Assuming MadNLP has a Python interface
    else:
        raise ValueError("Unsupported solver")

    return model