from pyomo.environ import Objective, minimize, value
from src.ModelBuilder_py.Hyperparameters import *
from src.helperFunctions import *

def define_objective_function_t_in_Tset(model, data, Tset=None, tSOC_hard=False):
    if Tset is None:
        Tset = data['Tset']

    objfun0 = data['objfun0']
    objfun2 = data['objfun2']

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow":
        # Set the objective function to zero for powerflow
        objfun = 0
        func_obj_est = None  # no objective function (to minimize) for powerflow

    elif objfun0 == "subsPowerCostMin":
        # Define the base objective function (generation cost minimization)
        LoadShapeCost = data['LoadShapeCost']
        delta_t = data['delta_t']
        kVA_B = data['kVA_B']
        C = LoadShapeCost
        dollars_per_pu = kVA_B
        P_Subs = model.P_Subs
        objfun = sum(
            dollars_per_pu * C[t-1] * P_Subs[t] * delta_t
            for t in Tset
        )
        func_obj_est = estimate_substation_power_cost

    elif objfun0 == "lineLossMin":
        l = model.l
        Lset = data['Lset']
        rdict_pu = data['rdict_pu']
        kVA_B = data['kVA_B']

        objfun = kVA_B * sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset for t in Tset
        )
        func_obj_est = estimate_line_losses

    # Equivalent of @pack! data = func_obj_est
    data['func_obj_est'] = func_obj_est

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd":
        Bset = data['Bset']
        eta_C = data['eta_C']
        eta_D = data['eta_D']
        η_C = eta_C
        η_D = eta_D
        alpha = estimate_alpha(data)
        print(f"alpha = {alpha}")
        alphaAppendix = trim_number_for_printing(alpha)
        # Equivalent of @pack! data = alpha, alphaAppendix
        data['alpha'] = alpha
        data['alphaAppendix'] = alphaAppendix
        α = alpha
        P_c = model.P_c
        P_d = model.P_d

        objfun += sum(
            α * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset for t in Tset
        )

    # Append the gamma term only if tSOC_hard is false
    if not tSOC_hard:
        T = data['T']
        Bset = data['Bset']
        Bref_pu = data['Bref_pu']
        gamma = estimate_gamma(data)
        print(f"gamma = {gamma}")
        gammaAppendix = trim_number_for_printing(gamma)
        # Equivalent of @pack! data = gamma, gammaAppendix
        data['gamma'] = gamma
        data['gammaAppendix'] = gammaAppendix
        B = model.B
        γ = gamma

        objfun += sum(
            γ * (B[j, T] - Bref_pu[j])**2
            for j in Bset
        )

    # Equivalent of @objective(model, Min, objfun)
    model.obj = Objective(expr=objfun, sense=minimize)

    modelDict = {
        'model': model,
        'data': data
    }
    return modelDict
