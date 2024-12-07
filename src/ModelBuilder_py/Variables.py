from pyomo.environ import Var, NonNegativeReals

def define_model_variables_1ph_NL_t_in_Tset(model, data, Tset=None):
    if Tset is None:
        Tset = data['Tset']

    # Unpack the data
    Nset = data['Nset']
    Lset = data['Lset']
    Dset = data['Dset']
    Bset = data['Bset']
    PSubsMax_kW = data['PSubsMax_kW']
    kVA_B = data['kVA_B']
    PSubsMax_pu = PSubsMax_kW / kVA_B

    # Define all variables as before, using the data parsed
    model.P_Subs = Var(Tset, bounds=(0, PSubsMax_pu))

    # Define variables over the set of branches Lset and time periods Tset
    model.P = Var([(i, j, t) for (i, j) in Lset for t in Tset], bounds=(None, PSubsMax_pu), domain=NonNegativeReals)
    model.Q = Var([(i, j, t) for (i, j) in Lset for t in Tset], domain=NonNegativeReals)
    model.l = Var([(i, j, t) for (i, j) in Lset for t in Tset], domain=NonNegativeReals)

    model.v = Var([(j, t) for j in Nset for t in Tset], domain=NonNegativeReals)
    model.q_D = Var([(j, t) for j in Dset for t in Tset], domain=NonNegativeReals)
    model.q_B = Var([(j, t) for j in Bset for t in Tset], domain=NonNegativeReals)
    model.P_c = Var([(j, t) for j in Bset for t in Tset], domain=NonNegativeReals)
    model.P_d = Var([(j, t) for j in Bset for t in Tset], domain=NonNegativeReals)

    model.B = Var([(j, t) for j in Bset for t in Tset], domain=NonNegativeReals)

    # Handle SOC constraints for previous time steps
    if 1 not in Tset:
        t0m1 = min(Tset) - 1  # Assuming Tset is sorted
        model.B_t0m1 = Var([(j, t0m1) for j in Bset], domain=NonNegativeReals)

    return model
