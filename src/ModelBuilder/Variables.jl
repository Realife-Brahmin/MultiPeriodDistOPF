module Variables

export 
    define_model_variables_1ph_NL_t_in_Tset,
    define_model_variables_1ph_NL_t_in_Tset

using Parameters: @unpack, @pack!
using JuMP

#region define_model_variables_1ph_NL_t_in_Tset
"""
    define_model_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)

Define the optimization model variables for a single-phase network with nonlinear loads over a given time set.

This function defines the decision variables for the optimization model stored in `modelDict`. 
It handles different sets of variables, including power flows, voltages, and state of charge (SOC) for batteries.
"""
function define_model_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict;
    if Tset === nothing
        Tset = data[:Tset]
    end

    # @show Tset

    @unpack Nset, Lset, Dset, Bset, PSubsMax_kW, kVA_B = data
    PSubsMax_pu = PSubsMax_kW / kVA_B

    # Define all variables as before, using the data parsed
    @variable(model, PSubsMax_pu >= P_Subs[t in Tset] >= 0)
    # Define variables over the set of branches Lset and time periods Tset
    @variable(model, P[(i, j) in Lset, t in Tset] <= PSubsMax_pu, base_name = "P")
    @variable(model, Q[(i, j) in Lset, t in Tset], base_name = "Q")
    @variable(model, l[(i, j) in Lset, t in Tset] >= 0, base_name = "l")

    @variable(model, v[j in Nset, t in Tset], base_name = "v")
    @variable(model, q_D[j in Dset, t in Tset], base_name = "q_D")
    @variable(model, q_B[j in Bset, t in Tset], base_name = "q_B")
    @variable(model, P_c[j in Bset, t in Tset], base_name = "P_c")
    @variable(model, P_d[j in Bset, t in Tset], base_name = "P_d")

    # Since SOC constraints for the t-th hour require SOCs of the previous time-steps i.e. B[j, t-1], we need to ensure that they too (in addition to B[j, t]) are defined. This is important in the case where Tset = [t0] (i.e. only one time-step) and t0 != 1 (for t0=1, B[j, t-1] is not required as it is separately taken from B0 stored in data)
    if !(1 ∈ Tset)
        t0m1 = Tset[1] - 1 # assuming Tset is sorted, which it should be
        Tset_appended = [t0m1; Tset]
        @variable(model, B[j ∈ Bset, t ∈ Tset_appended], base_name = "B")
    elseif 1 ∈ Tset
        @variable(model, B[j ∈ Bset, t ∈ Tset], base_name = "B")
    else
        error("Tset must contain at least one element")
    end

    # @variable(model, B[j in Bset, t in Tset], base_name = "B")

    @pack! modelDict = model
    return modelDict
end
#endregion

#region define_model_variables_1ph_NL_t_in_Tset
"""
    define_model_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)

Define the optimization model variables for a single-phase network with nonlinear loads over a given time set.

This function defines the decision variables for the optimization model stored in `modelDict`. 
It handles different sets of variables, including power flows, voltages, and state of charge (SOC) for batteries.
"""
function define_model_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
    if Tset === nothing
        Tset = data[:Tset]
    end

    # @show Tset

    @unpack Nset, Lset, Dset, Bset, PSubsMax_kW, kVA_B = data
    PSubsMax_pu = PSubsMax_kW / kVA_B

    # Define all variables as before, using the data parsed
    @variable(model, PSubsMax_pu >= P_Subs[t in Tset] >= 0)
    # Define variables over the set of branches Lset and time periods Tset
    @variable(model, P[(i, j) in Lset, t in Tset] <= PSubsMax_pu, base_name = "P")
    @variable(model, Q[(i, j) in Lset, t in Tset], base_name = "Q")

    @variable(model, v[j in Nset, t in Tset], base_name = "v")
    @variable(model, q_D[j in Dset, t in Tset], base_name = "q_D")
    @variable(model, q_B[j in Bset, t in Tset], base_name = "q_B")
    @variable(model, P_c[j in Bset, t in Tset], base_name = "P_c")
    @variable(model, P_d[j in Bset, t in Tset], base_name = "P_d")

    # Since SOC constraints for the t-th hour require SOCs of the previous time-steps i.e. B[j, t-1], we need to ensure that they too (in addition to B[j, t]) are defined. This is important in the case where Tset = [t0] (i.e. only one time-step) and t0 != 1 (for t0=1, B[j, t-1] is not required as it is separately taken from B0 stored in data)
    if !(1 ∈ Tset)
        t0m1 = Tset[1] - 1 # assuming Tset is sorted, which it should be
        Tset_appended = [t0m1; Tset]
        @variable(model, B[j ∈ Bset, t ∈ Tset_appended], base_name = "B")
    elseif 1 ∈ Tset
        @variable(model, B[j ∈ Bset, t ∈ Tset], base_name = "B")
    else
        error("Tset must contain at least one element")
    end

    # @variable(model, B[j in Bset, t in Tset], base_name = "B")

    @pack! modelDict = model
    return modelDict
end
#endregion

end