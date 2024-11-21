module Variables

export define_model_variables_1ph_NL_t_in_Tset

using Parameters: @unpack
using JuMP

function define_model_variables_1ph_NL_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

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
    @variable(model, B[j in Bset, t in Tset], base_name = "B")

    return model
end

end