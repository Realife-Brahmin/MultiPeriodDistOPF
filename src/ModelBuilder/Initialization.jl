module Initialization

export initialize_variables_1ph_NL_t_in_Tset

using JuMP
using Parameters: @unpack

function initialize_variables_1ph_NL_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    # Initialize substation power flow variable
    P_Subs = model[:P_Subs]
    for t in Tset
        set_start_value(P_Subs[t], 0.0)
    end

    # Initialize power flow variables
    @unpack Lset = data
    P = model[:P]
    Q = model[:Q]
    l = model[:l]
    for (i, j) in Lset, t in Tset
        set_start_value(P[(i, j), t], 0.0)
        set_start_value(Q[(i, j), t], 0.0)
        set_start_value(l[(i, j), t], 0.0)
    end

    # Initialize voltage variables
    @unpack Compset, V_Subs = data
    v = model[:v]
    for j in Compset, t in Tset
        set_start_value(v[j, t], (V_Subs)^2)
    end

    # Initialize PV inverter reactive dispatch variables
    @unpack Dset = data
    q_D = model[:q_D]
    for j in Dset, t in Tset
        set_start_value(q_D[j, t], 0.0)
    end

    # Initialize battery real and reactive dispatch variables
    @unpack Bset = data
    q_B = model[:q_B]
    P_c = model[:P_c]
    P_d = model[:P_d]
    for j in Bset, t in Tset
        set_start_value(q_B[j, t], 0.0)
        set_start_value(P_c[j, t], 0.0)
        set_start_value(P_d[j, t], 0.0)
    end

    # Initialize battery state of charge (SOC) variables
    @unpack B0_pu = data
    B = model[:B]
    for j in Bset, t in Tset
        set_start_value(B[j, t], B0_pu[j])
    end

    return model
end

end