module Initialization

export 
    initialize_variables_1ph_NL_t_in_Tset,
    initialize_variables_1ph_L_t_in_Tset

using JuMP
using Parameters: @unpack, @pack!

#region initialize_variables_1ph_NL_t_in_Tset
"""
    initialize_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)

Initialize the optimization variables for a single-phase network with nonlinear loads over a given time set.

This function sets the initial values for the optimization variables in the model stored in `modelDict`. 
It handles the initialization of various variables, including substation power flow, power flow variables, voltage variables, PV inverter reactive dispatch, battery real and reactive dispatch, and battery state of charge (SOC).

# Arguments
- `modelDict::Dict`: A dictionary containing the model and data.
- `Tset::Union{Nothing, Vector{Int}}`: An optional vector specifying the time steps to consider. If not provided, it defaults to the full time set in `data`.

# Returns
- `modelDict::Dict`: The updated dictionary with initialized variables.
"""
function initialize_variables_1ph_NL_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict

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
    @unpack Nset, V_Subs_pu = data
    v = model[:v]
    for j in Nset, t in Tset
        set_start_value(v[j, t], (V_Subs_pu)^2)
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

    @pack! modelDict = model
    return modelDict
end
#endregion

#region initialize_variables_1ph_L_t_in_Tset
"""
    initialize_variables_1ph_L_t_in_Tset(modelDict; Tset=nothing)

Initialize the optimization variables for a single-phase network with linearized loads over a given time set.

This function sets the initial values for the optimization variables in the model stored in `modelDict`. 
It handles the initialization of various variables, including substation power flow, power flow variables, voltage variables, PV inverter reactive dispatch, battery real and reactive dispatch, and battery state of charge (SOC).

# Arguments
- `modelDict::Dict`: A dictionary containing the model and data.
- `Tset::Union{Nothing, Vector{Int}}`: An optional vector specifying the time steps to consider. If not provided, it defaults to the full time set in `data`.

# Returns
- `modelDict::Dict`: The updated dictionary with initialized variables.
"""
function initialize_variables_1ph_L_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict

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
    for (i, j) in Lset, t in Tset
        set_start_value(P[(i, j), t], 0.0)
        set_start_value(Q[(i, j), t], 0.0)
    end

    # Initialize voltage variables
    @unpack Nset, V_Subs_pu = data
    v = model[:v]
    for j in Nset, t in Tset
        set_start_value(v[j, t], (V_Subs_pu)^2)
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

    @pack! modelDict = model
    return modelDict
end
#endregion

end # module Initialization