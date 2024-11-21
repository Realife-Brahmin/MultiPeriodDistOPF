# optimizer.jl
module Playbook_of_MPOPF

export optimize_MPOPF_1ph_NL_TemporallyBruteforced

using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function define_model_variables_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Nset, Lset, Dset, Bset, PSubsMax_kW, kVA_B = data
    PSubsMax_pu = PSubsMax_kW / kVA_B

    # Define all variables as before, using the data parsed
    @variable(model, PSubsMax_pu >= P_Subs[t in Tset] >= 0)
    # @variable(model, P_Subs[t in Tset] <= PSubsMax_pu, base_name = "P_Subs")
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

function nodalRealPowerBalance_substation_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack substationBus, Tset, L1set = data
    P_Subs = model[:P_Subs]
    P = model[:P]
    # Substation node
    j1 = substationBus

    # Constraint h_1b_j: Nodal real power balance at non-substation nodes
    for t in Tset
        @constraint(
            model,
            base_name = "SubstationRealPowerBalance_t_$(t)",
            P_Subs[t] - sum(P[(j1, j), t] for (j1, j) in L1set) == 0
        )
    end

    return model
end

function nodalRealPowerBalance_non_substation_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack Nm1set = data
    for t in Tset, j in Nm1set
        # Sum of real powers flowing from node j to its children
        @unpack children = data
        P = model[:P]
        sum_Pjk = isempty(children[j]) ? 0 : sum(P[(j, k), t] for k in children[j])

        # parent node i of node j
        @unpack parent = data
        i = parent[j]

        # Real power flow from parent i to node j
        P = model[:P]
        P_ij_t = P[(i, j), t]

        # Line losses on branch (i, j)
        @unpack rdict_pu = data
        l = model[:l]
        r_ij = rdict_pu[(i, j)]
        l_ij_t = l[(i, j), t]
        line_loss = r_ij * l_ij_t

        # Load at node j and time t
        @unpack NLset, p_L_pu = data
        p_L_j_t = (j in NLset) ? p_L_pu[(j, t)] : 0.0  # Check if node j has a load

        # PV generation at node j and time t
        @unpack Dset, p_D_pu = data
        p_D_j_t = (j in Dset) ? p_D_pu[(j, t)] : 0.0  # Check if node j has PV

        # Battery variables at node j and time t
        @unpack Bset = data
        P_c = model[:P_c]
        P_d = model[:P_d]
        P_d_j_t = (j in Bset) ? P_d[j, t] : 0.0
        P_c_j_t = (j in Bset) ? P_c[j, t] : 0.0

        @constraint(model,
            sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0,
            base_name = "NodeRealPowerBalance_Node_j_$(j)_t_$(t)")
    end

    return model
end

function nodalReactivePowerBalance_non_substation_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack Nm1set = data
    for t in Tset, j in Nm1set
        # Sum of reactive powers flowing from node j to its children
        @unpack children = data
        Q = model[:Q]
        sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])

        # parent node i of node j
        @unpack parent = data
        i = parent[j]

        # Reactive power flow from parent i to node j
        Q_ij_t = Q[(i, j), t]

        # Line reactive losses on branch (i, j)
        @unpack xdict_pu = data
        l = model[:l]
        x_ij = xdict_pu[(i, j)]
        l_ij_t = l[(i, j), t]
        line_reactive_loss = x_ij * l_ij_t

        # Reactive load at node j and time t
        @unpack NLset, q_L_pu = data
        q_L_j_t = (j in NLset) ? q_L_pu[(j, t)] : 0.0  # Assign 0.0 if j is not in Nset

        # Reactive power from PV inverter at node j and time t
        @unpack Dset = data
        q_D = model[:q_D]
        q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0  # Assign 0.0 if j is not in Dset

        # Reactive power from battery inverter at node j and time t
        @unpack Bset = data
        q_B = model[:q_B]
        q_B_j_t = (j in Bset) ? q_B[j, t] : 0.0  # Assign 0.0 if j is not in BCPF_NonSubstationBranch_set

        # h_2_j^t: Nodal Reactive Power Balance Constraint #
        @constraint(model,
            base_name = "h_2_j^t_NodeReactivePowerBalance_Node_j_$(j)_t_$(t)",
            sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0,
        )
    end

    return model
end

function KVL_substation_branches_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack L1set = data

    # Constraint h_3a: KVL for branches connected directly to the substation
    for t in Tset, (i, j) in L1set
        @unpack rdict_pu, xdict_pu = data
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P = model[:P]
        Q = model[:Q]
        l = model[:l]
        v = model[:v]
        P_ij_t = P[(i, j), t]
        Q_ij_t = Q[(i, j), t]
        l_ij_t = l[(i, j), t]
        v_i_t = v[i, t]
        v_j_t = v[j, t]
        @constraint(model,
            base_name = "KVL_SubstationBranch_i_$(i)_j_$(j)_t_$(t)",
            v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij^2 + x_ij^2) * l_ij_t == 0,
        )
    end

    return model
end

function KVL_non_substation_branches_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack Lm1set = data

    # Constraint h_3b: KVL for branches not connected directly to the substation
    for t in Tset, (i, j) in Lm1set
        @unpack rdict_pu, xdict_pu = data
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
        P = model[:P]
        Q = model[:Q]
        l = model[:l]
        v = model[:v]
        P_ij_t = P[(i, j), t]
        Q_ij_t = Q[(i, j), t]
        l_ij_t = l[(i, j), t]
        v_i_t = v[i, t]
        v_j_t = v[j, t]
        @constraint(model,
            base_name = "KVL_NonSubstationBranch_i_$(i)_j_$(j)_t_$(t)",
            v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij^2 + x_ij^2) * l_ij_t == 0,
        )
    end

    return model
end

function BCPF_substation_branches_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack L1set = data

    # Constraint h_4a_j^t: For branches connected directly to the substation
    for t in Tset, (i, j) in L1set
        P = model[:P]
        Q = model[:Q]
        v = model[:v]
        l = model[:l]
        P_ij_t = P[(i, j), t]
        Q_ij_t = Q[(i, j), t]
        v_i_t = v[i, t]
        l_ij_t = l[(i, j), t]
        @constraint(model,
            base_name = "BCPF_SubstationBranch_i_$(i)_j_$(j)_t_$(t)",
            (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
        )
    end

    return model
end

function BCPF_non_substation_branches_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack Lm1set = data

    # Constraint h_4b_j^t: For branches not connected directly to the substation
    for t in Tset, (i, j) in Lm1set
        P = model[:P]
        Q = model[:Q]
        v = model[:v]
        l = model[:l]
        P_ij_t = P[(i, j), t]
        Q_ij_t = Q[(i, j), t]
        v_i_t = v[i, t]
        l_ij_t = l[(i, j), t]
        @constraint(model,
            base_name = "BCPF_NonSubstationBranch_i_$(i)_j_$(j)_t_$(t)",
            (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
        )
    end

    return model
end

function battery_SOC_constraints_t_in_Tset(model, data; Tset=nothing, tSOC_hard=false)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Bset, delta_t, eta_C, eta_D, B0_pu, Bref_pu = data
    Δt = delta_t
    η_C = eta_C
    η_D = eta_D
    P_c = model[:P_c]
    P_d = model[:P_d]
    B = model[:B]

    # Constraint h_SOC_j^{t=1}: Initial SOC constraint
    if 1 in Tset
        for j in Bset
            @constraint(model,
                base_name = "h_SOC_j^{t=1}_Initial_SOC_Node_j_$(j)_t1",
                B[j, 1] - (B0_pu[j] + Δt * η_C[j] * P_c[j, 1] - Δt * (1 / η_D[j]) * P_d[j, 1]) == 0,
            )
        end
    end

    # Constraint h_SOC_j^{t=2 to T}: SOC trajectory for middle and final time periods
    @unpack T = data
    for j in Bset, t in Tset
        if t > 1
            @constraint(model,
                base_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t)",
                B[j, t] - (B[j, t-1] + Δt * η_C[j] * P_c[j, t] - Δt * (1 / η_D[j]) * P_d[j, t]) == 0,
            )
        end
    end

    # Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = Bref_j)
    if tSOC_hard && maximum(Tset) == data[:T]
        for j in Bset
            @constraint(model,
                base_name = "h_SOC_j^{T}_Terminal_SOC_Node_j_$(j)_t_$(T)",
                B[j, data[:T]] == Bref_pu[j],
            )
        end
    end

    return model
end

function fixed_substation_voltage_constraints_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack substationBus, V_Subs = data
    for t in Tset
        v = model[:v]
        @constraint(model,
            base_name = "fixed_substation_node_j1_voltage_t_$(t)",
            v[substationBus, t] == (V_Subs)^2,
        )
    end

    return model
end

function voltage_limits_constraints_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Compset, Vminpu_Comp, Vmaxpu_Comp, substationBus, V_Subs = data
    for t in Tset, j in Compset
        if j == substationBus # not expected as per my parsing as substation should normally not have any components
            # Fix substation voltage
            @constraint(model,
                base_name = "fixed_voltage_substation_component_node_j1_t_$(t)",
                v[substationBus, t] == (V_Subs)^2,
            )
        end

        # Use per-node voltage limits if available
        V_min_j_sq = Vminpu_Comp[j]^2
        V_max_j_sq = Vmaxpu_Comp[j]^2
        v = model[:v]

        @constraint(model,
            base_name = "g1_j^t_lower_voltage_bound_component_node_j_$(j)_t_$(t)",
            V_min_j_sq - v[j, t] <= 0
        )
        @constraint(model,
            base_name = "g2_j^t_upper_voltage_bound_component_node_j_$(j)_t_$(t)",
            v[j, t] - V_max_j_sq <= 0
        )
    end

    return model
end

function reactive_power_limits_PV_inverters_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack p_D_R_pu, Dset = data
    for t in Tset, j in Dset
        # Rated active power of the PV inverter at node j
        p_D_R_j = p_D_R_pu[j]

        # Active power output of PV at node j and time t
        @unpack p_D_pu = data
        p_D_j_t = p_D_pu[(j, t)]

        # Compute q_D_Max_j^t
        q_D_Max_j_t = sqrt((1.2 * p_D_R_j)^2 - (p_D_j_t)^2)

        q_D = model[:q_D]

        # g_3_j^t: Lower Limit of Reactive Power from PV Inverter #
        @constraint(model,
            base_name = "g_3_j^t_LowerReactivePowerLimit_PV_Node_j_$(j)_t_$(t)",
            -q_D_Max_j_t - q_D[j, t] <= 0,
        )

        # g_4_j^t: Upper Limit of Reactive Power from PV Inverter #
        @constraint(model,
            base_name = "g_4_j^t_UpperReactivePowerLimit_PV_Node_j_$(j)_t_$(t)",
            q_D[j, t] - q_D_Max_j_t <= 0,
        )
    end

    return model
end

function reactive_power_limits_battery_inverters_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Bset, P_B_R_pu = data

    # Precompute q_B_Max_j for each battery inverter
    q_B_Max = Dict{Int,Float64}()

    for j in Bset
        # Rated active power of the battery inverter at node j
        P_B_R_j = P_B_R_pu[j]  # P_B_R_pu[j] should be provided (we'll handle data parsing later)

        # Compute q_B_Max_j
        q_B_Max_j = sqrt((1.2 * P_B_R_j)^2 - (1.0 * P_B_R_j)^2)

        # Store q_B_Max_j in the dictionary
        q_B_Max[j] = q_B_Max_j
    end

    for t in Tset, j in Bset
        q_B = model[:q_B]

        ## g_5_j^t: Lower Limit of Reactive Power from Battery Inverter ##
        @constraint(model,
            base_name = "g_5_j^t_LowerReactivePowerLimit_Battery_Node_j_$(j)_t_$(t)",
            -q_B_Max[j] - q_B[j, t] <= 0
        )

        ## g_6_j^t: Upper Limit of Reactive Power from Battery Inverter ##
        @constraint(model,
            base_name = "g_6_j^t_UpperReactivePowerLimit_Battery_Node_j_$(j)_t_$(t)",
            q_B[j, t] - q_B_Max[j] <= 0,
        )
    end

    return model
end

function charging_power_limits_batteries_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack P_B_R_pu, Bset = data
    for t in Tset, j in Bset
        P_c = model[:P_c]

        # g_7_j^t: Non-negativity of Charging Power #
        @constraint(model,
            base_name = "g_7_j^t_NonNegativity_ChargingPower_Node_j_$(j)_t_$(t)",
            -P_c[j, t] <= 0,
        )

        # g_8_j^t: Maximum Charging Power Limit #
        @constraint(model,
            base_name = "g_8_j^t_MaxChargingPowerLimit_Node_j_$(j)_t_$(t)",
            P_c[j, t] - P_B_R_pu[j] <= 0,
        )
    end

    return model
end

function discharging_power_limits_batteries_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack P_B_R_pu, Bset = data
    for t in Tset, j in Bset
        P_d = model[:P_d]

        # g_9_j^t: Non-negativity of Discharging Power #
        @constraint(model,
            base_name = "g_9_j^t_NonNegativity_DischargingPower_Node_j_$(j)_t_$(t)",
            -P_d[j, t] <= 0,
        )

        # g_10_j^t: Maximum Discharging Power Limit #
        @constraint(model,
            base_name = "g_10_j^t_MaxDischargingPowerLimit_Node_j_$(j)_t_$(t)",
            P_d[j, t] - P_B_R_pu[j] <= 0,
        )
    end

    return model
end

function SOC_limits_batteries_t_in_Tset(model, data; Tset=nothing)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Bset, B_R_pu, soc_min, soc_max = data
    for j in Bset, t in Tset
        B = model[:B]

        # g_11_j^t: Minimum SOC Constraint #
        @constraint(model,
            base_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t)",
            soc_min[j] * B_R_pu[j] - B[j, t] <= 0,
        )

        # g_12_j^t: Maximum SOC Constraint #
        @constraint(model,
            base_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t)",
            B[j, t] - soc_max[j] * B_R_pu[j] <= 0,
        )
    end

    return model
end

function define_objective_function_t_in_Tset(model, data; Tset=nothing, tSOC_hard=false)
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack objfun0, objfun2 = data;

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow"
        # Set the objective function to zero for powerflow
        objfun = 0
    elseif objfun0 == "subsPowerCostMin"
        # Define the base objective function (generation cost minimization)
        @unpack LoadShapeCost, delta_t, kVA_B, delta_t = data;
        P_Subs = model[:P_Subs]
        C = LoadShapeCost
        dollars_per_kWh = C
        dollars_per_puh = dollars_per_kWh * kVA_B
        objfun = sum(
            dollars_per_puh[t] * P_Subs[t] * delta_t
            for t in Tset
        )
    elseif objfun0 == "lineLossMin"
        @unpack rdict_pu, Lset = data;
        l = model[:l]
        objfun = sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset, t in Tset
        )
    end

    # Append the alpha term only if objfun2 == "scd"
    if objfun2 == "scd"
        @unpack eta_C, eta_D, alpha, Bset = data;
        P_c = model[:P_c]
        P_d = model[:P_d]
        α = alpha
        η_C = eta_C
        η_D = eta_D
        objfun += sum(
            α * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset, t in Tset
        )
    end

    # Append the gamma term only if tSOC_hard is false
    if !tSOC_hard
        @unpack T, Bset, Bref_pu, gamma = data;
        B = model[:B]
        γ = gamma
        objfun += sum(
            γ * (B[j, T] - Bref_pu[j])^2
            for j in Bset
        )
    end

    @objective(model, Min, objfun)

    return model
end

function initialize_variables_t_in_Tset(model, data; Tset=nothing)
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

function build_MPOPF_1ph_NL_model_t_1toT(data)
    @unpack solver = data

    # Define the optimization model including any specific solver settings
    model = configure_solver(solver)


    @unpack Tset = data;
    model = define_model_variables_t_in_Tset(model, data, Tset=Tset)

    # ===========================
    # Constraints
    # ===========================

    # Implement all constraints as before, using the data and variables

    #---Real Power Balance Constraints---#

    # Constraint h_1a_j: Nodal real power balance at substation node
    model = nodalRealPowerBalance_substation_t_in_Tset(model, data, Tset=Tset)

    # Constraint h_1b_j: Nodal real power balance at non-substation nodes
    model = nodalRealPowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)

    #--- Nodal Reactive Power Balance Constraints ---#

    # Constraint h_2_j: Nodal real power balance at non-substation nodes
    model = nodalReactivePowerBalance_non_substation_t_in_Tset(model, data, Tset=Tset)

    #---KVL Constraints---#

    # Constraint h_3a: KVL for branches connected directly to the substation
    model = KVL_substation_branches_t_in_Tset(model, data, Tset=Tset)

    # Constraint h_3b: KVL for branches not connected directly to the substation
    model = KVL_non_substation_branches_t_in_Tset(model, data)

    #---Branch Complex Power Flow Equations (BCPF)---#

    # Constraint h_4a_j^t: For branches connected directly to the substation
    model = BCPF_substation_branches_t_in_Tset(model, data)

    # Constraint h_4b_j^t: For branches not connected directly to the substation
    model = BCPF_non_substation_branches_t_in_Tset(model, data)

    #---Battery SOC Trajectory Equality Constraints---#

    # Constraints h_SOC_j^{t=1}, h_SOC_j^{t=2 to T} and h_SOC_j^{T} (if tSOC_hard is true)

    tSOC_hard = true # for brute force MPOPF as is the case in this function
    
    model = battery_SOC_constraints_t_in_Tset(model, data, tSOC_hard=tSOC_hard, Tset=Tset)

    #---Voltage Constraints---#

    # Fixed substation voltage constraint
    model = fixed_substation_voltage_constraints_t_in_Tset(model, data, Tset=Tset)

    # Voltage limits constraints
    model = voltage_limits_constraints_t_in_Tset(model, data, Tset=Tset)

    #---Reactive Power Limits for PV Inverters---#
    model = reactive_power_limits_PV_inverters_t_in_Tset(model, data, Tset=Tset)

    #---Reactive Power Limits for Battery Inverters---#
    model = reactive_power_limits_battery_inverters_t_in_Tset(model, data, Tset=Tset)

    #---Charging Power Limits for Batteries---#
    model = charging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    #---Discharging Power Limits for Batteries---#
    model = discharging_power_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    #---SOC Limits for Batteries---#
    model = SOC_limits_batteries_t_in_Tset(model, data, Tset=Tset)

    # ===========================
    # Objective Function
    # ===========================

    model = define_objective_function_t_in_Tset(model, data, Tset=Tset, tSOC_hard=tSOC_hard)

    # ===========================
    # Initializing Variables
    # ===========================
    model = initialize_variables_t_in_Tset(model, data)

    return model
end

function optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
    model = build_MPOPF_1ph_NL_model_t_1toT(data)

    optimize!(model)

    # Check solver status and retrieve results
    if termination_status(model) == LOCALLY_SOLVED
        println("Optimal solution found.")
    else
        println("Optimization did not find an optimal solution.")
    end

    optimal_obj_value = objective_value(model)
    println("Optimal objective function value: ", optimal_obj_value)

    return model

end

function configure_solver(solver_name)
    if solver_name == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "max_iter", 10000)
        set_optimizer_attribute(model, "print_level", 5)
    elseif solver_name == "Gurobi"
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", 300)        # Limit time (in seconds)
    elseif solver == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
        model = Model(optimizer)
    elseif solver == "EAGO"
        model = Model(EAGO.Optimizer)
    elseif solver == "MadNLP"
        model = Model(MadNLP.Optimizer)
    else
        error("Unsupported solver")
    end

    return model
end

function update_with_forwardStep_solutions!(model::Model, model_ddp_t::Model;
    verbose::Bool=false)
    # Iterate over all variables in the surrogate model
    for var_ddp in all_variables(model_ddp_t)
        # Extract the variable's value
        var_value = value(var_ddp)

        if isnothing(var_value)
            @warn "Variable $(name(var_ddp)) in surrogate model has no value. Skipping."
            continue
        end

        # Get the variable name
        var_name = name(var_ddp)

        # Check if the variable exists in the overarching model
        if haskey(model, var_name)
            # Retrieve the corresponding variable from the overarching model
            var_in_main_model = model[:$(var_name)]

            # Update the value in the main model
            set_start_value(var_in_main_model, var_value)  # Set the value in the main model
        else
            @warn "Variable $(var_name) does not exist in the overarching model. Skipping."
        end
    end

    return model
end

function build_ddpMPOPF_1ph_NL_model_t_is_T(ddpModel, data;
    verbose::Bool=false)

    @unpack k_ddp = ddpModel;
    if k_ddp == 1
        # cold start
        # copy the model locs here (maybe carve them into functions)
        # save the model as
        model_ddp_t0 = Model() # placeholder
    elseif k_ddp >= 2
        @unpack models_ddp_vs_t_vs_k = ddpModel;
        model_ddp_t0_km1 = models_ddp_vs_t_vs_k[(t0, k_ddp-1)]
        model_ddp_t0 = deepcopy(model_ddp_t0_km1)
        @unpack model = ddpModel; # because it has previous iteration's model's values saved
        B = model[:B]
        # modify hsoc equations (how to index them correctly?)
        # no need for using μ for terminal time-step, only modify objective function with f0, fscd, ftsoc
        model_ddp_t0 = Model() # placeholder
    else
        @error "Invalid value of k_ddp: $k_ddp"
    end

    models_ddp_vs_t_vs_k[t0, k_ddp] = model_ddp_t0 # this loc assumes that models_ddp_vs_t_vs_k is at least an already defined dictionary (even if empty) in ddpModel

    @pack! ddpModel = models_ddp_vs_t_vs_k;

    return ddpModel
end

end # module Playbook_of_MPOPF
