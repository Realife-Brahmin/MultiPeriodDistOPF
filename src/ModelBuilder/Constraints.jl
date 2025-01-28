module Constraints

export 
    battery_SOC_constraints_t_in_Tset,
    BCPF_non_substation_branches_t_in_Tset,
    BCPF_substation_branches_t_in_Tset,
    charging_power_limits_batteries_t_in_Tset,
    discharging_power_limits_batteries_t_in_Tset,
    fixed_substation_voltage_constraints_t_in_Tset,
    KVL_non_substation_branches_t_in_Tset,
    KVL_substation_branches_t_in_Tset,
    nodalReactivePowerBalance_non_substation_t_in_Tset,
    nodalRealPowerBalance_non_substation_t_in_Tset,
    nodalRealPowerBalance_substation_t_in_Tset,
    reactive_power_limits_PV_inverters_t_in_Tset,
    reactive_power_limits_battery_inverters_t_in_Tset,
    SOC_limits_batteries_t_in_Tset,
    voltage_limits_constraints_t_in_Tset

using JuMP
using Parameters: @unpack, @pack!

# Define all constraint functions here...

#region nodalRealPowerBalance_substation_t_in_Tset
"""
    nodalRealPowerBalance_substation_t_in_Tset(modelDict; Tset=nothing)

Define the nodal real power balance constraints for the substation node over a given time set.

This function sets the nodal real power balance constraints for the optimization model stored in `modelDict`.
"""
function nodalRealPowerBalance_substation_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack substationBus, L1set = data
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

    return modelDict
end
#endregion

#region nodalRealPowerBalance_non_substation_t_in_Tset
"""
    nodalRealPowerBalance_non_substation_t_in_Tset(modelDict; Tset=nothing)

Define the nodal real power balance constraints for non-substation nodes over a given time set.

This function sets the nodal real power balance constraints for the optimization model stored in `modelDict`.
"""
function nodalRealPowerBalance_non_substation_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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
        p_L_j_t = (j in NLset) ? p_L_pu[j, t] : 0.0  # Check if node j has a load

        # PV generation at node j and time t
        @unpack Dset, p_D_pu = data
        p_D_j_t = (j in Dset) ? p_D_pu[j, t] : 0.0  # Check if node j has PV

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

    return modelDict
end
#endregion

#region nodalReactivePowerBalance_non_substation_t_in_Tset
"""
    nodalReactivePowerBalance_non_substation_t_in_Tset(modelDict; Tset=nothing)

Define the nodal reactive power balance constraints for non-substation nodes over a given time set.

This function sets the nodal reactive power balance constraints for the optimization model stored in `modelDict`.
"""
function nodalReactivePowerBalance_non_substation_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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
        q_L_j_t = (j in NLset) ? q_L_pu[j, t] : 0.0  # Assign 0.0 if j is not in Nset

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

    return modelDict
end
#endregion

#region KVL_substation_branches_t_in_Tset
"""
    KVL_substation_branches_t_in_Tset(modelDict; Tset=nothing)

Define the Kirchhoff's Voltage Law (KVL) constraints for substation branches over a given time set.

This function sets the KVL constraints for the optimization model stored in `modelDict`.
"""
function KVL_substation_branches_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region KVL_non_substation_branches_t_in_Tset
"""
    KVL_non_substation_branches_t_in_Tset(modelDict; Tset=nothing)

Define the Kirchhoff's Voltage Law (KVL) constraints for non-substation branches over a given time set.

This function sets the KVL constraints for the optimization model stored in `modelDict`.
"""
function KVL_non_substation_branches_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region BCPF_substation_branches_t_in_Tset
"""
    BCPF_substation_branches_t_in_Tset(modelDict; Tset=nothing)

Define the Branch Complex Power Flow (BCPF) constraints for substation branches over a given time set.

This function sets the BCPF constraints for the optimization model stored in `modelDict`.
"""
function BCPF_substation_branches_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region BCPF_non_substation_branches_t_in_Tset
"""
    BCPF_non_substation_branches_t_in_Tset(modelDict; Tset=nothing)

Define the branch complex power flow (BCPF) constraints for non-substation branches over a given time set.

This function sets the BCPF constraints for the optimization model stored in `modelDict`.
"""
function BCPF_non_substation_branches_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region battery_SOC_constraints_t_in_Tset
"""
    battery_SOC_constraints_t_in_Tset(modelDict; Tset=nothing, tSOC_hard=false)

Define the state of charge (SOC) constraints for batteries over a given time set.

This function sets the SOC constraints for the optimization model stored in `modelDict`. 
It includes initial SOC constraints, SOC trajectory constraints, and optionally terminal SOC constraints if `tSOC_hard` is true.
"""
function battery_SOC_constraints_t_in_Tset(modelDict; Tset=nothing, tSOC_hard=false)
    @unpack model, data = modelDict;
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
    @unpack T = data # Note that Tset now no longer necessarily contains T i.e. Tset could be [2], with no knowledge of T=24
    for j in Bset, t in Tset
        if t > 1
            @constraint(model,
                base_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t)",
                B[j, t] - (B[j, t-1] + Δt * η_C[j] * P_c[j, t] - Δt * (1 / η_D[j]) * P_d[j, t]) == 0,
            )
        end
    end

    # Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = Bref_j)
    @unpack T = data
    if tSOC_hard && maximum(Tset) == T
        for j in Bset
            @constraint(model,
                base_name = "h_SOC_j^{T}_Terminal_SOC_Node_j_$(j)_t_$(T)",
                B[j, T] == Bref_pu[j],
            )
        end
    end

    return modelDict
end
#endregion

#region fixed_substation_voltage_constraints_t_in_Tset
"""
    fixed_substation_voltage_constraints_t_in_Tset(modelDict; Tset=nothing)

Define the fixed voltage constraints for the substation over a given time set.

This function sets the fixed voltage constraints for the optimization model stored in `modelDict`.
"""
function fixed_substation_voltage_constraints_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region voltage_limits_constraints_t_in_Tset
"""
    voltage_limits_constraints_t_in_Tset(modelDict; Tset=nothing)

Define the voltage limits constraints over a given time set.

This function sets the voltage limits constraints for the optimization model stored in `modelDict`.
"""
function voltage_limits_constraints_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region reactive_power_limits_PV_inverters_t_in_Tset
"""
    reactive_power_limits_PV_inverters_t_in_Tset(modelDict; Tset=nothing)

Define the reactive power limits for PV inverters over a given time set.

This function sets the reactive power limits for the optimization model stored in `modelDict`.
"""
function reactive_power_limits_PV_inverters_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack p_D_R_pu, Dset = data
    for t in Tset, j in Dset
        # Rated active power of the PV inverter at node j
        p_D_R_j = p_D_R_pu[j]

        # Active power output of PV at node j and time t
        @unpack p_D_pu = data
        p_D_j_t = p_D_pu[j, t]

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

    return modelDict
end
#endregion

#region reactive_power_limits_battery_inverters_1ph_NL_t_in_Tset
"""
    reactive_power_limits_battery_inverters_1ph_NL_t_in_Tset(modelDict; Tset=nothing)

Define the reactive power limits for battery inverters over a given time set using the constraint (P_d_j - P_c_j)^2 + q_B_j^2 <= S_B_R_j^2.

This function sets the reactive power limits for the optimization model stored in `modelDict`.
"""
function reactive_power_limits_battery_inverters_1ph_NL_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
    if Tset === nothing
        Tset = data[:Tset]
    end

    @unpack Bset, S_B_R_pu = data

    for t in Tset, j in Bset
        q_B = model[:q_B]
        P_d = model[:P_d]
        P_c = model[:P_c]

        ## g_qB_j^T: Reactive Power Limit for Battery Inverter ##
        @constraint(model,
            base_name = "g_qB_j^T_ReactivePowerLimit_Battery_Node_j_$(j)_t_$(t)",
            (P_d[j, t] - P_c[j, t])^2 + q_B[j, t]^2 <= S_B_R_pu[j]^2
        )
    end

    return modelDict
end
#endregion

#region reactive_power_limits_battery_inverters_t_in_Tset
"""
    reactive_power_limits_battery_inverters_t_in_Tset(modelDict; Tset=nothing)

Define the reactive power limits for battery inverters over a given time set.

This function sets the reactive power limits for the optimization model stored in `modelDict`.
"""
function reactive_power_limits_battery_inverters_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region charging_power_limits_batteries_t_in_Tset
"""
    charging_power_limits_batteries_t_in_Tset(modelDict; Tset=nothing)

Define the charging power limits for batteries over a given time set.

This function sets the charging power limits for the optimization model stored in `modelDict`.
"""
function charging_power_limits_batteries_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region discharging_power_limits_batteries_t_in_Tset
"""
    discharging_power_limits_batteries_t_in_Tset(modelDict; Tset=nothing)

Define the discharging power limits for batteries over a given time set.

This function sets the discharging power limits for the optimization model stored in `modelDict`.
"""
function discharging_power_limits_batteries_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

#region SOC_limits_batteries_t_in_Tset
"""
    SOC_limits_batteries_t_in_Tset(modelDict; Tset=nothing)

Define the state of charge (SOC) limits for batteries over a given time set.

This function sets the SOC limits for the optimization model stored in `modelDict`.
"""
function SOC_limits_batteries_t_in_Tset(modelDict; Tset=nothing)
    @unpack model, data = modelDict
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

    return modelDict
end
#endregion

end