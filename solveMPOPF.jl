# solveMPOPF.jl
using Revise
using MultiPeriodDistOPF
using Parameters: @unpack
using Debugger

Revise.revise()
systemName = "ads10_1ph"
T = 5

# Parse all data
data = parse_all_data(systemName, T)

# Import necessary packages
using JuMP
using Ipopt

# Define the optimization model
model = Model(Ipopt.Optimizer)

# ===========================
# Variables
# ===========================

@unpack Tset, Nset, Lset, Dset, Bset = data;

# Define all variables as before, using the data parsed
@variable(model, P_Subs[t in Tset])
# @variable(model, P[i in Nset, j in Nset, t in Tset], base_name = "P")
# @variable(model, Q[i in Nset, j in Nset, t in Tset], base_name = "Q")
# @variable(model, l[i in Nset, j in Nset, t in Tset] >= 0, base_name = "l")

# Define variables over the set of branches Lset and time periods Tset
@variable(model, P[(i, j) in Lset, t in Tset], base_name = "P")
@variable(model, Q[(i, j) in Lset, t in Tset], base_name = "Q")
# Todo How to model nonnegative l? Is it always automatically taken care of? Should I just sneakily let it be as a JuMP constraint? Or should I expliclty define lower limit of l as an 'official' inequality constraint?
@variable(model, l[(i, j) in Lset, t in Tset] >= 0, base_name = "l")

@variable(model, v[j in Nset, t in Tset], base_name = "v")
@variable(model, q_D[j in Dset, t in Tset], base_name = "q_D")
@variable(model, q_B[j in Bset, t in Tset], base_name = "q_B")
@variable(model, P_c[j in Bset, t in Tset], base_name = "P_c")
@variable(model, P_d[j in Bset, t in Tset], base_name = "P_d")
@variable(model, B[j in Bset, t in 1:T], base_name = "B")



# ===========================
# Constraints
# ===========================

# Implement all constraints as before, using the data and variables

# Substation node
SubstationNode = 1

## Real Power Balance Constraints ##

@unpack L1set = data;
# Constraint h_1a: Nodal real power balance at the substation node
# for t in Tset
#     # @constraint(model,
#     #     P_Subs[t] - sum(P[SubstationNode, j, t] for (i, j) in L1set) == 0,
#     #     "SubstationRealPowerBalance_t$(t)")
#     @constraint(model, P_Subs[t] - sum(P[SubstationNode, j, t] for (i, j) in L1set) == 0, "SubstationRealPowerBalance_$(t)")

# end
for t in Tset
@constraint(
    model,
    base_name = "SubstationRealPowerBalance_$(t)",
    P_Subs[t] - sum(P[SubstationNode, j, t] for (i, j) in L1set) == 0
)
end

@unpack NLset, Nm1set, children, parent, rdict, xdict, p_L, p_D,  = data;
# Constraint h_1b_j: Nodal real power balance at non-substation nodes
for t in Tset, j in Nm1set
    # Sum of real powers flowing from node j to its children
    # sum_Pjk = sum(P[j, k, t] for k in children[j])
    sum_Pjk = isempty(children[j]) ? 0 : sum(P[j, k, t] for k in children[j])

    # parent node i of node j
    i = parent[j]

    # Real power flow from parent i to node j
    P_ij_t = P[i, j, t]

    # Line losses on branch (i, j)
    r_ij = rdict[(i, j)]
    l_ij_t = l[i, j, t]
    line_loss = r_ij * l_ij_t

    # # Load and PV generation at node j and time t
    # p_L_j_t = p_L[j][t]
    # p_D_j_t = haskey(p_D, j) ? p_D[j][t] : 0.0  # Check if node j has PV

    # Load at node j and time t
    p_L_j_t = (j in NLset) ? p_L[j][t] : 0.0  # Check if node j has a load

    # PV generation at node j and time t
    p_D_j_t = (j in Dset) ? p_D[j][t] : 0.0  # Check if node j has PV

    # # Battery variables at node j and time t
    # P_d_j_t = haskey(P_d, (j, t)) ? P_d[j, t] : 0.0
    # P_c_j_t = haskey(P_c, (j, t)) ? P_c[j, t] : 0.0
    P_d_j_t = (j in Bset && t in Tset) ? P_d[j, t] : 0.0
    P_c_j_t = (j in Bset && t in Tset) ? P_c[j, t] : 0.0

    @constraint(model,
        sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0,
        base_name = "NodeRealPowerBalance_Node$(j)_t$(t)")
end

@unpack q_L = data;
## Nodal Reactive Power Balance Constraints ##
for t in Tset, j in Nm1set
    # # 
    # sum_Qjk = sum(Q[j, k, t] for k in children[j])

    # Sum of reactive powers flowing from node j to its children
    sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[j, k, t] for k in children[j])

    # parent node i of node j
    i = parent[j]

    # Reactive power flow from parent i to node j
    Q_ij_t = Q[i, j, t]

    # Line reactive losses on branch (i, j)
    x_ij = xdict[(i, j)]
    l_ij_t = l[i, j, t]
    line_reactive_loss = x_ij * l_ij_t

    # # Reactive load at node j and time t
    # q_L_j_t = q_L[j][t]
    # Reactive load at node j and time t
    q_L_j_t = (j in NLset) ? q_L[j][t] : 0.0  # Assign 0.0 if j is not in Nset

    # # Reactive power from PV inverter at node j and time t
    # q_D_j_t = q_D[j, t]
    # Reactive power from PV inverter at node j and time t
    q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0  # Assign 0.0 if j is not in Dset

    # # Reactive power from battery inverter at node j and time t
    # q_B_j_t = q_B[j, t]
    # Reactive power from battery inverter at node j and time t
    q_B_j_t = (j in Bset) ? q_B[j, t] : 0.0  # Assign 0.0 if j is not in BCPF_NonSubstationBranch_set

    # # Static reactive power injection from capacitor bank at node j and time t
    # q_C_j_t = q_C[j][t]  # q_C_j_t is a parameter (from data)

    # ## h_2_j^t: Nodal Reactive Power Balance Constraint ##
    # @constraint(model,
    #     sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t - q_C_j_t == 0,
    #     "h_2_j^t_NodeReactivePowerBalance_Node$(j)_t$(t)")

    ## h_2_j^t: Nodal Reactive Power Balance Constraint ##
    @constraint(model,
        sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0,
        base_name = "h_2_j^t_NodeReactivePowerBalance_Node$(j)_t$(t)")
end

## KVL Constraints ##

# Constraint h_3a: KVL for branches connected directly to the substation
for t in Tset, (i, j) in L1set
    r_ij = rdict[(i, j)]
    x_ij = xdict[(i, j)]
    P_ij_t = P[i, j, t]
    Q_ij_t = Q[i, j, t]
    l_ij_t = l[i, j, t]
    v_i_t = v[i, t]
    v_j_t = v[j, t]
    @constraint(model,
        v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij^2 + x_ij^2) * l_ij_t == 0,
        "KVL_SubstationBranch_$(i)_$(j)_t$(t)")
end

# Constraint h_3b: KVL for branches not connected directly to the substation
for t in Tset, (i, j) in Lm1set
    r_ij = rdict[(i, j)]
    x_ij = xdict[(i, j)]
    P_ij_t = P[i, j, t]
    Q_ij_t = Q[i, j, t]
    l_ij_t = l[i, j, t]
    v_i_t = v[i, t]
    v_j_t = v[j, t]
    @constraint(model,
        v_i_t - v_j_t - 2 * (r_ij * P_ij_t + x_ij * Q_ij_t) + (r_ij^2 + x_ij^2) * l_ij_t == 0,
        "KVL_NonSubstationBranch_$(i)_$(j)_t$(t)")
end

## Branch Complex Power Flow Equations (BCPF) ##

# Constraint h_4a_j^t: For branches connected directly to the substation
for t in Tset, (i, j) in L1set
    P_ij_t = P[i, j, t]
    Q_ij_t = Q[i, j, t]
    v_i_t = v[i, t]
    l_ij_t = l[i, j, t]
    @constraint(model,
        (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
        "BCPF_SubstationBranch_$(i)_$(j)_t$(t)")
end

# Constraint h_4b_j^t: For branches not connected directly to the substation
for t in Tset, (i, j) in Lm1set
    P_ij_t = P[i, j, t]
    Q_ij_t = Q[i, j, t]
    v_i_t = v[i, t]
    l_ij_t = l[i, j, t]
    @constraint(model,
        (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
        "BCPF_NonSubstationBranch_$(i)_$(j)_t$(t)")
end

## Battery SOC Trajectory Equality Constraints ##

# Constraint h_SOC_j^{t=1}: Initial SOC constraint
for j in Bset
    @constraint(model,
        B[j, 1] - (B0[j] + Δt * η_C * P_c[j, 1] - Δt * (1 / η_D) * P_d[j, 1]) == 0,
        "h_SOC_j^{t=1}_Initial_SOC_Node$(j)_t1")
end

# Constraint h_SOC_j^{t=2 to T}: SOC trajectory for middle and final time periods
for j in Bset, t in 2:T
    @constraint(model,
        B[j, t] - (B[j, t-1] + Δt * η_C * P_c[j, t] - Δt * (1 / η_D) * P_d[j, t]) == 0,
        "h_SOC_j^{t=$(t)}_SOC_Trajectory_Node$(j)_t$(t)")
end

# Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = B_ref_j)
for j in Bset
    @constraint(model,
        B[j, T] == B_ref[j],
        "h_SOC_j^{T}_Final_SOC_Node$(j)_t$(T)")
end


# Voltage variables
@variable(model, v[j in Nset, t in Tset] >= 0, base_name = "v")

# Voltage limits (using per-node limits if available)
# Todo: what about voltage limits for buses with no loads?
for t in Tset, j in Nset
    if j == substationBus
        # Fix substation voltage
        @constraint(model, v[j, t] == (V_Subs)^2, "fixed_substation_voltage_t_$(t)")
    else
        # Use per-node voltage limits if available
        V_min_j_sq = (V_minpu[j] * kv_bus[j])^2
        V_max_j_sq = (V_maxpu[j] * kv_bus[j])^2

        @constraint(model, V_min_j_sq - v[j, t] <= 0, "g1_j^t_lower_voltage_bound_node_$(j)_time_$(t)")
        @constraint(model, v[j, t] - V_max_j_sq <= 0, "g2_j^t_upper_voltage_bound_node_$(j)_time_$(t)")
    end
end


## Reactive Power Limits for PV Inverters ##
for t in Tset, j in Dset
    # Rated active power of the PV inverter at node j
    p_D_R_j = p_D_R[j]

    # Active power output of PV at node j and time t
    p_D_j_t = p_D[j][t]

    # Compute q_D_Max_j^t
    q_D_Max_j_t = sqrt((1.2 * p_D_R_j)^2 - (p_D_j_t)^2)

    ## g_3_j^t: Lower Limit of Reactive Power from PV Inverter ##
    @constraint(model,
        -q_D_Max_j_t - q_D[j, t] <= 0,
        "g_3_j^t_LowerReactivePowerLimit_PV_Node$(j)_t$(t)")

    ## g_4_j^t: Upper Limit of Reactive Power from PV Inverter ##
    @constraint(model,
        q_D[j, t] - q_D_Max_j_t <= 0,
        "g_4_j^t_UpperReactivePowerLimit_PV_Node$(j)_t$(t)")
end

## Reactive Power Limits for Battery Inverters ##

# Precompute q_B_Max_j for each battery inverter
q_B_Max = Dict{Int,Float64}()

for j in Bset
    # Rated active power of the battery inverter at node j
    P_B_R_j = P_B_R[j]  # P_B_R[j] should be provided (we'll handle data parsing later)

    # Compute q_B_Max_j
    q_B_Max_j = sqrt((1.2 * P_B_R_j)^2 - (1.0 * P_B_R_j)^2)

    # Store q_B_Max_j in the dictionary
    q_B_Max[j] = q_B_Max_j
end

# Define constraints for each time period and battery inverter
for t in Tset, j in Bset
    ## g_5_j^t: Lower Limit of Reactive Power from Battery Inverter ##
    @constraint(model,
        -q_B_Max[j] - q_B[j, t] <= 0,
        "g_5_j^t_LowerReactivePowerLimit_Battery_Node$(j)_t$(t)")

    ## g_6_j^t: Upper Limit of Reactive Power from Battery Inverter ##
    @constraint(model,
        q_B[j, t] - q_B_Max[j] <= 0,
        "g_6_j^t_UpperReactivePowerLimit_Battery_Node$(j)_t$(t)")
end

## Charging Power Limits for Batteries ##
for t in Tset, j in Bset
    ## g_7_j^t: Non-negativity of Charging Power ##
    @constraint(model,
        -P_c[j, t] <= 0,
        "g_7_j^t_NonNegativity_ChargingPower_Node$(j)_t$(t)")

    ## g_8_j^t: Maximum Charging Power Limit ##
    @constraint(model,
        P_c[j, t] - P_B_R[j] <= 0,
        "g_8_j^t_MaxChargingPowerLimit_Node$(j)_t$(t)")
end

## Discharging Power Limits for Batteries ##
for t in Tset, j in Bset
    ## g_9_j^t: Non-negativity of Discharging Power ##
    @constraint(model,
        -P_d[j, t] <= 0,
        "g_9_j^t_NonNegativity_DischargingPower_Node$(j)_t$(t)")

    ## g_10_j^t: Maximum Discharging Power Limit ##
    @constraint(model,
        P_d[j, t] - P_B_R[j] <= 0,
        "g_10_j^t_MaxDischargingPowerLimit_Node$(j)_t$(t)")
end

## SOC Limits for Batteries ##

# Parameters:
# soc_min: Minimum SOC fraction (e.g., 0.2 for 20%)
# soc_max: Maximum SOC fraction (e.g., 0.9 for 90%)
# B_R[j]: Rated capacity of battery j (from data or assumptions)

# Constraints:
for j in Bset, t in 1:T
    ## g_11_j^t: Minimum SOC Constraint ##
    @constraint(model,
        soc_min * B_R[j] - B[j, t] <= 0,
        "g_11_j^t_MinSOC_Node$(j)_t$(t)")

    ## g_12_j^t: Maximum SOC Constraint ##
    @constraint(model,
        B[j, t] - soc_max * B_R[j] <= 0,
        "g_12_j^t_MaxSOC_Node$(j)_t$(t)")
end

# ===========================
# Objective Function
# ===========================

alpha = 1e-3  # Adjust based on your problem requirements

@objective(model, Min,
    sum(
        C[t] * P_Subs[t] +
        alpha * sum(
            (1 - η_C) * P_c[j, t] + (1 / η_D - 1) * P_d[j, t]
            for j in Bset
        )
        for t in Tset
    )
)

# ===========================
# Initializing Variables
# ===========================

# Initialize voltage variables
for j in Nset, t in Tset
    set_start_value(v[j, t], (V_base)^2)
end

# Initialize power flow variables
for (i, j) in Lset, t in Tset
    set_start_value(P[i, j, t], 0.0)
    set_start_value(Q[i, j, t], 0.0)
    set_start_value(l[i, j, t], 0.0)
end

# Initialize battery variables
for j in Bset, t in Tset
    set_start_value(P_d[j, t], 0.0)
    set_start_value(P_c[j, t], 0.0)
    set_start_value(q_B[j, t], 0.0)
end

# Initialize PV inverter variables
for j in Dset, t in Tset
    set_start_value(q_D[j, t], 0.0)
end

# ===========================
# Solver Settings
# ===========================

# Set IPOPT options
set_optimizer_attribute(model, "tol", 1e-6)
set_optimizer_attribute(model, "max_iter", 10000)
set_optimizer_attribute(model, "print_level", 5)

# ===========================
# Solve the Model
# ===========================

optimize!(model)

# ===========================
# Retrieve Results
# ===========================

# Check solver status and retrieve results
if termination_status(model) == MOI.OPTIMAL
    println("Optimal solution found.")
    # Retrieve and process results as needed
else
    println("Optimization did not find an optimal solution.")
end
