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
# Define variables over the set of branches Lset and time periods Tset
@variable(model, P[(i, j) in Lset, t in Tset], base_name = "P")
@variable(model, Q[(i, j) in Lset, t in Tset], base_name = "Q")
# Todo How to model nonnegative l? Is it always automatically taken care of (speaking for a practical real world problem) even without the constraint? Should I just sneakily let it be as a JuMP constraint? Or should I expliclty define lower limit of l as an 'official' inequality constraint?
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

@unpack substationBus = data;
# Substation node
j1 = substationBus

## Real Power Balance Constraints ##

@unpack L1set = data;

for t in Tset
    @constraint(
        model,
        base_name = "SubstationRealPowerBalance_t_$(t)",
        P_Subs[t] - sum(P[(j1, j), t] for (j1, j) in L1set) == 0
    )
end

@unpack NLset, Nm1set, children, parent, rdict, xdict, p_L, p_D,  = data;
# Constraint h_1b_j: Nodal real power balance at non-substation nodes
for t in Tset, j in Nm1set
    # Sum of real powers flowing from node j to its children
    # sum_Pjk = sum(P[j, k, t] for k in children[j])
    sum_Pjk = isempty(children[j]) ? 0 : sum(P[(j, k), t] for k in children[j]) 

    # parent node i of node j
    i = parent[j]

    # Real power flow from parent i to node j
    P_ij_t = P[(i, j), t]

    # Line losses on branch (i, j)
    r_ij = rdict[(i, j)]
    l_ij_t = l[(i, j), t]
    line_loss = r_ij * l_ij_t

    # Load at node j and time t
    p_L_j_t = (j in NLset) ? p_L[j][t] : 0.0  # Check if node j has a load

    # PV generation at node j and time t
    p_D_j_t = (j in Dset) ? p_D[j][t] : 0.0  # Check if node j has PV

    # # Battery variables at node j and time t
    P_d_j_t = (j in Bset && t in Tset) ? P_d[j, t] : 0.0
    P_c_j_t = (j in Bset && t in Tset) ? P_c[j, t] : 0.0

    @constraint(model,
        sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0,
        base_name = "NodeRealPowerBalance_Node_j_$(j)_t_$(t)")
end

@unpack q_L = data;
## Nodal Reactive Power Balance Constraints ##
for t in Tset, j in Nm1set
    # Sum of reactive powers flowing from node j to its children
    sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])

    # parent node i of node j
    i = parent[j]

    # Reactive power flow from parent i to node j
    Q_ij_t = Q[(i, j), t]

    # Line reactive losses on branch (i, j)
    x_ij = xdict[(i, j)]
    l_ij_t = l[(i, j), t]
    line_reactive_loss = x_ij * l_ij_t

    # Todo: Figure out whether p_L[i][j] be allowed to exist or should it be made the same as p_L[i, j].
    # Reactive load at node j and time t
    q_L_j_t = (j in NLset) ? q_L[j][t] : 0.0  # Assign 0.0 if j is not in Nset

    # Reactive power from PV inverter at node j and time t
    q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0  # Assign 0.0 if j is not in Dset

    # Reactive power from battery inverter at node j and time t
    q_B_j_t = (j in Bset) ? q_B[j, t] : 0.0  # Assign 0.0 if j is not in BCPF_NonSubstationBranch_set

    ## h_2_j^t: Nodal Reactive Power Balance Constraint ##
    @constraint(model,
        base_name =     "h_2_j^t_NodeReactivePowerBalance_Node_j_$(j)_t_$(t)",
        sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0,
    )
end

## KVL Constraints ##

# Constraint h_3a: KVL for branches connected directly to the substation
for t in Tset, (i, j) in L1set
    r_ij = rdict[(i, j)]
    x_ij = xdict[(i, j)]
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

@unpack Lm1set = data;
# Constraint h_3b: KVL for branches not connected directly to the substation
for t in Tset, (i, j) in Lm1set
    r_ij = rdict[(i, j)]
    x_ij = xdict[(i, j)]
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

## Branch Complex Power Flow Equations (BCPF) ##

# Constraint h_4a_j^t: For branches connected directly to the substation
for t in Tset, (i, j) in L1set
    P_ij_t = P[(i, j), t]
    Q_ij_t = Q[(i, j), t]
    v_i_t = v[i, t]
    l_ij_t = l[(i, j), t]
    @constraint(model,
        base_name = "BCPF_SubstationBranch_i_$(i)_j_$(j)_t_$(t)",
        (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
    )
end

# Constraint h_4b_j^t: For branches not connected directly to the substation
for t in Tset, (i, j) in Lm1set
    P_ij_t = P[(i, j), t]
    Q_ij_t = Q[(i, j), t]
    v_i_t = v[i, t]
    l_ij_t = l[(i, j), t]
    @constraint(model,
        base_name = "BCPF_NonSubstationBranch_i_$(i)_j_$(j)_t_$(t)",
        (P_ij_t)^2 + (Q_ij_t)^2 - v_i_t * l_ij_t == 0,
    )
end

## Battery SOC Trajectory Equality Constraints ##

# Constraint h_SOC_j^{t=1}: Initial SOC constraint
@unpack delta_t, eta_C, eta_D, B0 = data;
Δt = delta_t
η_C = eta_C
η_D = eta_D

for j in Bset
    @constraint(model,
        base_name = "h_SOC_j^{t=1}_Initial_SOC_Node_j_$(j)_t1",
        B[j, 1] - (B0[j] + Δt * η_C[j] * P_c[j, 1] - Δt * (1 / η_D[j]) * P_d[j, 1]) == 0,
    )
end

# Constraint h_SOC_j^{t=2 to T}: SOC trajectory for middle and final time periods
for j in Bset, t in 2:T
    @constraint(model,
        base_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t)",
        B[j, t] - (B[j, t-1] + Δt * η_C[j] * P_c[j, t] - Δt * (1 / η_D[j]) * P_d[j, t]) == 0,
    )
end

@unpack Bref = data;
# Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = Bref_j)
for j in Bset
    @constraint(model,
        base_name = "h_SOC_j^{T}_Final_SOC_Node_j_$(j)_t_$(T)",
        B[j, T] == Bref[j],
    )
end

@unpack substationBus, V_Subs = data;
for t in Tset
    @constraint(model,
        base_name = "fixed_substation_node_j1_voltage_t_$(t)", 
        v[substationBus, t] == (V_Subs)^2,
    )
end

@unpack Compset, Vminpu_Comp, Vmaxpu_Comp = data;
# Voltage limits (using per-node limits if available)
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

    @constraint(model, 
        base_name = "g1_j^t_lower_voltage_bound_component_node_j_$(j)_t_$(t)",
        V_min_j_sq - v[j, t] <= 0 
    )
    @constraint(model,
        base_name = "g2_j^t_upper_voltage_bound_component_node_j_$(j)_t_$(t)",
        v[j, t] - V_max_j_sq <= 0
    )
    
end

@unpack p_D_R = data;
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
        base_name = "g_3_j^t_LowerReactivePowerLimit_PV_Node_j_$(j)_t_$(t)",
        -q_D_Max_j_t - q_D[j, t] <= 0,
    )

    ## g_4_j^t: Upper Limit of Reactive Power from PV Inverter ##
    @constraint(model,
        base_name = "g_4_j^t_UpperReactivePowerLimit_PV_Node_j_$(j)_t_$(t)",
        q_D[j, t] - q_D_Max_j_t <= 0,
    )
end

## Reactive Power Limits for Battery Inverters ##

@unpack P_B_R = data;
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
        base_name = "g_5_j^t_LowerReactivePowerLimit_Battery_Node_j_$(j)_t_$(t)",
        -q_B_Max[j] - q_B[j, t] <= 0
    )

    ## g_6_j^t: Upper Limit of Reactive Power from Battery Inverter ##
    @constraint(model,
        base_name = "g_6_j^t_UpperReactivePowerLimit_Battery_Node_j_$(j)_t_$(t)",
        q_B[j, t] - q_B_Max[j] <= 0,
    )
end

@unpack P_B_R = data;
## Charging Power Limits for Batteries ##
for t in Tset, j in Bset
    ## g_7_j^t: Non-negativity of Charging Power ##
    @constraint(model,
        base_name = "g_7_j^t_NonNegativity_ChargingPower_Node_j_$(j)_t_$(t)",
        - P_c[j, t] <= 0,
    )

    ## g_8_j^t: Maximum Charging Power Limit ##
    @constraint(model,
        base_name = "g_8_j^t_MaxChargingPowerLimit_Node_j_$(j)_t_$(t)",
        P_c[j, t] - P_B_R[j] <= 0,
    )
end

@unpack P_B_R = data;
## Discharging Power Limits for Batteries ##
for t in Tset, j in Bset
    ## g_9_j^t: Non-negativity of Discharging Power ##
    @constraint(model,
        base_name = "g_9_j^t_NonNegativity_DischargingPower_Node_j_$(j)_t_$(t)",
        - P_d[j, t] <= 0,
    )

    ## g_10_j^t: Maximum Discharging Power Limit ##
    @constraint(model,
        base_name = "g_10_j^t_MaxDischargingPowerLimit_Node_j_$(j)_t_$(t)",
        P_d[j, t] - P_B_R[j] <= 0,
    )
end

@unpack Bset, Tset, B_R, soc_min, soc_max = data;
## SOC Limits for Batteries ##

# Constraints:
for j in Bset, t in Tset
    ## g_11_j^t: Minimum SOC Constraint ##
    @constraint(model,
        base_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t)",
        soc_min[j] * B_R[j] - B[j, t] <= 0,
    )

    ## g_12_j^t: Maximum SOC Constraint ##
    @constraint(model,
        base_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t)",
        B[j, t] - soc_max[j] * B_R[j] <= 0,
    )
end

# ===========================
# Objective Function
# ===========================

alpha = 1e-3  # Adjust based on your problem requirements
@unpack Tset, Bset, eta_C, eta_D, LoadShapeCost = data;
C, η_C, η_D = LoadShapeCost, eta_C, eta_D
@objective(model, Min,
    sum(
        C[t] * P_Subs[t] +
        alpha * sum(
            (1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t]
            for j in Bset
        )
        for t in Tset
    )
)

# ===========================
# Initializing Variables
# ===========================

@unpack Tset, Nset, V_Subs = data;
# Initialize voltage variables
for j in Nset, t in Tset
    set_start_value(v[j, t], (V_Subs)^2)
end

@unpack Tset, Lset = data;
# Initialize power flow variables
for (i, j) in Lset, t in Tset
    set_start_value(P[(i, j), t], 0.0)
    set_start_value(Q[(i, j), t], 0.0)
    set_start_value(l[(i, j), t], 0.0)
end

@unpack Tset, Bset = data;
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
