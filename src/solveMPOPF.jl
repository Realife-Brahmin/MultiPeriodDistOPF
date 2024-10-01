# main.jl

# Include the parseOpenDSSFiles.jl script
include("./parseOpenDSSFiles.jl")

# Call the parsing function to get the data
N, Nset, Nm1set, Lset, L1set, Lm1set, r, x, Parent, Children,
T, Tset, C, η_C, η_D, V_base, v_min, v_max, p_L, q_L, Dset, p_D, Bset, battery_params = parseOpenDSSFiles()

# Import necessary packages
using JuMP
using Ipopt

# Define the optimization model
model = Model(Ipopt.Optimizer)

# ===========================
# Variables
# ===========================

# Define all variables as before, using the data parsed
@variable(model, P_Subs[t in Tset])

@variable(model, P[i in Nset, j in Nset, t in Tset], base_name = "P")
@variable(model, Q[i in Nset, j in Nset, t in Tset], base_name = "Q")
@variable(model, l[i in Nset, j in Nset, t in Tset] >= 0, base_name = "l")
@variable(model, v[j in Nset, t in Tset] >= 0, base_name = "v")
@variable(model, P_d[j in Bset, t in Tset] >= 0, base_name = "P_d")
@variable(model, P_c[j in Bset, t in Tset] >= 0, base_name = "P_c")
@variable(model, q_D[j in Dset, t in Tset], base_name = "q_D")
@variable(model, q_B[j in Bset, t in Tset], base_name = "q_B")

# B[j, t]: SOC of battery j at time t (for t = 1 to T-1)
@variable(model, B[j in Bset, t in 1:(T-1)] >= 0, base_name = "B")

# For now, set Bref[j] = B0[j] (desired final SOC equals initial SOC)
for j in Bset
    Bref[j] = B0[j]
end


# ===========================
# Constraints
# ===========================

# Implement all constraints as before, using the data and variables

# Substation node
SubstationNode = 1

## Real Power Balance Constraints ##

# Constraint h_1a: Nodal real power balance at the substation node
for t in Tset
    @constraint(model,
        P_Subs[t] - sum(P[SubstationNode, j, t] for (i, j) in L1set) == 0,
        "SubstationRealPowerBalance_t$(t)")
end

# Constraint h_1b_j: Nodal real power balance at non-substation nodes
for t in Tset, j in Nm1set
    # Sum of real powers flowing from node j to its children
    sum_Pjk = sum(P[j, k, t] for k in Children[j])

    # Parent node i of node j
    i = Parent[j]

    # Real power flow from parent i to node j
    P_ij_t = P[i, j, t]

    # Line losses on branch (i, j)
    r_ij = r[(i, j)]
    l_ij_t = l[i, j, t]
    line_loss = r_ij * l_ij_t

    # Load and PV generation at node j and time t
    p_L_j_t = p_L[j][t]
    p_D_j_t = haskey(p_D, j) ? p_D[j][t] : 0.0  # Check if node j has PV

    # Battery variables at node j and time t
    P_d_j_t = haskey(P_d, (j, t)) ? P_d[j, t] : 0.0
    P_c_j_t = haskey(P_c, (j, t)) ? P_c[j, t] : 0.0

    @constraint(model,
        sum_Pjk - (P_ij_t - line_loss) + p_L_j_t - p_D_j_t - (P_d_j_t - P_c_j_t) == 0,
        "NodeRealPowerBalance_Node$(j)_t$(t)")
end

## Reactive Power Balance Constraints ##

# Constraint h_2_j^t: Nodal reactive power balance at non-substation nodes
for t in Tset, j in Nm1set
    # Sum of reactive powers flowing from node j to its children
    sum_Qjk = sum(Q[j, k, t] for k in Children[j])

    # Parent node i of node j
    i = Parent[j]

    # Reactive power flow from parent i to node j
    Q_ij_t = Q[i, j, t]

    # Line reactive losses on branch (i, j)
    x_ij = x[(i, j)]
    l_ij_t = l[i, j, t]
    line_reactive_loss = x_ij * l_ij_t

    # Reactive load at node j and time t
    q_L_j_t = q_L[j][t]

    # Reactive power from PV inverter at node j and time t
    q_D_j_t = haskey(q_D, (j, t)) ? q_D[j, t] : 0.0

    # Reactive power from battery inverter at node j and time t
    q_B_j_t = haskey(q_B, (j, t)) ? q_B[j, t] : 0.0

    @constraint(model,
        sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0,
        "NodeReactivePowerBalance_Node$(j)_t$(t)")
end

## KVL Constraints ##

# Constraint h_3a: KVL for branches connected directly to the substation
for t in Tset, (i, j) in L1set
    r_ij = r[(i, j)]
    x_ij = x[(i, j)]
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
    r_ij = r[(i, j)]
    x_ij = x[(i, j)]
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
        "SOC_Initial_Node$(j)_t1")
end

# Constraint h_SOC_j^{t=2 to T-1}: SOC trajectory for middle time periods
for j in Bset, t in 2:(T-1)
    @constraint(model,
        B[j, t] - (B[j, t-1] + Δt * η_C * P_c[j, t] - Δt * (1 / η_D) * P_d[j, t]) == 0,
        "SOC_Trajectory_Node$(j)_t$(t)")
end

# Constraint h_SOC_j^{t=T}: Final SOC constraint
for j in Bset
    @constraint(model,
        Bref[j] - (B[j, T-1] + Δt * η_C * P_c[j, T] - Δt * (1 / η_D) * P_d[j, T]) == 0,
        "SOC_Final_Node$(j)_t$(T)")
end

# Voltage limits at each node
for t in Tset, j in Nset
    @constraint(model,
        v_min <= v[j, t] <= v_max,
        "VoltageLimits_Node$(j)_t$(t)")
end

# Fix substation voltage
for t in Tset
    @constraint(model, v[SubstationNode, t] == (V_base)^2)
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
