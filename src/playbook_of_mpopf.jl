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

function nodalRealPowerBalance_substation_t_in_Tset(model, data)
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

function nodalRealPowerBalance_non_substation_t_in_Tset(model, data)
    @unpack Tset, Nm1set = data
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

    # ## Real Power Balance Constraints ##

    # Constraint h_1a_j: Nodal real power balance at substation node

    model = nodalRealPowerBalance_substation_t_in_Tset(model, data)

    # Constraint h_1b_j: Nodal real power balance at non-substation nodes

    model = nodalRealPowerBalance_non_substation_t_in_Tset(model, data)

    # ## Nodal Reactive Power Balance Constraints ##

    # @unpack Tset, Nm1set, NLset, Dset, Bset, children, parent, xdict_pu, q_L_pu = data
    # for t in Tset, j in Nm1set
    #     # Sum of reactive powers flowing from node j to its children
    #     sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q[(j, k), t] for k in children[j])

    #     # parent node i of node j
    #     i = parent[j]

    #     # Reactive power flow from paren t i to node j
    #     Q_ij_t = Q[(i, j), t]

    #     # Line reactive losses on branch (i, j)
    #     x_ij = xdict_pu[(i, j)]
    #     l_ij_t = l[(i, j), t]
    #     line_reactive_loss = x_ij * l_ij_t

    #     # Reactive load at node j and time t
    #     q_L_j_t = (j in NLset) ? q_L_pu[(j, t)] : 0.0  # Assign 0.0 if j is not in Nset

    #     # Reactive power from PV inverter at node j and time t
    #     q_D_j_t = (j in Dset) ? q_D[j, t] : 0.0  # Assign 0.0 if j is not in Dset

    #     # Reactive power from battery inverter at node j and time t
    #     q_B_j_t = (j in Bset) ? q_B[j, t] : 0.0  # Assign 0.0 if j is not in BCPF_NonSubstationBranch_set

    #     ## h_2_j^t: Nodal Reactive Power Balance Constraint ##
    #     @constraint(model,
    #         base_name = "h_2_j^t_NodeReactivePowerBalance_Node_j_$(j)_t_$(t)",
    #         sum_Qjk - (Q_ij_t - line_reactive_loss) + q_L_j_t - q_D_j_t - q_B_j_t == 0,
    #     )
    # end

    ## KVL Constraints ##

    @unpack Tset, L1set, rdict_pu, xdict_pu = data
    # Constraint h_3a: KVL for branches connected directly to the substation
    for t in Tset, (i, j) in L1set
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
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

    @unpack Tset, Lm1set, rdict_pu, xdict_pu = data
    # Constraint h_3b: KVL for branches not connected directly to the substation
    for t in Tset, (i, j) in Lm1set
        r_ij = rdict_pu[(i, j)]
        x_ij = xdict_pu[(i, j)]
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

    @unpack Tset, L1set = data
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

    @unpack Tset, Lm1set = data
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
    @unpack Bset, delta_t, eta_C, eta_D, B0_pu = data
    Δt = delta_t
    η_C = eta_C
    η_D = eta_D

    for j in Bset
        @constraint(model,
            base_name = "h_SOC_j^{t=1}_Initial_SOC_Node_j_$(j)_t1",
            B[j, 1] - (B0_pu[j] + Δt * η_C[j] * P_c[j, 1] - Δt * (1 / η_D[j]) * P_d[j, 1]) == 0,
        )
    end

    # Constraint h_SOC_j^{t=2 to T}: SOC trajectory for middle and final time periods
    @unpack T = data
    for j in Bset, t in 2:T
        @constraint(model,
            base_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t)",
            B[j, t] - (B[j, t-1] + Δt * η_C[j] * P_c[j, t] - Δt * (1 / η_D[j]) * P_d[j, t]) == 0,
        )
    end

    @unpack Bref_pu = data
    # Constraint h_SOC_j^{T}: Final SOC constraint (B_j^T = Bref_j)
    for j in Bset
        @constraint(model,
            base_name = "h_SOC_j^{T}_Final_SOC_Node_j_$(j)_t_$(T)",
            B[j, T] == Bref_pu[j],
        )
    end

    @unpack substationBus, V_Subs = data
    for t in Tset
        @constraint(model,
            base_name = "fixed_substation_node_j1_voltage_t_$(t)",
            v[substationBus, t] == (V_Subs)^2,
        )
    end

    @unpack Tset, Compset, Vminpu_Comp, Vmaxpu_Comp = data
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

    @unpack p_D_R_pu = data
    ## Reactive Power Limits for PV Inverters ##
    for t in Tset, j in Dset
        # Rated active power of the PV inverter at node j
        p_D_R_j = p_D_R_pu[j]

        # Active power output of PV at node j and time t
        p_D_j_t = p_D_pu[(j, t)]

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

    @unpack Bset, P_B_R_pu = data
    ## Reactive Power Limits for Battery Inverters ##

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

    @unpack P_B_R_pu = data
    ## Charging Power Limits for Batteries ##
    for t in Tset, j in Bset
        ## g_7_j^t: Non-negativity of Charging Power ##
        @constraint(model,
            base_name = "g_7_j^t_NonNegativity_ChargingPower_Node_j_$(j)_t_$(t)",
            -P_c[j, t] <= 0,
        )

        ## g_8_j^t: Maximum Charging Power Limit ##
        @constraint(model,
            base_name = "g_8_j^t_MaxChargingPowerLimit_Node_j_$(j)_t_$(t)",
            P_c[j, t] - P_B_R_pu[j] <= 0,
        )
    end

    @unpack P_B_R_pu = data
    ## Discharging Power Limits for Batteries ##
    for t in Tset, j in Bset
        ## g_9_j^t: Non-negativity of Discharging Power ##
        @constraint(model,
            base_name = "g_9_j^t_NonNegativity_DischargingPower_Node_j_$(j)_t_$(t)",
            -P_d[j, t] <= 0,
        )

        ## g_10_j^t: Maximum Discharging Power Limit ##
        @constraint(model,
            base_name = "g_10_j^t_MaxDischargingPowerLimit_Node_j_$(j)_t_$(t)",
            P_d[j, t] - P_B_R_pu[j] <= 0,
        )
    end

    @unpack Bset, Tset, B_R_pu, soc_min, soc_max = data
    ## SOC Limits for Batteries ##

    # Constraints:
    for j in Bset, t in Tset
        ## g_11_j^t: Minimum SOC Constraint ##
        @constraint(model,
            base_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t)",
            soc_min[j] * B_R_pu[j] - B[j, t] <= 0,
        )

        ## g_12_j^t: Maximum SOC Constraint ##
        @constraint(model,
            base_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t)",
            B[j, t] - soc_max[j] * B_R_pu[j] <= 0,
        )
    end

    # ===========================
    # Objective Function
    # ===========================

    @unpack objfun0, objfun2, Tset, Bset, eta_C, eta_D, LoadShapeCost, kVA_B = data
    C, η_C, η_D = LoadShapeCost, eta_C, eta_D

    # Assume objfun0 and objfun2 are passed to the function that defines the model.
    if objfun0 == "powerflow"
        # Set the objective function to zero for powerflow
        objfun = 0
    elseif objfun0 == "subsPowerCostMin"
        # Define the base objective function (generation cost minimization)
        dollars_per_pu = kVA_B
        objfun = sum(
            dollars_per_pu * C[t] * P_Subs[t] * delta_t
            for t in Tset
        )
    elseif objfun0 == "lineLossMin"
        objfun = sum(
            rdict_pu[(i, j)] * l[(i, j), t]
            for (i, j) in Lset, t in Tset
        )
    end

    # Append the alpha term only if objfun2 == "scd"
    @unpack alpha = data
    if objfun2 == "scd"
        objfun += sum(
            alpha * ((1 - η_C[j]) * P_c[j, t] + (1 / η_D[j] - 1) * P_d[j, t])
            for j in Bset, t in Tset
        )
    end

    # objfun represents our actual objective function which will be optimized for
    @objective(model, Min, objfun)

    # ===========================
    # Initializing Variables
    # ===========================

    @unpack Tset = data
    # Initialize substation power flow variable
    for t in Tset
        set_start_value(P_Subs[t], 0.0)
    end

    @unpack Tset, Lset = data
    # Initialize power flow variables
    for (i, j) in Lset, t in Tset
        set_start_value(P[(i, j), t], 0.0)
        set_start_value(Q[(i, j), t], 0.0)
        set_start_value(l[(i, j), t], 0.0)
    end

    @unpack Tset, Compset, V_Subs = data
    # Initialize voltage variables
    for j in Compset, t in Tset
        set_start_value(v[j, t], (V_Subs)^2)
    end

    @unpack Tset, Dset = data
    # Initialize PV inverter reactive dispatch variables
    for j in Dset, t in Tset
        set_start_value(q_D[j, t], 0.0)
    end

    @unpack Tset, Bset = data
    # Initialize battery real and reactive dispatch variables
    for j in Bset, t in Tset
        set_start_value(q_B[j, t], 0.0)
        set_start_value(P_c[j, t], 0.0)
        set_start_value(P_d[j, t], 0.0)
    end

    @unpack Tset, Bset, B0_pu = data
    # Initialize battery state of charge (SOC) variables
    for j in Bset, t in Tset
        set_start_value(B[j, t], B0_pu[j])
    end

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
