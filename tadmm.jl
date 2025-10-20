#! Activate project environment

import Pkg; Pkg.activate(".")
using Revise
using JuMP
import OpenDSSDirect as dss
using OpenDSSDirect
using Dates
using Printf
using Random
using Parameters: @unpack
using Gurobi
using Ipopt

# Import the ODD parser
Revise.includet("src/Parser/parseFromDSS.jl")
import .parseFromDSS as Parser


# --- USER DEFINED META VARIABLES ---
# System and simulation parameters (match copper plate example)
systemName = "ads10A_1ph"   # Change as needed
T = 24                      # Number of time steps
delta_t_h = 1.0               # Time step duration in hours

# --- COPPER PLATE VALUES (from admm_temporal_copper_plate.jl) ---
LoadShape = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2  # Normalized load shape [0, 1]
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2   # $/kWh
C_B = 1e-6 * minimum(LoadShapeCost)  # Quadratic cost coefficient for battery power: C_B * P_B^2

# If you want to keep PVShape, you can set it to zeros or as needed
# PVShape = zeros(T)

# --- Helper: update data dict with user-supplied kwargs ---
function update_data_with_kwargs!(data::Dict; kwargs...)
    for (k, v) in kwargs
        data[k] = v
    end
    return data
end

# --- PARSE SYSTEM AND APPLY USER OVERRIDES ---
data = Parser.get_system_config_from_dss(systemName, T)

user_overrides = Dict(
    :LoadShapeLoad => LoadShape,
    :LoadShapeCost => LoadShapeCost,
    :C_B => C_B,
    :delta_t_h => delta_t_h,
    # :LoadShapePV => PVShape,
)
update_data_with_kwargs!(data; user_overrides...)



# --- Print basic OpenDSS powerflow stats for quick validation ---
include("src/openDSSValidator.jl")
using .openDSSValidator
print_basic_powerflow_stats(data)


# --- SAVE DATA DICT KEYS (GROUPED BY CONTEXT WITH HEADERS) ---
const context_map = Dict(
    "Substation" => ["PSubsMax_kW", "V_Subs_pu", "substationBus"],
    "Simulation_Problem" => ["T", "Tset", "gedAppendix", "gedString", "inputForecastDescription", "objfun0", "objfun2", "objfunSense", "objfunString", "objfunPrefix", "objfunUnit", "objfunAppendix", "objfunConciseDescription", "tSOC_hard", "relax_terminal_soc_constraint", "LoadShapeLoad", "LoadShapePV", "LoadShapeCost", "C_B", "offPeakCost", "peakCost", "peakHoursFraction", "delta_t", "gedDict_ud", "systemName"],
    "Simulation_Algorithm_Params" => ["solver", "linearizedModel", "linearizedModelAppendix", "linearizedModelString", "algo_temporal_decmp", "alpha_fpi", "gamma_fpi", "numAreas", "threshold_conv_iters", "warmStart_mu", "temporal_decmp", "temporalDecmpAppendix", "temporalDecmpString", "simNatureAppendix", "simNatureString", "spatialDecAppendix", "spatialDecString"],
    "Simulation_Machine" => ["machine_ID"],
    "Simulation_Run" => ["macroItrsCompleted", "solution_time"],
    "Directories" => ["processedDataFolderPath", "rawDataFolderPath", "rootFolderPath", "srcFolderPath"],
    "System_Topology" => ["kVA_B", "kVA_B_dict", "kV_B", "kV_B_dict", "MVA_B", "MVA_B_dict", "Z_B", "Z_B_dict", "rdict", "rdict_pu", "xdict", "xdict_pu", "L1set", "LTset", "Lm1set", "LnotTset", "Lset", "N", "N1", "N1set", "NLset", "N_L", "Nc1", "Nc1set", "Nm1", "Nm1set", "Nnc1", "Nnc1set", "Nset", "children", "m", "m1", "mm1", "parent"],
    "System_Storage_Data" => ["B0", "B0_pu", "B_R", "B_R_pu", "Batt_percent", "Bref", "Bref_percent", "Bref_pu", "Bset", "P_B_R", "P_B_R_pu", "S_B_R", "S_B_R_pu", "Vmaxpu_B", "Vminpu_B", "eta_C", "eta_D", "n_B", "soc_0", "soc_max", "soc_min"],
    "System_DER_Data" => ["Dset", "n_D", "DER_percent", "p_D", "p_D_R", "p_D_R_pu", "p_D_pu", "S_D_R", "S_D_R_pu", "Vmaxpu_D", "Vminpu_D", "irrad"],
    "System_Load_Data" => ["NLset", "N_L", "p_L", "p_L_R", "p_L_R_pu", "p_L_pu", "q_L", "q_L_R", "q_L_R_pu", "q_L_pu", "Vmaxpu_L", "Vminpu_L"],
    "Voltage_Limits" => ["Vmaxpu", "Vminpu"],
    "Other_Misc" => []
)


all_keys_str = Set(string(k) for k in keys(data))
open("tadmm_data_keys.txt", "w") do io
    for (section, keyslist) in context_map
        println(io, "==== $section ====")
        for k in sort(keyslist)
            if k in all_keys_str
                println(io, k)
            end
        end
        println(io)
    end
    # Print any keys not matched
    matched_keys = reduce(vcat, values(context_map))
    unmatched = setdiff(all_keys_str, Set(matched_keys))
    if !isempty(unmatched)
        println(io, "==== Unmatched ====")
        for k in sort(collect(unmatched))
            println(io, k)
        end
    end
end

# --- SAVE ALL DATA VALUES (GROUPED BY CONTEXT WITH HEADERS) ---
open("tadmm_data_values.txt", "w") do io
    for (section, keyslist) in context_map
        println(io, "==== $section ====")
        for k in sort(keyslist)
            if k in all_keys_str
                val = get(data, Symbol(k), "<missing>")
                println(io, "[", k, "] => ", val)
            end
        end
        println(io)
    end
    # Print any keys not matched
    matched_keys = reduce(vcat, values(context_map))
    unmatched = setdiff(all_keys_str, Set(matched_keys))
    if !isempty(unmatched)
        println(io, "==== Unmatched ====")
        for k in sort(collect(unmatched))
            val = get(data, Symbol(k), "<missing>")
            println(io, "[", k, "] => ", val)
        end
    end
end

# --- MPOPF LinDistFlow Brute-Force Builder ---

# --- MPOPF LinDistFlow Brute-Force Builder (ModelBuilder-style variable convention) ---

# --- MPOPF LinDistFlow Brute-Force Builder (P_B only, copper plate style) ---
function build_MPOPF_with_LinDistFlow_BruteForced(data; solver=:ipopt)

    @unpack Nset, Lset, Dset, Bset, Tset, NLset = data
    @unpack p_L, p_D, P_B_R, B_R, kVA_B, kV_B, rdict_pu, xdict_pu, Vminpu, Vmaxpu, LoadShapeCost, C_B, delta_t_h, soc_min, soc_max = data
    Δt = delta_t_h  # Time step duration in hours
    model = Model()
    if solver == :gurobi
        set_optimizer(model, Gurobi.Optimizer)
    elseif solver == :ipopt
        set_optimizer(model, Ipopt.Optimizer)
    end
    # Variable definitions (P_B only)
    @variable(model, P_Subs[t in Tset] >= 0)
    @variable(model, P[(i, j) in Lset, t in Tset])
    @variable(model, Q[(i, j) in Lset, t in Tset])
    @variable(model, v[n in Nset, t in Tset])
    @variable(model, P_B[b in Bset, t in Tset])
    @variable(model, B[b in Bset, t in Tset])
    @variable(model, P_Subs[t in Tset])

    # Objective: substation cost + battery cost
    @expression(model, sub_cost, sum(LoadShapeCost[t] * P_Subs[t] for t in Tset))
    @expression(model, batt_cost, sum(C_B * P_B[b, t]^2 for b in Bset, t in Tset))
    @objective(model, Min, sub_cost + batt_cost)

    # Constraints
    for t in Tset
        # Substation node real power balance (bus 1)
        @constraint(model, P_Subs[t] + sum(P[(1, j), t] for (i, j) in Lset if i == 1) == 0)

        # Non-substation node real power balance
        for n in NLset
            incoming = sum(P[(i, n), t] for (i, j) in Lset if j == n)
            outgoing = sum(P[(n, j), t] for (i, j) in Lset if i == n)
            p_load = get(p_L, (n, t), 0.0)
            p_pv = get(p_D, (n, t), 0.0)
            p_batt = sum(P_B[b, t] for b in Bset if b == n)
            @constraint(model, incoming - outgoing + p_load - p_pv + p_batt == 0)
        end

        # Substation node reactive power balance (bus 1)
        @constraint(model, sum(Q[(1, j), t] for (i, j) in Lset if i == 1) == 0)

        # Non-substation node reactive power balance
        for n in NLset
            incoming = sum(Q[(i, n), t] for (i, j) in Lset if j == n)
            outgoing = sum(Q[(n, j), t] for (i, j) in Lset if i == n)
            # For now, assume no q_D or q_B (can be extended)
            @constraint(model, incoming - outgoing == 0)
        end

        # KVL constraints (LinDistFlow, using pu values)
        for (i, j) in Lset
            r = get(rdict_pu, (i, j), 0.0)
            x = get(xdict_pu, (i, j), 0.0)
            @constraint(model, v[j, t] == v[i, t] - 2 * (r * P[(i, j), t] + x * Q[(i, j), t]))
        end

        # Substation voltage (squared, 1.07^2)
        @constraint(model, v[1, t] == 1.07^2)

        # Voltage limits (squared)
        for n in Nset
            @constraint(model, Vminpu[n]^2 <= v[n, t] <= Vmaxpu[n]^2)
        end

        # Battery SOC update and limits (no efficiency)
        for b in Bset
            if t == 1
                @constraint(model, B[b, t] == 0.5 * B_R[b])
            else
                @constraint(model, B[b, t] == B[b, t-1] - P_B[b, t-1] * Δt)
            end
            @constraint(model, soc_min * B_R[b] <= B[b, t] <= soc_max * B_R[b])
            @constraint(model, -P_B_R[b] <= P_B[b, t] <= P_B_R[b])
        end
    end

    optimize!(model)
    status = termination_status(model)
    obj = objective_value(model)
    return Dict(:model => model, :status => status, :objective => obj, :P_Subs => P_Subs, :B => B, :P_B => P_B)
end

# --- CALL LIN DISTFLOW MPOPF BUILDER ---
println("\nSolving MPOPF with LinDistFlow Brute-Forced approach...")
modelDict = build_MPOPF_with_LinDistFlow_BruteForced(data; solver=:ipopt)
println("LinDistFlow solution status: ", modelDict[:status])
println("LinDistFlow objective value: ", modelDict[:objective])
