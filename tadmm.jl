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

# Import the ODD parser
Revise.includet("src/Parser/parseFromDSS.jl")
using .parseFromDSS


# --- USER DEFINED META VARIABLES ---
# System and simulation parameters (match copper plate example)
systemName = "ads10A_1ph"   # Change as needed
T = 24                      # Number of time steps


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
data = get_system_config_from_dss(systemName, T)

user_overrides = Dict(
    :LoadShapeLoad => LoadShape,
    :LoadShapeCost => LoadShapeCost,
    :C_B => C_B,
    # :LoadShapePV => PVShape,
)
update_data_with_kwargs!(data; user_overrides...)



# --- Print basic OpenDSS powerflow stats for quick validation ---
include("src/openDSSValidator.jl")
using .openDSSValidator
print_basic_powerflow_stats(data)


# --- BRUTE-FORCE LDF MPOPF SOLUTION (at end) ---
# include("solve_MPOPF_with_LDF_BruteForced.jl")
#
# println("\nSolving MPOPF with LDF Brute-Forced approach...")
# brute_result = solve_MPOPF_with_LDF_BruteForced(data; solver=:gurobi)
# println("Brute-force solution status: ", brute_result[:status])
# println("Brute-force objective value: ", brute_result[:objective])
#
# Save battery actions (SOC and P_B) for plotting
# using Serialization
# Serialization.serialize("brute_battery_SOC.jls", brute_result[:soc])
# Serialization.serialize("brute_battery_P_B.jls", brute_result[:P_B])

# (Later: add plotting and tadmm implementation here)


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
