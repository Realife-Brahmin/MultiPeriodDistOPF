#! Activate project environment
import Pkg; Pkg.activate(".")
using Revise
using OpenDSSDirect
using Dates
using Printf
using Random

# Import the ODD parser
Revise.includet("src/Parser/parseFromDSS.jl")
using .parseFromDSS

# --- USER DEFINED META VARIABLES ---
# System and simulation parameters
systemName = "ads10A_1ph"   # Change as needed
T = 24                      # Number of time steps

# Example: user can define LoadShape, PVShape, etc. here if needed
# For now, just placeholders
LoadShape = ones(T)         # Placeholder, replace with actual shape if needed
PVShape = ones(T)           # Placeholder

# --- PARSE SYSTEM ---
data = get_system_config_from_dss(systemName, T)

# --- PATCH OpenDSSDirect.Text.Command to suppress warnings unless 'Unknown Property' ---
let orig_cmd = OpenDSSDirect.Text.Command
    function OpenDSSDirect.Text.Command(cmd::AbstractString)
        io = IOBuffer()
        redirect_stderr(io) do
            res = orig_cmd(cmd)
            msg = String(take!(io))
            if occursin("Unknown Property", res)
                @warn "Result of running OpenDSS Command $cmd is: $res"
            end
            return res
        end
    end
end


# --- SAVE DATA DICT KEYS (GROUPED BY CONTEXT WITH HEADERS) ---
const context_map = Dict(
    "Substation" => ["PSubsMax_kW", "V_Subs_pu", "substationBus"],
    "Simulation_Problem" => ["T", "Tset", "gedAppendix", "gedString", "inputForecastDescription", "objfun0", "objfun2", "objfunSense", "objfunString", "objfunPrefix", "objfunUnit", "objfunAppendix", "objfunConciseDescription", "tSOC_hard", "relax_terminal_soc_constraint", "LoadShapeLoad", "LoadShapePV", "LoadShapeCost", "offPeakCost", "peakCost", "peakHoursFraction", "delta_t", "gedDict_ud", "systemName"],
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
        for k in sort(unmatched)
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
        for k in sort(unmatched)
            val = get(data, Symbol(k), "<missing>")
            println(io, "[", k, "] => ", val)
        end
    end
end
