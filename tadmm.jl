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

# --- SAVE DATA DICT KEYS (CONTEXTUAL + LEX SORT) ---
function contextual_sort(keysarr)
    context_order = [
        :meta, :network, :loads, :pvs, :batteries, :voltages, :impedances, :costs, :shapes, :solver, :misc
    ]
    context_map = Dict(
        :meta => r"^(systemName|T(set)?|numAreas|objfun|temporal|algo|PSubsMax|inputForecast|solver|tSOC|relax|linearized|ged|alpha|gamma|warmStart|threshold|machine_ID|macroItrs|solution_time|simNature|spatial|root|src|processed|rawData)",
        :network => r"^(N(set|1set|m1set|c1set|nc1set)?|L(set|1set|m1set|Tset|notTset)?|parent|children|m|N1|Nm1|Nc1|Nnc1|m1|mm1)",
        :loads => r"^(NLset|N_L|p_L(_R|_pu)?|q_L(_R|_pu)?|Vminpu_L|Vmaxpu_L|LoadShapeLoad)",
        :pvs => r"^(Dset|n_D|DER_percent|p_D(_R|_pu)?|S_D(_R|_pu)?|irrad|Vminpu_D|Vmaxpu_D|LoadShapePV)",
        :batteries => r"^(Bset|n_B|Batt_percent|B0(_pu)?|Bref(_pu|_percent)?|B_R(_pu)?|P_B_R(_pu)?|S_B_R(_pu)?|eta_[CD]|soc_(min|max|0)|Vminpu_B|Vmaxpu_B)",
        :voltages => r"^(Vminpu|Vmaxpu|V_Subs_pu|substationBus|delta_t)",
        :impedances => r"^(rdict(_pu)?|xdict(_pu)?|Z_B(_dict)?|kVA_B(_dict)?|kV_B(_dict)?|MVA_B(_dict)?|Z_B|kVA_B|kV_B|MVA_B)",
        :costs => r"^(peakCost|offPeakCost|peakHoursFraction|LoadShapeCost)",
        :shapes => r"Shape",
        :solver => r"^(solver|warmStart_mu|threshold_conv_iters)",
        :misc => r"."
    )
    grouped = Dict(c => String[] for c in context_order)
    for k in keysarr
        for c in context_order
            if occursin(context_map[c], String(k))
                push!(grouped[c], String(k))
                break
            end
        end
    end
    out = String[]
    for c in context_order
        append!(out, sort(grouped[c]))
    end
    return out
end

sorted_keys = contextual_sort(collect(keys(data)))
open("tadmm_data_keys.txt", "w") do io
    for k in sorted_keys
        println(io, k)
    end
end

# --- SAVE ALL DATA VALUES (CONTEXTUAL + LEX SORT) ---
open("tadmm_data_values.txt", "w") do io
    for k in sorted_keys
        println(io, "[", k, "] => ", data[Symbol(k)])
    end
end
