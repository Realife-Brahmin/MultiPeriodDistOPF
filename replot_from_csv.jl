# Reconstruct solution dict from CSV and call Plotter.jl
using Plots, Printf, LaTeXStrings, Statistics

# Load Plotter functions
include("envs/tadmm/Plotter.jl")

# Manually parse CSV
function parse_csv(filename)
    lines = readlines(filename)
    header = split(lines[1], ",")
    data = Dict(col => [] for col in header)

    for line in lines[2:end]
        vals = split(line, ",")
        for (i, col) in enumerate(header)
            try
                push!(data[col], parse(Float64, vals[i]))
            catch
                push!(data[col], NaN)
            end
        end
    end
    return data
end

# Load convergence data
println("Loading convergence_data.csv...")
conv_data = parse_csv("envs/tadmm/processedData/large10kC_1ph_T6/convergence_data.csv")

# Reconstruct convergence_history dict
convergence_history = Dict(
    :obj_history => conv_data["objective"],
    :r_norm_history => conv_data["r_norm"],
    :s_norm_history => conv_data["s_norm"],
    :ρ_history => conv_data["rho"],
    :α_history => Float64[],
    :restart_history => Int[],
    :use_faadmm => false
)

# Reconstruct sol_socp_tadmm dict
sol_socp_tadmm = Dict(
    :convergence_history => convergence_history,
    :timing => Dict(
        :iteration_effective_times => Float64[],
        :total_effective_time => NaN,
        :subproblem_times_matrix => Matrix{Float64}(undef, 0, 0),
        :subproblem_iters_matrix => Matrix{Int}(undef, 0, 0)
    )
)

# Dummy BF solution
sol_socp_bf = Dict(
    :objective => NaN,
    :wallclock_time => NaN
)

# Call Plotter
eps_pri = 1e-4
eps_dual = 1e-2

println("Generating convergence plot...")
plot_tadmm_ldf_convergence(sol_socp_tadmm, sol_socp_bf, eps_pri, eps_dual;
                           showPlots=false, savePlots=true,
                           filename="envs/tadmm/processedData/large10kC_1ph_T6/convergence/tadmm_convergence_socp.png")
