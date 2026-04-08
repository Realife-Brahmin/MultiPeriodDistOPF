#!/usr/bin/env julia
# Standalone timing plotter: reads saved CSV and generates exact Plotter.jl layout

using Plots, Printf, Statistics

# Manual CSV parser
function parse_csv(filename)
    lines = readlines(filename)
    header = split(lines[1], ",")
    data = Dict(col => [] for col in header)

    for line in lines[2:end]
        if isempty(strip(line))
            continue
        end
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

# Configuration
timing_csv = "envs/tadmm/processedData/large10kC_1ph_T6/subproblem_timing_details.csv"
convergence_csv = "envs/tadmm/processedData/large10kC_1ph_T6/convergence_data.csv"
output_path = "envs/tadmm/processedData/large10kC_1ph_T6/convergence/tadmm_timing_socp.png"

println("Loading timing data from: $timing_csv")
timing_data = parse_csv(timing_csv)

println("Loading convergence data from: $convergence_csv")
conv_data = parse_csv(convergence_csv)

# Extract convergence data
iterations = Int.(conv_data["iteration"])
cum_eff_time = conv_data["cum_eff_time"]
eff_time_this_k = conv_data["eff_time_this_k"]

# Get max iteration
n_iter = length(iterations)

# Extract timing data by iteration
iteration_times = timing_data["iteration"]
subproblem_times = timing_data["solve_time_sec"]

# Group by iteration to compute per-iteration statistics
max_sp_times = Float64[]
mean_sp_times = Float64[]
median_sp_times = Float64[]

for k in 1:n_iter
    mask = iteration_times .== k
    times_k = subproblem_times[mask]
    if !isempty(times_k)
        push!(max_sp_times, maximum(times_k))
        push!(mean_sp_times, mean(times_k))
        push!(median_sp_times, median(times_k))
    else
        push!(max_sp_times, NaN)
        push!(mean_sp_times, NaN)
        push!(median_sp_times, NaN)
    end
end

# Setup
gr()
theme(:mute)

# X-ticks (same logic as convergence plot)
markerstrokewidth = if n_iter <= 50
    1.5
elseif n_iter <= 100
    1.0
elseif n_iter <= 200
    0.5
else
    0.3
end

xtick_vals = if n_iter <= 10
    collect(iterations)
elseif n_iter <= 50
    ticks = collect(10:10:n_iter)
    if isempty(ticks) || (n_iter - last(ticks)) > 3
        vcat(1, ticks, n_iter)
    else
        vcat(1, ticks[1:end-1], n_iter)
    end
else
    step = max(1, div(n_iter, 5))
    ticks = collect(step:step:n_iter)
    if isempty(ticks) || (n_iter - last(ticks)) > step * 0.3
        vcat(1, ticks, n_iter)
    else
        vcat(1, ticks[1:end-1], n_iter)
    end
end
xtick_vals = sort(unique(xtick_vals))

# Subplot 1: Maximum subproblem time per iteration (bottleneck for parallelization)
p1 = plot(
    iterations, max_sp_times,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Time [seconds]",
    title="Max Subproblem Time (Parallelization Bottleneck)",
    lw=3,
    color=:crimson,
    markershape=:circle,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    label="Max SP time",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern",
    top_margin=2Plots.mm
)

# Subplot 2: Median subproblem time per iteration (typical case)
p2 = plot(
    iterations, median_sp_times,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Time [seconds]",
    title="Median Subproblem Time (Typical Case)",
    lw=3,
    color=:dodgerblue,
    markershape=:square,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    label="Median SP time",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern"
)

# Subplot 3: Total iteration time (sum of all subproblems, showing sequential cost)
total_times = mean_sp_times .* 12  # Rough estimate: mean * number of subproblems
p3 = plot(
    iterations, total_times,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Time [seconds]",
    title="Total Subproblem Time (Sequential Cost)",
    lw=3,
    color=:darkgreen,
    markershape=:diamond,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    label="Total SP time",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern"
)

# Subplot 4: Cumulative wall-clock time
p4 = plot(
    iterations, cum_eff_time,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Cumulative Time [seconds]",
    title="Cumulative Wall-Clock Time (Effective)",
    lw=3,
    color=:purple,
    markershape=:hexagon,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    label="Cumulative eff time",
    legend=:topleft,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern",
    bottom_margin=2Plots.mm
)

# Combine (4 subplots, 1 column)
p_combined = plot(p1, p2, p3, p4,
                 layout=(4, 1),
                 size=(900, 1300),
                 plot_title="tADMM Timing Analysis",
                 plot_titlefontsize=14,
                 plot_titlefontfamily="Computer Modern",
                 left_margin=8Plots.mm,
                 right_margin=5Plots.mm,
                 top_margin=8Plots.mm,
                 bottom_margin=5Plots.mm)

# Print summary
total_eff_time = last(cum_eff_time)
total_seq_time = sum(total_times)
avg_eff_time_per_iter = mean(eff_time_this_k)
max_eff_time = maximum(eff_time_this_k)

println("\ntADMM Timing Analysis:")
println("  Total iterations:           $n_iter")
@printf "  Total effective time:       %.2f seconds\n" total_eff_time
@printf "  Total sequential time:      %.2f seconds\n" total_seq_time
@printf "  Avg eff time per iteration: %.2f seconds\n" avg_eff_time_per_iter
@printf "  Max eff time (bottleneck):  %.2f seconds\n" max_eff_time
@printf "  Parallel efficiency:        %.1f%%\n" (total_seq_time / total_eff_time * 100)

# Save
@printf "Saving tADMM timing plot to: %s\n" output_path
savefig(p_combined, output_path)
println("✓ Done")
