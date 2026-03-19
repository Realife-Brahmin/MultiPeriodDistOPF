using Pkg

# Install required packages if needed
required_packages = ["Plots", "LaTeXStrings", "CSV", "DataFrames"]
for pkg in required_packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    end
end

using Plots
using LaTeXStrings
using CSV
using DataFrames
using Statistics
using Printf

# Define color scheme
primary = "#eb6f92"         # Rose/pink
secondary = "#31748f"       # Teal/blue
tertiary = "#9ccfd8"        # Light teal
accent = "#f6c177"          # Orange/gold
bg_color = "#FFFFFF"        # Pure white
text_color = "#000000"      # Black
grid_major = "#E8C4D4"      # Soft pink
grid_minor = "#D4D4E8"      # Soft lavender

# Command line argument or default
T = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 6

# Load data
csv_path = raw"c:\Users\arjha\OneDrive - Tesla\Documents\documents_general_addendum\MultiPeriodDistOPF\envs\tadmm\processedData\ieee2552_1ph_T" * string(T) * raw"\subproblem_timing_details.csv"

if !isfile(csv_path)
    error("CSV file not found: $csv_path")
end

println("Analyzing T=$T data from: $csv_path")

df = CSV.read(csv_path, DataFrame)

# For each iteration, find the worst (slowest) subproblem
worst_per_iter = combine(groupby(df, :iteration)) do group
    idx_max = argmax(group.solve_time_sec)
    group[idx_max, :]
end

# Count how often each time period is the bottleneck
bottleneck_counts = combine(groupby(worst_per_iter, :subproblem), nrow => :count)
sort!(bottleneck_counts, :subproblem)

# Get the distribution of worst times
worst_times = worst_per_iter.solve_time_sec

# Expected frequency (uniform distribution)
total_iters = nrow(worst_per_iter)
expected_freq = total_iters / T

println("\nCreating bottleneck analysis plots for T=$T...")

# =============================================================================
# PLOT 1: Bar chart of bottleneck frequency
# =============================================================================

p1 = bar(
    bottleneck_counts.subproblem,
    bottleneck_counts.count,
    label="",  # No legend needed for single series
    color=primary,
    xlabel=L"Time Period $t_0$ (hour)",
    ylabel=L"Count (times this period was slowest)",
    title="Which Time Period is the Bottleneck?\n" * L"ieee2552\_1ph, $T=%$T$, %$total_iters iterations",
    titlefontsize=14,
    titlefontcolor=text_color,
    labelfontsize=13,
    guidefontcolor=text_color,
    tickfontcolor=text_color,
    tickfontsize=11,
    background_color=bg_color,
    foreground_color=text_color,
    grid=true,
    minorgrid=true,
    gridcolor=grid_major,
    minorgridcolor=grid_minor,
    gridlinewidth=1.5,
    minorgridlinewidth=0.8,
    gridalpha=0.5,
    minorgridalpha=0.3,
    gridstyle=:solid,
    framestyle=:box,
    size=(800, 600),
    dpi=300,
    margin=5Plots.mm,
    xticks=1:T,
    alpha=0.85
)

# Add horizontal line for expected frequency
hline!(
    [expected_freq],
    label=@sprintf("Expected (uniform): %.2f", expected_freq),
    linewidth=3,
    linestyle=:dash,
    color=secondary,
    legend=:topright,
    legendfontcolor=text_color,
    legendfontsize=11
)

# Add horizontal line for 2x threshold
hline!(
    [2*expected_freq],
    label=L"Outlier threshold ($2 \times$ expected)",
    linewidth=3,
    linestyle=:dot,
    color=accent,
    legendfontcolor=text_color,
    legendfontsize=11
)

# Add count labels on top of bars
for row in eachrow(bottleneck_counts)
    annotate!(
        row.subproblem, row.count + 0.5,
        text(string(row.count), text_color, 10, :center)
    )
end

# Highlight outliers
outlier_mask = bottleneck_counts.count .> 2*expected_freq
if any(outlier_mask)
    outlier_data = bottleneck_counts[outlier_mask, :]
    bar!(
        p1,
        outlier_data.subproblem,
        outlier_data.count,
        label=L"Persistent bottlenecks ($>2 \times$ expected)",
        color=accent,
        alpha=0.85,
        legendfontcolor=text_color,
        legendfontsize=11
    )
end

output_path1 = joinpath(dirname(csv_path), "bottleneck_frequency.png")
savefig(p1, output_path1)
println("✓ Saved: $output_path1")

# =============================================================================
# PLOT 2: Distribution of worst subproblem times
# =============================================================================

p2 = histogram(
    worst_times,
    bins=20,
    label="",  # No legend for histogram itself
    color=primary,
    xlabel=L"Worst Subproblem Time (seconds)",
    ylabel=L"Frequency (iterations)",
    title="Distribution of Iteration Makespans\n" * L"(worst subproblem time per iteration)",
    titlefontsize=14,
    titlefontcolor=text_color,
    labelfontsize=13,
    guidefontcolor=text_color,
    tickfontcolor=text_color,
    tickfontsize=11,
    background_color=bg_color,
    foreground_color=text_color,
    grid=true,
    minorgrid=true,
    gridcolor=grid_major,
    minorgridcolor=grid_minor,
    gridlinewidth=1.5,
    minorgridlinewidth=0.8,
    gridalpha=0.5,
    minorgridalpha=0.3,
    gridstyle=:solid,
    framestyle=:box,
    size=(800, 600),
    dpi=300,
    margin=5Plots.mm,
    alpha=0.85
)

# Add statistics as vertical lines
mean_time = mean(worst_times)
median_time = median(worst_times)

vline!(
    [mean_time],
    label=@sprintf("Mean: %.2f s", mean_time),
    linewidth=4,
    linestyle=:solid,
    color=secondary,
    legend=:topright,
    legendfontcolor=text_color,
    legendfontsize=11
)

vline!(
    [median_time],
    label=@sprintf("Median: %.2f s", median_time),
    linewidth=4,
    linestyle=:dash,
    color=tertiary,
    legendfontcolor=text_color,
    legendfontsize=11
)

# Add text annotation with summary statistics
stats_text = @sprintf("Min: %.2f s\nMax: %.2f s\nStd: %.2f s",
                      minimum(worst_times), maximum(worst_times), std(worst_times))
annotate!(
    maximum(worst_times) * 0.75, maximum(ylims()) * 0.85,
    text(stats_text, text_color, 11, :left)
)

output_path2 = joinpath(dirname(csv_path), "makespan_distribution.png")
savefig(p2, output_path2)
println("✓ Saved: $output_path2")

# =============================================================================
# PLOT 3: Time series showing which period was slowest per iteration
# =============================================================================

p3 = scatter(
    worst_per_iter.iteration,
    worst_per_iter.subproblem,
    label="",
    color=primary,
    markersize=6,
    markerstrokewidth=0,
    xlabel=L"ADMM Iteration $k$",
    ylabel=L"Time Period $t_0$ (bottleneck)",
    title="Bottleneck Evolution Across Iterations\n" * L"Which period was slowest at each iteration?",
    titlefontsize=14,
    titlefontcolor=text_color,
    labelfontsize=13,
    guidefontcolor=text_color,
    tickfontcolor=text_color,
    tickfontsize=11,
    background_color=bg_color,
    foreground_color=text_color,
    grid=true,
    minorgrid=true,
    gridcolor=grid_major,
    minorgridcolor=grid_minor,
    gridlinewidth=1.5,
    minorgridlinewidth=0.8,
    gridalpha=0.5,
    minorgridalpha=0.3,
    gridstyle=:solid,
    framestyle=:box,
    size=(800, 600),
    dpi=300,
    margin=5Plots.mm,
    yticks=1:T,
    alpha=0.7
)

# Highlight the persistent outliers (if any)
if any(outlier_mask)
    outlier_periods = bottleneck_counts[outlier_mask, :subproblem]
    for t0 in outlier_periods
        outlier_iters = worst_per_iter[worst_per_iter.subproblem .== t0, :]
        scatter!(
            p3,
            outlier_iters.iteration,
            outlier_iters.subproblem,
            label=@sprintf("t₀=%d (persistent bottleneck)", t0),
            color=accent,
            markersize=8,
            markerstrokewidth=1,
            markerstrokecolor=text_color,
            alpha=0.9,
            legendfontcolor=text_color,
            legendfontsize=11,
            legend=:topright
        )
    end
end

output_path3 = joinpath(dirname(csv_path), "bottleneck_timeseries.png")
savefig(p3, output_path3)
println("✓ Saved: $output_path3")

# =============================================================================
# Summary statistics
# =============================================================================

println("\n" * "="^70)
println("BOTTLENECK ANALYSIS SUMMARY: T=$T")
println("="^70)
println("\nWhich time periods are slowest most often?")
println("-"^70)
for row in eachrow(bottleneck_counts)
    pct = 100 * row.count / total_iters
    is_outlier = row.count > 2*expected_freq
    marker = is_outlier ? " ⚠" : ""
    println(@sprintf("  t₀=%2d: %2d times (%5.1f%%)%s", row.subproblem, row.count, pct, marker))
end
println("-"^70)

if any(outlier_mask)
    println("\nPersistent bottlenecks detected:")
    for row in eachrow(bottleneck_counts[outlier_mask, :])
        ratio = row.count / expected_freq
        println(@sprintf("  t₀=%d: %.2f× expected frequency", row.subproblem, ratio))
    end
else
    println("\n✓ No persistent bottlenecks (all within 2× expected frequency)")
end

println("\n" * "="^70)
println("All plots generated successfully!")
println("="^70)
println("Location: $(dirname(csv_path))")
println("  1. bottleneck_frequency.png")
println("  2. makespan_distribution.png")
println("  3. bottleneck_timeseries.png")
println("="^70)
