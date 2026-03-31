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

# Use GR backend for 3D plotting
gr()

# Define color scheme
primary = "#eb6f92"         # Rose/pink
secondary = "#31748f"       # Teal/blue
tertiary = "#9ccfd8"        # Light teal
accent = "#f6c177"          # Orange/gold
bg_color = "#FFFFFF"        # Pure white
text_color = "#000000"      # Black
grid_major = "#E8C4D4"      # Soft pink
grid_minor = "#D4D4E8"      # Soft lavender

# Command line arguments or defaults
# Usage: julia plot_3d_bar_iterations.jl [T] [systemName] [num_samples]
T = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 24
systemName = length(ARGS) > 1 ? ARGS[2] : "ieee123A_1ph"
num_samples = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 4  # Excluding iteration 1

# Load data
base_path = raw"c:\Users\arjha\OneDrive - Tesla\Documents\documents_general_addendum\MultiPeriodDistOPF\envs\tadmm\processedData"
csv_path = joinpath(base_path, "$(systemName)_T$(T)", "subproblem_timing_details.csv")

if !isfile(csv_path)
    error("CSV file not found: $csv_path")
end

println("Analyzing T=$T data from: $csv_path")
df = CSV.read(csv_path, DataFrame)

# Get total iterations
total_iters = maximum(df.iteration)
println("Total ADMM iterations: $total_iters")

# Select sample iterations (excluding iteration 1, which we'll plot separately)
if num_samples >= (total_iters - 1)
    sample_iters = collect(2:total_iters)
else
    # Strategic sampling from iterations 2 onwards
    indices = round.(Int, range(2, total_iters, length=num_samples))
    sample_iters = unique(indices)
end

println("Iterations for warm-start plot: $sample_iters")

# Color palette for different iterations
colors = [secondary, tertiary, accent, "#f6c177", "#e0def4"]
while length(colors) < length(sample_iters)
    push!(colors, "#$(rand(UInt8))$(rand(UInt8))$(rand(UInt8))")
end

# =============================================================================
# 3D Bar Plot (Cityscape View) - Warm Iterations Only
# =============================================================================

println("\nCreating 3D cityscape visualization...")

# Helper function to create a 3D rectangular bar (box mesh)
function create_3d_box(x_center, y_center, z_height, width_x, width_y, color_val)
    # Box dimensions
    x_min = x_center - width_x/2
    x_max = x_center + width_x/2
    y_min = y_center - width_y/2
    y_max = y_center + width_y/2
    z_min = 0.0
    z_max = z_height

    # 8 vertices of the box
    vertices = [
        # Bottom face (z=0)
        [x_min, y_min, z_min],
        [x_max, y_min, z_min],
        [x_max, y_max, z_min],
        [x_min, y_max, z_min],
        # Top face (z=height)
        [x_min, y_min, z_max],
        [x_max, y_min, z_max],
        [x_max, y_max, z_max],
        [x_min, y_max, z_max]
    ]

    # Extract x, y, z coordinates
    xs = [v[1] for v in vertices]
    ys = [v[2] for v in vertices]
    zs = [v[3] for v in vertices]

    # Define 6 faces (each face is 2 triangles = 6 vertices in total, but mesh3d uses indices)
    # We'll use i, j, k for triangle indices (0-indexed in Plots)
    # Each face needs to be split into 2 triangles
    faces_i = Int[]
    faces_j = Int[]
    faces_k = Int[]

    # Bottom face (vertices 1,2,3,4) - split into triangles (1,2,3) and (1,3,4)
    append!(faces_i, [0, 0])
    append!(faces_j, [1, 2])
    append!(faces_k, [2, 3])

    # Top face (vertices 5,6,7,8) - split into triangles (5,6,7) and (5,7,8)
    append!(faces_i, [4, 4])
    append!(faces_j, [5, 6])
    append!(faces_k, [6, 7])

    # Front face (vertices 1,2,6,5)
    append!(faces_i, [0, 0])
    append!(faces_j, [1, 5])
    append!(faces_k, [5, 4])

    # Back face (vertices 4,3,7,8)
    append!(faces_i, [3, 3])
    append!(faces_j, [2, 6])
    append!(faces_k, [6, 7])

    # Left face (vertices 1,4,8,5)
    append!(faces_i, [0, 0])
    append!(faces_j, [3, 7])
    append!(faces_k, [7, 4])

    # Right face (vertices 2,3,7,6)
    append!(faces_i, [1, 1])
    append!(faces_j, [2, 6])
    append!(faces_k, [6, 5])

    return (xs, ys, zs, faces_i, faces_j, faces_k, color_val)
end

# Create the base plot
p2 = plot(
    camera=(45, 30),
    xlabel=L"Time Period $t_0$ (hour)",
    ylabel=L"Iteration Sample",
    zlabel=L"Ipopt Iterations",
    title="Warm-Start Iterations: Ipopt Iteration Count (3D Bar View)\n" *
          L"%$(systemName), $T=%$T$ periods, iterations %$(minimum(sample_iters))-%$(maximum(sample_iters))",
    titlefontsize=14,
    titlefontcolor=text_color,
    labelfontsize=13,
    guidefontcolor=text_color,
    tickfontcolor=text_color,
    tickfontsize=11,
    background_color=bg_color,
    foreground_color=text_color,
    grid=true,
    gridcolor=grid_major,
    gridlinewidth=1.5,
    gridalpha=0.5,
    framestyle=:box,
    size=(1200, 800),
    dpi=300,
    margin=8Plots.mm,
    yticks=(1:length(sample_iters), ["k=$(iter)" for iter in sample_iters]),
    legend=:outertopright,
    legendfontcolor=text_color,
    legendfontsize=13,
    legend_background_color=:white,
    legend_foreground_color=text_color,
    legendtitle="ADMM Iteration",
    legendtitlefontsize=13,
    xlims=(0, T+1),
    ylims=(0, length(sample_iters)+1)
)

# Add 3D bars for each iteration
bar_width_x = 0.6
bar_width_y = 0.6

for (idx, iter) in enumerate(sample_iters)
    df_iter = filter(row -> row.iteration == iter, df)
    sort!(df_iter, :subproblem)

    # Add all bars for this iteration
    for row in eachrow(df_iter)
        x_center = Float64(row.subproblem)
        y_center = Float64(idx)
        z_height = Float64(row.ipopt_iters)

        # Create the 3D box
        xs, ys, zs, i_faces, j_faces, k_faces, _ = create_3d_box(
            x_center, y_center, z_height, bar_width_x, bar_width_y, colors[idx]
        )

        # Add mesh to plot (only add label for first bar of each iteration)
        if row.subproblem == 1
            mesh3d!(
                p2,
                xs, ys, zs,
                connections=(i_faces, j_faces, k_faces),
                color=colors[idx],
                alpha=0.85,
                linecolor=:black,
                linewidth=0.5,
                label=@sprintf("k=%d", iter)
            )
        else
            mesh3d!(
                p2,
                xs, ys, zs,
                connections=(i_faces, j_faces, k_faces),
                color=colors[idx],
                alpha=0.85,
                linecolor=:black,
                linewidth=0.5,
                label=""
            )
        end
    end
end

output_path = joinpath(dirname(csv_path), "3d_bars_warm_iterations.png")
savefig(p2, output_path)
println("✓ Saved: $output_path")

# =============================================================================
# Summary statistics
# =============================================================================

println("\n" * "="^70)
println("3D CITYSCAPE VISUALIZATION - WARM ITERATIONS")
println("="^70)
println("System: $systemName, T=$T periods")
println("Total ADMM iterations: $total_iters")
println("Sampled iterations: $sample_iters")
println()

for (idx, iter) in enumerate(sample_iters)
    df_iter = filter(row -> row.iteration == iter, df)
    mean_val = mean(df_iter.ipopt_iters)
    median_val = round(Int, median(df_iter.ipopt_iters))
    min_val = minimum(df_iter.ipopt_iters)
    max_val = maximum(df_iter.ipopt_iters)
    hardest_row = df_iter[argmax(df_iter.ipopt_iters), :]

    println(@sprintf("Iteration k=%2d: mean=%.1f, median=%d, range=%d-%d, hardest=t₀=%d",
                     iter, mean_val, median_val, min_val, max_val, hardest_row.subproblem))
end

println()
println("="^70)
println("Plot saved: $(dirname(csv_path))")
println("  → 3d_bars_warm_iterations.png")
println("="^70)
