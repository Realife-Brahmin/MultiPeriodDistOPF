#!/usr/bin/env julia
# Standalone convergence plotter: reads saved CSV and generates exact Plotter.jl layout

using Plots, Printf, Statistics

# Manual CSV parser (no dependencies)
function parse_csv(filename)
    lines = readlines(filename)
    header = split(lines[1], ",")
    data = Dict(col => Float64[] for col in header)

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
csv_path = "envs/tadmm/processedData/large10kC_1ph_T6/convergence_data.csv"
output_path = "envs/tadmm/processedData/large10kC_1ph_T6/convergence/tadmm_convergence_socp.png"
eps_pri = 1e-4
eps_dual = 1e-2

println("Loading convergence data from: $csv_path")
conv_data = parse_csv(csv_path)

# Extract columns
obj_history = conv_data["objective"]
r_norm_history = conv_data["r_norm"]
s_norm_history = conv_data["s_norm"]
ρ_history = conv_data["rho"]

# Filter NaN values (exact replication of Plotter.jl line 1221)
valid_idx = .!isnan.(obj_history) .& .!isnan.(r_norm_history) .& .!isnan.(s_norm_history)
obj_history = obj_history[valid_idx]
r_norm_history = r_norm_history[valid_idx]
s_norm_history = s_norm_history[valid_idx]
ρ_history = ρ_history[valid_idx]

# Setup (exact replication of Plotter.jl lines 1233-1288)
gr()
theme(:mute)
line_colour_obj = :dodgerblue
line_colour_primal = :darkgreen
line_colour_dual = :darkorange2
line_colour_rho = :purple

n_iter = length(obj_history)
iterations = 1:n_iter
println("Plotting $n_iter iterations")

markerstrokewidth = if n_iter <= 50
    1.5
elseif n_iter <= 100
    1.0
elseif n_iter <= 200
    0.5
else
    0.3
end

residual_markersize = if n_iter <= 50
    4
elseif n_iter <= 100
    3
elseif n_iter <= 200
    2
else
    1.5
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

# Subplot 1: Objective (Plotter.jl lines 1291-1329)
p1 = plot(
    iterations, obj_history,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Objective Function [\$]",
    title="Objective",
    lw=3,
    color=line_colour_obj,
    markershape=:circle,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    label="tADMM Objective",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    minorgridstyle=:solid,
    minorgridalpha=0.15,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern",
    top_margin=2Plots.mm
)

final_obj_tadmm = last(obj_history)
plot!(p1, [n_iter], [final_obj_tadmm],
      seriestype=:scatter, markersize=0, label="tADMM Final = \$$(round(final_obj_tadmm, digits=2))")

# Subplot 2: Primal residual (Plotter.jl lines 1367-1400)
p2 = plot(
    iterations, r_norm_history,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Primal Residual ‖r‖ [log scale]",
    title="Primal Residual",
    lw=3,
    color=line_colour_primal,
    markershape=:square,
    markersize=residual_markersize,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    yscale=:log10,
    label="Primal Residual (‖r‖)",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    minorgridstyle=:solid,
    minorgridalpha=0.15,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern"
)

hline!(p2, [eps_pri],
       color=:red, lw=2, linestyle=:dash, alpha=0.7,
       label="Threshold ε_pri = $(eps_pri)")

# Subplot 3: Dual residual (Plotter.jl lines 1402-1436)
p3 = plot(
    iterations, s_norm_history,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Dual Residual ‖s‖ [log scale]",
    title="Dual Residual",
    lw=3,
    color=line_colour_dual,
    markershape=:diamond,
    markersize=residual_markersize,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    yscale=:log10,
    label="Dual Residual (‖s‖)",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    minorgridstyle=:solid,
    minorgridalpha=0.15,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern",
    bottom_margin=2Plots.mm
)

hline!(p3, [eps_dual],
       color=:red, lw=2, linestyle=:dash, alpha=0.7,
       label="Threshold ε_dual = $(eps_dual)")

# Subplot 4: Adaptive ρ (Plotter.jl lines 1438-1481)
p4 = plot(
    iterations, ρ_history,
    dpi=600,
    xlabel="Iteration (k)",
    ylabel="Penalty Parameter ρ [log scale]",
    title="Adaptive ρ Schedule",
    lw=3,
    color=line_colour_rho,
    markershape=:hexagon,
    markersize=4,
    markerstrokecolor=:black,
    markerstrokewidth=markerstrokewidth,
    yscale=:log10,
    label="ρ value",
    legend=:topright,
    legendfontsize=9,
    grid=true,
    gridstyle=:solid,
    gridalpha=0.3,
    minorgrid=false,
    minorgridstyle=:solid,
    minorgridalpha=0.15,
    xlims=(0.5, n_iter + 0.5),
    xticks=xtick_vals,
    titlefont=font(12, "Computer Modern"),
    guidefont=font(12, "Computer Modern"),
    tickfontfamily="Computer Modern",
    bottom_margin=2Plots.mm
)

# Combine (Plotter.jl lines 1559-1569)
p_combined = plot(p1, p2, p3, p4,
                 layout=(4, 1),
                 size=(900, 1300),
                 plot_title="tADMM Convergence Summary",
                 plot_titlefontsize=14,
                 plot_titlefontfamily="Computer Modern",
                 left_margin=8Plots.mm,
                 right_margin=5Plots.mm,
                 top_margin=8Plots.mm,
                 bottom_margin=5Plots.mm)

# Print summary (Plotter.jl lines 1572-1594)
final_obj = last(obj_history)
final_r_norm = last(r_norm_history)
final_s_norm = last(s_norm_history)
converged = (final_r_norm ≤ eps_pri) && (final_s_norm ≤ eps_dual)

println("\ntADMM Convergence Summary:")
println("  Total iterations: $(length(obj_history))")
@printf "  Final objective:  \$%.6f\n" final_obj
@printf "  Final ‖r‖:        %.2e (threshold: %.1e) %s\n" final_r_norm eps_pri (final_r_norm ≤ eps_pri ? "✓" : "✗")
@printf "  Final ‖s‖:        %.2e (threshold: %.1e) %s\n" final_s_norm eps_dual (final_s_norm ≤ eps_dual ? "✓" : "✗")
println("  Converged:        $(converged ? "✅ YES" : "❌ NO")")

# Save (Plotter.jl lines 1602-1605)
@printf "Saving tADMM convergence plot to: %s\n" output_path
savefig(p_combined, output_path)
println("✓ Done")
