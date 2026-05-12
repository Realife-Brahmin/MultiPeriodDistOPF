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

# Configuration — pass directories as args, or defaults to single path
# Usage:
#   julia plot_convergence_from_csv.jl                     # default single dir
#   julia plot_convergence_from_csv.jl --sweep-dir <path>  # all rho_* subdirs
#   julia plot_convergence_from_csv.jl <csv> <output>      # explicit paths
const DEFAULT_EPS_PRI  = 1e-3
const DEFAULT_EPS_DUAL = 1e-2

# Parse eps_pri / eps_dual from the run's results_socp_tadmm.txt so plotted
# threshold lines reflect the actual config (not a hardcoded default).
function parse_eps_thresholds(csv_path)
    results_path = joinpath(dirname(csv_path), "results_socp_tadmm.txt")
    eps_pri  = DEFAULT_EPS_PRI
    eps_dual = DEFAULT_EPS_DUAL
    isfile(results_path) || return (eps_pri, eps_dual)
    for line in eachline(results_path)
        m = match(r"Primal tolerance:\s*([\d.eE+-]+)", line)
        if !isnothing(m); eps_pri = parse(Float64, m.captures[1]); end
        m = match(r"Dual tolerance:\s*([\d.eE+-]+)", line)
        if !isnothing(m); eps_dual = parse(Float64, m.captures[1]); end
    end
    return (eps_pri, eps_dual)
end

function find_bf_objective(csv_path)
    dir = dirname(csv_path)
    for depth in 0:3
        bf_file = joinpath(dir, "results_socp_bf.txt")
        if isfile(bf_file)
            for line in eachline(bf_file)
                m = match(r"Total Cost: \$([0-9.]+)", line)
                if !isnothing(m)
                    return parse(Float64, m.captures[1])
                end
            end
        end
        dir = dirname(dir)
    end
    return NaN
end

function collect_jobs()
    if length(ARGS) >= 1 && ARGS[1] == "--sweep-dir"
        sweep_dir = ARGS[2]
        jobs = Tuple{String,String}[]
        for d in sort(readdir(sweep_dir))
            if startswith(d, "rho_")
                csv = joinpath(sweep_dir, d, "convergence_data.csv")
                out = joinpath(sweep_dir, d, "convergence", "tadmm_convergence_socp.png")
                if isfile(csv)
                    push!(jobs, (csv, out))
                end
            end
        end
        return jobs
    elseif length(ARGS) >= 2
        return [(ARGS[1], ARGS[2])]
    else
        return [("envs/tadmm/processedData/large10kC_1ph_T6/convergence_data.csv",
                 "envs/tadmm/processedData/large10kC_1ph_T6/convergence/tadmm_convergence_socp.png")]
    end
end

function generate_plot(csv_path, output_path)
    println("\nLoading convergence data from: $csv_path")
    conv_data = parse_csv(csv_path)

    obj_history = conv_data["objective"]
    r_norm_history = conv_data["r_norm"]
    s_norm_history = conv_data["s_norm"]
    ρ_history = conv_data["rho"]

    valid_idx = .!isnan.(obj_history) .& .!isnan.(r_norm_history) .& .!isnan.(s_norm_history)
    obj_history = obj_history[valid_idx]
    r_norm_history = r_norm_history[valid_idx]
    s_norm_history = s_norm_history[valid_idx]
    ρ_history = ρ_history[valid_idx]

    bf_obj = find_bf_objective(csv_path)
    rho_init = first(ρ_history)
    eps_pri, eps_dual = parse_eps_thresholds(csv_path)
    @printf "  eps_pri=%.1e eps_dual=%.1e (parsed from results file)\n" eps_pri eps_dual

    gr()
    theme(:default)
    # Colorblind-friendly palette (Wong 2011)
    line_colour_obj    = RGB(0/255, 114/255, 178/255)   # blue
    line_colour_primal = RGB(0/255, 158/255, 115/255)   # bluish green
    line_colour_dual   = RGB(230/255, 159/255, 0/255)   # orange
    line_colour_rho    = RGB(204/255, 121/255, 167/255) # reddish purple
    colour_bf_ref      = RGB(213/255, 94/255, 0/255)    # vermillion
    colour_threshold   = RGB(80/255, 80/255, 80/255)    # neutral gray

    n_iter = length(obj_history)
    iterations = 1:n_iter
    println("Plotting $n_iter iterations")

    # Subsample marker indices: target ~20 markers regardless of iteration count
    n_markers_target = 20
    mark_step = max(1, div(n_iter, n_markers_target))
    mark_idx = sort(unique(vcat(1, mark_step:mark_step:n_iter, n_iter)))

    markersize = 5.5
    markerstrokewidth = 0.8

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

    p1 = plot(
        iterations, obj_history,
        dpi=600,
        xlabel="Iteration (k)",
        ylabel="Objective Function [\$]",
        title="Objective",
        lw=2.5,
        color=line_colour_obj,
        label="tADMM Objective",
        legend=:topright,
        legendfontsize=9,
        grid=true, gridalpha=0.25,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5),
        xticks=xtick_vals,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=2Plots.mm,
        background_color_inside=:white,
    )
    scatter!(p1, iterations[mark_idx], obj_history[mark_idx],
             markershape=:circle, markersize=markersize,
             color=line_colour_obj,
             markerstrokecolor=:white, markerstrokewidth=markerstrokewidth,
             label="")

    final_obj_tadmm = last(obj_history)
    plot!(p1, [n_iter], [final_obj_tadmm],
          seriestype=:scatter, markersize=0, label="tADMM Final = \$$(round(final_obj_tadmm, digits=2))")

    if !isnan(bf_obj)
        hline!(p1, [bf_obj],
               color=colour_bf_ref, lw=1.8, linestyle=:dash, alpha=0.85,
               label="BF Optimal = \$$(round(bf_obj, digits=2))")
    end

    p2 = plot(
        iterations, r_norm_history,
        dpi=600,
        xlabel="Iteration (k)",
        ylabel="Primal Residual ‖r‖ [log scale]",
        title="Primal Residual",
        lw=2.5,
        color=line_colour_primal,
        yscale=:log10,
        label="Primal Residual (‖r‖)",
        legend=:topright,
        legendfontsize=9,
        grid=true, gridalpha=0.25,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5),
        xticks=xtick_vals,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        background_color_inside=:white,
    )
    scatter!(p2, iterations[mark_idx], r_norm_history[mark_idx],
             markershape=:square, markersize=markersize,
             color=line_colour_primal,
             markerstrokecolor=:white, markerstrokewidth=markerstrokewidth,
             label="")

    hline!(p2, [eps_pri],
           color=colour_threshold, lw=1.5, linestyle=:dash, alpha=0.75,
           label="Threshold ε_pri = $(eps_pri)")

    p3 = plot(
        iterations, s_norm_history,
        dpi=600,
        xlabel="Iteration (k)",
        ylabel="Dual Residual ‖s‖ [log scale]",
        title="Dual Residual",
        lw=2.5,
        color=line_colour_dual,
        yscale=:log10,
        label="Dual Residual (‖s‖)",
        legend=:topright,
        legendfontsize=9,
        grid=true, gridalpha=0.25,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5),
        xticks=xtick_vals,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        bottom_margin=2Plots.mm,
        background_color_inside=:white,
    )
    scatter!(p3, iterations[mark_idx], s_norm_history[mark_idx],
             markershape=:diamond, markersize=markersize,
             color=line_colour_dual,
             markerstrokecolor=:white, markerstrokewidth=markerstrokewidth,
             label="")

    hline!(p3, [eps_dual],
           color=colour_threshold, lw=1.5, linestyle=:dash, alpha=0.75,
           label="Threshold ε_dual = $(eps_dual)")

    p4 = plot(
        iterations, ρ_history,
        dpi=600,
        xlabel="Iteration (k)",
        ylabel="Penalty Parameter ρ [log scale]",
        title="Adaptive ρ Schedule",
        lw=2.5,
        color=line_colour_rho,
        yscale=:log10,
        label="ρ value",
        legend=:topright,
        legendfontsize=9,
        grid=true, gridalpha=0.25,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5),
        xticks=xtick_vals,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        bottom_margin=2Plots.mm,
        background_color_inside=:white,
    )
    scatter!(p4, iterations[mark_idx], ρ_history[mark_idx],
             markershape=:hexagon, markersize=markersize,
             color=line_colour_rho,
             markerstrokecolor=:white, markerstrokewidth=markerstrokewidth,
             label="")

    rho_str = rho_init == round(rho_init) ? @sprintf("%.0f", rho_init) : @sprintf("%.1f", rho_init)
    p_combined = plot(p1, p2, p3, p4,
                     layout=(4, 1),
                     size=(900, 1300),
                     plot_title="tADMM Convergence Summary  (rho_init = $rho_str)",
                     plot_titlefontsize=14,
                     plot_titlefontfamily="Computer Modern",
                     left_margin=8Plots.mm,
                     right_margin=5Plots.mm,
                     top_margin=8Plots.mm,
                     bottom_margin=5Plots.mm)

    final_obj = last(obj_history)
    final_r_norm = last(r_norm_history)
    final_s_norm = last(s_norm_history)
    converged = (final_r_norm ≤ eps_pri) && (final_s_norm ≤ eps_dual)

    println("  Total iterations: $(length(obj_history))")
    @printf "  Final objective:  \$%.6f\n" final_obj
    @printf "  Final ‖r‖:        %.2e (threshold: %.1e) %s\n" final_r_norm eps_pri (final_r_norm ≤ eps_pri ? "✓" : "✗")
    @printf "  Final ‖s‖:        %.2e (threshold: %.1e) %s\n" final_s_norm eps_dual (final_s_norm ≤ eps_dual ? "✓" : "✗")
    println("  Converged:        $(converged ? "YES" : "NO")")

    mkpath(dirname(output_path))
    savefig(p_combined, output_path)
    @printf "  Saved: %s\n" output_path
end

# Run all jobs
jobs = collect_jobs()
println("Generating $(length(jobs)) convergence plot(s)...")
for (csv, out) in jobs
    generate_plot(csv, out)
end
println("\nDone — $(length(jobs)) plot(s) generated.")
