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
eps_pri = 1e-3
eps_dual = 1e-2

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

    if !isnan(bf_obj)
        hline!(p1, [bf_obj],
               color=:darkorange, lw=2, linestyle=:dash, alpha=0.8,
               label="BF Optimal = \$$(round(bf_obj, digits=2))")
    end

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
