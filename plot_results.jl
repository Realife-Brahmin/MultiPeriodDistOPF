# ============================================================================
# plot_results.jl — Standalone plotting from CSV/TXT files
# ============================================================================
# Usage: edit the PLOT_CONFIG below, then run this script.
# No running Julia session or solution objects needed — reads files only.
# ============================================================================

ENV["GKSwstype"] = "png"

# Activate the tadmm environment (for Plots dependency)
import Pkg
Pkg.activate(joinpath(@__DIR__, "envs", "tadmm"))

using Plots
using Printf
using Statistics
using DelimitedFiles

gr()
theme(:mute)

# ============================================================================
# PLOT CONFIGURATION — edit these
# ============================================================================
const PLOT_CONFIG = (
    system = "large10kC_1ph",
    T = 24,
    eps_pri = 1e-3,
    eps_dual = 1e-2,
    show_plots = false,
    save_plots = true,
    dpi = 600,
)

const RESULTS_DIR = joinpath(@__DIR__, "envs", "tadmm", "processedData",
                             "$(PLOT_CONFIG.system)_T$(PLOT_CONFIG.T)")
const CONV_DIR = joinpath(RESULTS_DIR, "convergence")

# ============================================================================
# HELPERS
# ============================================================================

function read_csv_as_dict(filepath)
    if !isfile(filepath)
        println("File not found: $filepath")
        return nothing
    end
    lines = readlines(filepath)
    header = split(lines[1], ',')
    data = Dict{String, Vector{Float64}}()
    for col in header
        data[col] = Float64[]
    end
    for line in lines[2:end]
        vals = split(line, ',')
        for (i, col) in enumerate(header)
            push!(data[col], parse(Float64, vals[i]))
        end
    end
    return data
end

function load_bf_objective(results_dir)
    bf_file = joinpath(results_dir, "results_socp_bf.txt")
    if !isfile(bf_file)
        return NaN
    end
    for line in eachline(bf_file)
        m = match(r"Total Cost: \$([0-9.]+)", line)
        if !isnothing(m)
            return parse(Float64, m.captures[1])
        end
    end
    return NaN
end

function load_bf_wallclock(results_dir)
    bf_file = joinpath(results_dir, "results_socp_bf.txt")
    if !isfile(bf_file)
        return NaN
    end
    for line in eachline(bf_file)
        m = match(r"Wall-clock time: ([0-9.]+) seconds", line)
        if !isnothing(m)
            return parse(Float64, m.captures[1])
        end
    end
    return NaN
end

function adaptive_xticks(n_iter)
    if n_iter <= 10
        return collect(1:n_iter)
    elseif n_iter <= 50
        ticks = collect(10:10:n_iter)
        if isempty(ticks) || (n_iter - last(ticks)) > 3
            return sort(unique(vcat(1, ticks, n_iter)))
        else
            return sort(unique(vcat(1, ticks[1:end-1], n_iter)))
        end
    else
        step = max(1, div(n_iter, 5))
        ticks = collect(step:step:n_iter)
        if isempty(ticks) || (n_iter - last(ticks)) > step * 0.3
            return sort(unique(vcat(1, ticks, n_iter)))
        else
            return sort(unique(vcat(1, ticks[1:end-1], n_iter)))
        end
    end
end

function adaptive_markerstroke(n_iter)
    n_iter <= 50 ? 1.5 : n_iter <= 100 ? 1.0 : n_iter <= 200 ? 0.5 : 0.3
end

function adaptive_markersize(n_iter)
    n_iter <= 50 ? 4 : n_iter <= 100 ? 3 : n_iter <= 200 ? 2 : 1.5
end

# ============================================================================
# CONVERGENCE PLOT (4 panels: objective, primal residual, dual residual, rho)
# ============================================================================

function plot_convergence(; config=PLOT_CONFIG)
    conv_data = read_csv_as_dict(joinpath(RESULTS_DIR, "convergence_data.csv"))
    if isnothing(conv_data)
        println("No convergence_data.csv found in $RESULTS_DIR")
        return
    end

    iterations = Int.(conv_data["iteration"])
    obj = conv_data["objective"]
    r_norm = conv_data["r_norm"]
    s_norm = conv_data["s_norm"]
    rho = conv_data["rho"]
    n_iter = length(iterations)

    # Clamp zeros to a sensible floor (2 decades below smallest nonzero value)
    # so log10(0) doesn't kill the GR series or blow up the y-axis range
    function safe_floor(v)
        nz = filter(x -> x > 0, v)
        floor_val = isempty(nz) ? 1e-16 : minimum(nz) * 1e-2
        return max.(v, floor_val)
    end
    r_norm = safe_floor(r_norm)
    s_norm = safe_floor(s_norm)
    rho = safe_floor(rho)

    bf_obj = load_bf_objective(RESULTS_DIR)

    xticks = adaptive_xticks(n_iter)
    msw = adaptive_markerstroke(n_iter)
    ms = adaptive_markersize(n_iter)

    # Panel 1: Objective
    p1 = plot(iterations, obj,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Objective Function [\$]",
        title="Objective", lw=3, color=:dodgerblue,
        markershape=:circle, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        label="tADMM Objective", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern", top_margin=2Plots.mm)

    if isfinite(bf_obj)
        hline!(p1, [bf_obj], color=:darkorange3, lw=3, linestyle=:dash, alpha=0.9,
               label="BF Objective = \$$(round(bf_obj, digits=2))")
    end
    final_obj = last(obj)
    plot!(p1, [n_iter], [final_obj], seriestype=:scatter, markersize=0,
          label="tADMM Final = \$$(round(final_obj, digits=2))")

    # Panel 2: Primal residual
    p2 = plot(iterations, r_norm,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Primal Residual ||r|| [log]",
        title="Primal Residual", lw=3, color=:darkgreen,
        markershape=:square, markersize=ms, markerstrokecolor=:black, markerstrokewidth=msw,
        yscale=:log10, label="Primal (||r||)", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=true, minorgridstyle=:solid, minorgridalpha=0.15,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern")
    hline!(p2, [config.eps_pri], color=:red, lw=2, linestyle=:dash, alpha=0.7,
           label="eps_pri = $(config.eps_pri)")

    # Panel 3: Dual residual
    p3 = plot(iterations, s_norm,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Dual Residual ||s|| [log]",
        title="Dual Residual", lw=3, color=:darkorange2,
        markershape=:diamond, markersize=ms, markerstrokecolor=:black, markerstrokewidth=msw,
        yscale=:log10, label="Dual (||s||)", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=true, minorgridstyle=:solid, minorgridalpha=0.15,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern", bottom_margin=2Plots.mm)
    hline!(p3, [config.eps_dual], color=:red, lw=2, linestyle=:dash, alpha=0.7,
           label="eps_dual = $(config.eps_dual)")

    # Panel 4: Adaptive rho
    p4 = plot(iterations, rho,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Penalty Parameter rho [log]",
        title="Adaptive rho Schedule", lw=3, color=:purple,
        markershape=:hexagon, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        yscale=:log10, label="rho", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=true, minorgridstyle=:solid, minorgridalpha=0.15,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern", bottom_margin=2Plots.mm)

    # Combine
    combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(1400, 1100),
                    plot_title="tADMM Convergence — $(config.system) T=$(config.T)",
                    margin=5Plots.mm)

    if config.save_plots
        mkpath(CONV_DIR)
        outpath = joinpath(CONV_DIR, "tadmm_convergence_socp.png")
        savefig(combined, outpath)
        println("Convergence plot saved: $outpath")
    end
    config.show_plots && display(combined)
    return combined
end

# ============================================================================
# TIMING PLOT (4 panels: max SP, median SP, total, cumulative)
# ============================================================================

function plot_timing(; config=PLOT_CONFIG)
    conv_data = read_csv_as_dict(joinpath(RESULTS_DIR, "convergence_data.csv"))
    if isnothing(conv_data)
        println("No convergence_data.csv found in $RESULTS_DIR")
        return
    end

    iterations = Int.(conv_data["iteration"])
    max_sp = conv_data["max_sp_time_this_k"]
    eff_time = conv_data["eff_time_this_k"]
    cum_eff = conv_data["cum_eff_time"]
    n_iter = length(iterations)

    # Try to get median/total from subproblem_timing_details.csv
    sp_csv = joinpath(RESULTS_DIR, "subproblem_timing_details.csv")
    median_times = Float64[]
    total_times = Float64[]
    if isfile(sp_csv)
        sp_data = read_csv_as_dict(sp_csv)
        if !isnothing(sp_data)
            sp_iters = Int.(sp_data["iteration"])
            sp_times = sp_data["solve_time_sec"]
            for k in 1:n_iter
                mask = sp_iters .== k
                times_k = sp_times[mask]
                push!(median_times, isempty(times_k) ? 0.0 : median(times_k))
                push!(total_times, isempty(times_k) ? 0.0 : sum(times_k))
            end
        end
    end
    if isempty(median_times)
        median_times = max_sp  # fallback
        total_times = max_sp
    end

    bf_wallclock = load_bf_wallclock(RESULTS_DIR)

    xticks = adaptive_xticks(n_iter)
    msw = adaptive_markerstroke(n_iter)

    # Panel 1: Max subproblem time
    p1 = plot(iterations, max_sp,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Max SP Time (s)",
        title="Max Subproblem Time (Parallel Bottleneck)", lw=3, color=:crimson,
        markershape=:circle, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        label="Max subproblem", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern")
    mean_max = mean(max_sp)
    hline!(p1, [mean_max], color=:darkgreen, lw=1.5, linestyle=:dot, alpha=0.6,
           label="Mean = $(round(mean_max, digits=2))s")

    # Panel 2: Median subproblem time
    p2 = plot(iterations, median_times,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Median SP Time (s)",
        title="Median Subproblem Time (Typical Case)", lw=3, color=:darkorange2,
        markershape=:square, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        label="Median subproblem", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern")
    hline!(p2, [mean(median_times)], color=:darkgreen, lw=1.5, linestyle=:dot, alpha=0.6,
           label="Mean = $(round(mean(median_times), digits=2))s")

    # Panel 3: Total iteration time
    p3 = plot(iterations, total_times,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Total Iter Time (s)",
        title="Total Iteration Time (Sum of All Subproblems)", lw=3, color=:purple,
        markershape=:diamond, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        label="Sum all subproblems", legend=:topright, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern")
    hline!(p3, [mean(total_times)], color=:darkgreen, lw=1.5, linestyle=:dot, alpha=0.6,
           label="Mean = $(round(mean(total_times), digits=2))s")

    # Panel 4: Cumulative wall-clock
    p4 = plot(iterations, cum_eff,
        dpi=config.dpi, xlabel="Iteration (k)", ylabel="Cumulative Time (s)",
        title="Cumulative Effective Time", lw=3, color=:dodgerblue,
        markershape=:hexagon, markersize=4, markerstrokecolor=:black, markerstrokewidth=msw,
        label="Cumulative", legend=:topleft, legendfontsize=9,
        grid=true, gridstyle=:solid, gridalpha=0.3,
        minorgrid=false,
        xlims=(0.5, n_iter + 0.5), xticks=xticks,
        titlefont=font(12, "Computer Modern"), guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern", bottom_margin=2Plots.mm)
    if isfinite(bf_wallclock)
        hline!(p4, [bf_wallclock], color=:red, lw=2, linestyle=:dash, alpha=0.7,
               label="BF wall-clock = $(round(bf_wallclock, digits=1))s")
    end

    combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(1400, 1000),
                    plot_title="tADMM Timing Analysis — $(config.system) T=$(config.T)")

    if config.save_plots
        mkpath(CONV_DIR)
        outpath = joinpath(CONV_DIR, "tadmm_timing_socp.png")
        savefig(combined, outpath)
        println("Timing plot saved: $outpath")
    end
    config.show_plots && display(combined)
    return combined
end

# ============================================================================
# RUN BOTH PLOTS
# ============================================================================

println("\n" * "="^80)
println("GENERATING PLOTS FROM CSV/TXT FILES")
println("="^80)
println("Results dir: $RESULTS_DIR")

plot_convergence()
plot_timing()

println("\nDone.")
