#!/usr/bin/env julia
# plot_bf_vs_tadmm_trajectory.jl
#
# Beautiful overlay of BF (Ipopt interior-point) vs tADMM convergence
# trajectories on the same wall-clock time axis. Color-codes BF by
# feasibility (red dashed = constraint-violating intermediate iterates,
# green solid = feasible search) to highlight that BF cannot be early-stopped
# while tADMM produces a feasible solution at every iteration.
#
# Usage:
#   julia plot_bf_vs_tadmm_trajectory.jl <bf_dir> <tadmm_csv> [output_png]
#
# Example:
#   julia plot_bf_vs_tadmm_trajectory.jl \
#       envs/tadmm/processedData/large10kC_1ph_T48 \
#       envs/tadmm/processedData/large10kC_1ph_T48/rho_sweep/rho_30000/convergence_data.csv

using Plots, Printf, Statistics

const FEAS_TOL = 1e-4   # primal infeasibility threshold for "feasible"
const GAP_TOL  = 0.005  # 0.5% gap band

# --- BF Ipopt log parser ------------------------------------------------------

"""
Parse Ipopt log file. Returns vector of (k, obj, inf_pr, inf_du) tuples.
Lines look like: "  30  2.8197835e+03 4.44e-16 6.77e-07  -2.5 ..."
"""
function parse_ipopt_log(path::String)
    iters = Tuple{Int,Float64,Float64,Float64}[]
    rx = r"^\s*(\d+)r?\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
    for line in eachline(path)
        m = match(rx, line)
        isnothing(m) && continue
        k = parse(Int, m.captures[1])
        obj = parse(Float64, m.captures[2])
        inf_pr = parse(Float64, m.captures[3])
        inf_du = parse(Float64, m.captures[4])
        push!(iters, (k, obj, inf_pr, inf_du))
    end
    iters
end

function parse_bf_total_time(results_path::String)
    isfile(results_path) || return NaN
    for line in eachline(results_path)
        m = match(r"Wall-clock time:\s*([\d.]+)\s*seconds", line)
        if !isnothing(m)
            return parse(Float64, m.captures[1])
        end
    end
    return NaN
end

# --- tADMM CSV parser ---------------------------------------------------------

function parse_tadmm_csv(path::String)
    lines = readlines(path)
    header = split(lines[1], ",")
    cols = Dict(c => Float64[] for c in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, ",")
        for (i, c) in enumerate(header)
            try
                push!(cols[c], parse(Float64, vals[i]))
            catch
                push!(cols[c], NaN)
            end
        end
    end
    return cols
end

# --- main plotting -----------------------------------------------------------

function generate_overlay(bf_dir::String, tadmm_csv::String, output::String)
    println("\nLoading BF Ipopt log from: $bf_dir/ipopt_bf.log")
    bf_iters = parse_ipopt_log(joinpath(bf_dir, "ipopt_bf.log"))
    bf_total_time = parse_bf_total_time(joinpath(bf_dir, "results_socp_bf.txt"))

    println("Loading tADMM trajectory from: $tadmm_csv")
    tadmm = parse_tadmm_csv(tadmm_csv)

    # Final objectives
    bf_final = bf_iters[end][2]
    tadmm_final = last(tadmm["objective"])

    # Use BF as the reference for "optimal"; tADMM gap is computed against it.
    ref_obj = bf_final
    band_lo = ref_obj * (1 - GAP_TOL)
    band_hi = ref_obj * (1 + GAP_TOL)

    # Estimate per-iter time for BF (Ipopt doesn't log it; use uniform mean).
    bf_n = length(bf_iters)
    dt_bf = bf_total_time / bf_n
    bf_t = [(i-1) * dt_bf for i in 1:bf_n]
    bf_obj  = [it[2] for it in bf_iters]
    bf_inf  = [it[3] for it in bf_iters]
    bf_feas = bf_inf .<= FEAS_TOL

    # Find first feasible and first feasible+near-optimal
    first_feas_idx = findfirst(bf_feas)
    first_feas_t = isnothing(first_feas_idx) ? NaN : bf_t[first_feas_idx]
    first_near_idx = findfirst(i -> bf_feas[i] && abs(bf_obj[i] - ref_obj) <= GAP_TOL * ref_obj, 1:bf_n)
    first_near_t = isnothing(first_near_idx) ? NaN : bf_t[first_near_idx]

    # tADMM trajectory on cumulative effective time
    t_tadmm = tadmm["cum_eff_time"]
    obj_tadmm = tadmm["objective"]
    r_tadmm = tadmm["r_norm"]

    # Find tADMM near-optimal point (first k where r<=eps_pri and gap<=0.5%)
    # Parse eps_pri from the run's results file so we use the actual tolerance.
    tadmm_dir = dirname(tadmm_csv)
    results_path = joinpath(tadmm_dir, "results_socp_tadmm.txt")
    eps_pri = 1e-3  # fallback default
    if isfile(results_path)
        for line in eachline(results_path)
            m = match(r"Primal tolerance:\s*([\d.eE+-]+)", line)
            if !isnothing(m); eps_pri = parse(Float64, m.captures[1]); break; end
        end
    end
    @printf "Using tADMM eps_pri = %.2e (parsed from results file)\n" eps_pri
    near_idx_t = findfirst(i -> r_tadmm[i] <= eps_pri && abs(obj_tadmm[i] - ref_obj) <= GAP_TOL * ref_obj, 1:length(obj_tadmm))
    near_t_tadmm = isnothing(near_idx_t) ? NaN : t_tadmm[near_idx_t]
    near_obj_tadmm = isnothing(near_idx_t) ? NaN : obj_tadmm[near_idx_t]

    # Split BF trajectory into infeasible / feasible segments for color-coding
    bf_t_inf = [bf_feas[i] ? NaN : bf_t[i] for i in 1:bf_n]
    bf_obj_inf = [bf_feas[i] ? NaN : bf_obj[i] for i in 1:bf_n]
    bf_t_fea = [bf_feas[i] ? bf_t[i] : NaN for i in 1:bf_n]
    bf_obj_fea = [bf_feas[i] ? bf_obj[i] : NaN for i in 1:bf_n]

    # Time unit autoscale
    max_t = max(maximum(bf_t), last(t_tadmm))
    if max_t > 3600
        scale = 1/3600; tunit = "hours"
    elseif max_t > 60
        scale = 1/60; tunit = "minutes"
    else
        scale = 1.0; tunit = "seconds"
    end

    # ---- Plot ----
    gr()
    theme(:mute)

    title_str = @sprintf("BF vs tADMM Convergence Trajectory\n%s", basename(bf_dir))

    p = plot(
        dpi=600, size=(1100, 700),
        title=title_str,
        xlabel="Wall-clock time [$tunit]",
        ylabel="Objective value [\$]",
        titlefont=font(13, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        legendfontsize=10,
        legend=:bottomright,
        grid=true, gridalpha=0.3,
        left_margin=10Plots.mm, right_margin=8Plots.mm,
        top_margin=4Plots.mm, bottom_margin=8Plots.mm,
    )

    # ±0.5% gap band around BF optimal
    plot!(p, [0, max_t*scale], [band_lo, band_lo],
          fillrange=[band_hi, band_hi], fillalpha=0.15, fillcolor=:limegreen,
          linealpha=0, label="±0.5% optimality band")

    # BF infeasible portion (red dashed)
    plot!(p, bf_t_inf .* scale, bf_obj_inf,
          lw=1.8, lc=:firebrick, ls=:dash, label="BF (infeasible search trajectory)")

    # BF feasible portion (green solid)
    plot!(p, bf_t_fea .* scale, bf_obj_fea,
          lw=2.5, lc=:darkgreen, label="BF (feasible iterates)")

    # tADMM trajectory (blue solid; feasible at every iteration)
    plot!(p, t_tadmm .* scale, obj_tadmm,
          lw=2.5, lc=:dodgerblue, marker=:circle, ms=3, msw=0.4,
          label="tADMM (feasible per subproblem at every iter)")

    # Final objective horizontal line
    hline!(p, [ref_obj], lw=1, lc=:black, ls=:dot, alpha=0.6,
           label=@sprintf("BF optimal = \$%.2f", ref_obj))

    # Annotate the milestones
    if !isnan(first_feas_t)
        scatter!(p, [first_feas_t * scale], [bf_obj[first_feas_idx]],
                 mc=:darkgreen, ms=8, msw=1.5, msc=:black,
                 label=@sprintf("BF first feasible (k=%d, %.1f%% of run)", first_feas_idx, 100*first_feas_idx/bf_n))
    end
    if !isnan(first_near_t)
        scatter!(p, [first_near_t * scale], [bf_obj[first_near_idx]],
                 mc=:gold, ms=10, msw=1.5, msc=:black, marker=:star5,
                 label=@sprintf("BF first feasible+0.5%% (k=%d)", first_near_idx))
    end
    if !isnan(near_t_tadmm)
        scatter!(p, [near_t_tadmm * scale], [near_obj_tadmm],
                 mc=:dodgerblue, ms=10, msw=1.5, msc=:black, marker=:diamond,
                 label=@sprintf("tADMM near-optimal (k=%d)", near_idx_t))
    end

    mkpath(dirname(output))
    savefig(p, output)
    println("✓ Saved: $output")

    # Print summary
    println()
    println("="^70)
    println("SUMMARY: $(basename(bf_dir))")
    println("="^70)
    @printf("BF: %d iters, %.1f s total wall-clock, final obj \$%.2f\n", bf_n, bf_total_time, bf_final)
    if !isnan(first_feas_t)
        @printf("  First FEASIBLE (inf_pr<=%g): k=%d (%.1f%% of run, ~%.1f s)\n",
                FEAS_TOL, first_feas_idx, 100*first_feas_idx/bf_n, first_feas_t)
    end
    if !isnan(first_near_t)
        @printf("  First feasible & ≤0.5%% gap: k=%d (%.1f%% of run, ~%.1f s)\n",
                first_near_idx, 100*first_near_idx/bf_n, first_near_t)
    end
    @printf("tADMM: %d iters, %.1f s effective time, final obj \$%.2f\n",
            length(obj_tadmm), last(t_tadmm), tadmm_final)
    if !isnan(near_t_tadmm)
        @printf("  Near-optimal: k=%d at %.1f s\n", near_idx_t, near_t_tadmm)
        if !isnan(first_near_t)
            @printf("  Speedup vs BF first feasible+0.5%%: %.2fx\n", first_near_t / near_t_tadmm)
        end
    end
    println("="^70)
end

# --- entry point --------------------------------------------------------------

if length(ARGS) < 2
    println("Usage: julia plot_bf_vs_tadmm_trajectory.jl <bf_dir> <tadmm_csv> [output_png]")
    exit(1)
end

bf_dir = ARGS[1]
tadmm_csv = ARGS[2]
output = length(ARGS) >= 3 ? ARGS[3] : joinpath(bf_dir, "convergence", "bf_vs_tadmm_trajectory.png")

generate_overlay(bf_dir, tadmm_csv, output)
