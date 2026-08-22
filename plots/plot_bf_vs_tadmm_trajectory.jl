#!/usr/bin/env julia
# plot_bf_vs_tadmm_trajectory.jl
#
# Publication-quality overlay of BF (Ipopt interior-point) vs tADMM convergence
# trajectories on a common wall-clock time axis.  Both algorithms are TRUNCATED
# at their first feasible-and-near-optimal iterate ("first NO" milestone) —
# beyond that point either algorithm would be terminated in practice.
#
# Color & style language:
#   tADMM line                  = dodgerblue (matches convergence plot obj color)
#   tADMM NO terminus marker    = darkorange diamond (matches BF-ref color)
#   BF infeasible portion       = maroon dashed (constraint-violating search)
#   BF feasible-but-pre-NO      = forestgreen solid (feasible, not yet ≤0.5%)
#   BF NO terminus marker       = forestgreen square
#   BF optimal horizontal ref   = black dotted thin line
#   Grid                        = minor only, light gray
#   Per-iter subsampled markers = ~20 per line, with notable iters preserved
#
# Usage:
#   julia --project=envs/tadmm plot_bf_vs_tadmm_trajectory.jl <bf_dir> <tadmm_csv> [output_png]

using Plots, Printf, Statistics, LaTeXStrings

const FEAS_TOL = 1e-4   # primal infeasibility threshold for "feasible"
const GAP_TOL  = 0.005  # 0.5% gap band

# Colors --------------------------------------------------------------------
const COL_TADMM_LINE = :dodgerblue          # tADMM trusted portion (solid)
const COL_TADMM_BAD  = :firebrick           # tADMM untrusted portion (brick red dashdotdot)
const COL_TADMM_END  = :dodgerblue           # tADMM NO marker — matches the trusted line color
const COL_BF_INFEAS  = :maroon              # BF infeasible (wine red dashed)
const COL_BF_FEAS    = :forestgreen         # BF feasible (forest green solid)
const COL_BF_END     = :forestgreen         # BF NO marker
const COL_TEXT       = :black
const COL_BG         = :white
const COL_GRID       = RGB(180/255, 180/255, 180/255)  # neutral soft gray

# Marker target count -------------------------------------------------------
const N_MARKERS = 20

# --- BF Ipopt log parser ---------------------------------------------------

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

# --- tADMM CSV parser ------------------------------------------------------

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

# --- Marker subsampling ----------------------------------------------------

"""
Return ~`target` indices in [1, n] that include all `notable` indices.
"""
function subsample_markers(n::Int, notable::Vector{Int}; target::Int=N_MARKERS)
    n <= 0 && return Int[]
    if n <= target
        return collect(1:n)
    end
    sampled = round.(Int, range(1, n, length=target))
    combined = sort(unique(vcat(sampled, filter(i -> 1 <= i <= n, notable))))
    return combined
end

# --- system/horizon label from the BF directory name -----------------------

"""
Derive a human-friendly "<system>, T=<T>" title from a processedData dir name
like "large10kC_1ph_T48" or "ieee2522C_1ph_T144".  Falls back to the raw
basename if the pattern is unrecognized.
"""
function system_horizon_label(bf_dir::String)
    base = basename(rstrip(bf_dir, ['/', '\\']))
    m = match(r"^(.*?)_1ph_T(\d+)$", base)
    isnothing(m) && return base
    raw_sys, T = m.captures[1], m.captures[2]
    sysname = if occursin("large10k", raw_sys)
        "large10k"
    elseif occursin("2522", raw_sys)
        "medium2552"
    elseif occursin("123", raw_sys)
        "ieee123"
    else
        raw_sys
    end
    return "$sysname, T=$T"
end

# --- main plotting ---------------------------------------------------------

function generate_overlay(bf_dir::String, tadmm_csv::String, output::String)
    println("\nLoading BF Ipopt log from: $bf_dir/ipopt_bf.log")
    bf_iters = parse_ipopt_log(joinpath(bf_dir, "ipopt_bf.log"))
    bf_total_time = parse_bf_total_time(joinpath(bf_dir, "results_socp_bf.txt"))

    println("Loading tADMM trajectory from: $tadmm_csv")
    tadmm = parse_tadmm_csv(tadmm_csv)

    bf_final = bf_iters[end][2]
    ref_obj = bf_final

    # BF per-iter timing (uniform mean since Ipopt doesn't log per-iter wall-clock)
    bf_n = length(bf_iters)
    dt_bf = bf_total_time / bf_n
    bf_t = [(i-1) * dt_bf for i in 1:bf_n]
    bf_obj  = [it[2] for it in bf_iters]
    bf_inf  = [it[3] for it in bf_iters]
    bf_feas = bf_inf .<= FEAS_TOL

    first_feas_idx = findfirst(bf_feas)
    first_near_idx = findfirst(i -> bf_feas[i] && abs(bf_obj[i] - ref_obj) <= GAP_TOL * ref_obj, 1:bf_n)

    # tADMM cumulative effective time + parse eps_pri
    t_tadmm_full = tadmm["cum_eff_time"]
    obj_tadmm_full = tadmm["objective"]
    r_tadmm_full = tadmm["r_norm"]

    tadmm_dir = dirname(tadmm_csv)
    results_path = joinpath(tadmm_dir, "results_socp_tadmm.txt")
    eps_pri = 1e-3
    if isfile(results_path)
        for line in eachline(results_path)
            m = match(r"Primal tolerance:\s*([\d.eE+-]+)", line)
            if !isnothing(m); eps_pri = parse(Float64, m.captures[1]); break; end
        end
    end
    @printf "Using tADMM eps_pri = %.2e (parsed from results file)\n" eps_pri
    near_idx_t = findfirst(i -> r_tadmm_full[i] <= eps_pri && abs(obj_tadmm_full[i] - ref_obj) <= GAP_TOL * ref_obj,
                            1:length(obj_tadmm_full))

    # TRUNCATE both algorithms at their first NO iterate
    tadmm_end = isnothing(near_idx_t) ? length(obj_tadmm_full) : near_idx_t
    bf_end    = isnothing(first_near_idx) ? bf_n : first_near_idx

    # ---- tADMM: prepend (t=0, obj=obj[1]) so the trajectory starts at the origin,
    # matching BF's k=0 initial-point convention.  At t=0, no consensus exists
    # yet, so the leading point is treated as consensus-infeasible.
    t_tadmm   = vcat([0.0],            t_tadmm_full[1:tadmm_end])
    obj_tadmm = vcat([obj_tadmm_full[1]], obj_tadmm_full[1:tadmm_end])
    r_tadmm   = vcat([Inf],             r_tadmm_full[1:tadmm_end])  # Inf at t=0 ⇒ "infeasible" leading point

    # Split tADMM by consensus feasibility (r vs eps_pri) AND objective
    # sanity (obj >= f_BF).  Any iterate with f < f_BF cannot be truly
    # feasible (else it would be the global optimum), so we conservatively
    # mark those as infeasible regardless of r_norm.
    # Consensus feasibility per iterate (r vs eps_pri) AND objective sanity
    # (obj >= f_BF): any iterate with f < f_BF cannot be truly feasible (else
    # it would be the global optimum), so we conservatively mark those
    # infeasible regardless of r_norm.  Used to locate first feasibility; the
    # redesigned plot draws tADMM as one continuous hero curve rather than
    # separate feasible/infeasible segments.
    tadmm_feas_mask = (r_tadmm .<= eps_pri) .& (obj_tadmm .>= ref_obj)
    n_t = length(t_tadmm)

    # BF infeasible vs feasible portions (within the truncated [1, bf_end] range)
    if isnothing(first_feas_idx) || first_feas_idx > bf_end
        bf_t_inf = bf_t[1:bf_end]; bf_obj_inf = bf_obj[1:bf_end]
        bf_t_fea = Float64[]; bf_obj_fea = Float64[]
    else
        bf_t_inf = bf_t[1:first_feas_idx]
        bf_obj_inf = bf_obj[1:first_feas_idx]
        bf_t_fea = bf_t[first_feas_idx:bf_end]
        bf_obj_fea = bf_obj[first_feas_idx:bf_end]
    end

    near_t_tadmm   = isnothing(near_idx_t) ? NaN : t_tadmm_full[near_idx_t]
    near_obj_tadmm = isnothing(near_idx_t) ? NaN : obj_tadmm_full[near_idx_t]
    bf_no_t        = isnothing(first_near_idx) ? NaN : bf_t[first_near_idx]
    bf_no_obj      = isnothing(first_near_idx) ? NaN : bf_obj[first_near_idx]

    # Time-unit autoscale (driven by the TRUNCATED endpoints)
    max_t = max(isempty(bf_t_fea) ? maximum(bf_t_inf) : maximum(bf_t_fea), maximum(t_tadmm))
    if max_t > 3600
        scale = 1/3600; tunit = "hours"
    elseif max_t > 60
        scale = 1/60; tunit = "minutes"
    else
        scale = 1.0; tunit = "seconds"
    end

    # y-axis bounds: cap the enormous initial interior-point objective spike,
    # but leave headroom above the converged optimum for the speedup callout.
    obj_all  = vcat(filter(!isnan, bf_obj[1:bf_end]), filter(!isnan, obj_tadmm))
    obj_conv = vcat(filter(!isnan, bf_obj_fea), filter(!isnan, obj_tadmm))
    ycap     = quantile(obj_all, 0.95)
    ypeak    = isempty(obj_conv) ? ycap : maximum(obj_conv)
    yhi      = max(ycap, ypeak) * 1.14        # headroom for the "N× faster" arrow
    ylo      = minimum(obj_all) * 0.96

    # Mask any point above the cap so the steep initial spikes don't "shoot
    # out" the top of the frame as vertical streaks.
    capy(y) = [(!isnan(v) && v > yhi) ? NaN : v for v in y]
    bf_obj_inf = capy(bf_obj_inf)
    bf_obj_fea = capy(bf_obj_fea)

    # Money-unit y tick formatter (no scientific notation).
    yref = max(abs(yhi), abs(ylo))
    yfmt = if yref >= 1e6
        v -> @sprintf("\$%.1fM", v/1e6)
    elseif yref >= 1e3
        v -> @sprintf("\$%.1fk", v/1e3)
    else
        v -> @sprintf("\$%.0f", v)
    end

    # Subsample markers (fewer than before for a cleaner, bolder look) and
    # drop any that fall above the y-cap.
    bf_marks_global = subsample_markers(bf_end, [1, bf_end]; target=14)
    bf_marks_inf    = filter(i -> (isnothing(first_feas_idx) || i < first_feas_idx) && bf_obj[i] <= yhi, bf_marks_global)
    bf_marks_fea    = filter(i -> (!isnothing(first_feas_idx) && i >= first_feas_idx) && bf_obj[i] <= yhi, bf_marks_global)
    tadmm_mark      = subsample_markers(n_t, [1, n_t]; target=14)

    # ---- Plot ----
    gr()
    # Short, croppable title; the caption carries the "BF vs. tADMM" framing.
    title_str = system_horizon_label(bf_dir)

    p = plot(
        background_color=COL_BG,
        foreground_color=COL_TEXT,
        size=(1100, 650),
        dpi=300,
        title=title_str,
        titlefont=font(16, "Computer Modern"),
        titlefontcolor=COL_TEXT,
        xlabel="Wall-clock time [$tunit]",
        ylabel="Objective f(x)",
        guidefont=font(14, "Computer Modern"),
        guidefontcolor=COL_TEXT,
        tickfont=font(12, "Computer Modern"),
        tickfontcolor=COL_TEXT,
        yformatter=yfmt,
        legend=:bottomright,
        legendfont=font(11, "Computer Modern"),
        legendfontcolor=COL_TEXT,
        # Minor grid only
        grid=false,
        minorgrid=true,
        minorgridcolor=COL_GRID,
        minorgridlinewidth=0.6,
        minorgridalpha=0.30,
        minorgridstyle=:solid,
        framestyle=:box,
        xlims=(0, max_t * scale * 1.04),
        ylims=(ylo, yhi),
        left_margin=10Plots.mm, right_margin=8Plots.mm,
        top_margin=6Plots.mm, bottom_margin=8Plots.mm,
    )

    # BF infeasible search (maroon dashed) — the long slow interior-point climb
    plot!(p, bf_t_inf .* scale, bf_obj_inf,
          lw=4, lc=COL_BF_INFEAS, ls=:dash, label="BF — infeasible search")
    # BF feasible (forest green solid) — the brief feasible tail before NO
    if !isempty(bf_t_fea)
        plot!(p, bf_t_fea .* scale, bf_obj_fea,
              lw=6, lc=COL_BF_FEAS, ls=:solid, label="BF — feasible")
    end
    # tADMM hero curve — one bold continuous line from first feasibility to the
    # NO terminus.  The long low-objective ramp-up before first feasibility is
    # omitted so the curve reads as a fast approach, not a vertical streak.
    tfeas_first = findfirst(tadmm_feas_mask)
    tadmm_lo    = isnothing(tfeas_first) ? 1 : tfeas_first
    plot!(p, t_tadmm[tadmm_lo:end] .* scale, capy(obj_tadmm[tadmm_lo:end]),
          lw=8, lc=COL_TADMM_LINE, ls=:solid, label="tADMM")

    # Per-iter markers — light texture; B&W-safe shapes (BF triangle, tADMM circle)
    !isempty(bf_marks_inf) && scatter!(p, (bf_t[bf_marks_inf]) .* scale, bf_obj[bf_marks_inf],
             marker=:utriangle, mc=:white, ms=6, msw=1.4, msc=COL_BF_INFEAS, label="")
    !isempty(bf_marks_fea) && scatter!(p, (bf_t[bf_marks_fea]) .* scale, bf_obj[bf_marks_fea],
             marker=:utriangle, mc=:white, ms=6, msw=1.4, msc=COL_BF_FEAS, label="")
    tadmm_mark_hero = filter(i -> i >= tadmm_lo && obj_tadmm[i] <= yhi, tadmm_mark)
    !isempty(tadmm_mark_hero) && scatter!(p, (t_tadmm[tadmm_mark_hero]) .* scale,
             obj_tadmm[tadmm_mark_hero],
             marker=:circle, mc=:white, ms=6, msw=1.4, msc=COL_TADMM_LINE, label="")

    # Optimum reference line
    hline!(p, [ref_obj], lw=1.6, lc=COL_TEXT, ls=:dot, alpha=0.7, label="optimum")

    # ---- Terminus "arrived" markers — bold stars, no legend clutter ----
    if !isnan(bf_no_t)
        scatter!(p, [bf_no_t * scale], [min(bf_no_obj, yhi)],
                 mc=COL_BF_END, ms=18, msw=2.2, msc=:black, marker=:star5, label="")
    end
    if !isnan(near_t_tadmm)
        scatter!(p, [near_t_tadmm * scale], [min(near_obj_tadmm, yhi)],
                 mc=COL_TADMM_END, ms=18, msw=2.2, msc=:black, marker=:star5, label="")
    end

    # ---- Speedup callout — the headline ----
    if !isnan(bf_no_t) && !isnan(near_t_tadmm) && near_t_tadmm < bf_no_t
        spd  = bf_no_t / near_t_tadmm
        xa   = near_t_tadmm * scale
        xb   = bf_no_t * scale
        yarr = yhi * 0.93
        plot!(p, [xa, xb], [yarr, yarr], lc=:black, lw=2.2,
              arrow=arrow(:closed, :both), label="")
        annotate!(p, (xa + xb)/2, yhi*0.975,
                  text(@sprintf("%.2f× faster", spd), 16, :black, "Computer Modern"))
    end

    mkpath(dirname(output))
    savefig(p, output)
    println("✓ Saved: $output")

    # Print summary
    println()
    println("="^70)
    println("SUMMARY: $(basename(bf_dir))")
    println("="^70)
    @printf("BF: %d iters total, %.1f s wall-clock, final obj \$%.2f\n", bf_n, bf_total_time, bf_final)
    if !isnothing(first_feas_idx)
        @printf("  First FEASIBLE: k=%d (%.1f%% of run, ~%.1f s)\n",
                first_feas_idx, 100*first_feas_idx/bf_n, bf_t[first_feas_idx])
    end
    if !isnothing(first_near_idx)
        @printf("  First near-opt: k=%d (%.1f%% of run, ~%.1f s)  ← PLOT TRUNCATES HERE\n",
                first_near_idx, 100*first_near_idx/bf_n, bf_t[first_near_idx])
    end
    @printf("tADMM: %d iters total, plot truncates at k=%d (eff time %.1f s)\n",
            length(obj_tadmm_full), tadmm_end, last(t_tadmm))
    n_trusted = count(tadmm_feas_mask)
    n_untrusted = length(tadmm_feas_mask) - n_trusted
    n_obj_below = count(obj_tadmm .< ref_obj)
    n_r_violate = count(r_tadmm .> eps_pri)
    @printf("  tADMM split (plotted): %d trusted, %d untrusted (%d with r>eps_pri, %d with f<f_BF)\n",
            n_trusted, n_untrusted, n_r_violate, n_obj_below)
    if !isnan(near_t_tadmm) && !isnothing(first_near_idx)
        @printf("  Speedup (tADMM NO vs BF first feas+0.5%%): %.2fx\n",
                bf_t[first_near_idx] / near_t_tadmm)
    end
    println("="^70)
end

# --- entry point -----------------------------------------------------------

if length(ARGS) < 2
    println("Usage: julia plot_bf_vs_tadmm_trajectory.jl <bf_dir> <tadmm_csv> [output_png]")
    exit(1)
end

bf_dir = ARGS[1]
tadmm_csv = ARGS[2]
output = length(ARGS) >= 3 ? ARGS[3] : joinpath(bf_dir, "convergence", "bf_vs_tadmm_trajectory.png")

generate_overlay(bf_dir, tadmm_csv, output)
