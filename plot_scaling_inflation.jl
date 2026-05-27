#!/usr/bin/env julia
# Three-panel scaling inflation figure (one panel per system).
# Within each panel: grouped bar chart with BF and tADMM bars per T.
# y = inflation factor = time / (T * single-period BF time)
# y = 1 line marks perfectly linear scaling in T.
# Output: <paper>/figures/scaling_inflation.png

using Plots, Printf

# Use serif fonts to match the paper aesthetic.
Plots.default(fontfamily = "Computer Modern")

# Single-period BF wall-clock per system (T=1 reference).
# Numbers match tab:scaling_inflation in results.tex.
const BF_1 = Dict(
    "ieee123"     => 1.26,
    "medium2552"  => 5.17,
    "large10k"    => 38.57,
)

# Bar palette: darker = BF, lighter = tADMM.
# Hues match the row-colour vocabulary used in the paper tables:
#   ieee123    = soft blue   (rowieee123)     -> steelblue family
#   medium2552 = soft cream  (rowmedium2552)  -> goldenrod family
#   large10k   = soft mint   (rowlarge10k)    -> forestgreen family
const PALETTE = Dict(
    "ieee123"    => (bf = RGB(31/255,  81/255, 134/255),  ta = RGB(120/255, 180/255, 220/255)),
    "medium2552" => (bf = RGB(184/255, 134/255,  11/255), ta = RGB(238/255, 203/255, 110/255)),
    "large10k"   => (bf = RGB(34/255,  120/255,  34/255), ta = RGB(120/255, 200/255, 130/255)),
)

# ----- ieee123: read live data from rho_sweep_eps1e-4 dirs -----
const IEEE123_T_VALUES   = [6, 12, 24, 48, 96, 144]
const IEEE123_DATA_ROOT  = joinpath(@__DIR__, "envs", "tadmm", "processedData")
const IEEE123_SYSTEM     = "ieee123C_1ph"
const SWEEP_DIR_NAME     = "rho_sweep_eps1e-4"

function parse_bf_walltime(bf_file::AbstractString)
    isfile(bf_file) || return NaN
    for line in eachline(bf_file)
        m = match(r"Wall-clock time:\s+([0-9.]+)\s+seconds", line)
        isnothing(m) || return parse(Float64, m.captures[1])
    end
    return NaN
end

function pick_winner_no_time(sweep_dir::AbstractString)
    isdir(sweep_dir) || return NaN
    best = Inf
    for entry in readdir(sweep_dir; join=true)
        isdir(entry) && startswith(basename(entry), "rho_") || continue
        csv = joinpath(entry, "near_optimal_summary.csv")
        isfile(csv) || continue
        lines = readlines(csv)
        length(lines) >= 2 || continue
        header = split(lines[1], ",")
        vals   = split(lines[2], ",")
        idx    = findfirst(==("near_opt_eff_time"), header)
        isnothing(idx) && continue
        t = try parse(Float64, vals[idx]); catch; NaN; end
        isfinite(t) && t < best && (best = t)
    end
    return isfinite(best) ? best : NaN
end

ieee123_T  = Float64[]
ieee123_bf = Float64[]
ieee123_ta = Float64[]
for T in IEEE123_T_VALUES
    sysdir = joinpath(IEEE123_DATA_ROOT, "$(IEEE123_SYSTEM)_T$(T)")
    bf     = parse_bf_walltime(joinpath(sysdir, "results_socp_bf.txt"))
    ta     = pick_winner_no_time(joinpath(sysdir, SWEEP_DIR_NAME))
    push!(ieee123_T,  Float64(T))
    push!(ieee123_bf, bf)
    push!(ieee123_ta, ta)
end

# ----- medium2552 + large10k: hardcoded from tab:scaling_inflation -----
# medium2552 T=24 BF uses the worst-case of 6 reps (296.32 s) per Table I footnote.
# T=24 tADMM = rho_14000 / tau_decr=5 winner (251.70 s, gap 0.45%).
medium_T  = [6.0,    12.0,    24.0,    48.0,    72.0,    96.0,    144.0]
medium_bf = [37.62,  127.01,  296.32,  719.16,  1155.32, 1622.20, 12780.14]
medium_ta = [29.71,  94.12,   251.70,  545.96,  736.46,  738.30,  1372.60]

large_T   = [6.0,    12.0,    24.0,    48.0]
large_bf  = [664.3,  1534.9,  4332.0,  17479.0]
large_ta  = [594.3,  1113.9,  2264.6,  5709.3]

# Build one grouped-bar panel for a given system.
function panel(Ts, bf_secs, ta_secs, sysname, palette;
               panel_title, show_ylabel = false, show_legend = false)
    n     = length(Ts)
    xs    = collect(1.0:Float64(n))
    width = 0.38
    bf_inf = bf_secs ./ (Ts .* BF_1[sysname])
    ta_inf = ta_secs ./ (Ts .* BF_1[sysname])

    # Drop any NaN tADMM entries (e.g. ieee123 T=96 pending sweep)
    bf_xs, bf_ys = filter_finite(xs .- width/2, bf_inf)
    ta_xs, ta_ys = filter_finite(xs .+ width/2, ta_inf)

    p = bar(bf_xs, bf_ys;
        bar_width = width, color = palette.bf, linecolor = :black, linewidth = 0.4,
        label = show_legend ? "BF" : "",
        legend = show_legend ? :topleft : false,
        xticks = (xs, string.(Int.(Ts))),
        xlabel = "T",
        ylabel = show_ylabel ? "Inflation factor" : "",
        title  = panel_title,
        titlefontsize = 11,
        grid = true,
        ylim = (0, max(maximum(skipnan(bf_inf)), maximum(skipnan(ta_inf))) * 1.15),
    )
    bar!(p, ta_xs, ta_ys;
        bar_width = width, color = palette.ta, linecolor = :black, linewidth = 0.4,
        label = show_legend ? "tADMM" : "",
    )
    hline!(p, [1.0]; line = (:dash, 1, :gray40), label = show_legend ? "linear" : "")

    # Annotate the numeric value above each bar (small text).
    for (x, y) in zip(bf_xs, bf_ys)
        annotate!(p, x, y * 1.04, text(@sprintf("%.2f", y), 7, :center, :bottom))
    end
    for (x, y) in zip(ta_xs, ta_ys)
        annotate!(p, x, y * 1.04, text(@sprintf("%.2f", y), 7, :center, :bottom))
    end
    return p
end

function filter_finite(xs, ys)
    keep = isfinite.(ys)
    return xs[keep], ys[keep]
end

skipnan(v) = filter(isfinite, v)

p_ieee  = panel(ieee123_T, ieee123_bf, ieee123_ta, "ieee123",    PALETTE["ieee123"];
                panel_title = "ieee123",    show_ylabel = true,  show_legend = true)
p_med   = panel(medium_T,  medium_bf,  medium_ta,  "medium2552", PALETTE["medium2552"];
                panel_title = "medium2552", show_ylabel = false, show_legend = false)
p_large = panel(large_T,   large_bf,   large_ta,   "large10k",   PALETTE["large10k"];
                panel_title = "large10k",   show_ylabel = false, show_legend = false)

fig = plot(p_ieee, p_med, p_large;
           layout = (1, 3), size = (1200, 380), dpi = 150,
           left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
           top_margin = 3Plots.mm, right_margin = 2Plots.mm)

const PAPER_FIG_DIR = joinpath(@__DIR__, "..",
    "IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition",
    "figures")
const OUTPUT_PATH = joinpath(PAPER_FIG_DIR, "scaling_inflation.png")

mkpath(PAPER_FIG_DIR)
savefig(fig, OUTPUT_PATH)

# Echo numeric values for cross-checking against tab:scaling_inflation.
function dump_panel(Ts, bf_secs, ta_secs, sysname)
    @printf("\n=== %s (BF1 = %.2f s) ===\n", sysname, BF_1[sysname])
    for i in eachindex(Ts)
        bfi = bf_secs[i] / (Ts[i] * BF_1[sysname])
        tai = ta_secs[i] / (Ts[i] * BF_1[sysname])
        @printf("T=%3d  BF_inf=%.3f  tADMM_inf=%s\n",
            Int(Ts[i]), bfi, isnan(tai) ? "NaN" : @sprintf("%.3f", tai))
    end
end
dump_panel(ieee123_T, ieee123_bf, ieee123_ta, "ieee123")
dump_panel(medium_T,  medium_bf,  medium_ta,  "medium2552")
dump_panel(large_T,   large_bf,   large_ta,   "large10k")
@printf("\nSaved: %s\n", OUTPUT_PATH)
