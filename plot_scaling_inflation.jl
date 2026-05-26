#!/usr/bin/env julia
# Two-panel scaling inflation figure across all 3 systems.
# Both panels share log-x (T) and log-y (inflation factor) axes.
# Left:  BF wall-clock / (T * BF_1) vs T
# Right: tADMM NO eff time / (T * BF_1) vs T
# Inflation = 1.0 means perfectly linear scaling in T.
# Output: <paper>/figures/scaling_inflation.png

using Plots, Printf

# BF wall-clock time at T=1 (single-period reference) per system.
# Numbers match tab:scaling_inflation in results.tex.
const BF_1 = Dict(
    "ieee123"     => 1.26,
    "medium2552"  => 5.17,
    "large10k"    => 38.57,
)

# ----- ieee123: read live data from rho_sweep_eps1e-4 dirs -----
const IEEE123_T_VALUES = [6, 12, 24, 48, 96, 144]
const IEEE123_DATA_ROOT = joinpath(@__DIR__, "envs", "tadmm", "processedData")
const IEEE123_SYSTEM = "ieee123C_1ph"
const SWEEP_DIR_NAME = "rho_sweep_eps1e-4"

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
        vals = split(lines[2], ",")
        idx = findfirst(==("near_opt_eff_time"), header)
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
    if isfinite(bf) && isfinite(ta)
        push!(ieee123_T, Float64(T))
        push!(ieee123_bf, bf)
        push!(ieee123_ta, ta)
    elseif isfinite(bf)  # BF only (e.g. T=96 tADMM pending)
        push!(ieee123_T, Float64(T))
        push!(ieee123_bf, bf)
        push!(ieee123_ta, NaN)
    end
end

# ----- medium2552 + large10k: hardcoded from tab:scaling_inflation -----
medium_T  = [48.0,    72.0,    96.0,    144.0]
medium_bf = [719.16,  1155.32, 1622.20, 12780.14]
medium_ta = [545.96,  736.46,  738.30,  1372.60]

large_T   = [6.0,    12.0,    24.0,    48.0]
large_bf  = [664.3,  1534.9,  4332.0,  17479.0]
large_ta  = [594.3,  1113.9,  2264.6,  5709.3]

# Compute inflation factors
ieee123_bf_inf = ieee123_bf ./ (ieee123_T .* BF_1["ieee123"])
ieee123_ta_inf = ieee123_ta ./ (ieee123_T .* BF_1["ieee123"])
medium_bf_inf  = medium_bf  ./ (medium_T  .* BF_1["medium2552"])
medium_ta_inf  = medium_ta  ./ (medium_T  .* BF_1["medium2552"])
large_bf_inf   = large_bf   ./ (large_T   .* BF_1["large10k"])
large_ta_inf   = large_ta   ./ (large_T   .* BF_1["large10k"])

# Drop NaNs per-line so missing tADMM points don't break plotting
function strip_nan(xs, ys)
    keep = isfinite.(ys)
    return xs[keep], ys[keep]
end

# Colors picked to match the row-colour vocabulary in the paper tables
const COL_IEEE123 = :tomato
const COL_MEDIUM  = :goldenrod
const COL_LARGE   = :steelblue
const MARK_BF     = :circle
const MARK_TA     = :diamond

# Left panel — BF inflation
p_left = plot(
    xscale = :log10, yscale = :log10,
    xlabel = "T (periods)", ylabel = "BF inflation: BF / (T · BF₁)",
    title  = "Centralised BF scaling",
    titlefontsize = 10,
    legend = :topleft,
    grid = true,
    minorgrid = true,
)
plot!(p_left, ieee123_T, ieee123_bf_inf;
    label = "ieee123",    marker = (MARK_BF, 7, COL_IEEE123), line = (:solid, 2, COL_IEEE123))
plot!(p_left, medium_T, medium_bf_inf;
    label = "medium2552", marker = (MARK_BF, 7, COL_MEDIUM),  line = (:solid, 2, COL_MEDIUM))
plot!(p_left, large_T, large_bf_inf;
    label = "large10k",   marker = (MARK_BF, 7, COL_LARGE),   line = (:solid, 2, COL_LARGE))
hline!(p_left, [1.0]; line = (:dash, 1, :gray), label = "linear (1×)")

# Right panel — tADMM inflation
p_right = plot(
    xscale = :log10, yscale = :log10,
    xlabel = "T (periods)", ylabel = "tADMM inflation: tADMM_NO / (T · BF₁)",
    title  = "Temporal-ADMM scaling",
    titlefontsize = 10,
    legend = :topleft,
    grid = true,
    minorgrid = true,
)
xs, ys = strip_nan(ieee123_T, ieee123_ta_inf)
plot!(p_right, xs, ys;
    label = "ieee123",    marker = (MARK_TA, 7, COL_IEEE123), line = (:solid, 2, COL_IEEE123))
plot!(p_right, medium_T, medium_ta_inf;
    label = "medium2552", marker = (MARK_TA, 7, COL_MEDIUM),  line = (:solid, 2, COL_MEDIUM))
plot!(p_right, large_T, large_ta_inf;
    label = "large10k",   marker = (MARK_TA, 7, COL_LARGE),   line = (:solid, 2, COL_LARGE))
hline!(p_right, [1.0]; line = (:dash, 1, :gray), label = "linear (1×)")

fig = plot(p_left, p_right; layout = (1, 2), size = (1000, 400), dpi = 150,
           left_margin = 6Plots.mm, bottom_margin = 5Plots.mm, top_margin = 3Plots.mm)

const PAPER_FIG_DIR = joinpath(@__DIR__, "..",
    "IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition",
    "figures")
const OUTPUT_PATH = joinpath(PAPER_FIG_DIR, "scaling_inflation.png")

mkpath(PAPER_FIG_DIR)
savefig(fig, OUTPUT_PATH)

# Echo the inflation values for sanity-checking against tab:scaling_inflation
@printf("\n=== ieee123 (BF_1 = %.2f s) ===\n", BF_1["ieee123"])
for i in eachindex(ieee123_T)
    @printf("T=%3d  BF_inf=%.3f  tADMM_inf=%s\n",
        Int(ieee123_T[i]), ieee123_bf_inf[i],
        isnan(ieee123_ta_inf[i]) ? "NaN" : @sprintf("%.3f", ieee123_ta_inf[i]))
end
@printf("\n=== medium2552 (BF_1 = %.2f s) ===\n", BF_1["medium2552"])
for i in eachindex(medium_T)
    @printf("T=%3d  BF_inf=%.3f  tADMM_inf=%.3f\n",
        Int(medium_T[i]), medium_bf_inf[i], medium_ta_inf[i])
end
@printf("\n=== large10k (BF_1 = %.2f s) ===\n", BF_1["large10k"])
for i in eachindex(large_T)
    @printf("T=%3d  BF_inf=%.3f  tADMM_inf=%.3f\n",
        Int(large_T[i]), large_bf_inf[i], large_ta_inf[i])
end
@printf("\nSaved: %s\n", OUTPUT_PATH)
