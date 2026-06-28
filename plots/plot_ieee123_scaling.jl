#!/usr/bin/env julia
# Two-panel scaling figure for ieee123 across T at eps_pri=1e-4.
# Left:  log-y BF wall-clock + tADMM NO eff time vs T
# Right: linear-y tADMM/BF wall-clock multiplier vs T
# Output: <paper>/figures/ieee123_scaling.png

using Plots, Printf

const SYSTEM = "ieee123C_1ph"
const T_VALUES = [6, 12, 24, 48, 96, 144]
const DATA_ROOT = joinpath(@__DIR__, "envs", "tadmm", "processedData")
const SWEEP_DIR_NAME = "rho_sweep_eps1e-4"
const PAPER_FIG_DIR = joinpath(@__DIR__, "..",
    "IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition",
    "figures")
const OUTPUT_PATH = joinpath(PAPER_FIG_DIR, "ieee123_scaling.png")

function parse_bf_walltime(bf_file::AbstractString)
    isfile(bf_file) || return NaN
    for line in eachline(bf_file)
        m = match(r"Wall-clock time:\s+([0-9.]+)\s+seconds", line)
        isnothing(m) || return parse(Float64, m.captures[1])
    end
    return NaN
end

# Read each rho_*/near_optimal_summary.csv under sweep_dir, return the row
# with the smallest near_opt_eff_time as the (rho, NO eff time, gap) winner.
function pick_winner(sweep_dir::AbstractString)
    isdir(sweep_dir) || return (NaN, NaN, NaN)
    best_rho, best_time, best_gap = NaN, Inf, NaN
    for entry in readdir(sweep_dir; join=true)
        isdir(entry) || continue
        startswith(basename(entry), "rho_") || continue
        rho = parse(Float64, replace(basename(entry), "rho_" => ""))
        csv = joinpath(entry, "near_optimal_summary.csv")
        isfile(csv) || continue
        lines = readlines(csv)
        length(lines) >= 2 || continue
        header = split(lines[1], ",")
        vals = split(lines[2], ",")
        idx_time = findfirst(==("near_opt_eff_time"), header)
        idx_gap  = findfirst(==("near_opt_gap_pct"), header)
        isnothing(idx_time) && continue
        t = try parse(Float64, vals[idx_time]); catch; NaN; end
        g = isnothing(idx_gap) ? NaN : (try parse(Float64, vals[idx_gap]); catch; NaN; end)
        if !isnan(t) && t < best_time
            best_time = t; best_rho = rho; best_gap = g
        end
    end
    return (best_rho, isfinite(best_time) ? best_time : NaN, best_gap)
end

# Gather (T, BF, tADMM_NO, rho_win, gap)
points = NamedTuple{(:T, :bf, :tadmm, :rho, :gap), NTuple{5, Float64}}[]
for T in T_VALUES
    sysdir   = joinpath(DATA_ROOT, "$(SYSTEM)_T$(T)")
    bf       = parse_bf_walltime(joinpath(sysdir, "results_socp_bf.txt"))
    sweep    = joinpath(sysdir, SWEEP_DIR_NAME)
    rho, no_time, gap = pick_winner(sweep)
    push!(points, (T=Float64(T), bf=bf, tadmm=no_time, rho=rho, gap=gap))
    @printf("T=%3d  BF=%7.2fs  tADMM_NO=%8.2fs  rho_win=%g  gap=%.2f%%\n",
            T, bf, no_time, rho, gap)
end

# Filter complete points for plotting
complete = filter(p -> isfinite(p.bf) && isfinite(p.tadmm), points)
isempty(complete) && error("No complete (BF + tADMM) points found.")

Ts        = [p.T     for p in complete]
bf_times  = [p.bf    for p in complete]
no_times  = [p.tadmm for p in complete]
mults     = no_times ./ bf_times

# Left panel: log-y times
p_left = plot(
    Ts, bf_times;
    label = "BF wall-clock",
    marker = (:circle, 7, :forestgreen),
    line = (:solid, 2, :forestgreen),
    yscale = :log10,
    xlabel = "T (periods)",
    ylabel = "Time (s)",
    legend = :topleft,
    grid = true,
    title = "BF vs tADMM (NO) time",
    titlefontsize = 10,
)
plot!(p_left, Ts, no_times;
    label = "tADMM near-optimal",
    marker = (:square, 7, :dodgerblue),
    line = (:solid, 2, :dodgerblue),
)

# Right panel: multiplier
p_right = plot(
    Ts, mults;
    label = "tADMM / BF",
    marker = (:diamond, 7, :firebrick),
    line = (:solid, 2, :firebrick),
    xlabel = "T (periods)",
    ylabel = "Multiplier (×)",
    legend = :topleft,
    grid = true,
    title = "tADMM overhead vs BF",
    titlefontsize = 10,
)
hline!(p_right, [1.0]; line = (:dash, 1, :gray), label = "parity")

fig = plot(p_left, p_right; layout = (1, 2), size = (900, 360), dpi = 150,
           left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)

mkpath(PAPER_FIG_DIR)
savefig(fig, OUTPUT_PATH)
@printf("\nSaved: %s\n", OUTPUT_PATH)
@printf("Plotted %d / %d T values (%d missing data).\n",
        length(complete), length(T_VALUES), length(T_VALUES) - length(complete))
