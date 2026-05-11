#!/usr/bin/env julia
# build_sweep_summary.jl --- aggregate per-rho results into a comparison table
#
# Usage:
#   julia build_sweep_summary.jl <system_T_dirname>
#   e.g.  julia build_sweep_summary.jl ieee123C_1ph_T96
#
# Reads:  envs/tadmm/processedData/<dir>/results_socp_bf.txt
#         envs/tadmm/processedData/<dir>/rho_sweep/rho_*/results_socp_tadmm.txt
# Writes: envs/tadmm/processedData/<dir>/rho_sweep/sweep_summary.csv
# Prints: human-readable comparison table to stdout

using Printf

const PROCESSED_DATA_DIR = joinpath(@__DIR__, "envs", "tadmm", "processedData")

if length(ARGS) < 1
    println("Usage: julia build_sweep_summary.jl <system_T_dirname>")
    println("  e.g.  julia build_sweep_summary.jl ieee123C_1ph_T96")
    exit(1)
end

const TARGET_DIR = ARGS[1]
const T_DIR = joinpath(PROCESSED_DATA_DIR, TARGET_DIR)
const SWEEP_DIR = joinpath(T_DIR, "rho_sweep")

isdir(SWEEP_DIR) || error("Sweep directory not found: $SWEEP_DIR")

# --- helpers ------------------------------------------------------------------

function extract_field(text, pattern; default=nothing, T=Float64)
    m = match(pattern, text)
    isnothing(m) && return default
    try
        return T == String ? m.captures[1] : parse(T, m.captures[1])
    catch
        return default
    end
end

function parse_bf(bf_path)
    isfile(bf_path) || return nothing
    txt = read(bf_path, String)
    return (
        objective = extract_field(txt, r"Total Cost:\s*\$([0-9.eE+-]+)"; default=NaN),
        wall_time = extract_field(txt, r"Wall-clock time:\s*([0-9.]+)\s*seconds"; default=NaN),
        solver_time = extract_field(txt, r"Solver time:\s*([0-9.]+)\s*seconds"; default=NaN),
    )
end

function parse_tadmm(rho_dir)
    res_path = joinpath(rho_dir, "results_socp_tadmm.txt")
    isfile(res_path) || return nothing
    txt = read(res_path, String)

    # Near-optimal line is special; parse it carefully
    near_opt_line = match(r"Effective time \(near-optimal\):\s*([^\n]+)", txt)
    near_opt_str = isnothing(near_opt_line) ? "" : near_opt_line.captures[1]
    near_opt_time   = extract_field(near_opt_str, r"^\s*([0-9.]+)\s*seconds"; default=NaN)
    near_opt_k      = extract_field(near_opt_str, r"\(k=(\d+)"; default=-1, T=Int)
    near_opt_gap    = extract_field(near_opt_str, r"gap=([0-9.eE+-]+)%"; default=NaN)
    near_opt_r      = extract_field(near_opt_str, r"r=([0-9.eE+-]+)"; default=NaN)
    near_opt_s      = extract_field(near_opt_str, r"s=([0-9.eE+-]+)"; default=NaN)

    return (
        rho_init   = extract_field(txt, r"Initial rho:\s*([0-9.]+)"; default=NaN),
        adaptive   = extract_field(txt, r"Adaptive rho:\s*(true|false)"; default="?", T=String),
        max_iter   = extract_field(txt, r"Max iterations:\s*(\d+)"; default=0, T=Int),
        eps_pri    = extract_field(txt, r"Primal tolerance:\s*([0-9.eE+-]+)"; default=NaN),
        eps_dual   = extract_field(txt, r"Dual tolerance:\s*([0-9.eE+-]+)"; default=NaN),
        iters      = extract_field(txt, r"Iterations:\s*(\d+)"; default=0, T=Int),
        final_r    = extract_field(txt, r"Final primal residual:\s*([0-9.eE+-]+)"; default=NaN),
        final_s    = extract_field(txt, r"Final dual residual:\s*([0-9.eE+-]+)"; default=NaN),
        objective  = extract_field(txt, r"Total Cost:\s*\$([0-9.eE+-]+)"; default=NaN),
        wall_time  = extract_field(txt, r"Wall-clock time:\s*([0-9.]+)\s*seconds"; default=NaN),
        true_eff   = extract_field(txt, r"Effective time \(true\):\s*([0-9.]+)\s*seconds"; default=NaN),
        near_eff   = near_opt_time,
        near_k     = near_opt_k,
        near_gap   = near_opt_gap,
        near_r     = near_opt_r,
        near_s     = near_opt_s,
        spd_true   = extract_field(txt, r"tADMM eff\(true\) vs BF:\s*([0-9.]+)x"; default=NaN),
        spd_near   = extract_field(txt, r"tADMM eff\(near-optimal\) vs BF:\s*([0-9.]+)x"; default=NaN),
    )
end

# --- main ---------------------------------------------------------------------

bf = parse_bf(joinpath(T_DIR, "results_socp_bf.txt"))

println("="^120)
println("RHO SWEEP SUMMARY: $TARGET_DIR")
println("="^120)
if !isnothing(bf)
    @printf("BF baseline: obj=\$%.4f, wall=%.2fs, solver=%.2fs\n", bf.objective, bf.wall_time, bf.solver_time)
else
    println("BF baseline: NOT FOUND (results_socp_bf.txt missing)")
end
println("-"^120)

# Discover rho_* dirs and sort by rho_init numerically
rho_dirs = filter(d -> startswith(d, "rho_") && isdir(joinpath(SWEEP_DIR, d)),
                  readdir(SWEEP_DIR))
sort!(rho_dirs, by = d -> begin
    s = replace(d, "rho_" => "")
    try parse(Float64, s) catch; Inf end
end)

rows = NamedTuple[]
for d in rho_dirs
    rho_dir = joinpath(SWEEP_DIR, d)
    r = parse_tadmm(rho_dir)
    if isnothing(r)
        @printf("  [%s] no results file — skipped\n", d)
        continue
    end
    push!(rows, merge((dirname=d,), r))
end

if isempty(rows)
    println("No tADMM results found in any rho_* directory.")
    exit(0)
end

# Print table
println()
@printf("%-12s %-9s %-7s %-10s %-10s %-12s %-7s %-9s %-12s %-10s %-10s %-10s\n",
        "rho_init", "adaptive", "iters", "final_r", "final_s", "obj", "near_k", "near_gap%", "near_eff_s", "true_eff_s", "spd_near", "spd_true")
println("-"^120)
for r in rows
    converged_str = (r.final_r <= r.eps_pri && r.final_s <= r.eps_dual) ? "Y" : "N"
    @printf("%-12.0f %-9s %-7d %-10.2e %-10.2e %-12.2f %-7d %-9s %-12s %-10.2f %-10s %-10.3fx\n",
        r.rho_init,
        r.adaptive,
        r.iters,
        r.final_r,
        r.final_s,
        r.objective,
        r.near_k,
        isnan(r.near_gap)  ? "—" : @sprintf("%.3f", r.near_gap),
        isnan(r.near_eff)  ? "—" : @sprintf("%.2f", r.near_eff),
        r.true_eff,
        isnan(r.spd_near)  ? "—" : @sprintf("%.3fx", r.spd_near),
        r.spd_true,
    )
end
println("-"^120)

# Identify winners
converged = filter(r -> r.final_r <= r.eps_pri && r.final_s <= r.eps_dual, rows)
near_opt_valid = filter(r -> !isnan(r.near_eff) && r.near_k > 0, rows)

println("\n--- WINNERS ---")
if isempty(converged)
    println("True convergence: NONE (no run satisfied both r<=eps_pri and s<=eps_dual)")
else
    fastest_true = argmin(r -> r.true_eff, converged)
    @printf("Fastest true convergence: rho=%.0f, k=%d, eff=%.2fs (%.3fx vs BF)\n",
            fastest_true.rho_init, fastest_true.iters, fastest_true.true_eff, fastest_true.spd_true)
end
if isempty(near_opt_valid)
    println("Near-optimal: NONE")
else
    fastest_near = argmin(r -> r.near_eff, near_opt_valid)
    @printf("Best near-optimal:        rho=%.0f, k=%d, gap=%.3f%%, eff=%.2fs (%.3fx vs BF)\n",
            fastest_near.rho_init, fastest_near.near_k, fastest_near.near_gap,
            fastest_near.near_eff, fastest_near.spd_near)
end
println("="^120)

# Write CSV
csv_path = joinpath(SWEEP_DIR, "sweep_summary.csv")
open(csv_path, "w") do io
    println(io, "rho_init,adaptive,iters,final_r,final_s,objective,near_k,near_gap_pct,near_eff_s,true_eff_s,spd_near,spd_true,converged")
    for r in rows
        converged = (r.final_r <= r.eps_pri && r.final_s <= r.eps_dual) ? "Y" : "N"
        @printf(io, "%.0f,%s,%d,%.6e,%.6e,%.4f,%d,%.6f,%.4f,%.4f,%.6f,%.6f,%s\n",
            r.rho_init, r.adaptive, r.iters, r.final_r, r.final_s, r.objective,
            r.near_k, r.near_gap, r.near_eff, r.true_eff, r.spd_near, r.spd_true, converged)
    end
end
@printf("CSV summary saved: %s\n", csv_path)
