# ============================================================================
# run_rho_sweep.jl — T=24 Rho sweep (vanilla ADMM, no adaptive rho)
# ============================================================================
# Usage:  julia run_rho_sweep.jl
#   (threads are forwarded to subprocesses automatically)
#
# Sweeps rho = 500, 1000, 1500, 2000 at T=24 with:
#   - adaptive_rho = false (fixed rho throughout)
#   - stagnation_window = 10 (hardcoded in run_tadmm.jl)
#   - stagnation_threshold = 1e-3
# BF results already exist in large10kC_1ph_T24/
# ============================================================================

using Printf

# Load config just to know directory paths
include("config.jl")

# Determine thread count for subprocesses
const N_THREADS = max(Threads.nthreads(), parse(Int, get(ENV, "JULIA_NUM_THREADS", "16")))
const JULIA = Base.julia_cmd()
const SCRIPT_DIR = @__DIR__

println("="^80)
println("T=24 RHO SWEEP — VANILLA ADMM (no adaptive rho)")
println("="^80)
println("Julia:   $JULIA")
println("Threads: $N_THREADS")
println("System:  $SYSTEM_NAME")
println("="^80)

function run_julia_script(script; env_overrides=Dict{String,String}())
    env = copy(ENV)
    for (k, v) in env_overrides
        env[k] = v
    end
    cmd = setenv(`$JULIA --threads=$N_THREADS $(joinpath(SCRIPT_DIR, script))`, env)
    println("\n>>> Running: $script with overrides: $env_overrides")
    run(cmd)
end

function copy_results(src_dir, dst_dir; files=nothing)
    mkpath(dst_dir)
    result_files = isnothing(files) ? [
        "results_socp_tadmm.txt", "convergence_data.csv",
        "near_optimal_summary.csv", "subproblem_timing_details.csv",
        "tadmm_run.log", "sol_socp_tadmm.jls",
    ] : files
    for fname in result_files
        src = joinpath(src_dir, fname)
        if isfile(src)
            cp(src, joinpath(dst_dir, fname); force=true)
        end
    end
    # Copy convergence plots subdir
    conv_src = joinpath(src_dir, "convergence")
    if isdir(conv_src)
        conv_dst = joinpath(dst_dir, "convergence")
        mkpath(conv_dst)
        for f in readdir(conv_src)
            cp(joinpath(conv_src, f), joinpath(conv_dst, f); force=true)
        end
    end
end

# ============================================================================
# RHO SWEEP: T=24, vanilla ADMM
# ============================================================================

const RHO_VALUES = [25000.0, 10000.0, 12000.0]
const T24_DIR = joinpath(PROCESSED_DATA_DIR, "$(SYSTEM_NAME)_T24")
const SWEEP_DIR = joinpath(T24_DIR, "rho_sweep")

println("\n" * "="^80)
println("RHO SWEEP (T=24, vanilla ADMM)")
println("  Values: ", join([@sprintf("%.0f", r) for r in RHO_VALUES], ", "))
println("  adaptive_rho = false")
println("  stagnation_window = 10")
println("  stagnation_threshold = 1e-3")
println("  Estimated time: ~$(length(RHO_VALUES) * 60) minutes")
println("="^80)

for (i, rho_val) in enumerate(RHO_VALUES)
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    rho_dir = joinpath(SWEEP_DIR, rho_label)

    println("\n>>> [$i/$(length(RHO_VALUES))] rho = $rho_val")

    t_start = time()
    run_julia_script("run_tadmm.jl"; env_overrides=Dict(
        "T_OVERRIDE" => "24",
        "RHO_OVERRIDE" => string(rho_val),
        "ADAPTIVE_RHO_OVERRIDE" => "false",
    ))
    elapsed = time() - t_start

    copy_results(T24_DIR, rho_dir)
    @printf(">>> rho=%d done in %.1f minutes. Results saved to %s\n",
            round(Int, rho_val), elapsed / 60, rho_dir)
end

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^80)
println("T=24 RHO SWEEP COMPLETE")
println("="^80)
println("Results locations:")
for rho_val in RHO_VALUES
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    println("  rho=$rho_label: $(joinpath(SWEEP_DIR, rho_label))")
end
println("  BF reference:  $(joinpath(T24_DIR, "results_socp_bf.txt"))")
println("="^80)
println("Done.")
