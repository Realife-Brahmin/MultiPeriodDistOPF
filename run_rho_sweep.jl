# ============================================================================
# run_rho_sweep.jl — Automated experiment pipeline
# ============================================================================
# Usage:  julia run_rho_sweep.jl
#   (threads are forwarded to subprocesses automatically)
#
# Phase 1: Rho sweep for T=12 (rho = 500, 1000, 1500)
# Phase 2: T=24 BF solve
# Phase 3: T=24 tADMM solve
# ============================================================================

using Printf

# Load config just to know directory paths
include("config.jl")

# Determine thread count for subprocesses
const N_THREADS = max(Threads.nthreads(), parse(Int, get(ENV, "JULIA_NUM_THREADS", "16")))
const JULIA = Base.julia_cmd()
const SCRIPT_DIR = @__DIR__

println("="^80)
println("EXPERIMENT PIPELINE")
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
# PHASE 1: Rho sweep at T=12
# ============================================================================

const RHO_VALUES = [500.0, 1000.0, 1500.0]
const T12_DIR = joinpath(PROCESSED_DATA_DIR, "$(SYSTEM_NAME)_T12")
const SWEEP_DIR = joinpath(T12_DIR, "rho_sweep")

println("\n" * "="^80)
println("PHASE 1: RHO SWEEP (T=12)")
println("  Values: ", join([@sprintf("%.0f", r) for r in RHO_VALUES], ", "))
println("  Estimated time: ~$(length(RHO_VALUES) * 35) minutes")
println("="^80)

for (i, rho_val) in enumerate(RHO_VALUES)
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    rho_dir = joinpath(SWEEP_DIR, rho_label)

    println("\n>>> [$i/$(length(RHO_VALUES))] rho = $rho_val")

    t_start = time()
    run_julia_script("run_tadmm.jl"; env_overrides=Dict(
        "T_OVERRIDE" => "12",
        "RHO_OVERRIDE" => string(rho_val),
    ))
    elapsed = time() - t_start

    copy_results(T12_DIR, rho_dir)
    @printf(">>> rho=%d done in %.1f minutes. Results saved to %s\n",
            round(Int, rho_val), elapsed / 60, rho_dir)
end

# Also generate plots for each rho (reads from the copied convergence_data.csv)
println("\n>>> Generating plots for rho sweep...")
for rho_val in RHO_VALUES
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    rho_dir = joinpath(SWEEP_DIR, rho_label)
    # Plot generation would need convergence_data.csv in SYSTEM_DIR,
    # but we don't want to overwrite. Plots can be generated manually.
end

# ============================================================================
# PHASE 2: T=24 Brute Force
# ============================================================================

println("\n" * "="^80)
println("PHASE 2: BRUTE FORCE (T=24)")
println("  This may take 30-90+ minutes for large10kC")
println("="^80)

t_start = time()
run_julia_script("run_bf.jl"; env_overrides=Dict(
    "T_OVERRIDE" => "24",
))
elapsed = time() - t_start
@printf(">>> BF T=24 done in %.1f minutes\n", elapsed / 60)

# ============================================================================
# PHASE 3: T=24 tADMM
# ============================================================================

println("\n" * "="^80)
println("PHASE 3: tADMM (T=24)")
println("  Using default rho scaling: 3000 * sqrt(24/24) = 3000")
println("="^80)

t_start = time()
run_julia_script("run_tadmm.jl"; env_overrides=Dict(
    "T_OVERRIDE" => "24",
))
elapsed = time() - t_start
@printf(">>> tADMM T=24 done in %.1f minutes\n", elapsed / 60)

# ============================================================================
# PHASE 4: T=24 Plots
# ============================================================================

println("\n" * "="^80)
println("PHASE 4: GENERATING T=24 PLOTS")
println("="^80)

run_julia_script("plot_results.jl"; env_overrides=Dict(
    "T_OVERRIDE" => "24",
))

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^80)
println("ALL EXPERIMENTS COMPLETE")
println("="^80)
println("Results locations:")
for rho_val in RHO_VALUES
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    println("  T=12 rho=$rho_label: $(joinpath(SWEEP_DIR, rho_label))")
end
println("  T=24 BF+tADMM:    $(joinpath(PROCESSED_DATA_DIR, "$(SYSTEM_NAME)_T24"))")
println("="^80)
println("Done.")
