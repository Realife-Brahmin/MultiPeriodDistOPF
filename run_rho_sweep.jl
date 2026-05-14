# ============================================================================
# run_rho_sweep.jl — Rho sweep (adaptive rho, varying initial rho)
# ============================================================================
# Usage:  julia run_rho_sweep.jl
#   (threads are forwarded to subprocesses automatically)
#
# Sweeps starting rho values with:
#   - adaptive_rho = true (rho adjusts during solve)
#   - stagnation_window = 10 (hardcoded in run_tadmm.jl)
#   - stagnation_threshold = 1e-3
# ============================================================================

using Printf, Dates

# Load config just to know directory paths
include("config.jl")

# Determine thread count for subprocesses
const N_THREADS = max(Threads.nthreads(), parse(Int, get(ENV, "JULIA_NUM_THREADS", "16")))
const JULIA = Base.julia_cmd()
const SCRIPT_DIR = @__DIR__

const SWEEP_T = haskey(ENV, "SWEEP_T") ? parse(Int, ENV["SWEEP_T"]) : 24

println("="^80)
println("T=$(SWEEP_T) RHO SWEEP — adaptive rho, varying rho_init")
println("="^80)
println("Julia:   $JULIA")
println("Threads: $N_THREADS")
println("System:  $SYSTEM_NAME")
println("="^80)

function run_julia_script(script; env_overrides=Dict{String,String}(), stdout_log=nothing)
    env = copy(ENV)
    for (k, v) in env_overrides
        env[k] = v
    end
    cmd = setenv(`$JULIA --threads=$N_THREADS $(joinpath(SCRIPT_DIR, script))`, env)
    println("\n>>> Running: $script with overrides: $env_overrides")
    if isnothing(stdout_log)
        run(cmd)
    else
        mkpath(dirname(stdout_log))
        # Tee: stream subprocess stdout/stderr to both console AND log file
        open(stdout_log, "w") do io
            process = open(pipeline(cmd, stderr=Base.stderr), "r")
            for line in eachline(process)
                println(line)
                println(io, line)
                flush(io)
            end
            wait(process)
        end
    end
end

function write_params_file(rho_dir, rho_val)
    params_path = joinpath(rho_dir, "params.txt")
    open(params_path, "w") do io
        @printf(io, "tADMM Hyperparameters -- T=%s, rho=%.0f\n", get(ENV, "SWEEP_T", "?"), rho_val)
        @printf(io, "%s\n", "="^60)
        @printf(io, "system_override   = %s\n", get(ENV, "SYSTEM_OVERRIDE", "(default)"))
        @printf(io, "T                 = %s\n", get(ENV, "SWEEP_T", "?"))
        @printf(io, "rho_init          = %.1f\n", rho_val)
        @printf(io, "adaptive_rho      = %s\n", get(ENV, "ADAPTIVE_RHO_OVERRIDE", "(default true)"))
        @printf(io, "tau_incr          = %s\n", get(ENV, "TAU_INCR_OVERRIDE", "(default 2.0)"))
        @printf(io, "tau_decr          = %s\n", get(ENV, "TAU_DECR_OVERRIDE", "(default 5.0)"))
        @printf(io, "mu_balance        = %s\n", get(ENV, "MU_BALANCE_OVERRIDE", "(default 5.0)"))
        @printf(io, "update_interval   = %s\n", get(ENV, "UPDATE_INTERVAL_OVERRIDE", "(default 5)"))
        @printf(io, "eps_pri           = %s\n", get(ENV, "EPS_PRI_OVERRIDE", "(default 2e-3)"))
        @printf(io, "eps_dual          = %s\n", get(ENV, "EPS_DUAL_OVERRIDE", "(default 1e-2)"))
        @printf(io, "max_iter          = 100 (hardcoded)\n")
        @printf(io, "watchdog          = true (window=15, factor=2.5, trigger=r>2*eps_pri)\n")
        @printf(io, "two_phase         = true (Phase 1 unconditional increase; Phase 2 mu_balance-driven)\n")
    end
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
# RHO SWEEP: T=24, adaptive rho
# ============================================================================

const RHO_VALUES = haskey(ENV, "SWEEP_RHOS") ? parse.(Float64, split(ENV["SWEEP_RHOS"], ",")) : [25000.0, 10000.0, 12000.0]
const SWEEP_TDIR = joinpath(PROCESSED_DATA_DIR, "$(SYSTEM_NAME)_T$(SWEEP_T)")
const SWEEP_DIR = joinpath(SWEEP_TDIR, "rho_sweep")

println("\n" * "="^80)
println("RHO SWEEP (T=$(SWEEP_T), adaptive rho)")
println("  Values: ", join([@sprintf("%.0f", r) for r in RHO_VALUES], ", "))
println("  adaptive_rho = true")
println("  stagnation_window = 10")
println("  stagnation_threshold = 1e-3")
println("  Estimated time: ~$(length(RHO_VALUES) * 60) minutes")
println("="^80)

for (i, rho_val) in enumerate(RHO_VALUES)
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    rho_dir = joinpath(SWEEP_DIR, rho_label)

    println("\n>>> [$i/$(length(RHO_VALUES))] rho = $rho_val")

    t_start = time()
    mkpath(rho_dir)
    stdout_log = joinpath(rho_dir, "run_stdout.log")
    crashed = false
    try
        run_julia_script("run_tadmm.jl"; env_overrides=Dict(
            "T_OVERRIDE" => string(SWEEP_T),
            "RHO_OVERRIDE" => string(rho_val),
            "ADAPTIVE_RHO_OVERRIDE" => get(ENV, "ADAPTIVE_RHO_OVERRIDE", "true"),
        ), stdout_log=stdout_log)
    catch e
        crashed = true
        @printf(">>> rho=%d CRASHED: %s\n", round(Int, rho_val), sprint(showerror, e))
        # Mark the failure in the rho_dir so we know later
        open(joinpath(rho_dir, "CRASH.txt"), "w") do io
            println(io, "Subprocess crashed at $(now()).")
            println(io, "Error: ", sprint(showerror, e))
            println(io, "Partial files (whatever made it to disk before crash) copied below.")
        end
    end
    elapsed = time() - t_start

    # Always copy whatever made it to disk -- even on crash, the latest
    # tadmm_run.log etc. in the shared system_T dir is THIS rho's data,
    # and the next rho would overwrite it. Snapshot it now.
    try
        copy_results(SWEEP_TDIR, rho_dir)
        write_params_file(rho_dir, rho_val)
    catch e2
        @printf(">>> rho=%d copy_results failed: %s\n", round(Int, rho_val), sprint(showerror, e2))
    end

    status = crashed ? "CRASHED" : "done"
    @printf(">>> rho=%d %s in %.1f minutes. Results saved to %s\n",
            round(Int, rho_val), status, elapsed / 60, rho_dir)
end

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^80)
println("T=$(SWEEP_T) RHO SWEEP COMPLETE")
println("="^80)
println("Results locations:")
for rho_val in RHO_VALUES
    rho_label = @sprintf("rho_%d", round(Int, rho_val))
    println("  rho=$rho_label: $(joinpath(SWEEP_DIR, rho_label))")
end
println("  BF reference:  $(joinpath(SWEEP_TDIR, "results_socp_bf.txt"))")
println("="^80)
println("Done.")
