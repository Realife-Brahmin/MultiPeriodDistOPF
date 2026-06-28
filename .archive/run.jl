#!/usr/bin/env julia
"""
Launcher script that ensures Julia runs with sufficient threads.
Simply run: julia run.jl
It will automatically restart with 16 threads if needed.
"""

using Distributed

const TARGET_THREADS = 16

# Check current thread count
current_threads = Threads.nthreads()

if current_threads < TARGET_THREADS
    println("="^80)
    println("⚠️  Current Julia has $current_threads threads, but we need $TARGET_THREADS")
    println("="^80)
    println("\nRestarting Julia with $TARGET_THREADS threads...")
    println()

    # Relaunch Julia with correct thread count
    cmd = `julia --threads=$TARGET_THREADS $(PROGRAM_FILE)`
    exec(cmd)
else
    println("="^80)
    println("✓ Julia running with $current_threads threads (target: $TARGET_THREADS)")
    println("="^80)
    println()

    # Load and run the main script
    include("tadmm_socp.jl")
end
