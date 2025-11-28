#!/bin/bash
# Run tADMM SOCP with multi-threading enabled
# Uses all available cores (16 threads on your system)

export JULIA_NUM_THREADS=16
echo "Starting Julia with $JULIA_NUM_THREADS threads..."
julia tadmm_socp.jl 2>&1 | tee socp_watchdog_threaded.log
