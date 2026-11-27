#!/bin/bash
# usage, from the repository root: run_storage_coordinator.sh system T   (settings through SC_* and TADMM_* variables)
export JULIA_DEPOT_PATH='C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia_depot'
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TADMM_LINEAR_SOLVER=ma57
export TADMM_IPOPT_OPTIONS='hsllib=C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build/libma57.dll;linear_system_scaling=none'
unset NO_COLOR FORCE_COLOR
J="/c/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.tools/julia-1.12.6/bin/julia.exe"
"$J" --startup-file=no --project=problems/mpopf/tadmm/review_pilot problems/mpopf/tadmm/review_pilot/storage_coordinator.jl "$1" "$2" problems/mpopf/tadmm/review_pilot/results/storage_coordinator_nov15/"$1"_T"$2"
