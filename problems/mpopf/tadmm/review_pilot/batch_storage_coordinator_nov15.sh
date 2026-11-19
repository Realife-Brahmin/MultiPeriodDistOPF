#!/bin/bash
# Sequential verification batch: new scheme and baseline, ieee123 T=24/48/96, then med2522 T=24.
export JULIA_DEPOT_PATH='C:\Users\Aryan Ritwajeet Jha\Documents\documents_general\MultiPeriodDistOPF\.tools\julia_depot'
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TADMM_LINEAR_SOLVER=ma57
export PROFILE_PERIODIC=1 REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0
export TADMM_IPOPT_OPTIONS='hsllib=C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build/libma57.dll;linear_system_scaling=none'
unset NO_COLOR FORCE_COLOR
J="/c/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.tools/julia-1.12.6/bin/julia.exe"
P=problems/mpopf/tadmm/review_pilot; R=$P/results/storage_coordinator_nov15
run(){ "$J" --startup-file=no --project=$P "$@" 2>&1 | grep -E 'CASE|CENTRAL|SC_NO|SC_DONE|TADMM_NO|TADMM_DONE|BASELINE_DONE|ERROR|exported'; }
for T in 48 96; do run problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl ieee123C_1ph $T; done
for T in 24 48 96; do
  SC_TAG=horizon SC_RHO=1000,3000,10000 SC_MAX_ITER=150 run $P/storage_coordinator.jl ieee123C_1ph $T $R/ieee123C_1ph_T$T
  TADMM_MAX_ITER=1000 run $P/baseline_rounds.jl ieee123C_1ph $T $R/ieee123C_1ph_T$T
done
SC_TAG=horizon SC_RHO=3000,10000,30000 SC_MAX_ITER=150 run $P/storage_coordinator.jl ieee2522C_1ph 24 $R/ieee2522C_1ph_T24
echo BATCH_FINISHED
