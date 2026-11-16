#!/usr/bin/env bash
# Centralized conic baseline on Kamiak: Ipopt/MA57 and a conic solver on the same
# node in the same job. usage: kamiak_conic.sh C system T [solvers]
set -euo pipefail
C=$(realpath "$1"); SYS=$2; T=$3; SOLVERS=${4:-ipopt,clarabel}
P=problems/mpopf/tadmm/review_pilot
OUT="$C/results/conic_${SYS}_T${T}_${SLURM_JOB_ID}"
mkdir -p "$OUT"; cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
D=/home/aryan.jha/kamiak_jobs
export JULIA_DEPOT_PATH="$C/julia-depot:$D/tadmm_review_nov11_large10k_interim_ice/julia-depot:$D/tadmm_review_nov08_horizon_sweep/julia-depot:$D/tadmm_review_nov08_watchdog_rho/julia-depot:$D/tadmm_review_nov07_accuracy_fix/julia-depot:$D/tadmm_review_nov07_parallel/julia-depot:$D/tadmm_review_nov07_v2/julia-depot:$D/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR LD_PRELOAD
LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57 TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=july REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0
case $SYS in
 ieee2522C_1ph) export TADMM_SCREEN=substation,ell;;
 *) export TADMM_SCREEN=substation,vupper,ell;;
esac
{ hostname; date -Is; lscpu | grep -E 'Model name|^CPU\(s\)'; taskset -pc $$; } > "$OUT/allocation.txt"
julia --startup-file=no --project="$P/conic_env" "$P/conic_baseline.jl" "$SYS" "$T" "$OUT" "$SOLVERS" > "$OUT/solve.log" 2> "$OUT/solve.err"
grep -E '^CONIC' "$OUT/solve.log"
date -Is > "$OUT/COMPLETE"
