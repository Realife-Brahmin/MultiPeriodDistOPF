#!/usr/bin/env bash
# Kamiak: gurobi_stage_bench.jl on finished runs. usage: gurobi_stage_bench.sh C run_dir [run_dir ...]
set -euo pipefail
C=$(realpath "$1"); shift
P=problems/mpopf/tadmm/review_pilot
cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
D=/home/aryan.jha/kamiak_jobs
export JULIA_DEPOT_PATH="$C/julia-depot:$D/tadmm_review_nov11_large10k_interim_ice/julia-depot:$D/tadmm_review_nov08_horizon_sweep/julia-depot:$D/tadmm_review_nov08_watchdog_rho/julia-depot:$D/tadmm_review_nov07_accuracy_fix/julia-depot:$D/tadmm_review_nov07_parallel/julia-depot:$D/tadmm_review_nov07_v2/julia-depot:$D/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR LD_PRELOAD
LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57 TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=july REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0 SC_ETA=0.95
export GRB_LICENSE_FILE=/opt/apps/gurobi/13.0.1/licenses/gurobi.lic SC_STAGE_SOLVER=gurobi
for R in "$@"; do
 case $R in *ieee2522C_1ph*) export TADMM_SCREEN=substation,ell;; *) export TADMM_SCREEN=substation,vupper,ell;; esac
 echo "== $R"
 julia --startup-file=no --project="$P/conic_env" "$P/gurobi_stage_bench.jl" "$R/instance.jls" "$R/lean_trace_timed.csv" 2>&1 | grep -E "^BENCH|ERROR|rror" || true
done
