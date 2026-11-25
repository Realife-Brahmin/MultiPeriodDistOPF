#!/usr/bin/env bash
# Centralized Gurobi under four settings on one node, each against Ipopt/MA57:
# as in the solver table (all bounds, tolerance 1e-8), with the Ipopt model's
# screening, with Gurobi's default cone tolerance (1e-6), and both.
# usage: gurobi_fairness_kamiak.sh C system T
set -euo pipefail
C=$(realpath "$1"); SYS=$2; T=$3
P=problems/mpopf/tadmm/review_pilot
OUT="$C/results/fairness_${SYS}_T${T}_${SLURM_JOB_ID}"; mkdir -p "$OUT"; cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
D=/home/aryan.jha/kamiak_jobs
export JULIA_DEPOT_PATH="$C/julia-depot:$D/tadmm_review_nov11_large10k_interim_ice/julia-depot:$D/tadmm_review_nov08_horizon_sweep/julia-depot:$D/tadmm_review_nov08_watchdog_rho/julia-depot:$D/tadmm_review_nov07_accuracy_fix/julia-depot:$D/tadmm_review_nov07_parallel/julia-depot:$D/tadmm_review_nov07_v2/julia-depot:$D/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR LD_PRELOAD
LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57 TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=july REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0 SC_ETA=${SC_ETA:-0.95}
case $SYS in ieee2522C_1ph) export TADMM_SCREEN=substation,ell;; *) export TADMM_SCREEN=substation,vupper,ell;; esac
export GRB_LICENSE_FILE=/opt/apps/gurobi/13.0.1/licenses/gurobi.lic
run(){ tag=$1; shift; mkdir -p "$OUT/$tag"; env "$@" julia --startup-file=no --project="$P/conic_env" "$P/${CONIC_SCRIPT:-conic_baseline.jl}" "$SYS" "$T" "$OUT/$tag" ${SOLVERS:-gurobi} > "$OUT/$tag/solve.log" 2> "$OUT/$tag/solve.err"; grep -E "^CONIC solver" "$OUT/$tag/solve.log" | sed "s/^/$tag /"; }
SOLVERS=ipopt run ipopt X=1
run base CONIC_SCREEN=0
run screen CONIC_SCREEN=1
run deftol CONIC_SCREEN=0 CONIC_QCPTOL=1e-6
run screen_deftol CONIC_SCREEN=1 CONIC_QCPTOL=1e-6
date -Is > "$OUT/COMPLETE"
