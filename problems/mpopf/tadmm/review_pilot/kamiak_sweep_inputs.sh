#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); REF=$(realpath "$2")
cd "$C/source"
P=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_watchdog_rho/julia-depot:$REF/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PROFILE_PERIODIC=1 TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
cp "$C/input.sha256" "$C/campaign/base-input.sha256"
for T in 48 96; do
 julia --startup-file=no --project="$P" problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl large10kC_1ph "$T" > "$C/logs/export-T$T.log" 2>&1
done
sha256sum problems/mpopf/spar_ddp/results/network_filterddp/network_data_*.jls > "$C/input.sha256"
date -u +%FT%TZ > "$C/SWEEP_INPUTS_READY"
export TADMM_DIAGNOSTIC_REFERENCE=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_july_weights TADMM_DIAGNOSTICS=typed
bash "$P/kamiak_july_gate.sh" "$C" "$REF"
