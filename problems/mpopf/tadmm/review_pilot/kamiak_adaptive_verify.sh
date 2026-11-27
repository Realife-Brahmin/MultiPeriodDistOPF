#!/usr/bin/env bash
# Separately labeled calibration; this does not change the frozen benchmark queue.
set -euo pipefail
C=$(realpath "$1"); BASE=$(realpath "$2"); SYS=$3
cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH="$BASE/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export TADMM_NWORKERS=6
export TADMM_ADDITIONAL_IPOPT_OPTIONS='mu_strategy=adaptive;compl_inf_tol=1e-8'
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
for SOLVER in ma57 ma97; do
 bash "$C/harness/kamiak_case_adaptive.sh" "$C" "$SYS" 6 "$SOLVER" "adaptive_verify_$SLURM_JOB_ID"
done
sha256sum "$0" > "$C/campaign/adaptive_verify_$SLURM_JOB_ID.sha256"
date -u +%FT%TZ > "$C/campaign/${SYS}_adaptive_$SLURM_JOB_ID.complete"
