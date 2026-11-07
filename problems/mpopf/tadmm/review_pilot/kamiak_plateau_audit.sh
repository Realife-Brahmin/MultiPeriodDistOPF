#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); BASE=$(realpath "$2"); CASE=$(realpath "$3")
SCRIPT=${4:-audit_plateau.jl}
VARIANTS=${5:-default,centered,strict_comp,adaptive}; STEPS=${6:-8}
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=ma57
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
OUT="$C/diagnostics/plateau_audit_$SLURM_JOB_ID"
julia --startup-file=no --project=problems/mpopf/tadmm/review_pilot "$C/harness/$SCRIPT" "$C/source" "$CASE" "$OUT" "$VARIANTS" "$STEPS"
sha256sum "$C/harness/$SCRIPT" "$0" > "$OUT/orchestration.sha256"
date -u +%FT%TZ > "$OUT/AUDIT_FINISHED"
