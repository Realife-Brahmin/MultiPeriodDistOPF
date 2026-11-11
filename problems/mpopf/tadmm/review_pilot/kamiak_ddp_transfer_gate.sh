#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); REF=$(realpath "$2")
cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$REF/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TADMM_LINEAR_SOLVER=ma57
export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
julia --startup-file=no --project=problems/mpopf/tadmm/review_pilot problems/mpopf/tadmm/review_pilot/validate_ddp_transfer.jl "$C/source" "$REF" "$C/transfer-gate-$SLURM_JOB_ID"
date -u +%FT%TZ > "$C/DDP_TRANSFER_VALIDATED"
