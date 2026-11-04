#!/usr/bin/env bash
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
C="$SLURM_SUBMIT_DIR"
CAP=/home/aryan.jha/kamiak_jobs/compact_campaign_capture_oct23_v3
export VALUE_FOLLOWUP_BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3
export VALUE_REPLAY_NUMERICAL_ROOT="$CAP/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$CAP/julia-depot:$VALUE_FOLLOWUP_BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
while IFS= read -r K;do unset "$K";done < <(compgen -v | grep '^FILTERDDP_' || true)
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_BLOCKED_SOLVE=16
sha256sum -c source.sha256
sha256sum -c "$CAP/runtime.sha256"
julia --version
julia --startup-file=no --project="$VALUE_REPLAY_NUMERICAL_ROOT/envs/spar_ddp" late.jl "$CAP/results/large10kC_1ph_T24_campaign_capture/captures/pass0023_k0014_j008"
