#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1")
SOURCE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22/source
CAP=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22/results/ieee2522C_1ph_T24_dense_capture/captures/pass0061_k0052_j008
cd "$SOURCE"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
unset FILTERDDP_VALUE_GROUPING FILTERDDP_FROZEN_CAPTURE_DIR NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
timeout --signal=TERM --kill-after=20s 5m julia --startup-file=no --project=envs/spar_ddp "$C/diagnose_reference.jl" "$SOURCE" "$CAP" "$C/results/late" > "$C/logs/reference.log" 2>&1
date -u +%FT%TZ > "$C/reference.finished"
