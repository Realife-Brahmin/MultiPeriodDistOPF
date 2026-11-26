#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1")
CAP=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3
export VALUE_FOLLOWUP_BASE="$BASE"
export VALUE_REPLAY_NUMERICAL_ROOT=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2/source
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
unset FILTERDDP_VALUE_GROUPING FILTERDDP_FROZEN_CAPTURE_DIR FILTERDDP_PAIRED_INPLACE NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
sha256sum -c "$BASE/source.sha256" > "$C/logs/inherited-source-check.log"
(cd "$CAP/source";sha256sum -c "$CAP/input.sha256") > "$C/logs/input-check.log"
julia --startup-file=no "$C/direct_initial_unit.jl" > "$C/logs/unit.log" 2>&1
DATA="$CAP/results/large10kC_1ph_T24_dense_capture"
SNAP=$(awk -F, '$1=="late" {print $2}' "$DATA/captures/coverage.csv")
[[ $SNAP =~ ^pass[0-9]+_k[0-9]+_j[0-9]+$ ]]
INPUT="$CAP/source/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls"
timeout --signal=TERM --kill-after=20s 25m julia --startup-file=no --project="$VALUE_REPLAY_NUMERICAL_ROOT/envs/spar_ddp" "$C/direct_initial_probe.jl" "$DATA/captures/$SNAP" "$INPUT" "$C/results/large_late" > "$C/logs/large_late.log" 2>&1
grep -q '^DIRECT_INITIAL_COMPLETE' "$C/logs/large_late.log"
date -u +%FT%TZ > "$C/direct_initial.finished"
