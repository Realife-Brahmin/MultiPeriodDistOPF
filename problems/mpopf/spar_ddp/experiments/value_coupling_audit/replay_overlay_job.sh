#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2
CAP=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
export VALUE_REPLAY_NUMERICAL_ROOT="$BASE/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
unset FILTERDDP_VALUE_GROUPING FILTERDDP_FROZEN_CAPTURE_DIR NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
(cd "$CAP/source";sha256sum -c "$CAP/input.sha256") > "$C/logs/input-check.log"
julia --startup-file=no "$C/hierarchy_unit.jl" > "$C/logs/unit.log" 2>&1
run_case() {
  local SYS=$1 T=$2 W=$3
  python3 "$VALUE_REPLAY_NUMERICAL_ROOT/spar_ddp/experiments/value_coupling_audit/verify_coverage.py" "$CAP" "$SYS" "$T" > "$C/logs/${SYS}_T${T}_coverage.json"
  timeout --signal=TERM --kill-after=20s 12m julia --startup-file=no --project="$VALUE_REPLAY_NUMERICAL_ROOT/envs/spar_ddp" "$C/replay_series.jl" "$CAP" "$SYS" "$T" "$W" "$C/results/${SYS}_T${T}" > "$C/logs/${SYS}_T${T}_series.log" 2>&1
  grep -q '^FROZEN_SERIES_COMPLETE' "$C/logs/${SYS}_T${T}_series.log"
}
run_case ieee123C_1ph 3 0
date -u +%FT%TZ > "$C/smoke.complete"
run_case ieee2522C_1ph 24 0
run_case large10kC_1ph 24 4
date -u +%FT%TZ > "$C/replay.finished"
