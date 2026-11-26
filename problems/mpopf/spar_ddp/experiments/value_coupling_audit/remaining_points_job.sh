#!/usr/bin/env bash
# New independent outputs for ONLY the two points missing after the series cap.
# No optimizer resume and no repeats of the thirteen completed replay points.
set -euo pipefail
C=$(realpath "$1"); BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3
CAP=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
export VALUE_REPLAY_NUMERICAL_ROOT=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2/source
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
unset FILTERDDP_VALUE_GROUPING FILTERDDP_FROZEN_CAPTURE_DIR NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
test -s "$BASE/smoke.complete"
test -s "$BASE/results/ieee2522C_1ph_T24/late/metadata.toml"
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
sha256sum -c "$BASE/source.sha256" > "$C/logs/inherited-source-check.log"
(cd "$CAP/source";sha256sum -c "$CAP/input.sha256") > "$C/logs/input-check.log"
SYS=large10kC_1ph;T=24
DATA="$CAP/results/${SYS}_T${T}_dense_capture"
INPUT="$CAP/source/spar_ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls"
python3 "$VALUE_REPLAY_NUMERICAL_ROOT/spar_ddp/experiments/value_coupling_audit/verify_coverage.py" "$CAP" "$SYS" "$T" > "$C/logs/coverage.json"
mkdir -p "$C/results/${SYS}_T${T}"
for TAG in middle late; do
  # Refuse to duplicate any completed point from the interrupted attempt.
  test ! -s "$BASE/results/${SYS}_T${T}/$TAG/metadata.toml"
  SNAP=$(awk -F, -v tag="$TAG" '$1==tag {print $2}' "$DATA/captures/coverage.csv")
  [[ $SNAP =~ ^pass[0-9]+_k[0-9]+_j[0-9]+$ ]]
  OUT="$C/results/${SYS}_T${T}/$TAG"
  timeout --signal=TERM --kill-after=20s 8m julia --startup-file=no --project="$VALUE_REPLAY_NUMERICAL_ROOT/envs/spar_ddp" "$BASE/replay.jl" "$DATA/captures/$SNAP" "$INPUT" "$OUT" 4 > "$C/results/${SYS}_T${T}/${TAG}.log" 2>&1
  grep -q '^FROZEN_REPLAY_COMPLETE' "$C/results/${SYS}_T${T}/${TAG}.log"
done
date -u +%FT%TZ > "$C/remaining.finished"
