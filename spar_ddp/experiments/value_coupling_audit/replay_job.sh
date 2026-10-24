#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); CAP=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
unset FILTERDDP_VALUE_GROUPING FILTERDDP_FROZEN_CAPTURE_DIR NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
test -s "$CAP/capture.finished"
(cd "$CAP/source";sha256sum -c "$CAP/input.sha256") > "$C/logs/input-check.log"
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg;Pkg.develop(path="spar_ddp/DDP4OPF.jl");Pkg.instantiate()' > "$C/logs/setup.log" 2>&1
julia --startup-file=no spar_ddp/experiments/value_coupling_audit/hierarchy_unit.jl > "$C/logs/unit.log" 2>&1
run_case() {
  local SYS=$1 T=$2 W=$3 D="$CAP/results/${1}_T${2}_dense_capture"
  python3 spar_ddp/experiments/value_coupling_audit/verify_coverage.py "$CAP" "$SYS" "$T" > "$C/logs/${SYS}_T${T}_coverage.json"
  local INPUT="$CAP/source/spar_ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls"
  local OUT="$C/results/${SYS}_T${T}";mkdir -p "$OUT"
  # Exactly the five selected frozen points; no optimizer invocation.
  while IFS=, read -r TAG SNAP REST; do
    [[ $TAG != tag ]] || continue
    timeout --signal=TERM --kill-after=20s 8m julia --startup-file=no --project=envs/spar_ddp spar_ddp/experiments/value_coupling_audit/replay.jl "$D/captures/$SNAP" "$INPUT" "$OUT/$TAG" "$W" > "$OUT/${TAG}.log" 2>&1
    grep -q '^FROZEN_REPLAY_COMPLETE' "$OUT/${TAG}.log"
  done < "$D/captures/coverage.csv"
}
run_case ieee123C_1ph 3 0
date -u +%FT%TZ > "$C/smoke.complete"
run_case ieee2522C_1ph 24 0
run_case large10kC_1ph 24 4
date -u +%FT%TZ > "$C/replay.finished"
