#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1")
CAP=/home/aryan.jha/kamiak_jobs/compact_campaign_capture_oct23_v3
CASE="$CAP/results/large10kC_1ph_T24_campaign_capture"
test -s "$CAP/campaign_capture.finished";test -s "$CASE/COMPLETE"
grep -q 'CAMPAIGN_CAPTURE_PASS' "$CASE/capture-check.log"
grep -q 'validated = true' "$CASE/capture_validation.toml"
export VALUE_FOLLOWUP_BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3
export VALUE_REPLAY_NUMERICAL_ROOT="$CAP/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$CAP/julia-depot:$VALUE_FOLLOWUP_BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
while IFS= read -r K;do unset "$K";done < <(compgen -v | grep '^FILTERDDP_' || true)
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_BLOCKED_SOLVE=16
unset NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
sha256sum -c "$CAP/source.sha256" > "$C/logs/capture-source-check.log"
sha256sum -c "$CAP/runtime.sha256" > "$C/logs/runtime-check.log"
sha256sum -c "$CAP/input.sha256" > "$C/logs/input-check.log"
cp "$CASE/capture_validation.toml" "$C/capture_validation.toml"
julia --startup-file=no "$C/campaign_compact_unit.jl" > "$C/logs/unit.log" 2>&1
SNAP=$(awk -F, '$1=="late" {print $2}' "$CASE/captures/coverage.csv")
[[ $SNAP =~ ^pass[0-9]+_k[0-9]+_j[0-9]+$ ]]
INPUT="$CAP/source/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls"
timeout --signal=TERM --kill-after=20s 25m julia --startup-file=no --project="$VALUE_REPLAY_NUMERICAL_ROOT/envs/spar_ddp" "$C/campaign_probe.jl" "$CASE/captures/$SNAP" "$INPUT" "$C/results/large_late" > "$C/logs/large_late.log" 2>&1
grep -q '^CAMPAIGN_PROBE_COMPLETE' "$C/logs/large_late.log"
date -u +%FT%TZ > "$C/campaign_replay.finished"
