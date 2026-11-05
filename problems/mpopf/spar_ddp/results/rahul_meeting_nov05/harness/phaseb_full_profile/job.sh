#!/usr/bin/env bash
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
C=$PWD; V=$1; T=$2
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
M="$C/jobs/$SLURM_JOB_ID"; mkdir -p "$M"
{ hostname; julia --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q 'Gold 6438Y+' "$M/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c harness.sha256
S="$C/variants/$V/source"; cd "$S"; sha256sum -c "$C/$V.sha256" > "$M/source-check.log"
for r in 1 2 3; do
  O="$C/results/${V}_T${T}_ddp_r${r}"
  export PROFILE_RUN_DIR=
  if [[ $r == 3 ]]; then export PROFILE_RUN_DIR="$O"; fi
  bash spar_ddp/experiments/battery_count/run.sh "$C" "$V" "$T" ddp "$r"
  if [[ $r == 3 ]]; then
    test -s "$O/profile-coordinator.txt"
    test "$(find "$O" -name 'profile-worker-*.txt' -size +0c | wc -l)" -eq 8
  fi
done
date -u +%FT%TZ > "$M/COMPLETE"
