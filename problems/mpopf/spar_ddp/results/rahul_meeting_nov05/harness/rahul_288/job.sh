#!/usr/bin/env bash
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
C=$PWD
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
S="$C/variants/10k_1020/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
{ hostname; julia --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "hardware-$SLURM_JOB_ID.txt"
grep -q 'Gold 6438Y+' "hardware-$SLURM_JOB_ID.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c harness.sha256
(cd "$S"; sha256sum -c "$C/source.sha256" >/dev/null)
if [[ $1 == pair ]]; then
  T=96
  cd "$S"; sha256sum -c "$BASE/original-input.sha256"
  for mode in balanced grouped; do
    export SPAR_MEETING_GROUPED=0
    [[ $mode == balanced ]] || export SPAR_MEETING_GROUPED=1
    bash spar_ddp/experiments/battery_count/run.sh "$C" 10k_1020 "$T" ddp "$mode"
  done
else
  T=$1
  cd "$S"
  F="spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T${T}_periodic.jls"
  sha256sum -c "$BASE/original-input.sha256"
  test ! -e "$F"
  ln -s "$BASE/inputs/${F##*/}" "$F"
  sha256sum "$F" > "$C/input-T$T.sha256"
  export SPAR_MEETING_GROUPED=0
  for rep in 1 2; do
    sha256sum -c "$C/input-T$T.sha256"
    bash spar_ddp/experiments/battery_count/run.sh "$C" 10k_1020 "$T" ddp "$rep"
  done
fi
date -u +%FT%TZ > "$C/complete-$1"
