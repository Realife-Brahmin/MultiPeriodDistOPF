#!/usr/bin/env bash
set -euo pipefail
C=$(realpath -m "$1");ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
[[ ! -e "$C" && -z $(git -C "$ROOT" status --porcelain --untracked-files=no) ]] || exit 2
mkdir -p "$C/source" "$C/inputs" "$C/jobs" "$C/results"
git -C "$ROOT" rev-parse HEAD > "$C/source.commit"
git -C "$ROOT" archive HEAD AGENTS.md spar_ddp/DDP4OPF.jl spar_ddp/examples/power_system spar_ddp/experiments/battery_count spar_ddp/experiments/screening_ablation spar_ddp/experiments/tree_ablation envs/spar_ddp envs/tadmm rawData/large10kC_1ph rawData/LoadShapeDefault.dss rawData/LoadShapePVDefault.dss rawData/LoadShapePSubsCostDefault.dss | tar --no-same-owner -x -C "$C/source"
for T in 24 96 144 288; do cp "/home/aryan.jha/kamiak_jobs/proper_logging_campaign_oct16/repo/ddp/results/network_filterddp/network_data_large10kC_1ph_T${T}_periodic.jls" "$C/inputs/";done
sha256sum "$C"/inputs/*.jls > "$C/original-input.sha256"
cp /home/aryan.jha/kamiak_jobs/proper_logging_campaign_oct16/check_ipopt_logs.py "$C/check_ipopt_logs.py"
sha256sum "$C/check_ipopt_logs.py" > "$C/frozen-validator.sha256"
FLAGS=(--parsable --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=600G --time=04:00:00 --partition=kamiak --account=dubey '--constraint=sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64')
JOB=$(sbatch "${FLAGS[@]}" --job-name=battery-prepare --output="$C/jobs/prepare-%j.out" --error="$C/jobs/prepare-%j.err" "$C/source/spar_ddp/experiments/battery_count/job.sh" "$C" prepare)
printf 'job,variant,T,dependency\n%s,prepare,3,\n' "$JOB" > "$C/submitted.csv"
for T in 24 96 144 288; do
 for V in 10k_X 10k_Y; do
  PREV=$JOB
  JOB=$(sbatch "${FLAGS[@]}" --dependency="afterok:$PREV" --job-name="battery-${V}-T${T}" --output="$C/jobs/case-%j.out" --error="$C/jobs/case-%j.err" "$C/source/spar_ddp/experiments/battery_count/job.sh" "$C" case "$V" "$T")
  printf '%s,%s,%s,%s\n' "$JOB" "$V" "$T" "$PREV" >> "$C/submitted.csv"
 done
done
cat "$C/submitted.csv"
