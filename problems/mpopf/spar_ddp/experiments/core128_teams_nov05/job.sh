#!/usr/bin/env bash
set -euo pipefail
C=$SLURM_SUBMIT_DIR; cd "$C"
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1
export STUDY_CORE_BUDGET=128
T=96
M="$C/jobs/$SLURM_JOB_ID"; mkdir "$M"
exec > >(tee "$M/stdout.log") 2> >(tee "$M/stderr.log" >&2)
{ date -u +%FT%TZ; hostname; julia --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q '6760P' "$M/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c harness.sha256 > "$M/harness-check.log"
S="$C/variants/10k_1020/source"
(cd "$S"; sha256sum -c "$C/execution-source.sha256") > "$M/source-check.log"
run_one() {
 local arm=$1 rep=$2
 export STUDY_SCAN=grouped STUDY_KERNEL=right
 case "$arm" in
 w48t2) export STUDY_WORKERS=48 STUDY_BLAS_THREADS=2 ;;
 w24t4) export STUDY_WORKERS=24 STUDY_BLAS_THREADS=4 ;;
 *) exit 2 ;;
 esac
 local label="${arm}_${rep}"
 local mapdir="$M/${T}_${label}"
 mkdir "$mapdir"
 (cd "$mapdir"; python3 "$C/affinity.py" "$STUDY_WORKERS" "$STUDY_BLAS_THREADS")
 export STUDY_CPU_MAP="$mapdir/cpu_groups.csv"
 export STUDY_OUT="$C/results/10k_1020_T${T}_ddp_r${label}"
 echo "TEAMS_BEGIN arm=$arm T=$T repeat=$rep workers=$STUDY_WORKERS blas=$STUDY_BLAS_THREADS"
 bash "$S/spar_ddp/experiments/battery_count/run.sh" "$C" 10k_1020 "$T" ddp "$label"
 test -s "$STUDY_OUT/COMPLETE"
 sha256sum "$STUDY_OUT/solution.jls" >> "$M/solution.sha256"
 echo "TEAMS_COMPLETE arm=$arm T=$T repeat=$rep"
}
for arm in w48t2 w24t4; do run_one "$arm" 1; done
for arm in w24t4 w48t2; do run_one "$arm" 2; done
(cd "$S"; sha256sum -c "$C/execution-source.sha256") > "$M/source-check-after.log"
sha256sum -c harness.sha256 > "$M/harness-check-after.log"
date -u +%FT%TZ > T96.finished

