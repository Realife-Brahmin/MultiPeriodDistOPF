#!/usr/bin/env bash
set -euo pipefail
C="$SLURM_SUBMIT_DIR"; cd "$C"
T=${1:-96}
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export OPENBLAS_MAIN_FREE=1
M="$C/jobs/$SLURM_JOB_ID"; mkdir -p "$M"
exec > >(tee "$M/stdout.log") 2> >(tee "$M/stderr.log" >&2)
{ date -u +%FT%TZ; hostname; julia --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q 'Gold 6438Y+' "$M/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c harness.sha256 > "$M/harness-check.log"
S="$C/variants/10k_1020/source"
(cd "$S"; sha256sum -c "$C/execution-source.sha256") > "$M/source-check.log"
run_one() {
    local arm=$1 horizon=$2 rep=$3
    export STUDY_SCAN=grouped STUDY_KERNEL=current STUDY_WORKERS=63 STUDY_BLAS_THREADS=1
    case "$arm" in
        c1) export STUDY_SCAN=balanced ;;
        c2) ;;
        c3) export STUDY_KERNEL=right ;;
        c4) export STUDY_KERNEL=right STUDY_WORKERS=31 STUDY_BLAS_THREADS=2 ;;
        *) exit 2 ;;
    esac
    local actual=$((horizon<STUDY_WORKERS?horizon:STUDY_WORKERS))
    local label="${arm}_${rep}"
    local mapdir="$C/jobs/$SLURM_JOB_ID/${horizon}_${label}"
    mkdir "$mapdir"
    (cd "$mapdir"; python3 "$C/affinity.py" "$actual" "$STUDY_BLAS_THREADS")
    export STUDY_CPU_MAP="$mapdir/cpu_groups.csv"
    export STUDY_OUT="$C/results/10k_1020_T${horizon}_ddp_r${label}"
    printf 'FULL_STUDY_BEGIN arm=%s T=%s repeat=%s workers=%s blas=%s\n' "$arm" "$horizon" "$rep" "$actual" "$STUDY_BLAS_THREADS"
    bash "$S/spar_ddp/experiments/battery_count/run.sh" "$C" 10k_1020 "$horizon" ddp "$label"
    test -s "$STUDY_OUT/COMPLETE"
    printf 'FULL_STUDY_COMPLETE arm=%s T=%s repeat=%s\n' "$arm" "$horizon" "$rep"
}
if [[ $T == 96 ]]; then
    for arm in c1 c2 c3 c4; do run_one "$arm" 3 smoke; done
    date -u +%FT%TZ > smoke.finished
    for arm in c1 c2 c3 c4; do run_one "$arm" 96 1; done
    for arm in c4 c3 c2 c1; do run_one "$arm" 96 2; done
else
    test -s T96.validated
    test -s selected_candidate.txt
    candidate=$(cat selected_candidate.txt)
    [[ $candidate == c3 || $candidate == c4 ]] || exit 2
    for arm in c1 "$candidate"; do run_one "$arm" "$T" 1; done
    for arm in "$candidate" c1; do run_one "$arm" "$T" 2; done
fi
date -u +%FT%TZ > "T${T}.finished"
