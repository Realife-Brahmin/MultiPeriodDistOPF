#!/usr/bin/env bash
# First assessment: one chosen backend and one measured run, after warm-up.
set -euo pipefail
C=$(realpath "$1"); BASE=$(realpath "$2"); SYS=$3; T=$4; SOLVER=$5
case $SYS in ieee2522C_1ph|large10kC_1ph) ;; *) exit 2;; esac
[[ $T == 24 ]] || { echo 'This launch protocol is frozen for T=24'; exit 2; }
case $SOLVER in ma57|ma97) ;; *) exit 2;; esac
cd "$C/source"
PILOT=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH="$BASE/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export TADMM_NWORKERS=$T
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
unset TADMM_ADDITIONAL_IPOPT_OPTIONS TADMM_DIAG_CHECKPOINT
test -s "$C/campaign/ieee2522C_1ph.verified"
test -s "$C/campaign/large10kC_1ph.verified"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
TOKEN=${SLURM_JOB_ID}_a${SLURM_RESTART_COUNT:-0}
exec > >(tee "$C/logs/firstpass_$TOKEN.log") 2>&1
HW="$C/campaign/firstpass_${TOKEN}_hardware.txt"
{ hostname; lscpu; julia --version; scontrol show job "$SLURM_JOB_ID"; scontrol show node "$(hostname -s)" -o; } > "$HW"
grep -q '6438Y+' "$HW"
grep -q 'NumCPUs=64' "$HW"
grep -q 'CPUTot=64' "$HW"
grep -q 'ThreadsPerCore=1' "$HW"
sha256sum "$0" > "$C/campaign/firstpass_${TOKEN}_orchestration.sha256"
bash "$PILOT/kamiak_case.sh" "$C" "$SYS" "$T" "$SOLVER" "firstpass_$TOKEN"
test -s "$C/results/${SYS}_T${T}_${SOLVER}_firstpass_$TOKEN/COMPLETE"
date -u +%FT%TZ > "$C/campaign/${SYS}_T${T}_${SOLVER}_firstpass_$TOKEN.complete"
