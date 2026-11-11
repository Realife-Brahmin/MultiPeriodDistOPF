#!/usr/bin/env bash
# Orchestration only: run against the already frozen and verified numerical source.
set -euo pipefail
C=$(realpath "$1"); BASE=$(realpath "$2"); MODE=$3
cd "$C/source"
PILOT=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH="$BASE/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
test -s "$C/parallel-gate.complete"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
TOKEN=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-single}_a${SLURM_RESTART_COUNT:-0}
mkdir -p "$C/campaign"
exec > >(tee "$C/logs/${MODE}_${TOKEN}.log") 2>&1
case $MODE in
 verify)
  SYS=$4; T=6; REP=0
  case $SYS in ieee2522C_1ph|large10kC_1ph) ;; *) exit 2;; esac
  ;;
 benchmark)
  test -s "$C/campaign/ieee2522C_1ph.verified"
  test -s "$C/campaign/large10kC_1ph.verified"
  INDEX=$SLURM_ARRAY_TASK_ID
  SYSTEMS=(ieee123C_1ph ieee2522C_1ph large10kC_1ph)
  SYS=${SYSTEMS[$((INDEX / 4))]}
  if (( (INDEX / 2) % 2 == 0 )); then T=6; else T=24; fi
  REP=$((INDEX % 2 + 1))
  ;;
 *) exit 2;;
esac
export TADMM_NWORKERS=$T
HW="$C/campaign/${MODE}_${TOKEN}_hardware.txt"
{ hostname; lscpu; julia --version; scontrol show job "$SLURM_JOB_ID"; scontrol show node "$(hostname -s)" -o; } > "$HW"
grep -q '6438Y+' "$HW"
if [[ $MODE == benchmark ]]; then
 # Kamiak disallows --exclusive. Reserve all 64 CPUs on a 64-CPU node instead;
 # OverSubscribe=NO prevents CPU sharing. Numerical work uses only T+1 cores.
 # Never label shared verification as benchmark timing.
 grep -q 'NumCPUs=64' "$HW"
 grep -q 'CPUTot=64' "$HW"
 grep -q 'ThreadsPerCore=1' "$HW"
fi
sha256sum "$0" > "$C/campaign/${MODE}_${TOKEN}_orchestration.sha256"
SOLVERS=(ma57 ma97)
if (( REP == 2 )); then SOLVERS=(ma97 ma57); fi
LABEL=${MODE}_r${REP}_${TOKEN}
for SOLVER in "${SOLVERS[@]}"; do
 bash "$PILOT/kamiak_case.sh" "$C" "$SYS" "$T" "$SOLVER" "$LABEL"
done
# Completion requires both solvers' full optimization and independent validation.
for SOLVER in ma57 ma97; do
 test -s "$C/results/${SYS}_T${T}_${SOLVER}_${LABEL}/COMPLETE"
done
if [[ $MODE == verify ]]; then
 printf '%s %s\n' "$TOKEN" "$LABEL" > "$C/campaign/$SYS.verified"
else
 printf '%s %s %s %s\n' "$SYS" "$T" "$REP" "$LABEL" > "$C/campaign/benchmark_${INDEX}.complete"
fi
date -u +%FT%TZ
