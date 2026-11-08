#!/usr/bin/env bash
# Two controlled arms, not repeat trials or a July input/profile replay.
set -euo pipefail
C=$(realpath "$1"); REF=$(realpath "$2"); SYSTEM=$3; MAXITER=${4:-200}
PAIR_REFERENCE=${5:-}
case $SYSTEM in
 ieee2522C_1ph) RHO=3800; POLICY=relative; SCREEN=substation,ell;;
 large10kC_1ph) RHO=48; POLICY=fixed; SCREEN=substation,vupper,ell;;
 *) exit 2;;
esac
RHO=${TADMM_REVIEW_RHO:-$RHO}
POLICY=${TADMM_REVIEW_POLICY:-$POLICY}
cd "$C/source"
PILOT=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$REF/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=ma57 TADMM_NWORKERS=24 TADMM_RHO_INIT=$RHO TADMM_RHO_POLICY=$POLICY
export TADMM_CB_PROFILE=july TADMM_LOCAL_SOLVE_POLICY=exact TADMM_SCREEN=$SCREEN
export TADMM_LOCAL_PRIMAL_REPORT=1
if [[ -n "$PAIR_REFERENCE" ]]; then export TADMM_LOCAL_PRIMAL_REPORT=0; fi
export TADMM_RHO_BALANCE_RATIO=10 TADMM_PRIMAL_TOL=2.5e-7 TADMM_DUAL_TOL=1e-4
export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
export REDUCED_PROFILE=periodic TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
test -s "$C/JULY_PROFILE_VALIDATED"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
OUT="$C/results/${SYSTEM}_T24_july_pair_$SLURM_JOB_ID"
if [[ -n "$PAIR_REFERENCE" ]]; then OUT="$C/results/${SYSTEM}_T24_july_light_$SLURM_JOB_ID"; fi
if [[ -n "${TADMM_HISTORICAL_REFERENCE:-}" ]]; then OUT="$C/results/${SYSTEM}_T24_historical_$SLURM_JOB_ID"; fi
[[ ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
export TADMM_CPU_IDS=$(python3 -c 'import sys; a=[]
for p in sys.argv[1].split(","):
 b=p.split("-"); a.extend(range(int(b[0]),int(b[-1])+1))
assert len(a)>=25; print(",".join(map(str,a[:25])))' "$CPU_RANGE")
taskset -apc "$TADMM_CPU_IDS" "$$" > "$OUT/affinity.txt"
{ hostname; lscpu; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(TADMM_|OMP_|OPENBLAS_|JULIA_)'; } > "$OUT/allocation.txt"
grep -q '6438Y+' "$OUT/allocation.txt"
grep -q 'NumCPUs=64' "$OUT/allocation.txt"
python3 "$PILOT/capture.py" "$OUT" "$SYSTEM"
sha256sum "$0" > "$OUT/orchestration.sha256"
if [[ -n "${TADMM_HISTORICAL_REFERENCE:-}" ]]; then
 taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$PILOT" "$PILOT/assess_historical_stop.jl" "$SYSTEM" 24 "$OUT" "$TADMM_HISTORICAL_REFERENCE" > "$OUT/solve.log" 2> "$OUT/solve.err" &
elif [[ -n "$PAIR_REFERENCE" ]]; then
 taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$PILOT" "$PILOT/assess_cap_light.jl" "$SYSTEM" 24 "$OUT" "$PAIR_REFERENCE" > "$OUT/solve.log" 2> "$OUT/solve.err" &
else
 taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$PILOT" "$PILOT/compare_july.jl" "$SYSTEM" 24 "$OUT" "$MAXITER" > "$OUT/solve.log" 2> "$OUT/solve.err" &
fi
PID=$!
python3 "$PILOT/sample_linux.py" "$PID" "$OUT/process_samples.csv" > "$OUT/sampler.log" 2>&1 &
SAMPLER=$!
wait "$PID"
wait "$SAMPLER"
test -s "$OUT/ASSESSMENT_VALIDATED"
date -u +%FT%TZ > "$OUT/CASE_RECORDED"
