#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); T=$2; RHOS=$3
cd "$C/source"
P=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_watchdog_rho/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=ma57 TADMM_NWORKERS=24 TADMM_RHO_POLICY=legacy_phase
export TADMM_CB_PROFILE=july TADMM_LOCAL_SOLVE_POLICY=exact TADMM_SCREEN=substation,vupper,ell
export TADMM_DIAGNOSTICS=typed TADMM_RHO_TAU_DECREASE=5 TADMM_RHO_WATCHDOG=1 TADMM_RHO_PRIMAL_GATE=1e-4
export TADMM_LOCAL_PRIMAL_REPORT=1 TADMM_MAX_ITER=1000
export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
export REDUCED_PROFILE=periodic TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
test -s "$C/JULY_PROFILE_VALIDATED"
OUT="$C/results/large10k_T${T}_$SLURM_JOB_ID"
test ! -e "$OUT";mkdir -p "$OUT"
CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
export TADMM_MONITOR_CPU=$(python3 -c 'import sys; a=[]
for p in sys.argv[1].split(","):
 b=p.split("-"); a.extend(range(int(b[0]),int(b[-1])+1))
assert len(a)>=26; print(a[25])' "$CPU_RANGE")
export TADMM_CPU_IDS=$(python3 -c 'import sys; a=[]
for p in sys.argv[1].split(","):
 b=p.split("-"); a.extend(range(int(b[0]),int(b[-1])+1))
assert len(a)>=25; print(",".join(map(str,a[:25])))' "$CPU_RANGE")
taskset -apc "$TADMM_CPU_IDS" "$$" > "$OUT/affinity.txt"
{ hostname; lscpu; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(TADMM_|OMP_|OPENBLAS_|JULIA_)'; } > "$OUT/allocation.txt"
grep -q '6438Y+' "$OUT/allocation.txt";grep -q 'NumCPUs=64' "$OUT/allocation.txt"
python3 "$P/capture.py" "$OUT" large10kC_1ph
sha256sum "$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so" > "$OUT/hsl.sha256"
taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$P" "$P/run_horizon_sweep.jl" large10kC_1ph "$T" "$OUT" "$RHOS" > "$OUT/solve.log" 2> "$OUT/solve.err" &
PID=$!
taskset -c "$TADMM_MONITOR_CPU" python3 "$P/sample_linux.py" "$PID" "$OUT/process_samples.csv" > "$OUT/sampler.log" 2>&1 &
SAMPLER=$!
wait "$PID";wait "$SAMPLER"
test -s "$OUT/SWEEP_VALIDATED"
date -u +%FT%TZ > "$OUT/CASE_RECORDED"
CG=/sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_$SLURM_JOB_ID
for F in memory.max_usage_in_bytes memory.limit_in_bytes; do
 if [[ -f "$CG/$F" ]]; then printf '%s %s\n' "$F" "$(cat "$CG/$F")" >> "$OUT/cgroup_memory.txt"; fi
done
