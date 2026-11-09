#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1");T=$2;MODE=$3;WARM=${4:-primal_dual};RHO=${5:-120000};SOLVER=${6:-ma57}
SYSTEM=${7:-large10kC_1ph}
[[ $SOLVER == ma57 || $SOLVER == mumps ]]
[[ $SYSTEM == large10kC_1ph || ( $SYSTEM == ieee2522C_1ph && $T == 96 && $MODE != gate ) ]]
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_horizon_sweep/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_watchdog_rho/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
export JULIA_LOAD_PATH="$C/source/$P:$C/source/$P/slurm:@stdlib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=$SOLVER TADMM_CLUSTER_BACKEND=slurm TADMM_NWORKERS=$T TADMM_EXECUTION_POLICY=one_per_period
export TADMM_WARM_START=$WARM TADMM_OBJECTIVE_UPDATE=batched
export TADMM_RHO_POLICY=${TADMM_RHO_POLICY:-legacy_phase} TADMM_CB_PROFILE=july TADMM_LOCAL_SOLVE_POLICY=${TADMM_LOCAL_SOLVE_POLICY:-exact}
export TADMM_SCREEN=substation,vupper,ell TADMM_DIAGNOSTICS=typed TADMM_RHO_TAU_DECREASE=5
export TADMM_RHO_WATCHDOG=1 TADMM_RHO_PRIMAL_GATE=1e-4 TADMM_LOCAL_PRIMAL_REPORT=1 TADMM_MAX_ITER=1000
if [[ $SOLVER == ma57 ]];then
 export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
else
 export TADMM_IPOPT_OPTIONS="linear_system_scaling=none"
fi
export TADMM_JULIA_EXE=/home/aryan.jha/opt/julia-1.12.4/bin/julia
export TADMM_COORDINATOR_HOST=$(hostname)
export REDUCED_PROFILE=periodic TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
OUT="$C/results/${MODE}_T${T}_${WARM}_$SLURM_JOB_ID";test ! -e "$OUT";mkdir -p "$OUT/pids"
export TADMM_PID_REGISTRY="$OUT/pids"
export TADMM_CPU_IDS=$(python3 -c 'import os
s=next(x.split()[1] for x in open("/proc/self/status") if x.startswith("Cpus_allowed_list:"));a=[]
for p in s.split(","):
 b=p.split("-");a.extend(range(int(b[0]),int(b[-1])+1))
assert len(a)>=50;print(",".join(map(str,a)))')
{ hostname;lscpu;scontrol show job "$SLURM_JOB_ID";env | grep -E '^(TADMM_|SLURM_|JULIA_|OMP_|OPENBLAS_)'; } > "$OUT/allocation.txt"
grep -q '6438Y+' "$OUT/allocation.txt"
python3 "$P/capture.py" "$OUT" "$SYSTEM"
if [[ $SOLVER == ma57 ]];then
 sha256sum "$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so" > "$OUT/hsl.sha256"
fi
NODES=$SLURM_JOB_NUM_NODES
# Reserve 64 CPUs per node in sbatch, then launch exactly T one-thread workers.
# Slurm gives the step the allocated node CPUs; wrapper checks distinct pins.
export SLURM_NTASKS=$T SLURM_NTASKS_PER_NODE=$((T/NODES)) SLURM_CPUS_PER_TASK=1 SLURM_CPU_BIND=none
export TADMM_NODES=$NODES TADMM_WORKERS_PER_NODE=$((T/NODES))
export PATH="$C/source/$P/slurm/bin:$PATH"
finish() { touch "$OUT/SAMPLING_DONE"; }
trap finish EXIT
if [[ $MODE == gate ]]; then
 taskset -c "${TADMM_CPU_IDS%%,*}" julia --startup-file=no --project="$P" "$P/validate_stage_reuse.jl" > "$OUT/stage-reuse.log" 2>&1
 taskset -c "${TADMM_CPU_IDS%%,*}" julia --startup-file=no --project="$P" "$P/validate_full_period.jl" "$T" "$OUT" > "$OUT/solve.log" 2> "$OUT/solve.err"
 test -s "$OUT/FULL_PERIOD_GATE_VALIDATED"
else
 test -s "$C/FULL_PERIOD_READY"
 taskset -c "${TADMM_CPU_IDS%%,*}" julia --startup-file=no --project="$P" "$P/run_horizon_sweep.jl" "$SYSTEM" "$T" "$OUT" "$RHO" > "$OUT/solve.log" 2> "$OUT/solve.err"
 test -s "$OUT/SWEEP_VALIDATED"
fi
date -u +%FT%TZ > "$OUT/CASE_RECORDED"
