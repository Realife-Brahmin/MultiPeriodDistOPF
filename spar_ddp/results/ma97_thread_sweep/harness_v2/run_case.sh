#!/usr/bin/env bash
# One independently restartable, fully warmed benchmark case per allocation.
set -euo pipefail
ROOT=${SLURM_SUBMIT_DIR:?}
SYS=$1 T=$2 WORKERS=$3 SOLVER=$4 SCHEDULE=$5 REPEAT=$6
test -s "$ROOT/smoke_verified.txt"
cd "$ROOT/repo"
export PATH=$HOME/bin:$PATH JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${OMP_SWEEP_THREADS:?} OMP_DYNAMIC=FALSE
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
OUT=$ROOT/logs/${SLURM_JOB_ID}_${SYS}_T${T}_${SOLVER}_${SCHEDULE}_r${REPEAT}_attempt${SLURM_RESTART_COUNT:-0}
mkdir -p "$OUT"
exec > >(tee "$OUT/stdout.log") 2> >(tee "$OUT/stderr.log" >&2)
trap 'c=/sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_${SLURM_JOB_ID}; for f in memory.max_usage_in_bytes memory.limit_in_bytes; do [[ ! -f "$c/$f" ]] || echo "$f $(cat "$c/$f")"; done > "$OUT/node_memory.txt"' EXIT
(cd "$ROOT" && sha256sum -c source.sha256) > "$OUT/source_check.txt"
sha256sum -c "$ROOT/input_${SYS}_T${T}.sha256" > "$OUT/input_check.txt"
{ date -Is; hostname; lscpu; scontrol show job "$SLURM_JOB_ID"; env | grep -E 'OMP_|BLAS|JULIA_NUM_THREADS'; } > "$OUT/metadata.txt"
grep -q '6438Y+' "$OUT/metadata.txt"
for ((waited=0; waited<=180; waited+=5)); do
 squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%i %T %j' > "$OUT/node_jobs.txt"
 if awk -v mine="$SLURM_JOB_ID" '$1!=mine {bad=1} END {exit bad}' "$OUT/node_jobs.txt"; then break; fi
 ((waited<180)) || { echo NODE_GUARD_FAILED; exit 1; }
 sleep 5
done
tag=properlog_${SLURM_JOB_ID}_a${SLURM_RESTART_COUNT:-0}_${SOLVER}_${SCHEDULE}_r${REPEAT}
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
if [[ "$SOLVER" == ddp ]]; then
 export FILTERDDP_SLURM_WORKERS=0 FILTERDDP_PROCS_REMOTE_PAIR=0 FILTERDDP_PROCS_GC=pass
 export FILTERDDP_PROCS_WORKERS=$WORKERS FILTERDDP_PROCS_PAIR_SCHEDULE=$SCHEDULE
 export FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_PROCS_ROUND_DIAGNOSTIC=1
 export PAIRED_WARMUP=full PAIRED_ROUNDS=1
 if [[ "$SYS" == ieee123C_1ph ]]; then export PAIRED_TOL=1e-6; else export PAIRED_TOL=1e-5; fi
 bash ddp/results/paired_pass/run_paired.sh "$SYS" "$T" procs "$tag" > "$OUT/launcher.log" 2>&1
 cp "ddp/results/parallel_in_time/logs/fddp_${SYS}_T${T}_paired_${tag}.log" "$OUT/solve.log"
 cp "ddp/results/opendss_check/parsim_${SYS}_T${T}_paired_${tag}.log" "$OUT/opendss.log"
 cp "ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CBsystem.jls" "$OUT/solution.jls"
 python3 "$ROOT/check_logs.py" "$OUT/solve.log" | tee "$OUT/logging_check.json"
else
 [[ "$SOLVER" == ma97 || "$SOLVER" == ma57 ]] || exit 2
 export IPOPT_SCREEN=substation,vupper,ell,psubs IPOPT_LOADFLOW_START=1 IPOPT_MATCH_DDP_START=1 IPOPT_TRACE=1 IPOPT_WARMUP=full
 if [[ "$SOLVER" == ma97 ]]; then
  lib=$HOME/opt/hsl-private/install/hsl_ma97-2.8.1-metis-openmp/lib/libhsl_ma97.so
  extra=ma97_order=amd
 else
  lib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
  extra=ma57_pivot_order=2
 fi
 sha256sum "$lib" > "$OUT/private_library.sha256"
 export IPOPT_EXTRA_OPTIONS="linear_solver=$SOLVER;hsllib=$lib;linear_system_scaling=none;tol=1e-6;mu_strategy=monotone;inf_pr_output=original;$extra"
 export IPOPT_SAVE_DISPATCH=$OUT/solution.jls
 export MA97_SWEEP_LIBRARY=$lib MA97_SWEEP_OUT=$OUT
 julia --startup-file=no --project=envs/ddp2026 "$ROOT/thread_entry.jl" "$SYS" "$T" "$OUT/ipopt_internal.log" > "$OUT/solve.log" 2>&1
 grep -q 'status=LOCALLY_SOLVED' "$OUT/solve.log"
 python3 "$ROOT/check_ipopt_logs.py" "$OUT/solve.log" "$OUT/ipopt_internal.log" | tee "$OUT/logging_check.json"
 bash ddp/examples/power_system/run_opendss_check.sh "$SYS" "$T" "$OUT/solution.jls" "$tag" > "$OUT/opendss_launcher.log" 2>&1
 cp "ddp/results/opendss_check/$tag.log" "$OUT/opendss.log"
fi
grep -q 'OPENDSS_SUMMARY.*voltage.*below limit: 0.*above limit: 0' "$OUT/opendss.log"
if grep -E 'ERROR:|Unhandled Task|UNHANDLED TASK|ProcessExitedException' "$OUT/solve.log"; then exit 1; fi
echo "CASE_COMPLETE system=$SYS T=$T solver=$SOLVER schedule=$SCHEDULE repeat=$REPEAT" | tee "$OUT/COMPLETE"
