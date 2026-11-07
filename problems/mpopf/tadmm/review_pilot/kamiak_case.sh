#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); SYS=$2; T=$3; SOLVER=$4; LABEL=$5
cd "$C/source"
PILOT=problems/mpopf/tadmm/review_pilot
OUT="$C/results/${SYS}_T${T}_${SOLVER}_${LABEL}"
[[ ! -e "$OUT" ]] || { echo "Refusing overwrite: $OUT"; exit 2; }
mkdir -p "$OUT"
export TADMM_LINEAR_SOLVER=$SOLVER TADMM_RHO_POLICY=relative
case $SOLVER in
 ma57) LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so; EXTRA=ma57_pivot_order=2;;
 ma97) LIB=$HOME/opt/hsl-private/install/hsl_ma97-2.8.1-metis-openmp/lib/libhsl_ma97.so; EXTRA=ma97_order=amd;;
 *) echo 'This campaign permits MA57/MA97 only'; exit 2;;
esac
export TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;$EXTRA${TADMM_ADDITIONAL_IPOPT_OPTIONS:+;$TADMM_ADDITIONAL_IPOPT_OPTIONS}"
cp "$0" "$OUT/case-wrapper.sh"
sha256sum "$OUT/case-wrapper.sh" > "$OUT/case-wrapper.sha256"
case $SYS in ieee123C_1ph) CB=.001;; ieee2522C_1ph) CB=.0019;; large10kC_1ph) CB=.000024;; *) exit 2;; esac
export TADMM_RHO_INIT=$(python3 -c 'import sys;print(2*float(sys.argv[1])*1000**2*int(sys.argv[2])/24)' "$CB" "$T")
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
W=${TADMM_NWORKERS:-0}
export TADMM_NWORKERS=$W
if (( W>0 )); then
 BUDGET=$((W+1))
 CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
 export TADMM_CPU_IDS=$(python3 -c 'import sys; a=[]
for p in sys.argv[1].split(","):
 b=p.split("-"); a.extend(range(int(b[0]),int(b[-1])+1))
n=int(sys.argv[2]); assert len(a)>=n; print(",".join(map(str,a[:n])))' "$CPU_RANGE" "$BUDGET")
 taskset -apc "$TADMM_CPU_IDS" "$$" > "$OUT/case-affinity.txt"
 if [[ $SOLVER == ma97 ]]; then export OMP_NUM_THREADS=$BUDGET; else export OMP_NUM_THREADS=1; fi
 export OMP_WAIT_POLICY=PASSIVE GOMP_SPINCOUNT=0
else
 export OMP_NUM_THREADS=1
fi
sha256sum "$LIB" > "$OUT/hsl.sha256"
{ hostname; date -Is; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(OMP_|OPENBLAS_|JULIA_|TADMM_)'; } > "$OUT/allocation.txt"
python3 "$PILOT/capture.py" "$OUT" "$SYS"
sha256sum "problems/mpopf/spar_ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls" > "$OUT/input.sha256"
if (( W>0 )); then
 taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$PILOT" "$PILOT/run.jl" "$SYS" "$T" "$OUT" > "$OUT/solve.log" 2> "$OUT/solve.err" &
else
 julia --startup-file=no --project="$PILOT" "$PILOT/run.jl" "$SYS" "$T" "$OUT" > "$OUT/solve.log" 2> "$OUT/solve.err" &
fi
OPT_PID=$!
python3 "$PILOT/sample_linux.py" "$OPT_PID" "$OUT/process_samples.csv" > "$OUT/sampler.log" 2>&1 &
SAMPLE_PID=$!
wait "$OPT_PID"
wait "$SAMPLE_PID"
for METHOD in central tadmm; do
 julia --startup-file=no --project="$PILOT" problems/mpopf/spar_ddp/examples/power_system/opendss_check.jl "$SYS" "$T" "$OUT/${METHOD}_dispatch.jls" "$METHOD" > "$OUT/${METHOD}_opendss.log" 2>&1
done
python3 "$PILOT/validate.py" "$OUT" > "$OUT/validation.log"
date -u +%FT%TZ > "$OUT/COMPLETE"
CG=/sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_$SLURM_JOB_ID
for F in memory.max_usage_in_bytes memory.limit_in_bytes; do
 [[ ! -f "$CG/$F" ]] || printf '%s %s\n' "$F" "$(cat "$CG/$F")" >> "$OUT/cgroup_memory.txt"
done
echo "CASE_COMPLETE $OUT"
