#!/usr/bin/env bash
# One Kamiak case of the lean storage-coordinated tADMM with memory sampling.
# usage: kamiak_storage_coordinator_lean.sh C system T label
# SC_DRIVER=parallel runs the diagnostic driver instead, with the same sampling.
# SC_STAGE_SOLVER=mosek|clarabel runs tADMM's own solves with that conic solver
# (project conic_env) and first times the centralized problem with it too.
set -euo pipefail
C=$(realpath "$1"); SYS=$2; T=$3; LABEL=$4
P=problems/mpopf/tadmm/review_pilot
OUT="$C/results/${SYS}_T${T}_${LABEL}_${SLURM_JOB_ID}"
[[ ! -e "$OUT" ]] || { echo "Refusing overwrite: $OUT"; exit 2; }
mkdir -p "$OUT"; cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
D=/home/aryan.jha/kamiak_jobs
export JULIA_DEPOT_PATH="$C/julia-depot:$D/tadmm_review_nov11_large10k_interim_ice/julia-depot:$D/tadmm_review_nov08_horizon_sweep/julia-depot:$D/tadmm_review_nov08_watchdog_rho/julia-depot:$D/tadmm_review_nov07_accuracy_fix/julia-depot:$D/tadmm_review_nov07_parallel/julia-depot:$D/tadmm_review_nov07_v2/julia-depot:$D/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57
export TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=${TADMM_CB_PROFILE:-july}
case $SYS in
 ieee2522C_1ph) SCREEN=substation,ell;;
 ieee9500xf_sub1_1ph|smartds_p1u_xf_sub1_1ph) SCREEN=substation,ell; export SC_REPAIR_MAX_ITER=${SC_REPAIR_MAX_ITER:-3000};;
 *) SCREEN=substation,vupper,ell;;
esac
export TADMM_SCREEN=${TADMM_SCREEN:-$SCREEN}
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0
export SC_NWORKERS=${SC_NWORKERS:-$T} SC_RHO=${SC_RHO:-2:0.0625:20} SC_BALANCE=${SC_BALANCE:-2}
DRIVER=storage_coordinator_${SC_DRIVER:-lean}.jl
PROJECT=$P; [[ ${SC_STAGE_SOLVER:-ipopt} == ipopt ]] || PROJECT=$P/conic_env
{ hostname; date -Is; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(OMP_|OPENBLAS_|JULIA_|TADMM_|SC_)'; } > "$OUT/allocation.txt"
sha256sum "$P"/storage_coordinator*.jl "$P/model.jl" "$LIB" "problems/mpopf/spar_ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls" > "$OUT/sources.sha256"
# The solver keeps the first allowed CPUs; the sampler sits on the last one.
LAST=$(python3 -c 'import os; print(max(os.sched_getaffinity(0)))')
# The threaded driver runs every worker as a thread of one process.
JT=1; [[ ${SC_DRIVER:-lean} == threaded ]] && JT=$SC_NWORKERS
if [[ ${SC_STAGE_SOLVER:-ipopt} != ipopt ]]; then
 mkdir -p "$OUT/central_conic"
 julia --startup-file=no --project="$PROJECT" "$P/conic_baseline.jl" "$SYS" "$T" "$OUT/central_conic" "$SC_STAGE_SOLVER" > "$OUT/central_conic/solve.log" 2> "$OUT/central_conic/solve.err"
 grep -E "^CONIC" "$OUT/central_conic/solve.log" || true
fi
julia --threads="$JT" --startup-file=no --project="$PROJECT" "$P/$DRIVER" "$SYS" "$T" "$OUT" > "$OUT/solve.log" 2> "$OUT/solve.err" &
PID=$!
taskset -c "$LAST" python3 "$P/sample_memory.py" sample "$PID" "$OUT/memory_samples.csv" "${SC_SAMPLE_S:-2}" &
wait "$PID"
grep -E 'CENTRAL|CASE|SC_LEAN|SC_THREADED|SC_CONVERGED|SC_REPAIR_NO|SC_SUMMARY|SC_WARNING' "$OUT/solve.log" || true
python3 "$P/sample_memory.py" summarize "$OUT/memory_samples.csv" "$OUT/solve.log" | tee "$OUT/memory_summary.txt"
for name in central converged; do
 [[ -s "$OUT/${name}_dispatch.jls" ]] || { echo "OPENDSS $name: no dispatch saved"; continue; }
 SEQUENTIAL_DSS_GZIP=1 julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$OUT/instance.jls" "$OUT/${name}_dispatch.jls" "$OUT/${name}_sequential" model > "$OUT/${name}_opendss.log" 2>&1 || echo "OPENDSS $name: replay failed, see ${name}_opendss.log"
 grep -h SEQUENTIAL_DSS "$OUT/${name}_opendss.log" | sed "s/^/OPENDSS $name: /" || true
done
date -Is > "$OUT/COMPLETE"
