#!/usr/bin/env bash
# One Kamiak case of the parallel storage-coordinated tADMM (and its centralized
# reference on the same node). usage: kamiak_storage_coordinator.sh C system T label
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
# Some nodes lack the system OpenBLAS that libma57.so links to; a private copy makes every node usable.
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57
export TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=${TADMM_CB_PROFILE:-july}
case $SYS in
 ieee2522C_1ph) SCREEN=substation,ell;;
 *) SCREEN=substation,vupper,ell;;
esac
export TADMM_SCREEN=${TADMM_SCREEN:-$SCREEN}
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0
export SC_NWORKERS=${SC_NWORKERS:-$T}
{ hostname; date -Is; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(OMP_|OPENBLAS_|JULIA_|TADMM_|SC_)'; } > "$OUT/allocation.txt"
sha256sum "$P"/storage_coordinator*.jl "$P/model.jl" "$LIB" "problems/mpopf/spar_ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls" > "$OUT/sources.sha256"
julia --startup-file=no --project="$P" "$P/storage_coordinator_parallel.jl" "$SYS" "$T" "$OUT" > "$OUT/solve.log" 2> "$OUT/solve.err"
grep -E 'CENTRAL|CASE|SC_NO|SC_REPAIR_NO|SC_CONVERGED|SC_DONE|SC_SUMMARY|SC_WARNING' "$OUT/solve.log"
# Chronological OpenDSS replay of every saved dispatch: live Storage elements,
# energy carried through the horizon, never reset. Outside all solve clocks.
for name in central converged repair_no no; do
 [[ -s "$OUT/${name}_dispatch.jls" ]] || { echo "OPENDSS $name: no dispatch saved"; continue; }
 SEQUENTIAL_DSS_GZIP=1 julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$OUT/instance.jls" "$OUT/${name}_dispatch.jls" "$OUT/${name}_sequential" model > "$OUT/${name}_opendss.log" 2>&1 || echo "OPENDSS $name: replay failed, see ${name}_opendss.log"
 grep -h SEQUENTIAL_DSS "$OUT/${name}_opendss.log" | sed "s/^/OPENDSS $name: /" || true
done
date -Is > "$OUT/COMPLETE"
