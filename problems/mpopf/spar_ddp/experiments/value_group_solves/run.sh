#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1");SYS=$2;T=$3;MODE=$4;cd "$C/source"
OUT="$C/results/${SYS}_T${T}_${MODE}";test ! -e "$OUT";mkdir -p "$OUT"
while IFS= read -r K; do unset "$K"; done < <(compgen -v | grep '^FILTERDDP_' || true)
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_PROCS_WORKERS=$((T==3?2:24)) FILTERDDP_WARMUP=full
export FILTERDDP_TREE_KKT=1 FILTERDDP_NETWORK_BACKEND=tree FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TYPED_EQUATIONS=1
export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_DIRECT_DIAG_HESSIAN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HALF_BLOCK=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_PROCS_GC=pass FILTERDDP_PROCS_BALANCED=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_MAX_ITERATIONS=150 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_SCREEN_ABLATION_ORIGINAL_BOUNDS=1
export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5;[[ $SYS != ieee123C_1ph ]] || export FILTERDDP_OPTIMALITY_TOLERANCE=1e-6
export FILTERDDP_VALUE_GROUPING=$MODE FILTERDDP_VALUE_REFERENCE_DIR="$OUT/captures"
SOL="spar_ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CBsystem.jls"
test ! -e "$SOL"
timeout --signal=TERM --kill-after=30s 15m julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/procs_main.jl "$SYS" "$T" solve > "$OUT/solve.log" 2>&1
test -s "$SOL";mv "$SOL" "$OUT/solution.jls"
python3 spar_ddp/experiments/screening_ablation/validate.py "$OUT/solve.log" "$SYS" > "$OUT/numerical-check.log"
julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/opendss_check.jl "$SYS" "$T" "$OUT/solution.jls" "value_group_$MODE" > "$OUT/opendss.log" 2>&1
python3 spar_ddp/experiments/tree_ablation/validate_opendss.py "$OUT/opendss.log" > "$OUT/opendss-check.log"
julia --startup-file=no spar_ddp/experiments/value_group_solves/validate_physical.jl "$PWD" "$T" "$OUT/solution.jls" "$SYS" > "$OUT/physical.log" 2>&1
env | grep -E '^(FILTERDDP_|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|JULIA_NUM_THREADS|REDUCED_|TERMINAL_SOC_)' | sort > "$OUT/settings.txt"
date -u +%FT%TZ > "$OUT/COMPLETE"
