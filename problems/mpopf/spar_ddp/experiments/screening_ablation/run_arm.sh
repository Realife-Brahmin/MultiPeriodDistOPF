#!/usr/bin/env bash
set -euo pipefail
SYS=$1; T=$2; REP=$3; MASK=$4; OUT=$(realpath -m "$5")
ROOT=$(cd "$(dirname "$0")/../../.." && pwd); cd "$ROOT"
[[ -f .screening-ablation-private && ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
WORKERS=${SCREEN_WORKERS:-63}
while IFS= read -r NAME; do unset "$NAME"; done < <(compgen -v FILTERDDP_)
export FILTERDDP_PROCS_WORKERS=$WORKERS
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_NETWORK_BACKEND=tree FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TYPED_EQUATIONS=1
export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_DIRECT_DIAG_HESSIAN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HALF_BLOCK=1
export FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1 FILTERDDP_PROCS_GC=pass
export FILTERDDP_PROCS_BALANCED=1 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_WARMUP=full
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_SCREEN_ABLATION_ORIGINAL_BOUNDS=1 FILTERDDP_SCREEN=substation
(( MASK & 1 )) && FILTERDDP_SCREEN+=,vupper
(( MASK & 2 )) && FILTERDDP_SCREEN+=,ell
(( MASK & 4 )) && FILTERDDP_SCREEN+=,psubs
case "$SYS" in ieee123C_1ph) export FILTERDDP_OPTIMALITY_TOLERANCE=1e-6;;
 ieee2522C_1ph|large10kC_1ph) export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5;; *) exit 2;; esac
printf 'SCREEN_ABLATION_RUN system=%s T=%s repeat=%s mask=%s screen=%s\n' "$SYS" "$T" "$REP" "$MASK" "$FILTERDDP_SCREEN" > "$OUT/solve.log"
env | grep -E '^(FILTERDDP_|REDUCED_|TERMINAL_SOC_|JULIA_NUM_THREADS=|OPENBLAS_NUM_THREADS=|OMP_NUM_THREADS=|MKL_NUM_THREADS=)' | LC_ALL=C sort > "$OUT/environment.txt"
SOL="spar_ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CBsystem.jls"
[[ ! -e "$SOL" ]] || mv "$SOL" "$OUT/previous-unclaimed-solution.jls"
timeout --signal=TERM --kill-after=30s 25m julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/procs_main.jl "$SYS" "$T" solve >> "$OUT/solve.log" 2>&1
[[ -s "$SOL" ]] || exit 3
mv "$SOL" "$OUT/solution.jls"
python3 spar_ddp/experiments/screening_ablation/validate.py "$OUT/solve.log" "$SYS"
julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/opendss_check.jl "$SYS" "$T" "$OUT/solution.jls" screening_ablation > "$OUT/opendss.log" 2>&1
python3 spar_ddp/experiments/tree_ablation/validate_opendss.py "$OUT/opendss.log"
date -u +%FT%TZ > "$OUT/arm.complete"
