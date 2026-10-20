#!/usr/bin/env bash
# Runs only in a private snapshot (created by submit.sh). Full solver convergence.
set -euo pipefail
[[ $# -ge 4 ]] || { echo 'usage: run_pair.sh SYSTEM T REPEAT OUTDIR [tree|umfpack first]' >&2; exit 2; }
SYS=$1; HORIZON=$2; REP=$3; OUT=$(realpath -m "$4"); FIRST=${5:-tree}
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT"
[[ -f .tree-ablation-private ]] || { echo 'Refusing to write solutions outside an isolated ablation snapshot' >&2; exit 2; }
mkdir -p "$OUT"
[[ ! -e "$OUT/pair.complete" ]] || { echo 'Pair already completed; refusing overwrite' >&2; exit 2; }
# Clear inherited experimental solver toggles, preserving only explicit run controls.
WANT_WORKERS=${FILTERDDP_PROCS_WORKERS:-63}; WANT_WARMUP=${FILTERDDP_WARMUP:-full}
while IFS= read -r NAME; do unset "$NAME"; done < <(compgen -v FILTERDDP_)
export FILTERDDP_PROCS_WORKERS=$WANT_WORKERS
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
export REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TYPED_EQUATIONS=1
export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_DIRECT_DIAG_HESSIAN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HALF_BLOCK=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_PROCS_GC=pass FILTERDDP_PROCS_BALANCED=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_WARMUP=$WANT_WARMUP
case "$SYS" in ieee123C_1ph) export FILTERDDP_OPTIMALITY_TOLERANCE=1e-6;;
 ieee2522C_1ph|large10kC_1ph) export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5;; *) exit 2;; esac
unset FILTERDDP_DIAG_HESSIAN FILTERDDP_NEAR_OPT_REFERENCE FILTERDDP_NEAR_OPT_GAP FILTERDDP_NEAR_OPT_PRIMAL
unset FILTERDDP_SKIP_SOLUTION_WRITE TREE_ABLATION_CAPTURE_DIR
PTAG=${REDUCED_PROFILE:+_$REDUCED_PROFILE}
SOL="spar_ddp/results/network_filterddp/filterddp_solution_${SYS}_T${HORIZON}${PTAG}_CBsystem.jls"
if [[ "$FIRST" == tree ]]; then ARMS=(tree umfpack); elif [[ "$FIRST" == umfpack ]]; then ARMS=(umfpack tree); else exit 2; fi
for ARM in "${ARMS[@]}"; do
    [[ ! -e "$OUT/$ARM.log" && ! -e "$OUT/$ARM.jls" ]] || { echo 'Existing arm output; choose a new OUTDIR' >&2; exit 2; }
    rm -f "$SOL"  # private snapshot only; never a shared campaign solution
    export FILTERDDP_NETWORK_BACKEND=$ARM
    printf 'TREE_ABLATION_RUN system=%s T=%s repeat=%s backend=%s source=%s\n' "$SYS" "$HORIZON" "$REP" "$ARM" "${TREE_ABLATION_SOURCE_SHA:-local}" > "$OUT/$ARM.log"
    julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/procs_main.jl "$SYS" "$HORIZON" solve >> "$OUT/$ARM.log" 2>&1
    [[ -s "$SOL" ]] || { echo "Missing $ARM solution" >&2; exit 1; }
    cp "$SOL" "$OUT/$ARM.jls"
    # Independent original-model power-flow check, outside optimizer timing.
    julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/opendss_check.jl \
        "$SYS" "$HORIZON" "$OUT/$ARM.jls" "tree_ablation_${ARM}" > "$OUT/$ARM.opendss.log" 2>&1
    python3 "spar_ddp/experiments/tree_ablation/validate_opendss.py" "$OUT/$ARM.opendss.log"
done
julia --startup-file=no --project=envs/spar_ddp spar_ddp/experiments/tree_ablation/validate_pair.jl \
 "$OUT/tree.jls" "$OUT/umfpack.jls" "$OUT/tree.log" "$OUT/umfpack.log" | tee "$OUT/pair.validation.log"
printf 'completed %s\n' "$(date -u +%FT%TZ)" > "$OUT/pair.complete"
