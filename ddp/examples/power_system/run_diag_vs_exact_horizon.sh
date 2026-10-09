#!/usr/bin/env bash
# Does the diagonal stage Hessian degrade as the horizon grows? ieee123 is
# small enough for the exact stage Hessian at very long horizons, so run both
# to FilterDDP's strict tolerance (no near-optimality stop) and compare
# iteration counts and outcomes. UMFPACK + blocked solve, typed residuals,
# structured dynamics, per-system C_B; timing is not the point.
#
#   bash ddp/examples/power_system/run_diag_vs_exact_horizon.sh [system] [horizons]

set -u
cd "$(dirname "$0")/../../.." || exit 1
SYS="${1:-ieee123C_1ph}"; HORIZONS="${2:-96 384 1536}"
OUT=ddp/results/diag_vs_exact_horizon
mkdir -p "$OUT"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 FILTERDDP_MAX_ITERATIONS=400
export FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_TIMING_DIAGNOSTIC=1
export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_FACTOR_BACKED_POLICY=1
export FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_SKIP_SOLUTION_WRITE=1
JL="julia --startup-file=no"
for T in $HORIZONS; do
  DATA="ddp/results/network_filterddp/network_data_${SYS}_T${T}_periodic.jls"
  if [ ! -s "$DATA" ]; then
    PROFILE_PERIODIC=1 $JL --project=envs/tadmm ddp/examples/power_system/export_ieee123c_data.jl "$SYS" "$T" \
      > "$OUT/export_${SYS}_T${T}.log" 2>&1
    [ -s "$DATA" ] || { echo "export failed for $SYS T=$T"; continue; }
  fi
  for ARM in diag exact; do
    LOG="$OUT/fddp_${SYS}_T${T}_${ARM}.log"
    grep -q "solve complete" "$LOG" 2>/dev/null && { echo "skip $SYS T=$T $ARM"; continue; }
    if [ "$ARM" = diag ]; then
      FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 FILTERDDP_DIRECT_DIAG_HESSIAN=1 \
        $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve > "$LOG" 2>&1
    else
      $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve > "$LOG" 2>&1
    fi
    echo "[$(date '+%H:%M:%S')] $SYS T=$T $ARM: $(grep -E 'solve complete' "$LOG" | cut -c1-70) | $(grep -oE 'FilterDDP objective=\S+' "$LOG") | $(grep -oE 'final residuals: primal=\S+ dual=\S+' "$LOG")"
  done
done
echo "DIAG_VS_EXACT_DONE"
