#!/usr/bin/env bash
# Self-contained experiment pipeline for the linear-algebra agenda items.
# Pure Julia underneath; no analysis step needs a human or an LLM. It writes CSVs
# and raw logs, and summarize_agenda_pipeline.jl turns them into one table.
#
#   bash ddp/examples/power_system/run_agenda_pipeline.sh [small|medium|full] [outdir]
#
#   small   ieee123 only                             ~12 min
#   medium  + ieee2522                               ~1 h
#   full    + large10k                               ~6 h   (run overnight)
#
# RESUMABLE: every step is skipped if its log/CSV already exists, so re-running
# after an interrupt picks up where it stopped. Delete a log to force a redo.
# TIMING-SENSITIVE: run on an otherwise idle machine, one pipeline at a time.
#
# Covers:
#   A  diagonal stage Hessian: convergence + time benefit   (agenda "Hessian")
#   B  stale factor as preconditioner: is A reusable?       (agenda "X = Ainv*B")
#   C  UMFPACK vs MUMPS, thread sweep, pivot tolerance,
#      symmetric strategy                                   (agenda "Ax=b solvers")

set -u
cd "$(dirname "$0")/../../.." || exit 1
REPO="$PWD"
SCALE="${1:-small}"
OUT="${2:-$REPO/ddp/results/agenda_pipeline}"
JL="julia --startup-file=no --project=envs/ddp2026"
DRIVER="ddp/examples/power_system/ieee123c_filterddp.jl"

mkdir -p "$OUT/logs" "$OUT/captures"

case "$SCALE" in
  small)  SYSTEMS="ieee123C_1ph" ;;
  medium) SYSTEMS="ieee123C_1ph ieee2522C_1ph" ;;
  full)   SYSTEMS="ieee123C_1ph ieee2522C_1ph large10kC_1ph" ;;
  *) echo "scale must be small|medium|full"; exit 1 ;;
esac

export REDUCED_PROFILE=periodic
export REDUCED_CB=1e-3
export FILTERDDP_SKIP_SOLUTION_WRITE=1

say() { echo "[$(date '+%H:%M:%S')] $*"; }
done_marker() { [ -s "$1" ] && grep -q "PIPELINE_STEP_OK" "$1"; }
finish() { echo "PIPELINE_STEP_OK" >> "$1"; }

# ---------------------------------------------------------------- A: Hessian --
# Floors only swept on the cheapest system; bigger systems use the value that
# worked there. Baseline (exact) is the control arm for each system.
for SYS in $SYSTEMS; do
  for ARM in exact diag1e-8 diag1e-4 diag1e-1; do
    case "$SYS:$ARM" in
      ieee123C_1ph:*) ;;                                  # all four arms
      *:exact|*:diag1e-8) ;;                              # only these two
      *) continue ;;
    esac
    LOG="$OUT/logs/A_${SYS}_T3_${ARM}.log"
    done_marker "$LOG" && { say "skip A $SYS $ARM"; continue; }
    say "A: $SYS T=3 $ARM"
    if [ "$ARM" = exact ]; then
      unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR
    else
      export FILTERDDP_DIAG_HESSIAN=1
      export FILTERDDP_DIAG_HESSIAN_FLOOR="${ARM#diag}"
    fi
    FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_NNZ_DIAGNOSTIC=1 \
      $JL "$DRIVER" "$SYS" 3 solve > "$LOG" 2>&1
    finish "$LOG"
  done
done
unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR

# ------------------------------------------------- B: stale-factor reuse -----
# Needs every iteration's KKT, which is GBs at large10k, so the periodic dump is
# capped at the two smaller systems. The mechanism (Sigma = z/s dominating K) is
# instance-independent; only the step counts are not.
for SYS in $SYSTEMS; do
  [ "$SYS" = large10kC_1ph ] && { say "skip B large10k (capture too large)"; continue; }
  CAPDIR="$OUT/captures/${SYS}_T3"
  LOG="$OUT/logs/B_capture_${SYS}.log"
  if ! done_marker "$LOG"; then
    say "B: capturing per-iteration KKT for $SYS"
    mkdir -p "$CAPDIR"
    FILTERDDP_PERIODIC_CAPTURE_DIR="$CAPDIR" FILTERDDP_PERIODIC_CAPTURE_STRIDE=1 \
      $JL "$DRIVER" "$SYS" 3 solve > "$LOG" 2>&1
    finish "$LOG"
  fi
  ALOG="$OUT/logs/B_analyze_${SYS}.log"
  done_marker "$ALOG" && { say "skip B analyze $SYS"; continue; }
  say "B: stale-factor preconditioner study, $SYS"
  $JL ddp/examples/power_system/stale_factor_preconditioner.jl \
      "$CAPDIR" 1 "$OUT/B_stale_factor_${SYS}.csv" > "$ALOG" 2>&1
  finish "$ALOG"
  rm -rf "$CAPDIR"          # hundreds of MB, and the CSV is the deliverable
done

# ------------------------------------------------------- C: linear solvers ---
# One stage-1 KKT per system, then every solver family at each thread count.
for SYS in $SYSTEMS; do
  CAP="$OUT/captures/${SYS}_T3_stage1.jls"
  LOG="$OUT/logs/C_capture_${SYS}.log"
  if [ ! -s "$CAP" ] && ! done_marker "$LOG"; then
    say "C: capturing stage-1 KKT for $SYS"
    FILTERDDP_CAPTURE_KKT="$CAP" FILTERDDP_CAPTURE_STAGE=1 \
      $JL "$DRIVER" "$SYS" 3 solve > "$LOG" 2>&1
    finish "$LOG"
  fi
  [ -s "$CAP" ] || { say "C: no capture for $SYS, skipping"; continue; }
  for NT in 1 2 4 8 16; do
    SLOG="$OUT/logs/C_${SYS}_threads${NT}.log"
    done_marker "$SLOG" && { say "skip C $SYS threads=$NT"; continue; }
    say "C: $SYS, $NT thread(s)"
    OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT \
    KKT_BENCH_MUMPS_ICNTL16=$NT \
      $JL ddp/examples/power_system/benchmark_kkt_solvers.jl \
          "$CAP" 5 "$OUT/C_solvers_${SYS}_threads${NT}.csv" > "$SLOG" 2>&1
    finish "$SLOG"
  done
done

# ------------------------------------------------------------------ summary --
say "summarising"
$JL ddp/examples/power_system/summarize_agenda_pipeline.jl "$OUT" \
    | tee "$OUT/AGENDA_SUMMARY.txt"
say "done -- results in $OUT"
