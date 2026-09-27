#!/usr/bin/env bash
# Centralized Ipopt with HSL MA57 / HSL_MA97 against its default MUMPS, on the
# nine Table II instances at the per-system C_B (agenda of 2026-10-02). The HSL
# libraries are built locally from the user's licensed sources
# (hsl/build_hsl_windows.sh) and loaded through Ipopt's hsllib option; they are
# not in this repository. MC19 is not in the licensed packages, so every run
# uses linear_system_scaling = none (MUMPS included, for a like-for-like
# comparison). MA97 runs on 1 and on 8 OpenMP threads; MUMPS is Ipopt's
# sequential build. One case at a time, background load logged.
#
# BLAS1=1 pins OpenBLAS (which MUMPS, MA57 and MA97 all call) to one thread
# and tags the logs _blas1, so the only parallelism left is MA97's OpenMP.
#
#   bash ddp/examples/power_system/run_ipopt_hsl.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
HSL="${HSL_LIB_DIR:-C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build}"
OUT=ddp/results/ipopt_hsl
mkdir -p "$OUT/logs"
TAG=""
if [ "${BLAS1:-0}" = 1 ]; then export OPENBLAS_NUM_THREADS=1; TAG=_blas1; fi
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export PATH="$HSL:$PATH"
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl

for cell in ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 \
            ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 \
            large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48; do
  S=${cell%%:*}; T=${cell##*:}
  for LS in mumps ma57 ma97 ma97t8; do
    case $LS in
      mumps)  OPTS="linear_system_scaling=none"; OMP=1;;
      ma57)   OPTS="linear_solver=ma57;hsllib=$HSL/libma57.dll;linear_system_scaling=none"; OMP=1;;
      ma97)   OPTS="linear_solver=ma97;hsllib=$HSL/libhsl_ma97.dll;linear_system_scaling=none"; OMP=1;;
      ma97t8) OPTS="linear_solver=ma97;hsllib=$HSL/libhsl_ma97.dll;linear_system_scaling=none"; OMP=8;;
    esac
    LOG="$OUT/logs/ipopt_${LS}${TAG}_${S}_T${T}.log"
    [ -s "$LOG" ] && grep -q "CENTRAL_IPOPT " "$LOG" && { echo "skip $LS $S T=$T"; continue; }
    $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
    STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_${LS}${TAG}_${S}_T${T}.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    IPOPT_EXTRA_OPTIONS="$OPTS" OMP_NUM_THREADS=$OMP \
      $JL --project=envs/ddp2026 ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/ipopt_${LS}${TAG}_${S}_T${T}_ipoptlog.txt" >> "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    echo "[$(date '+%H:%M:%S')] $LS $(grep -oE 'CENTRAL_IPOPT system=\S+ T=[0-9]+ .*status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG")"
  done
done
echo "IPOPT_HSL_DONE"
