#!/usr/bin/env bash
# OpenDSS feasibility check (opendss_check.jl) of one saved solution, on the
# matched instance (periodic profile, per-system C_B, soft terminal SOC).
# Writes ddp/results/opendss_check/<tag>.log and prints the summary lines.
# OpenDSSDirect's precompile warnings under Julia 1.12 go to the log only.
#
#   bash ddp/examples/power_system/run_opendss_check.sh <system> <T> <solution.jls> <tag> ["label"]
set -u
cd "$(dirname "$0")/../../.." || exit 1
SYS=$1; T=$2; SOL=$3; TAG=$4; LABEL="${5:-$TAG}"
OUT=ddp/results/opendss_check
mkdir -p "$OUT"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
julia --startup-file=no --project=envs/tadmm ddp/examples/power_system/opendss_check.jl \
    "$SYS" "$T" "$SOL" "$LABEL" 2> "$OUT/$TAG.stderr" | grep -E "^OPENDSS_|^  +(t|[0-9]+) " > "$OUT/$TAG.log"
grep -q "OPENDSS_SUMMARY" "$OUT/$TAG.log" && rm -f "$OUT/$TAG.stderr" || { echo "FAILED: see $OUT/$TAG.stderr"; tail -5 "$OUT/$TAG.stderr"; }
grep "OPENDSS_SUMMARY" "$OUT/$TAG.log" | cut -c1-330
