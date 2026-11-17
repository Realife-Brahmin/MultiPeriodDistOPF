#!/usr/bin/env bash
# The largest load scale lambda at which one substation's OPF solves: for each
# lambda (highest first) scale the deck (scale_deck.py), export T, and solve the
# full model (no screening) with Ipopt MA57 from the load-flow start. Stops at
# the first lambda Ipopt solves.
#   bash .../lambda_search.sh <variant deck> <source deck> <T> <lambda ...>
set -u
cd "$(dirname "$0")/../../../../.." || exit 1
VAR=$1; SRC=$2; T=$3; shift 3
E=problems/mpopf/spar_ddp/experiments/scaling_effort_nov16
HSL="${HSL_LIB_DIR:-C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build}"; export PATH="$HSL:$PATH"
OUT=$E/lambda_search; mkdir -p "$OUT"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 IPOPT_LOADFLOW_START=1
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.dll;linear_system_scaling=none;tol=1e-6;max_iter=600"
unset IPOPT_SCREEN IPOPT_SAVE_DISPATCH
for L in "$@"; do
  python "$E/scale_deck.py" "test_networks/$VAR" "test_networks/$SRC" "$L" > /dev/null || exit 1
  PROFILE_PERIODIC=1 julia --startup-file=no --project=problems/mpopf/spar_ddp problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl "$VAR" "$T" > "$OUT/export_${VAR}_T${T}_l$L.log" 2>&1 || { echo "$VAR lambda=$L: export failed"; continue; }
  LOG="$OUT/ipopt_${VAR}_T${T}_l$L.log"
  julia --startup-file=no --project=problems/mpopf/spar_ddp problems/mpopf/spar_ddp/examples/power_system/centralized_ipopt_matched.jl "$VAR" "$T" > "$LOG" 2>&1
  R=$(grep -oE "status=[A-Z_]+ iterations=[0-9]+ objective=[-0-9.eE+]+ solve_time_s=[0-9.]+" "$LOG" | head -1)
  echo "$VAR lambda=$L T=$T: ${R:-failed: $(grep -m1 -oE 'ERROR: .{0,120}' "$LOG")}"
  echo "$R" | grep -qE "status=(LOCALLY_SOLVED|ALMOST_LOCALLY_SOLVED)" && break
done
