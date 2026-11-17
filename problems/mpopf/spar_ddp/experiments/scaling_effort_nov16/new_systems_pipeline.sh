#!/usr/bin/env bash
# First runs on the single-substation med3980 and sf6937: export, model checks
# (typed equations bitwise, Jacobian by central differences), Ipopt MA57 at
# T=24 (full model, load-flow start, OpenDSS check) and a sequential DDP solve
# at T=6 (run_paired.sh count, OpenDSS check). Verification, not timing.
#   bash .../new_systems_pipeline.sh [system ...]
set -u
cd "$(dirname "$0")/../../../../.." || exit 1
unset NO_COLOR FORCE_COLOR
P=problems/mpopf/spar_ddp; E=$P/experiments/scaling_effort_nov16; L=$E/first_runs; mkdir -p "$L"
SYSTEMS=${*:-"ieee9500xf_sub1_1ph smartds_p1u_xf_sub1_1ph"}
for S in $SYSTEMS; do
  for T in 6 24; do
    PROFILE_PERIODIC=1 julia --startup-file=no --project=$P $P/examples/power_system/export_ieee123c_data.jl "$S" "$T" > "$L/export_${S}_T$T.log" 2>&1 || { echo "$S T=$T export failed"; continue 2; }
  done
  julia --startup-file=no --project=$P $P/examples/power_system/check_typed_equations.jl "$S" 6 2>&1 | grep -vE "Warning|@ " | tail -2
  julia --startup-file=no --project=$P $E/check_jacobian.jl "$S" 6 2>&1 | grep -E "JACOBIAN_CHECK|ERROR"
  CELLS="$S:24" SOLVERS=ma57 IPOPT_TAG_SUFFIX=_ref IPOPT_LOADFLOW_START=1 IPOPT_MORE_OPTIONS="tol=1e-6;max_iter=1000" \
    bash $P/examples/power_system/run_ipopt_hsl.sh 2>&1 | grep -E "OPENDSS_SUMMARY|no dispatch" | cut -c1-260
  grep -hoE "CENTRAL_IPOPT system=\S+ T=[0-9]+ .*status=\S+ iterations=[0-9]+ objective=[-0-9.eE+]+ solve_time_s=[0-9.]+" $P/results/ipopt_hsl/logs/ipopt_ma57_ref_${S}_T24.log
  bash $P/examples/power_system/run_paired.sh "$S" 6 count count_first 2>&1 | grep -E "OPENDSS_SUMMARY|solve complete|OPENDSS_CHECK" | cut -c1-260
done
