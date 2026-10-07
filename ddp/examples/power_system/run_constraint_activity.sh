#!/bin/bash
# Constraint-activity count for one case (constraint_activity.jl), on the
# matched instance used everywhere else: periodic profile, per-system C_B,
# soft terminal SOC. Ipopt runs with MA57 only to get an accurate solution;
# no timing is taken from these runs.
#
#   bash ddp/examples/power_system/run_constraint_activity.sh <system> <T>
HSL="${HSL_DIR:-C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build}"
export PATH="$HSL:$PATH" REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.dll;linear_system_scaling=none"
R=ddp/results/constraint_activity
mkdir -p $R
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/constraint_activity.jl \
    "$1" "$2" "$R/${1}_T${2}_ipoptlog.txt" > "$R/activity_${1}_T${2}.log" 2>&1
grep -E "CENTRAL_IPOPT |CONSTRAINT_ACTIVITY|^  [a-zA-Z]|SCREEN_RULES|ORACLE|VOLTAGE_|BATTERY_IMPLIED|SOC_LOOSE|LOWER_EST" \
    "$R/activity_${1}_T${2}.log" | cut -c1-270
grep -B2 -A12 "ERROR" "$R/activity_${1}_T${2}.log" | head -30 | cut -c1-200
