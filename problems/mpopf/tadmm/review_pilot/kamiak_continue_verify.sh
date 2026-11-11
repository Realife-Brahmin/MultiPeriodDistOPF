#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); CHECKER=$(realpath "$2"); cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH="$C/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_CHECK_BACKENDS=ma57,ma97 PROFILE_PERIODIC=1
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
PILOT=problems/mpopf/tadmm/review_pilot
exec > >(tee "$C/logs/continue-verify-$SLURM_JOB_ID.log") 2>&1
sha256sum "$0" "$CHECKER" > "$C/portability-fix.sha256"
{ hostname; julia --version; lscpu; scontrol show job "$SLURM_JOB_ID"; } > "$C/continue-hardware.txt"
grep -q '6438Y+' "$C/continue-hardware.txt"
julia --startup-file=no --project="$PILOT" "$CHECKER" "$C/backend_checks.toml" > "$C/logs/backend-checks-$SLURM_JOB_ID.log" 2>&1
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
 for T in 6 24; do
  julia --startup-file=no --project="$PILOT" problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl "$SYS" "$T" > "$C/logs/export_${SYS}_T${T}.log" 2>&1
 done
done
find problems test_networks -type f -not -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
sha256sum problems/mpopf/spar_ddp/results/network_filterddp/network_data_*.jls > "$C/input.sha256"
date -u +%FT%TZ > "$C/setup.complete"
for SOLVER in ma57 ma97; do
 bash "$PILOT/kamiak_case.sh" "$C" ieee123C_1ph 6 "$SOLVER" "serial_verify_$SLURM_JOB_ID"
done
date -u +%FT%TZ > "$C/serial-verification.complete"
