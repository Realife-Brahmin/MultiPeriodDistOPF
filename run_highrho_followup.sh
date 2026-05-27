#!/usr/bin/env bash
# T=24 sweet-spot probe at rho=24000 (between 16000 winner with ramp-up and 32000 decay-only-but-slow).
set -euo pipefail

cd /c/repos_addendum/MultiPeriodDistOPF

run_sweep() {
  local T="$1"
  local RHO="$2"
  echo "=================================================================="
  echo "[$(date +%H:%M:%S)] Starting T=${T} rho=${RHO}"
  echo "=================================================================="
  SYSTEM_OVERRIDE=ieee2552C_1ph \
  SWEEP_T="${T}" \
  SWEEP_RHOS="${RHO}" \
  ADAPTIVE_RHO_OVERRIDE=true \
  EPS_PRI_OVERRIDE=1e-4 \
  EPS_DUAL_OVERRIDE=5e-4 \
  TAU_INCR_OVERRIDE=2 \
  TAU_DECR_OVERRIDE=4 \
  MU_BALANCE_OVERRIDE=5 \
    julia --project=envs/tadmm --threads=16 run_rho_sweep.jl
  echo "[$(date +%H:%M:%S)] Done T=${T} rho=${RHO}"
}

run_sweep 24 24000
echo "[$(date +%H:%M:%S)] All runs complete."
