#!/usr/bin/env bash
# T=24 sweet-spot probe: hold rho=16000, drop tau_decr from 4 to 3.
# Hypothesis: smaller ρ-decay steps keep ρ in the "good" range longer,
# avoiding over-shoot into the cheap-but-imprecise regime.
# Writes to rho_sweep_taudecr/ to preserve the original rho_sweep/rho_16000 winner.
set -euo pipefail

cd /c/repos_addendum/MultiPeriodDistOPF

run_sweep() {
  local T="$1"
  local RHO="$2"
  local TAUDECR="$3"
  local DIRNAME="rho_sweep_taudecr${TAUDECR}"
  echo "=================================================================="
  echo "[$(date +%H:%M:%S)] Starting T=${T} rho=${RHO} tau_decr=${TAUDECR} (-> ${DIRNAME})"
  echo "=================================================================="
  SYSTEM_OVERRIDE=ieee2552C_1ph \
  SWEEP_T="${T}" \
  SWEEP_RHOS="${RHO}" \
  SWEEP_DIR_NAME="${DIRNAME}" \
  ADAPTIVE_RHO_OVERRIDE=true \
  EPS_PRI_OVERRIDE=1e-4 \
  EPS_DUAL_OVERRIDE=5e-4 \
  TAU_INCR_OVERRIDE=2 \
  TAU_DECR_OVERRIDE="${TAUDECR}" \
  MU_BALANCE_OVERRIDE=5 \
    julia --project=envs/tadmm --threads=16 run_rho_sweep.jl
  echo "[$(date +%H:%M:%S)] Done T=${T} rho=${RHO} tau_decr=${TAUDECR}"
}

run_sweep 24 16000 3

echo "=================================================================="
echo "[$(date +%H:%M:%S)] tau_decr probe complete."
echo "=================================================================="
echo
echo "Comparison vs original rho_sweep/rho_16000 (tau_decr=4, NO=283.87s, k=20):"
DIR=envs/tadmm/processedData/ieee2552C_1ph_T24
CSV="${DIR}/rho_sweep_taudecr3/rho_16000/near_optimal_summary.csv"
if [ -f "${CSV}" ]; then
  printf "  tau_decr=3  "
  awk -F, 'NR==2 {printf "NO_time=%-10s k=%-4s gap=%s%%\n", $2, $1, $4}' "${CSV}"
fi
