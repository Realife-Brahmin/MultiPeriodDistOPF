#!/usr/bin/env bash
# Sweep rho around the T=24 winner (16000) hoping to lower NO eff_time below 283.87s.
# Trying rho ∈ {12000, 14000, 18000, 20000}. Same hyperparams as prior T=24 sweep.
set -euo pipefail

cd /c/repos_addendum/MultiPeriodDistOPF

run_sweep() {
  local T="$1"
  local RHO="$2"
  echo "=================================================================="
  echo "[$(date +%H:%M:%S)] Starting T=${T} rho=${RHO}"
  echo "=================================================================="
  SYSTEM_OVERRIDE=ieee2522C_1ph \
  SWEEP_T="${T}" \
  SWEEP_RHOS="${RHO}" \
  ADAPTIVE_RHO_OVERRIDE=true \
  EPS_PRI_OVERRIDE=1e-4 \
  EPS_DUAL_OVERRIDE=5e-4 \
  TAU_INCR_OVERRIDE=2 \
  TAU_DECR_OVERRIDE=4 \
  MU_BALANCE_OVERRIDE=5 \
    julia --project=envs/tadmm --threads=16 envs/tadmm/root_level/run_rho_sweep.jl
  echo "[$(date +%H:%M:%S)] Done T=${T} rho=${RHO}"
}

for RHO in 12000 14000 18000 20000; do
  run_sweep 24 "${RHO}"
done

echo "=================================================================="
echo "[$(date +%H:%M:%S)] All T=24 sweet-spot probe runs complete."
echo "=================================================================="
echo
echo "Summary (NO eff_time per rho):"
DIR=envs/tadmm/processedData/ieee2522C_1ph_T24/rho_sweep
for RHO in 4000 8000 12000 14000 16000 18000 20000 24000 32000; do
  CSV="${DIR}/rho_${RHO}/near_optimal_summary.csv"
  if [ -f "${CSV}" ]; then
    T=$(awk -F, 'NR==2 {print $2}' "${CSV}")
    K=$(awk -F, 'NR==2 {print $1}' "${CSV}")
    G=$(awk -F, 'NR==2 {print $4}' "${CSV}")
    printf "  rho=%-6s  NO_time=%-10s k=%-4s gap=%s%%\n" "${RHO}" "${T}" "${K}" "${G}"
  fi
done
