#!/usr/bin/env bash
# T=24 probe: tau_decr=5 across rho ∈ {14000, 16000, 20000}.
# Hypothesis: bigger ρ drops in Phase 2 cut per-iter Ipopt cost faster than
# tau_decr=4, possibly shaving below 283.87s if ρ doesn't over-shoot.
# Writes to rho_sweep_taudecr5/.
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
  SYSTEM_OVERRIDE=ieee2522C_1ph \
  SWEEP_T="${T}" \
  SWEEP_RHOS="${RHO}" \
  SWEEP_DIR_NAME="${DIRNAME}" \
  ADAPTIVE_RHO_OVERRIDE=true \
  EPS_PRI_OVERRIDE=1e-4 \
  EPS_DUAL_OVERRIDE=5e-4 \
  TAU_INCR_OVERRIDE=2 \
  TAU_DECR_OVERRIDE="${TAUDECR}" \
  MU_BALANCE_OVERRIDE=5 \
    julia --project=envs/tadmm --threads=16 envs/tadmm/root_level/run_rho_sweep.jl
  echo "[$(date +%H:%M:%S)] Done T=${T} rho=${RHO} tau_decr=${TAUDECR}"
}

for RHO in 14000 16000 20000; do
  run_sweep 24 "${RHO}" 5
done

echo "=================================================================="
echo "[$(date +%H:%M:%S)] tau_decr=5 probe complete."
echo "=================================================================="
echo
echo "Comparison (baseline tau_decr=4):"
echo "  rho=14000 tau_decr=4 -> 307.79s k=21 gap 0.27%"
echo "  rho=16000 tau_decr=4 -> 283.87s k=20 gap 0.26% (winner)"
echo "  rho=20000 tau_decr=4 -> 379.82s k=20 gap 0.23%"
echo
echo "New results (tau_decr=5):"
DIR=envs/tadmm/processedData/ieee2522C_1ph_T24/rho_sweep_taudecr5
for RHO in 14000 16000 20000; do
  CSV="${DIR}/rho_${RHO}/near_optimal_summary.csv"
  if [ -f "${CSV}" ]; then
    printf "  rho=%-6s tau_decr=5  " "${RHO}"
    awk -F, 'NR==2 {printf "NO_time=%-10s k=%-4s gap=%s%%\n", $2, $1, $4}' "${CSV}"
  fi
done
