#!/usr/bin/env bash
# Re-run BF for ieee2552C T=24 (existing 287.8s may be lucky-fast outlier).
# Save to a side directory so original isn't clobbered, then we compare.
set -euo pipefail

cd /c/repos_addendum/MultiPeriodDistOPF

echo "[$(date +%H:%M:%S)] Starting T=24 BF re-run"

SYSTEM_OVERRIDE=ieee2552C_1ph \
T_OVERRIDE=24 \
  julia --project=envs/tadmm --threads=16 envs/tadmm/root_level/run_bf.jl

echo "[$(date +%H:%M:%S)] Done T=24 BF re-run"
