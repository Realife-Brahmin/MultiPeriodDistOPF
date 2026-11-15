#!/usr/bin/env bash
set -euo pipefail
C=$1
test -s "$C/ANCHOR_GATE_VALIDATED"
[[ $SLURM_JOB_NUM_NODES == 1 && $SLURM_CPUS_PER_TASK == 128 ]]
export TADMM_CPU_MODEL_MATCH=6760P
export TADMM_RHO_POLICY=legacy_phase TADMM_RHO_FLOOR=48 TADMM_LOCAL_SOLVE_POLICY=exact
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 96 run primal_dual 20000 ma57 ieee2522C_1ph
