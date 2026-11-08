#!/usr/bin/env bash
# Read completed snapshots only; this is not a numerical solve benchmark.
set -euo pipefail
C=$(realpath "$1"); A=$(realpath "$2"); shift 2
cd "$A"
sha256sum --check --quiet source.sha256
export TADMM_AUDIT_LAYOUT="$C/source/problems/mpopf/spar_ddp/examples/power_system/control_layout.jl"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export TADMM_AUDIT_VALIDATED_ONLY=1
sha256sum "$0" "$TADMM_AUDIT_LAYOUT" > "orchestration-$SLURM_JOB_ID.sha256"
for CASE in "$@"; do
 T=${CASE%:*}; J=${CASE#*:}
 INPUT="$C/results/large10k_T${T}_$J"
 if [[ ! -s "$INPUT/SWEEP_VALIDATED" ]]; then
  compgen -G "$INPUT/validation_rho*.json" > /dev/null
 fi
 /home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no "$A/audit_sweep_energy.jl" "$INPUT" "$A/T${T}_${J}-all-points.toml"
done
