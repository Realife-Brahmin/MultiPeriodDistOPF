#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); BASE=$(realpath "$2"); cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH="$BASE/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
PILOT=problems/mpopf/tadmm/review_pilot
exec > >(tee "$C/logs/parallel-gate-$SLURM_JOB_ID.log") 2>&1
test -s "$BASE/serial-verification.complete"
{ hostname; lscpu; julia --version; scontrol show job "$SLURM_JOB_ID"; } > "$C/gate-hardware.txt"
grep -q '6438Y+' "$C/gate-hardware.txt"
julia --startup-file=no --project="$PILOT" "$PILOT/setup.jl" > "$C/logs/setup.log" 2>&1
mkdir -p problems/mpopf/spar_ddp/results/network_filterddp
cp "$BASE/source/problems/mpopf/spar_ddp/results/network_filterddp/"network_data_*_periodic.jls problems/mpopf/spar_ddp/results/network_filterddp/
find problems test_networks -type f -not -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
sha256sum problems/mpopf/spar_ddp/results/network_filterddp/network_data_*.jls > "$C/input.sha256"
# First prove deterministic serial and process-parallel execution on the same
# frozen instance. Each run is fully warmed and independently OpenDSS checked.
export TADMM_NWORKERS=0
bash "$PILOT/kamiak_case.sh" "$C" ieee123C_1ph 6 ma57 serial_gate
export TADMM_NWORKERS=2
for SOLVER in ma57 ma97; do
 bash "$PILOT/kamiak_case.sh" "$C" ieee123C_1ph 6 "$SOLVER" parallel_gate
done
python3 "$PILOT/compare_gate.py" "$C" > "$C/gate-comparison.json"
date -u +%FT%TZ > "$C/parallel-gate.complete"
