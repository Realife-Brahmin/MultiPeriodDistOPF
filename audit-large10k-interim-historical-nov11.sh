#!/usr/bin/env bash
set -euo pipefail
JOBS=/home/aryan.jha/kamiak_jobs
C="$JOBS/tadmm_review_nov11_large10k_interim_ice"
P=problems/mpopf/tadmm/review_pilot
OUT="$C/historical-milestone-audit-$SLURM_JOB_ID"
mkdir "$OUT"
source "$JOBS/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh"
export JULIA_DEPOT_PATH="$C/julia-depot:$JULIA_DEPOT_PATH"
export PYTHONPATH="$JOBS/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}"
export SEQUENTIAL_DSS_GZIP=1
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
for CASE in 48:31680924 72:31680922; do
 T=${CASE%:*}
 J=${CASE#*:}
 R="$C/results/run_T${T}_primal_dual_$J"
 test -s "$R/SWEEP_VALIDATED"
 POINT=$(python3 - "$R" "$OUT/T$T-milestone.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
try:
    import tomllib
except ImportError:
    import tomli as tomllib
arm=Path(sys.argv[1])/'rho15000'
row=tomllib.loads((arm/'milestones.toml').read_text())['paper_cost_1e-3']
point=arm/row['snapshot']
assert hashlib.sha256(point.read_bytes()).hexdigest()==row['snapshot_sha256']
Path(sys.argv[2]).write_text(json.dumps(row,indent=2)+'\n')
print(point)
PY
 )
 julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$R/instance.jls" "$POINT" "$OUT/T$T" model
 python3 "$C/campaign/audit-socp-storage-nov09.py" "$OUT/T$T" test_networks/large10kC_1ph/Storage.dss
done
date -u +%FT%TZ > "$OUT/RECORDED"
