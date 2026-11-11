#!/usr/bin/env bash
set -euo pipefail
J=$1
C=${2:-/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_caps_sequence}
T=${3:-96}
RHO=${4:-20000}
SYSTEM=${5:-ieee2522C_1ph}
OUT="$C/results/run_T${T}_primal_dual_$J"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
test -s "$OUT/sweep.toml"
RATINGS="$C/source/test_networks/$SYSTEM/Storage.dss"
python3 "$C/campaign/audit-socp-storage-nov09.py" "$OUT/central_sequential" "$RATINGS"
IFS=, read -r -a RHO_ARMS <<< "$RHO"
for INITIAL in "${RHO_ARMS[@]}"; do
 python3 "$C/campaign/audit-socp-storage-nov09.py" "$OUT/rho${INITIAL}_endpoint_sequential" "$RATINGS"
 python3 "$C/campaign/summarize-tadmm-round-timing-nov10.py" "$OUT/rho$INITIAL"
done
python3 - "$OUT" <<'PY'
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1]);peaks={};samples=0
for f in p.glob('process_samples.*.csv'):
 rows=list(csv.DictReader(f.open()));samples+=len(rows)
 peaks[f.name]={key:max((float(r[key]) for r in rows),default=0) for key in ('working_set_bytes','private_bytes','pss_bytes','process_count')}
assert peaks and samples>0,'Missing process memory samples'
(p/'posthoc-memory-peaks.json').write_text(json.dumps({'node_samples':peaks,'samples':samples,'scope':'Registered coordinator/workers; sampled one second; RSS sums include shared pages, PSS apportions them; these are sampled peaks, not exact allocation maxima.'},indent=2)+'\n')
PY
sacct -j "$J" --units=M --format=JobID,State,Elapsed,MaxRSS,MaxVMSize,ReqMem,AllocCPUS,NodeList -P > "$OUT/posthoc-slurm-memory.txt"
date -u +%FT%TZ > "$OUT/EXPLICIT_SOCP_STORAGE_MEMORY_AUDIT_RECORDED"
