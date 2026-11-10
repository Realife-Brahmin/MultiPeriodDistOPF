#!/usr/bin/env bash
set -euo pipefail
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_caps_sequence
J=$1
OUT="$C/results/run_T96_primal_dual_$J"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
test -s "$OUT/sweep.toml"
for REPLAY in central_sequential rho20000_endpoint_sequential;do
 python3 "$C/campaign/audit-socp-storage-nov09.py" "$OUT/$REPLAY" "$C/source/test_networks/ieee2522C_1ph/Storage.dss"
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
