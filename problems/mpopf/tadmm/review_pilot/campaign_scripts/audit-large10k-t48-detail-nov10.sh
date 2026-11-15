#!/usr/bin/env bash
set -euo pipefail
JOBS=/home/aryan.jha/kamiak_jobs
AUDIT="$JOBS/large10k-t48-detail-nov10"
export PYTHONPATH="$JOBS/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}"
for CASE in cap10:31679860 exact:31679947; do
 POLICY=${CASE%:*}
 J=${CASE#*:}
 C="$JOBS/tadmm_review_nov10_large10k_t48_$POLICY"
 OUT="$C/results/run_T48_primal_dual_$J"
 test -s "$OUT/EXPLICIT_SOCP_STORAGE_MEMORY_AUDIT_RECORDED"
 for REPLAY in central_sequential rho120000_endpoint_sequential; do
  python3 "$AUDIT/audit-socp-storage-nov09.py" "$OUT/$REPLAY" "$C/source/test_networks/large10kC_1ph/Storage.dss" physical-error-audit.json
 done
 python3 "$AUDIT/summarize-tadmm-round-timing-nov10.py" "$OUT/rho120000"
 date -u +%FT%TZ > "$OUT/DETAILED_PHYSICAL_AND_TIMING_AUDIT_RECORDED"
done
