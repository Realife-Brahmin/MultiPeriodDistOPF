#!/usr/bin/env bash
set -euo pipefail
C=${TADMM_COLLECTION_CAMPAIGN:-/home/aryan.jha/kamiak_jobs/tadmm_review_nov10_large10k_t48_tuning_v3}
T=${TADMM_COLLECTION_HORIZON:-48}
cd "$C"
for J in "$@"; do
 [[ $J =~ ^[0-9]+$ ]]
 R="results/run_T${T}_primal_dual_$J"
 test -s "$R/EXPLICIT_SOCP_STORAGE_MEMORY_AUDIT_RECORDED"
 test -s "$R/SWEEP_VALIDATED"
 A=$(awk -F '\t' -v job="$J" '$(NF-1) == job { print $NF }' campaign/submissions.tsv)
 [[ $A =~ ^[0-9]+$ ]]
 FILES=(campaign/submissions.tsv campaign/run.sh campaign/posthoc.sh input.sha256 execution-source.sha256
        logs/solve-"$J".* logs/audit-"$A".*
        "$R"/{allocation.txt,base_commit.txt,source.sha256,sweep.toml,solve.log,solve.err,phase_events.csv,central_trace.csv}
        "$R"/central_sequential/{summary.toml,socp-storage-audit.json}
        "$R"/rho*_endpoint_sequential/{summary.toml,socp-storage-audit.json}
        "$R"/rho*/{assessment.toml,milestones.toml,trace.csv,trace_subproblems.csv,round-timing.csv}
        "$R"/rho*/warmup/trace.csv "$R"/process_samples.*.csv
        "$R"/{posthoc-memory-peaks.json,posthoc-slurm-memory.txt,validation_complete.json})
 if [[ -f logs/controller-gate.log ]]; then FILES+=(logs/controller-gate.log); fi
 if [[ -f campaign/profile-job.txt ]]; then FILES+=(campaign/profile-job.txt); fi
 tar -czf "campaign/tuning-$J.tgz" "${FILES[@]}"
 sha256sum "campaign/tuning-$J.tgz"
done
