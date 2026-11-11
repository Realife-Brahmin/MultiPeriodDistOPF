# Matched MUMPS full-period comparison

Frozen numerical source:386ecf283cf5ad0e85a12ae132d7a99d102f6871.
Native MUMPS gate31656680 passed all248 repeated-solve/multiplier/objective
checks. T48 job31656681 is running with48 period workers, rho120000,
primal-dual warm starts, batched updates, exact local solves and matched
centralized MUMPS reference. The data, model, controller, target, warmup/reset,
Gold hardware, thread policy and128-GiB allocation match the MA57 comparison.
No speedup is claimed before convergence, independent validation and analysis.
Passive energy-audit job31656682 follows the measured run. Inspect T48 before
extending the MUMPS comparison to T96.

## Completed outcome

The numerical solve and all saved network checks finished:322 outer iterations,
635.6532 s tADMM core versus179.6677 s full centralized MUMPS. Both use the
same native backend and quality protocol. The owned/all-copy1e-5 RMS and1e-2
dual goal is met; final physical error remains approximately0.35764 kWh of
integrated energy, as in the MA57 result. MUMPS does not reverse the timing
comparison. The owned1e-3 plus0.5% cost milestone is210.3594 s at95 iterations;
every-copy1e-3 plus0.1% cost milestone is313.7923 s at135 iterations, also slower
than the full central reference.

The original scheduler job is FAILED because the old validator compared the
lowercase backend setting to uppercase `MUMPS 5.9.0` in the completed Ipopt log.
Posthoc validation uses a case-insensitive exact backend-token check. Original
logs, source snapshots and failed job state are retained; numerical solves were
not rerun or altered. The fixed validator independently passes all source,
warmup/measured-path, milestone, worker/status and OpenDSS checks. The driver
did not reach its final metadata-writing block, so campaign_total_wall_s is
absent; use the preserved Slurm accounting, not an invented completed-driver
wall time. Passive audit was rerun as31657409 after posthoc collection.

See [matched MA57/MUMPS plot](analysis-matched/implementation-comparison.png)
and [metrics](analysis-matched/summary.json). Complete/posthoc archive
mumps-completed-posthoc.tgz SHA256:
eb6521443d32becd5539b1833929c5b7a7203ea6ae0621163b67d150764a8a5f.

See [backend protocol](../../MUMPS_NOV08.md). Initial setup, native unit output,
source/input hashes and submission scripts are preserved in
mumps-initial-protocol.tgz, SHA256:
64ad09e6bff64a0e3b28d00d4e0e971f5243f8234f1b3287e4f8d0dcf598a275.
