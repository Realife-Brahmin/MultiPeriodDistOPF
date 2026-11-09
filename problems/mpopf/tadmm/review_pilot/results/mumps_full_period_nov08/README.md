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

See [backend protocol](../../MUMPS_NOV08.md). Initial setup, native unit output,
source/input hashes and submission scripts are preserved in
mumps-initial-protocol.tgz, SHA256:
64ad09e6bff64a0e3b28d00d4e0e971f5243f8234f1b3287e4f8d0dcf598a275.
