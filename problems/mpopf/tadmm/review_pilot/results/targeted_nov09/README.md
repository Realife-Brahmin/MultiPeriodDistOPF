# Targeted T48 experiments

Source4a074c167403ba9decb31a578a8df4de1e47f0c0, gate31657404 passed425 tests:
174 forced-interruption/exact-recovery checks,3 fixed-penalty schedule checks,
and248 existing repeated-solve/objective/multiplier checks.

Two separate full-warm/reset MA57 comparisons are running:

-31657405: legacy rho controller with exact primal-only recovery after a warm
 attempt reaches20 Ipopt iterations. Both attempts' work is timed and counted.
-31657407: ordinary exact local solves with high initial rho120000 ramped
 down by5 every5 iterations to48, then frozen; no later penalty cycling.

Both use48 workers, same modified large10k T48 data/weights/physical model,
Gold node64-CPU/128-GiB allocation, same stopping/quality checks and max1000
outer iterations. Passive energy audits31657406/31657408 follow the jobs.
Compare each to the334.2166-s baseline before combining changes or extending
to T96. No performance gain is claimed from the gate or submission.

See [predeclared protocol](../../TARGETED_NOV09.md). Exact setup and native gate
output are preserved in targeted-initial-protocol.tgz, SHA256:
a3407f265ea5b16fbd095d1318d48259a71cf8fcd6bdfa090e92ea72b1bb6b0c.
