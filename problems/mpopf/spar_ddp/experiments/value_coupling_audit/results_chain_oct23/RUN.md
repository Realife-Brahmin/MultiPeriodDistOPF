# Compact-to-compact bounded test

Job31368962, overlay eb4f0439, snapshot
/home/aryan.jha/kamiak_jobs/value_compact_chain_oct23.
Branch ddp/oct21-battery-count-study, base ef43e34e.
Only the existing original1020 large-late four-period window; no optimizer.
Gold6438Y+64CPU300GiB,35min batch/25min process cap, no-requeue, occupancy guard.
Private first depot/read-only caches; no other jobs active at submission.
Synthetic compact/component and compact/chain unit gates passed locally and
on compute node before loading the immutable capture.

See COMPACT_CHAIN.md for semantics, charged costs, growing retained-interior
representation, local rank-growth fallback and direction/negative-control gates.
Job COMPLETED0:0 in1m28, cn164. Both unit gates, actual capture/source/input
checks, independent direction-consumer check and all3alternating warmed pairs
validated. Strict analyze_compact_chain.py collector also recomputes overlay
checksums and checks Slurm accounting, all value/direction errors, negative
control and unchanged scope. See FINDINGS.md and analysis/validation.json.

Raw evidence.tgz SHA256:
2010fcc23d3b96c5ebf6dbddd503194851a53c2d63b1d2805fe89c3b75024f06.
Archive contains actual overlay source, source/input/inherited check logs,
hardware/affinity/allocation, all raw CSV, stdout/stderr, units and accounting;
no private depot or duplicated full capture/numerical tree.
