# October23 bounded follow-up provenance

Continues25dc3e46. Read FOLLOWUP.md for the three-point scope, explicit gates,
timer boundaries and dense-expansion limitations. No production/paper/model
or benchmark result changes. Local synthetic tests use Julia1.12.7; compute
nodes use1.12.4. Every snapshot is immutable and previous attempts remain.

First snapshot value_block_followup_oct23, sourcee5e4497a, job31368638 FAILED1:0
in4s: unit setup sought hierarchy.jl beside the overlay instead of the read-only
base. No frozen point ran. Archive setup_failed_evidence.tgz SHA256
79e918658c229b14ec2fb31ec8b72352fe621904892c950b69237773879b0017.

Replacement value_block_followup_oct23_v2, source8ba23465, job31368651 FAILED1:0
in34s. Unit gate passed; medium-late dense saved probes exactly0 and ordinary
value error1.98218e-12 passed. New measurement function then failed before any
approximation: nested norm accumulator reused `b`, overwriting the enclosing
battery-path collection via Julia closure capture. Archive
measurement_failed_evidence.tgz SHA256
8e02798fd59d985772832c7905b971e86665f0d84d3d3d5f4396a28cf1e2036a.

Repair only renames norm accumulators and adds a mocked regression executing
the actual measurement function for battery paths and all four full direction
fields. Numerical compression, physical states and accuracy gates unchanged.
The separate original unit decimal equality used approximate equality for
1e-6/sqrt(1e8) before initial submission; this is not a scientific threshold.

All attempts: Gold6438Y+,64CPU300GiB, occupancy guard, private first depot,
BLAS/OMP/Julia1, no requeue,75minute batch /22minute point caps. Read-only binary
states and old numerical source are reused. No other jobs at submission and
home5.3GiB free; no copies of full sources or captures, no deleted old work.
Final collection must verify all three points and exactly three alternating
same-input warm pairs for each qualifying candidate; partial is not complete.

Current corrected snapshot value_block_followup_oct23_v3 uses sourcef8f60d39,
job31368668. Both synthetic representation and actual measurement-function
regressions passed on the compute node. It reuses identical frozen points;
earlier attempts never reached an approximate scan. Read-only inherited source
checks and overlay SHA manifests are recorded independently.

Job31368668 then hit its22minute medium-point PROCESS cap: FAILED9:0 at22m30,
not the75minute Slurm cap, preemption, OOM or an accuracy assertion. The stack
was in LLVM instruction scheduling. No point metadata or warmed measurements
completed. Medium-late block1e-6 and conditioned1e-6 controls finished with
scaled errors2.38152e-7 /2.23660e-7 and payload ratios0.19216 /0.20918; both
missed the unchanged1e-7 gate. Conditioned1e-9 was incomplete; neither large
nor early control started. No accuracy or runtime conclusion for tighter
settings is supported by this interruption. ALL partial outputs/source/logs
are in compiler_interrupted_evidence.tgz, SHA256
0a559f039eb1b51178e56060f5c33d8ac852b7f47e8881006f8c4d716f17f96b.

Numerically equivalent fixed recursive node payload code was prepared without
touching the running snapshot. It bounds compiler specialization on entire
NamedTuple tree/dense-choice types. Local synthetic comparison to pinnedf8f60d39
passes bit-identical matrices, ranks, error/fallback and coefficient counts for
all four budgets. New snapshot must also pass actual medium-late first-leaf
equivalence at0/1e-6 before the same three-point replay. No arithmetic order,
SVD budget, direction gate, physical input or case scope is changed. These
compiler-heavy exploratory observations remain outside warmed benchmark data.

Replacement value_block_followup_oct23_v4, overlaye2a86a50, job31368755,
was submitted ONLY after31368668 stopped. It keeps the same three points,
hardware/memory/thread policy and time caps. A separate legacy-source hash
check and synthetic/real-leaf bit-identical equivalence gate precede replay.
Only payload container implementation changes; no new optimizer is queued.

Final:31368755 COMPLETED0:0 in3m32, all three frozen points strictly validated.
Large late and medium early pass conditioned1e-6 and each complete3alternating
warmed pairs; medium late rejects every tested candidate and has no warmed
comparison. No optimizer was run because total-cost and general accuracy
criteria do not justify it. Full evidence, hashes and assessment: FINDINGS.md
and analysis/. All earlier failed attempts remain preserved as separate files.
