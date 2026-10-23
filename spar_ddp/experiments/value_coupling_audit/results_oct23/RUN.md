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
