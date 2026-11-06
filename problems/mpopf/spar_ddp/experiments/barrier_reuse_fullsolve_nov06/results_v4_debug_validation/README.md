# Complete integrated T24 diagnostic validation

The repaired process-parallel optimizer completed a full original1020-battery
T24 solve on cn179 with grouped b2/right kernel,24 workers x2 BLAS. Strict
posthoc validate_debug.py passes against the retained actual complete output.

23 passes,14 accepted steps,15 fresh populations and8 eligible hits. All1196
fresh shadow merge nodes and every saved/current recovery field, full-stage
policy and multiple-scale control/equality/lower/upper actions are covered.
All34878 recorded comparisons are bit exact, with zero absolute/scaled/relative
error. All168 actual stage negative controls and416 deliberately stale-gradient
node controls reject. There are no observed factorization fallbacks. Measured
first pass populates rather than reusing warmup payloads; retirement is warmed,
then charged within measured exit. Minimal live cache holds1,298,419,200 numeric
bytes (52 merge entries x3 dense1020-square matrices).

Original final numerical residuals pass1e-5 (primal5.085987878709e-7,
dual3.635565146322e-10, complementarity7.829303945997e-6); original bounds have
zero violation. Independent original physical primal6.821210263297e-13;
OpenDSS maximum voltage difference8.18e-9pu and zero voltage-limit violations.
All actual worker OS-task masks remain disjoint two-core teams on one socket
per team, with49 distinct physical cores including coordinator; actual BLAS2,
OpenBLAS backend/default and coordinatorBLAS1 audited before/after. Source,
harness, input and runtime binary/library checks pass before/after this stage.

This diagnostic retains fingerprints, fresh all-merge shadows, full policies,
negatives, extra allocations and probes. Its339.584s recorded optimizer wall is
NOT a clean performance observation. Clean identical-source baseline/cache
ABBA timing runs remain in progress, followed by T96 only after T24 validation.
There is no speedup conclusion here. Raw warmup/solve/physical/OpenDSS/audits,
CSV checks and provenance are retained. Diagnostic archive SHA256:
a32e001084494618ec6f2c1dd9a11b5a1b567c744c20355725969c8df640535c.
