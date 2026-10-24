# Direct initial factors: partial generation verified, coverage gap remains

A direct generator works for the available period1 elimination record: it
does not form a full boundary quadratic then discover factors by dense SVD.
However the chain consumes periods2:4, whose elimination records were NOT
saved. Their1.258s discovery and the incoming true-future boundary's0.136s
discovery therefore remain charged. This experiment cannot establish the
cost of replacing those missing inputs by direct generation.

The honest partial-pipeline median is2.59433s versus1.33727s dense (1.94x
slower, explicit cleanup included). No fresh-input win, full direct pipeline,
PhaseB wall improvement or optimizer speedup is demonstrated. No trajectory
recapture, medium case, longer horizon or optimization test was launched.

## What is available and what was implemented

Verified capture inventory/metadata and fc_cache_stage source: full stage S
exists only for1,12,24. The other stages have dense quadratic leaves, not
their local primal/dual Newton states or network separator factors. The
packet's nominal field is LOCAL, not a full trajectory. The final solution
writer stores x/u, not phi/zl/zu needed for identical Newton reconstruction.
One cannot recover unique network KKT factors from the reduced leaves or
substitute period1 numerical factors into2:4.

For period1, the generator recreates local elimination blocks from sparse K,
existing network responses Bd and energy-elimination PN/NP/Ninv. It combines
these sparse pieces with the existing root separator factors using the raw
block-local-minus-U*Mroot^-1*U' identity. Dynamics scaling and sparse Rhat/D
terms give all three boundary blocks and full vectors directly (DIRECT_INITIAL.md).
Left/right factors remain distinct; neither raw skew nor Mab is projected to
symmetry. Upper completion occurs only where production requires it.

Initial full-boundary dense assembly/discovery SVD count is zero for this
generated stage. Known factors still undergo QR/small-core SVD recompression:
1,884 core calls, maximum core dimension8; separator width6, maximum sparse
cross interface5 and retained rank3 on this point. There are972 charged dense
cross-block fallbacks. Local elimination dimension is at most20 here; it
could be much larger on another topology. No fixed-rank, small-group or
horizon-independent guarantee follows from this large10k instance.

The true future at period5 represents the remaining20period suffix. Neither
the required complete suffix factor records nor a validated compact propagation
over that suffix are available. Its existing dense input is retained and its
numerical factor discovery is charged, not assumed free or replaced by zero.

## Accuracy and consumption

Final job31369093/41b98bcb COMPLETED0:0 in1m34; all source/input/capture/stage
hashes, hardware/occupancy and synthetic gates verified. Saved dense full
probes reproduce exactly. The independent dense-operator consumer still passes
before any generated input is used.

All three generated leaf1 blocks differ from the saved raw leaf by only
3.204e-12 absolute Frobenius /1.340e-20 unit-floor scaled; both vectors match
exactly. Repeating generation after poisoning the saved dense Sh0 with NaN
produces unchanged output, proving that generator does not read Sh0. Raw skew
differences1.73e-14 are recorded separately, not silently symmetrized.

Leaves2:4 retain the existing discovered factors and compact two-merge/closure
path. All three blocks/two vectors, closed value actions and vectors pass.
Generated leaf1's Maa/Mab/Mbb feed additional closed-period1 value actions,
gradient and battery-path recovery without dense boundary export. The largest
scaled value-action discrepancy is6.911e-9 (absolute1.298e-5 amid strong Schur
cancellation); this passes1e-7 but is not bit-exactness or a general late-
iterate sensitivity guarantee.

Battery, selected full control/equality/bound directions pass unchanged1e-7
scaled AND1e-5 relative/floor1e-10 gates at8.69e-13/5.33e-12. Stage CG needs
at most13iterations, prefix/battery solves12. Corrupting the composed gradient
fails both aggregate and full-control checks. All three repeated comparisons
pass; validation-only matrix expansion is excluded from execution and never
feeds it. No new feasibility/objective/OpenDSS claims arise without optimization.

## Warmed accounting and limitations

Strictly validated THREE alternating same-input pairs in one process. Every
measured phase has zero Julia compile_time AND recompile_time. Per-phase
medians: direct period1 generation0.04078s; missing2:4 discovery1.25772s;
incoming future discovery0.13580s; chained updates0.17493s; condition0.00631s;
closed-period1 actions0.11022s; directions0.46448s; cleanup0.40632s.

Dense period1 leaf generation0.01623s uses its EXISTING dense Sh0. It is not
the same output representation as direct compact generation and is not a
measurement of dense-to-compact discovery for period1. Dense merges0.38272s,
closed-period1 actions0.02016s, directions0.51904s and cleanup0.39964s are
charged. The extra prefix workload means these totals should not be directly
treated as a speedup/regression against the earlier narrower chaining study.

Median total2.59433s partial-direct vs1.33727s dense; without explicit cleanup,
2.19272s vs0.93762s. Setup generation/discovery1.43635s vs0.01623s is separated
from recurring1.16194s vs1.32155s but NOT omitted from fresh totals. Fixed
structural topology and existing phaseA network factors are given inputs;
their cold construction is not a new optimizer benchmark. Numerical factors
are freshly generated/discovered each repetition, not reused across iterates.
Identical untimed premeasurement GC resets the benchmark; in-phase and explicit
postmeasurement cleanup GC are charged. Phase medians need not sum to total median.

Cumulative allocation is about1.526GB partial-direct vs539MB dense; direct
leaf1 generation61.0MB, missing2:4 discovery623MB, future70.6MB. This prototype
has not eliminated allocation overhead. Generated leaf1 stores580,704 coefficient
bytes (raw24,985,920), object graph1,131,264bytes; neither is resident memory.
Chain interior dimensions still1020/2040/3060, coefficient storage0.285/1.359/
6.072MB and maximum recompressed rank18. Network-factor and full-stage consumer
objects still contain dense Sh0; generator independence is NOT proof the whole
pipeline can omit that matrix. Sampled batch MaxRSS2,792,312KiB includes both
paths/references/validation, not isolated candidate memory savings.

## Decision

Direct factor generation is structurally and numerically feasible for the
available stage, not disproved. Its application to the ACTUAL expensive chain
inputs is untested because the matching2:4 elimination records are missing.
The next prerequisite is an explicitly scoped way to obtain those records at
the same Newton state, with leaf/hash equivalence checks; the incoming suffix
must remain charged until separately supported. That would require new capture
or sufficiently complete saved primal/dual states, neither supplied here.
No new trajectory/campaign was silently substituted for immutable data.

Production defaults, benchmark rows and paper remain unchanged. Failed logging
attempt and completed mixed-compilation pilot are preserved separately and
excluded from warmed results (HARNESS_FAILURE.md). Reproduce final strict
analysis with analyze_direct_initial.py EXTRACTED_FINAL_ROOT NEW_OUTPUT_DIR.
