# Compact-to-compact composition: bounded assessment

The compact output CAN feed another merge and a direction consumer without
intermediate dense boundary export on this immutable original1020 large10k
late four-period window. Accuracy passes. Initial factor discovery remains
the obstacle to a fresh-input pipeline benefit: median2.46198s compact versus
1.29874s dense (1.90x slower), including explicit cleanup GC. Conditional on
initial factors already being available for this SAME unchanged iterate,
the measured recurring workload is1.02396s versus1.29874s (about1.27x faster).
This is NOT a fully compact arbitrary scan, PhaseB wall or optimizer speedup.
No medium/full-optimization extension was launched.

## What was composed

The existing period3/4 Schur object feeds period2's merge directly. Rather
than export its three dense matrices or assume fixed-rank boundary closure,
the new representation retains its compact interior factors/history. It forms
and factors a hierarchical augmented two-interior system, then closes against
the true period5 future by retaining a third interior. The original capture
has24periods; this is ONLY its first-four-period window, not a full compressed
T24 trajectory. The true future remains an existing dense input whose initial
hierarchical discovery cost is charged.

Original-coordinate raw diagonal upper completion is explicit, before spatial
permutation. Cross matrices remain raw and unsymmetrized; full vectors and
strong local blocks are retained. No physical battery/location/dynamics/limit
or original constraint was modified. Algebraic Schur actions can differ from
production's Cholesky/copytri! rounding; that difference is measured against
both the matching dense chain and unchanged production values/directions.

Interior dimension grows1020 ->2040 ->3060. Their hierarchical coefficient
storage grows284,944 ->1,358,720 ->6,072,000bytes. Maximum stored low-rank
interface widths at these representations are2,5,18; the recompression routine
records a maximum retained rank18, Woodbury core36 and QR/SVD core44. These
are observed values, not rank bounds or horizon independence. Initial true
future factors already have maximum low-rank width11.

The prior Schur object's input/interface factors are reused, but its LU is
NOT reused to factor the larger retained-interior system. Both refactorizations
are charged. Cross-factor QR/core-SVD updates remain necessary; this is not
an SVD-free or exact closed separator algebra. It makes one bounded composition
economical with factors available, not unlimited temporal composition.

## Accuracy and downstream consumption

Job31368962/overlay eb4f0439 COMPLETED0:0 in1m28. Gold6438Y+,64CPU300GiB,
Julia1.12.4/BLAS1/OMP1, occupancy guard, no overlap. All synthetic and actual
source/input/leaf/stage hash gates passed. Dense saved control/equality probes
reproduce exactly; independent ordinary/production value discrepancy5.35e-11
passes its unchanged1e-7 gate.

All three quadratic blocks and both vectors pass at both merges, followed
by final closed value/vector checks. Maximum unit-floor scaled discrepancy
is1.551e-12. Validation-only full matrix exports do NOT feed execution and are
outside execution timers. Final operator symmetry defect is2.18e-41 here;
this is measured, not presumed. Actual original-upper-triangle first-stage
elimination inertia matches dense at the existing machine-scaled threshold.

The downstream consumer reuses the saved stage network elimination factors
and raw dense stage Sh0 input. It solves its halved battery Schur system with
bounded matrix-free preconditioned CG using composed future actions, then
recovers ALL selected full control, equality and bound-multiplier directions.
This is not a claim based only on battery slices. First, with the unchanged
dense future, this independent consumer passes at2.75e-14 scaled/6.25e-13
relative. Then the compact chain, selected zero/sine stage1 directions and
four-period battery trajectory pass at8.69e-13 scaled/5.33e-12 relative.
Unchanged gates are1e-7 scaled AND1e-5 relative with1e-10 floor.

Full-stage CG used at most13iterations; the first battery-state solve used12.
Its residual tolerance2e-13, true-residual check,300iteration cap and non-SPD
rejection are separate from, and do not relax, the direction gates. Explicit
corruption of the composed gradient produces scaled direction/path discrepancy
1.76e7 and fails: the direction probes depend on the composed computation.
Every warmed repeat passes the unchanged gates.

The compact update path records2,942 dense CROSS-block fallbacks; their
construction/factor/action costs are included. Actual local dense factor size
never exceeds60. A bounded local rank-growth fallback up to192 is supported
but was NOT needed on this real point. No full boundary conversion was needed
by the algorithm (algorithm_dense_exports=0). A larger representation breakdown
would be rejected explicitly, not hidden as successful compression.
Maximum local relative recompression tail is approximately9.89e-13 against
the unchanged1e-12 update budget; this is not a global direction-error bound.

## Warmed cost accounting

Three identical-input alternating pairs in ONE process, after warming the
exact discovery/update/condition/consumer paths. Medians of phases:
dense three merges0.38514s, downstream directions0.51544s, cleanup0.39803s.
Compact fresh input factor discovery1.41532s, chained updates and terminal
closure0.18281s, condition probes0.00642s, downstream directions0.43451s,
cleanup0.40300s. Phase medians need not sum to median totals.

Totals including explicit cleanup: dense1.28893--1.32504s, median1.29874s;
compact fresh2.42120--2.52580s, median2.46198s. Subtracting ONLY fresh numerical
discovery for valid within-iterate reuse gives compact median1.02396s.
Without explicit cleanup the corresponding medians are dense0.90058s,
compact fresh2.05897s, compact with factors available0.62831s. Neither inclusion
nor exclusion of common explicit cleanup changes the fresh-input conclusion.

All per-phase allocation and GC clocks are preserved. Cumulative allocation
per measured pair is about472MB dense versus1.301GB compact fresh; discovery
alone allocates694MB. With input factors available, compact's other phases
allocate about607MB, still above dense. Actual within-phase GC plus explicit
cleanup is recorded, not inferred from RSS. Identical untimed premeasurement
GC establishes each benchmark boundary; the separately recorded postmeasurement
cleanup is charged. These phase sums are NOT the entire diagnostic job wall.

Fresh totals assume an existing structural topology, as appropriate for a fixed
network/battery set. Initial topology setup0.88257s was measured cold separately
(includes compilation, NOT a warmed benchmark or free first-use setup).
Initial cold discovery2.12310s is also separately recorded, not mixed into warm
clocks. No numerical-factor reuse across changed Newton iterates is assumed.
Condition probes are charged uncertified estimates, not condition-number bounds.

The dense consumer builds its existing factor-backed policy and necessary
two remaining battery suffixes. The compact consumer computes only the selected
zero/sine full-stage probes plus the battery trajectory through actions. Thus
the measured downstream saving is specific to this consumption workload;
it does not establish the cost of every policy/scan operation in an optimizer.
The independent dense-operator consumer gate prevents attributing a consumer
implementation error to compression.

## Storage is not memory or solve performance

Initial factors/vectors for the three leaves plus true future occupy3,041,440
coefficient bytes versus83,289,120 raw coefficient bytes (about27.38x fewer).
The raw metadata's3,033,280 count omitted the terminal vector's8,160bytes;
the posthoc collector explicitly corrects this without changing raw evidence.
Those inputs do not include all subsequently constructed factor/work arrays.
The composed object graph occupies15,920,263bytes by summarysize, including
indices/tree/object overhead, not allocator/process RSS. Existing dense stage
Sh0 and true-future inputs each occupy8,323,200 coefficient bytes and remain
present. We have NOT implemented direct stage/network generation of those
initial factors or eliminated every dense input.

Slurm sampled batch MaxRSS2,229,660KiB covers loading, BOTH implementations,
references and validation-only dense exports. Sys.maxrss is separately saved;
neither is an isolated candidate/dense resident-memory comparison. No actual
optimizer RAM saving, communication-only improvement, parallel PhaseB wall
gain, original-model feasibility/quality, iteration counts, near-optimal time
or full-convergence time is demonstrated. No compressed optimizer has run.

## Decision and next bounded question

Chaining and compatible consumers pass on this point and remove the old
algorithm-required dense export. Growing retained interfaces are manageable
over this short window, but fresh discovery still exceeds the recurring
compute saving. The correct conclusion is conditional compact-computation
progress, NOT a fresh-input pipeline win and NOT failure of compression.

The next proposed extension is to assess direct initial factors from the
existing stage/network elimination representation on this SAME frozen point,
including the incoming true-future boundary. It must preserve raw triangle
semantics, full-stage directions and charge construction; no assumed free
terminal factor or reuse across Newton iterates. This is separate unresolved
work, not automatically launched. Longer chains would additionally need
rank/storage/factor-growth checks; this experiment does not justify general
compact closure or a full optimization campaign.

Production defaults, published benchmark rows and paper are unchanged. Raw
evidence and strict collector outputs are committed with this assessment.
Reproduce: analyze_compact_chain.py EXTRACTED_EVIDENCE_ROOT NEW_ANALYSIS_DIR.
