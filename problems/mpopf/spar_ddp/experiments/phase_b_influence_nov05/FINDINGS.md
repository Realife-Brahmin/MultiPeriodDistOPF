# Phase B influence audit: barrier curvature is promising; stale gradients are sensitive

The bounded audit supports investigating exact curvature reuse across a
barrier-only transition. It does not support general previous-Newton-iterate
reuse or a new low-rank implementation. At the sampled barrier transition,
H/B_L/B_R are unchanged and the five matrix-retention interventions reproduce
all tested outputs exactly. The gradient input h must remain current.
Early/middle ordinary-iterate retention is sensitive; late matrix retention
looks substantially safer for recovered directions, but value/action checks
still expose discrepancies. These are mixed-input sensitivity diagnostics,
not internally consistent or deployable Newton updates.

## Dataset, reference and coverage

Original1020-battery large10k, fullT24, original balanced/current in-place
kernel; every variable, constraint and raw coupling retained. Eight captured
states form consecutive accepted early(k0->1,j0,mu1), middle(k7->8,j3,mu.008),
late(k12->13,j7,mu9.27718e-6) pairs and barrier(k4,j0->1,mu1->.2).
The late pair is late optimization, not two converged iterates. Actual pass
numbers/accepted-step links, regularization-before/reported-after fields,
source semantics and every leaf/stage SHA are in the capture metadata.
Reported reg_last after successful procs is reset by the original source;
these are recorded fields, not a guessed applied-regularization claim.

Capture31580863 completed0:0 after the worker experiment,9m49,cn156.
Original numerical/bounds1e-5 pass; independent original physical residual
6.82e-13, OpenDSS discrepancy8.18e-9pu/no voltage violations. Replay31596954
completed0:0 in20m50 on Gold6438Y+,64CPU300GiB,Julia1.12.4/OpenBLAS ILP64/BLAS1.
No expensive overlap. Source/runtime/input checks pass before/after.

All2536 production-reference rows have zero numerical error: all24 stages,
zero/sine/actualSOC saved probes, full controls/equalities/lower+upper bounds,
actual SOC propagation and the unchanged-input kernel control. All16 deliberate
composed-gradient negative controls reject (minimum scaled error90.7).
The independent battery closure agrees with reference to max1.07e-14 scaled;
its errors are separately retained. Dummy reconstruction lambda is never
reported as an original dual residual. The observed intervention errors below
are therefore not unexplained dense-reference reconstruction error.

Four pairs x six interventions at INTERIOR merge(depth2,left5,right7,
interval5:8): previous wholeH/B_L/B_R/h, Hdiagonal-only or B_Loffdiagonal-only.
Every remaining global round is replayed; the actual dependency DAG changes
closed values1:5, consumed as futures by stages1:4. Full forward recovery
continues through all24 stages with CURRENT network/barrier/reconstruction
records and the true full-horizon terminal, not a truncated-window boundary.
Coverage:1632 merge dependencies,16128 full-direction probe rows,3456 full
trajectory rows,8064 value-action rows,1200 complete suffix block/vector rows.
No factorization failure or production dense fallback occurred.

## What changes, and where

At the intervention merge, relative H/B_L/B_R/h changes are respectively:
* early:6.04%,3.46%,3.57%,19.64%;
* middle:2.33%,0.410%,2.20%,94.05%;
* late:5.46e-8,3.60e-8,4.99e-8,99.91%;
* barrier:0,0,0,1889%.
Late h is small absolutely (change1.26e-5), but its previous norm is also
1.27e-5: it is not relatively stable. At deeper late merge9/13, B_L/B_R norms
are below1e-10; unit scaling, denominator floors and explicit near-zero flags
prevent interpreting relative values there as substantial coefficient changes.

Distribution matters. At merge5/7, H delta diagonal energy is99.10% early,
99.9926% middle and64.42% late. Middle B_L delta is99.9242% diagonal, yet
omitting ONLY that offdiagonal change (keeping prior entries) produces
6.94e-5 maximum scaled
recovered-direction error. A small energy fraction is not a safe omission.
The largest11battery rows (~1%) hold only2.74-3.65% of early/middle H delta
energy, rising to11.13% late: changes become more concentrated, not confined
to a handful of coordinates. Largest H row buses are1057/4001/19100 across
those phases. All row/column/topology-block quantiles, outliers and bus IDs
remain raw; topology groups are diagnostic partitions, not eliminated variables.

The DELTAS are not generally much lower rank than their full matrices.
At merge5/7 and relative singular-value cutoff1e-3, delta H/B_L/B_R ranks
are1020/1020/1020 early,1020/1017/1020 middle,1020/999/1020 late;
all are1020 at1e-6. Current matrices are1020 at1e-3 there. The deeper late
B_L delta has rank159 at1e-3 but absolute norm1.02e-11, already below the
near-zero scale. This is no fixed-rank or horizon-independent result.
All24 actual-original-UPPER inertia checks are positive definite; selected
merge5/7 condition numbers are1.08-1.58, while all sampled locations span
1.06-10.61. Raw asymmetry is recorded and never symmetrically averaged away.

## Influence and first amplification

All six early and all six middle retentions exceed recovered-direction gates.
Middle stale h has max0.6495 scaled upper-bound direction error despite
max closed-gradient error5.17e-4. At the barrier, all five matrix tests have
ZERO errors at every measured merge/value/action/direction/trajectory output;
stale h instead gives0.0188 max fixed-probe direction and0.0796 forward-
trajectory error. Stable matrix curvature does not imply stable right-hand sides.

Late matrix tests have max recovered fixed-probe errors2.16e-9 to7.56e-9,
within unchanged1e-7 scaled AND1e-5 relative/floor1e-10 direction gates.
However every one has an intermediate/value/action threshold exceedance;
for example stale B_L gives suffix1.08e-7 and action2.62e-7. These are promising
local direction observations, not full-output reuse certificates.

Late stale h is the clearest small-value-error counterexample. Every merge,
complete suffix and action passes, with maximum normalized gradient error
1.264e-9 (absolute7.589e-6 against norm6003.8). First relative threshold crossing
is stage2 upper-bound recovery: absolute8.66e-12 against norm1.99e-7,
relative4.34e-5; this tiny absolute discrepancy is explicitly distinguished.
The first SCALED crossing is stage3 equality recovery,3.046e-6. Stage4 equality
error reaches6.204e-6 (reference norm.08876, relative6.99e-5), with maximum
coordinate2.053e-6 at equality coordinate41283, mapped bus1. Full controls
reach3.539e-6 and upper-bound directions5.802e-6; lower bounds stay8.25e-11.
The recursively recovered trajectory confirms6.204e-6 equality error.
Thus the important amplification first appears in current-stage consumption,
not a failed merge factorization. Ratios of differently scaled output metrics
are not asserted to be operator condition numbers.

## Cost implications and supported next experiment

Ordinary h needs refreshing even late; updating a1020-vector is cheap relative
to factoring/combining dense1020x1020 matrices. Exact barrier matrix stability
is the useful expensive-operation candidate: cached curvature factorization
and matrix solves/products could plausibly survive a pure barrier transition,
while vectors, right-hand sides and policy/bound recovery still update.
Only three locations and one barrier transition were established here; this
is not proof that every global curvature factor is reusable.

Diagnostic costs:24 injected dense combinations total3.737s (~.156s each);
36delta SVDs total10.884s (~.302s each),36current SVDs11.128s. SVD-based
change discovery would itself be costly. Matrix subtract/norm detection is
~.012s per sampled matrix, versus ~20microseconds for a vector; complete
allocation/retention-copy, outlier and IO costs are not individually timed.
All clocks include diagnostic context/possible first-use compilation and
repeated reference/consumer work. They are not clean PhaseB benchmarks or an
optimizer speedup. Batch MaxRSS15,897,568KiB (~15.2GiB) includes captures,
references/JIT/validation, not a private-memory saving.

The supported NEXT experiment is a bounded exact barrier-transition cache
prototype/test: first verify unchanged curvature at ALL24 leaves and EVERY
global merge, define conservative invalidation on primal/dual/factor/
regularization changes, then charge storage, detection, RHS refresh and
complete current-stage consumption against a fresh baseline. Preserve the
same direction gates. No ordinary-iteration reuse, broad sweep, compression
optimizer, default, paper or main-benchmark change follows from this audit.

## Evidence and reproduction

Run collect_replay.py EXTRACTED_REPLAY_ROOT NEW_OUTPUT_DIR. Strict collection
passes with no missing coverage/failures; threshold-exceeding interventions
remain in raw CSVs and are not relabeled successful solver updates.
results_replay contains the complete code/log/hardware/accounting/hash/raw
CSV archive and analysis (including separate first scaled/relative crossings,
row concentration and independent reference errors). results_capture contains
source/runtime/textual capture provenance and192 packet hashes. Immutable
binary captures remain at /home/aryan.jha/kamiak_jobs/phase_b_influence_capture_nov05;
replay snapshot /home/aryan.jha/kamiak_jobs/phase_b_influence_replay_nov05.
Archive SHA25665277407dae8044f87f234825fe6769756ecb3a1c489e63c228d717b95e21247.
