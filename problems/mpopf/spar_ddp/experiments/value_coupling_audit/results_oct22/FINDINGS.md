# Coupling-retaining frozen replay: final assessment (October 22)

Retaining cross-group interactions is substantially more accurate than dropping
them, including on the medium system where natural-root grouping is trivial.
However, the present global-Frobenius compression budget does not control late
Newton-direction errors. Do not advance this prototype to a full-solve or timing
campaign yet. Production defaults, benchmark results and paper are unchanged.

## Evidence and scope

Strict collection validates 15 frozen points / 75 mode-point observations:
five small T3 gates, five medium T24 full suffixes, and five original
1,020-battery four-period windows with the true production-dense future boundary.
Each includes early, middle, converged late, and the first barrier-before/after
pair. All physical batteries, locations, dynamics and limits remain present.
All three quadratic blocks, both vectors, scan depths, and selected recovered
control/equality/bound directions are recorded. Dense reproduction of saved
zero/sine probes passes the unchanged 1e-7 gate before approximation.

This is not a compressed full-T24 large solve. The capture trajectories passed
original physical/objective/bound and independent OpenDSS checks; replay does
not establish feasibility or convergence of a new approximate optimizer.

Job 31358012 remains FAILED after its large-series 12-minute process cap,
with 13 completed points and partial middle evidence preserved. Continuation
31358119 completed 0:0 in 10m36 and computed only missing large middle/late
points using unchanged numerical code. Final validation explicitly combines
13+2 with point-level attempt provenance; it does not manufacture completion
of the interrupted attempt. Earlier reference/symmetry failures are preserved.

## What the implementation does today

A battery-weighted topology hierarchy below the root (maximum leaf size 8)
preserves within-group quadratic blocks and both full vectors. Adaptive SVD
compresses only off-diagonal blocks of the joint two-boundary Hessian, not the
whole matrix. Ranks are selected from the error budget, never fixed in advance.
Unprofitable blocks retain dense storage. Slight raw asymmetry is not silently
discarded: if symmetric projection/reconstruction cannot meet the same budget,
the full raw element is retained bit-identically and charged as dense payload.
Mab is not assumed symmetric.

The accuracy prototype stores BOTH compact payload and dense expansion, and
uses dense production combinations. Closed suffix PintVal objects remain
dense. Thus it is not a compact combination kernel and saves no demonstrated
resident RAM. Payload ratios below count coefficients only, excluding indices,
objects and dense expansion, and are medians over represented elements rather
than total process memory.

## Accuracy: useful coupling, consequential late amplification

Direction errors below use ||candidate-reference||F / max(1,||reference||F),
maximized over the selected full direction fields/probes. These are not always
pure relative percentages; strict relative and absolute errors remain in CSVs.

At the medium late point, dropping coupling gives 0.3284 direction error.
Retaining coupling with budget 1e-8 reduces it to 6.923e-5, whereas budget 1e-5
still gives 0.1059. Even the zero-truncation representation gives 5.685e-6
direction error despite quadratic-block scaled error 7.747e-13 and one raw
fallback. Zero truncation is therefore not a bit-exact policy guarantee.

On original 1,020-battery windows, budget 1e-8 gives direction errors
1.288e-10, 1.621e-9, 1.015e-10 and 7.912e-9 at early, barrier-before,
barrier-after and middle respectively. At late it deteriorates to 0.001716,
versus 0.008728 for dropping. This is only about a fivefold improvement there.
The loose 1e-5 representation selects rank zero for all compressed cross blocks
at late and matches dropping's error. Late degradation affects controls and
battery directions too, not just equality multipliers: retained 1e-8 has
control scaled error 0.001100 and battery strict relative error about 0.001157.

A local joint-Hessian global error budget does not bound every small subblock,
composed suffix or recovered direction. For example, 1e-8 late maximum scaled
quadratic-block errors are 0.004009 (medium) and 0.006486 (large window).
No inertia changes were observed, and all 2,000 recorded systems have the
conservative actual-upper-triangle sign certificate described in INERTIA_NOTE.
That diagnostic is not a convergence or feasibility guarantee.

## Exactness, rank and storage

Exact compactness is conditional on the full time-expanded separator, not a
small single-period interface. Temporal composition may concatenate interfaces
and grow rank. See SEPARATOR_ASSESSMENT.md: a closed exact compact algebra is
not established. Zero-truncation/dense fallback is an exactness control, not
proof of such an algebra, and floating-point policy amplification remains.

At budget 1e-8, median coefficient-payload ratios across the five points are
0.07752--0.13109 on medium (7.63--12.90x smaller) and 0.006751--0.008201 on
large windows (121.93--148.13x smaller). These are promising representation
sizes, NOT measured RAM savings. Large-window maximum retained ranks vary
5/6/6/5/2 across early/barrier-before/barrier-after/middle/late. Lower late rank
does not mean better direction accuracy. Zero-truncation median payload ratios
are about 0.6803 medium and 0.6694 large, largely symmetry packing rather than
evidence of a new fast exact separator method. Rank/spectrum growth by scan
depth, including adverse high exact ranks, remains in the archived CSVs.

## Cost and smallest next experiment

There is no demonstrated runtime win. Construction/materialization costs are
additional to dense combination. For example medium-late 1e-8 construction is
16.19s and materialization 20.84s, versus dense serial combination 0.237s;
large-window construction ranges 4.91--53.42s versus dense combination around
0.915--0.919s. These single observations mix JIT and warmer execution and are
NOT benchmarks or estimates of an optimized compact kernel. Phase B wall time
also includes shared-memory access, scheduling, synchronization and waiting;
none of these measurements establish a communication-only bottleneck.

Reject the fixed 1e-5 candidate for advancement. Keep 1e-8 as a research lead,
not a validated optimizer configuration. The smallest next controlled test is
medium late plus large late, with one early control, reusing these frozen
states. Compare global-budget compression against block/separator-scaled or
Schur/energy/direction-aware error control, with conditioning-aware tightening
and dense fallback. Preserve local blocks and every physical variable. Check
all recovered direction fields, not Hessian error alone, and keep the existing
reference gates. No assumption of fixed rank or whole-matrix truncation.

For a surviving candidate, warm construction, materialization and combination
separately, then alternate three identical-input repetitions in one process;
record ranks, storage, allocation/GC and fallback costs. A factor-aware compact
combination kernel is needed before claiming computational savings. Only after
accuracy/inertia/storage/total-cost checks justify it should a bounded dense
versus candidate convergence test run. Later clean repeated comparisons must
include extra iterations and representation setup, with feasible near-optimality
and full convergence reported separately. No new full solve is queued here.

## Reproduction and provenance

Use analyze_replay.py PRIMARY_ROOT OUTPUT --continuation REMAINING_ROOT.
final_analysis/validation.json records strict 15-point completion and preserves
the old timeout; point_sources.csv and replay_files.sha256 pin attempts and
inputs. Full raw primary evidence is replay_v3_interrupted_evidence.tgz;
remaining evidence is replay_remaining_evidence.tgz, SHA256
2d0d56acca579c2002088953441acef2743df6c49cc14f2d0481a2b0603984f0.
Original capture binaries remain immutable on Kamiak, as documented in CAPTURE.
The collector additionally accepts Julia's actual `signal 15: Terminated`
spelling alongside `signal (15): Terminated`; scientific gates are unchanged.
