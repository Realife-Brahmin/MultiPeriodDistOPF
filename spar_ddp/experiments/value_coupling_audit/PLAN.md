# Coupling-retaining value representation: staged investigation

October23 compact component follow-on completed: COMPACT.md and
results_compact_oct23/FINDINGS.md. A recursive factor/operator merge passes
large-late value and direction checks with fast recurring component cost, but
fresh factor discovery and dense legacy export prevent a total-cost benefit.
No medium/optimizer extension or production/paper/benchmark changes.

October23 bounded follow-up also completed: see FOLLOWUP.md and
results_oct23/FINDINGS.md. Block-scaled/conditioned compression passes large-late
and medium-early frozen direction gates; medium late remains unresolved.
Three warmed alternating pairs at each passing point show higher prototype
total cost than dense. No approximate optimizer or larger campaign was run.

Final October 22 status: all 15 frozen points strictly validated using 13
preserved primary points plus 2 completed continuation points. See
results_oct22/FINDINGS.md and final_analysis/ for the assessment and provenance.
Retaining coupling improves accuracy, but late direction amplification and
unproven compact-kernel cost prevent advancement to a full-solve campaign.
No further solve is queued; historical submission notes below are retained.

User-authorized October22, continuing bb69011f. Production defaults, benchmark
results and paper stay unchanged. No candidate full-solve campaign before a
frozen replay passes. Do not confuse dense trajectory acquisition below with
a performance benchmark or a test of a new approximate optimizer.

## 1. Close actual coverage gap

Prior j==0 gate captured only early passes. Existing leaf matrices are reusable
for those same iterates. Saved solution.jls files contain x/u but NOT phi,zl,zu,
barrier mu or local factors; reconstructing the identical Newton systems from
them alone would be unjustified. A new immutable dense-only capture snapshot
therefore runs one IEEE123T3 smoke, then small/medium/original1020 T24 dense
trajectories. Midpoints use prior dense iteration counts:11/26/7, not j==0.

New opt-in FILTERDDP_FROZEN_CAPTURE_DIR does not activate approximation or old
value-group audits. Every pass caches original two-boundary leaves in memory;
selected workers copy local S derivative/KKT data BEFORE phase C overwrites
scratch. Only five labels are persisted: early, middle, first barrier-before/
after pair, converged late. Duplicated labels can share one snapshot. The
barrier selector mirrors the solver condition without changing it; before/
after must have equal k, consecutive j/pass and decreasing mu. Late requires
actual original stopping tolerance, not a guessed iteration count.

All24 leaves are saved at each selected point. Full local records and exact
state/control/equality/bound multipliers are saved for periods1/12/24 (or1/2/3
in smoke), plus dense control/equality-direction probes at zero and a fixed
small sine state perturbation. This supports replay of actual Newton systems,
not reinitializing missing multipliers. Preserve f=fu_active diagonals when
replaying battery directions. S contains pre-phase-C scratch, cl, dimensions
and current nominal point. Capture adds overhead; clocks are not benchmarks.

No warm-up solve is needed for data acquisition (WARMUP0,GCauto avoids disabling
collection during cold compilation). Same physical start, costs, algorithms,
exact curvature, screening and stopping tolerances otherwise. Whole Gold6438Y+
64CPU300GiB allocation,Julia1.12.4, private writable depot plus read-only prior
cache;2h batch cap,25min per capture process, no expensive overlap. Expected
new binary data~4--6GB; home had12GB free before submission. Monitor quota,
preserve old evidence and do not delete other work or silently omit captures.

verify_coverage.py supports cluster Python3.9 and checks actual logged labels,
late k versus final convergence, early<middle<late pass order, nonzero late j,
paired barrier transition, file counts and all snapshot SHA. COMPLETE only
after numerical/original physical/OpenDSS AND coverage gates. A smoke failure
must block larger captures. Repairs require new immutable snapshots.

## 2. Identical-iterate replay and topology partition

First reuse prior early matrices to regression-test the new dense replay.
For each new frozen point, compute ordinary dense suffix values and production
pairwise scan; require all-stage Vx/Vxx agreement before approximations.
Compare ALL Maa/Mab/Mbb/ma/mb at matching interval/stage/scan-depth, with
absolute and reference-scaled errors. Do not symmetrize Mab. Gradients remain
full vectors, not grouped scalars. Check conditioning/inertia of eliminated
K=left.Mbb+right.Maa and terminal closings at representative depths and where
errors grow. Record raw extrema and thresholds, not just Cholesky success.

Build an automatic graph-based hierarchy BELOW the root: recursively separate
connected tree components by battery-weighted vertex separators until small
leaf groups (initial cap8). Keep separator batteries as individual variables.
Keep graph separator identity/size and all battery-row permutations explicit.
Unlike root-only or arbitrary chunks, this preserves strongly interacting
local blocks while exploring general coupling boundaries. Never hardcode102
feeders or pretend a small per-period separator gives fixed temporal rank.

Exact assessment FIRST: local network KKT tree separator elimination provides
group-to-separator factors; saved local S/lay/stat/fac/b0 allow checking this
at selected stages. Temporal composition can concatenate separator interfaces
and grow rank with interval length. If an exact compact bound cannot be
established, retain exact dense/hierarchical fallback and label compression
approximate. Zero-truncation hierarchical storage must reconstruct the whole
quadratic to rounding and includes all vectors. Symmetry packing alone can
save storage but not cubic combination cost; report that separately.

## 3. Smallest compression candidate / cost controls

Initial approximate candidate: hierarchical off-diagonal compression of the
JOINT symmetric two-boundary Hessian [Maa Mab; Mab' Mbb], using the same
battery hierarchy on both boundaries. Preserve full within-leaf blocks and
both vectors. Compress cross-group blocks adaptively (e.g.SVD as diagnostic),
with a global Frobenius budget distributed across off-diagonal pairs. Retain
dense blocks whenever rank-factor storage/work is worse. No whole-matrix
low-rank truncation, no fixed rank, no dropping cross coupling by default.
Log actual reconstruction error and ranks/storage at every depth, including
construction/recompression/materialization costs. Dense expansion followed
by dense combination is only an accuracy prototype, NOT evidence of a faster
compact kernel. Rank growth may reject the candidate without a full solve.

Start on medium T24 complete suffix replay at representative points; it had
strong stress-group errors. Start original1020 on a four-period window with
the true dense suffix at the window boundary, then extend only if worthwhile.
This bounds large SVD costs while testing late/barrier behavior. Gate each
stage by source/input SHA, interval consistency and finite errors; preserve
negative results and inertia changes rather than shifting tolerances.

## 4. Recovered directions and candidate decision

Recover battery state/power directions from the ORIGINAL stage quadratic with
each candidate's future values, not from a different physical problem:
b=-(Mbb+Vnext)^-1(Mab'*a+mb+vnext), dp=(b-a)./f. Compare to dense at the same
frozen iterate and same initial perturbation. For selected periods use saved
S with _pp_solve and _tp_post in an isolated replay solver/update record to
recover full controls, equality and bound-multiplier directions. First require
dense replay to reproduce the saved zero/sine direction probes<1e-7. Do not
claim complete network directions from battery-only checks or infer original
dual residual from a dummy lambda used for a replay policy construction.

Report all errors over early/middle/late and barrier pairs; decide whether
retaining adaptive cross-group coupling improves directions over dropping it
without prohibitive construction/rank/storage cost. Only a candidate passing
those checks can proceed to a bounded dense-versus-candidate convergence test.
Then, and separately, clean repeated full and feasible near-optimal timing,
including extra iterations/setup/representation construction. No paper edits
or new benchmark claims from these diagnostics.

## Implementation status / continuation pointers

Dense capture job31357503 completed in
/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22, pinnede6b8b37b.
The small T3 coverage gate has already verified all5tags including late j9.
Do not alter that immutable snapshot. Cold capture timing is diagnostic only.

hierarchy.jl now provides a graph/battery-weighted centroid hierarchy below
the root, retaining separator batteries, and an INITIAL approximate cross-block
representation of the joint Hessian. hierarchy_unit.jl passed a synthetic
strongly coupled rank3 example: all states retained, local blocks full,
nonzero cross-group coupling preserved, global Frobenius budget satisfied.
This is NOT yet real-data validation, exact compactness evidence, or a compact
composition kernel. It returns BOTH compact payload and dense materialization;
numeric_payload_bytes counts coefficients of the compact payload only, not
resident memory of the prototype. Index/tree/object overhead is additional.

Capture completed and all four coverage/numerical/physical/OpenDSS gates passed;
see results_oct22/CAPTURE.md and the raw evidence archive (pushed96154b8a).
Replay implementation is now in replay.jl/replay_job.sh; read REPLAY.md for
scope, output schemas and exactness/performance limits. Original31357649,
source96154b8a in value_coupling_replay_oct22, FAILED at the medium late dense
policy gate; only9of15 points completed. Full failure evidence is preserved.
Independent diagnosis31357837 exactly reproduced all12 saved late probes using
the production IN-PLACE arithmetic; ordinary recurrence missed one by2.17e-6.
See results_oct22/REFERENCE_FAILURE.md; original1e-7 gate is unchanged.

New immutable replacement31357859 is submitted in
/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v2, numerical base96154b8a
plus committed replay overlay4f85b789 (45minute cap). The source.sha256 manifest
pins actual files; inherited numerical source is unchanged. Reference and
matched combinations now use production in-place arithmetic, with ordinary
recursion still checked independently. Original capture binaries are reused;
no new trajectory acquisition or optimizer solve. Home has5.4GiB free: no
more full source copies/captures without checking quota; delete no old work.
All-five-point T3 dense full-direction gates precede mediumT24 full suffixes
and original1020 first-four-period windows with the true dense future boundary.
No approximate optimizer/full-solve test is active. Inspect both stderr and
each per-point log; failure must preserve this snapshot, never weaken gates.
Synthetic hierarchy unit and parsing passed locally; the full replay module
requires cluster private DDP dependencies unavailable in the local environment.
V2 also FAILED after9complete points: medium-late dense probes were EXACT0,
but zero-truncation symmetric hierarchy exceeded its reconstruction bound on
raw slightly asymmetric input. See the appended REFERENCE_FAILURE.md and
medium_late_symmetry.csv. No approximate late mode or optimizer had failed.

Current replacement31358012 is in
/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3, overlay7b75cd80 with
read-only numerical root value_coupling_replay_oct22_v2/source. No full source
copy/capture added. If symmetric hierarchy cannot meet the SAME budget, preserve
the raw element bit-for-bit and charge/log FULL dense payload and the reason.
Reserve symmetry's share of the budget before truncation. Synthetic and real
failing-leaf tests passed; no tolerance increased. Raw fallback is not compression.

Same15points, one Julia process per system to amortize JIT; medium late first
after the full five-point small gate. Twelve minutes per system process,
45minute batch. Every point reloads its own state; no optimizer progression.
replay_series.jl/replay_overlay_job.sh and collector support case-series logs.
Keep v1/v2/v3 evidence separate; defaults/results/paper unchanged.
V3 reached its12minute LARGE-series process cap: job31358012 FAILED9:0 at
20m22,13of15points completed and independently posthoc-validated. All small/
medium and original1020early/barrier-before/barrier-after passed. The partial
large middle and untouched large late are NOT complete. Raw interruption,
including stderr, partial outputs and source overlay, is pushed atd1305334.

ONLY missing middle/late now run as31358119 in
/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_remaining, harnessd1305334
calling UNCHANGED v3 overlay7b75cd80/read-only numerical source. Eight minutes
per independent point,18minute batch. No repeats of completed13, no optimizer,
new captures or tolerance changes. Use collector --continuation REMAINING_ROOT
for final13+2 coverage; it must preserve the old failed status/timeout and
require new clean completion/source/input gates and unique point provenance.
Never create a fake replay.finished for the interrupted batch. See
results_oct22/REPLAY_INTERRUPTION.md. Final strict validation awaits both points.
Never infer performance from materialize+combine alone. Real-data replay,
exact-structure assessment and candidate decision remain to be completed.
