# Coupling-retaining value representation: staged investigation

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

Dense capture job31357503 is running in
/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22, pinnede6b8b37b.
The small T3 coverage gate has already verified all5tags including late j9.
Do not alter that active snapshot. Cold capture timing is diagnostic only.

hierarchy.jl now provides a graph/battery-weighted centroid hierarchy below
the root, retaining separator batteries, and an INITIAL approximate cross-block
representation of the joint Hessian. hierarchy_unit.jl passed a synthetic
strongly coupled rank3 example: all states retained, local blocks full,
nonzero cross-group coupling preserved, global Frobenius budget satisfied.
This is NOT yet real-data validation, exact compactness evidence, or a compact
composition kernel. It returns BOTH compact payload and dense materialization;
numeric_payload_bytes counts coefficients of the compact payload only, not
resident memory of the prototype. Index/tree/object overhead is additional.

Replay runner still needs implementation after capture gates pass: ordinary
reference, matching full scan-depth records, inertia/conditioning, frozen
direction recovery from saved S, rank/cost/spectrum CSVs and candidate decision.
Use this committed library, regression-test on smoke, then bounded medium and
large-window replays. Never infer performance from materialize+combine alone.
