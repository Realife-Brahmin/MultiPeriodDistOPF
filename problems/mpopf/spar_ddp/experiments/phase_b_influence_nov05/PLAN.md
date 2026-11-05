# Phase B numerical change and influence audit

User request: attachment34e66b83-bad0-4dcc-8c20-131d6f769eac/Pasted text.txt.
Do not interrupt/change the ongoing24x4 job31580246. All expensive capture/replay
work must follow it sequentially. No approximate optimizer, compression
implementation, worker sweep, production-default or paper change.

## Existing evidence and dataset gap

Read value_coupling_audit/README.md and results_campaign_oct23/FINDINGS.md,
results_oct22/FINDINGS.md and results_oct23/FINDINGS.md. Compact coefficients
alone failed to control directions; late policy/equality/bound reconstruction
amplified very small errors. Accurate compact composition did not establish
a fresh total-cost win. This audit is not a repetition of that campaign.

Inspected immutable value_coupling_audit_oct22 large10kT24 coverage:
early(k0,j0,pass1); barrier(k4,j0/1,passes5/6); middle(k7,j3,pass11);
late(k14,j8,pass23). Only stages1/12/24 have full S. The subsequent
compact_campaign_capture_oct23_v3 has ONLY late(k14,j8), fullS1:4.
Thus the barrier leaves are consecutive, but ordinary consecutive accepted
pairs and all-stage consumers are unavailable. solution.jls lacks Newton
duals/factors and cannot manufacture the missing records.

One targeted ORIGINAL1020 T24 BASELINE capture is justified. Fixed ORIGINAL
balanced scan/current dense in-place kernel/BLAS1 to match the reference path,
not the new scheduling candidate. All source/input manifests are pinned.
Cold no-warmup capture is diagnostic, not a timing benchmark.
No new horizon/network sweep or extra optimizer candidate.

Chronology from the validated baseline gives these selected states:
* early:(k0,j0)->(k1,j0), same mu1;
* middle:(k7,j3)->(k8,j3), same mu0.008;
* late:(k12,j7)->(k13,j7), same mu~9.27718e-6;
* first barrier:(k4,j0)->(k4,j1), mu1->0.2.
Late pair is late optimization, not a falsely labeled converged pair.
Validate ACTUAL consecutive pass numbers, accepted-step linkage for ordinary
pairs, identical barrier index/mu there; identical k and decreasing mu at
barrier. Never infer stability from scattered states. If chronology differs,
preserve/diagnose; do not relabel missing pairs or loosen coverage gates.

Save all24 exact leaves plus pre-phase-C deep-copied S/nominal primal/dual/
bound/barrier data for ALL24 stages at eight selected states. Existing S
packets~85.9MB/stage and leaves~0.6GB/state imply~22GB total; require>=30GiB
free before capture, preserve all prior work. Retain selected states only,
not every merge intermediate. Store raw regularization-before-pass,
reported-after-pass and documented source semantics, actual step outcomes,
schedule/kernel/backend metadata and hashes. No guessed regularization0.
Replay reconstructs intervals/intermediates under the same fixed association.

Save live zero/sine and actual linearized SOC-trajectory full control/equality/
lower+upper-bound direction probes. The actual SOC trajectory is propagated
through CURRENT original stage feedback and fx/fu, starting at zero fixed
initial-state perturbation. Additional scaled cosine and topology block-sign
probes are constructed deterministically offline. Do not substitute next
accepted-state differences for a current Newton direction.

## Baseline reference gate first

At each selected state verify all leaves/packets/input SHA and production
in-place scan association. Replay CURRENT S with _pp_solve/_tp_post; all24
live saved battery/control/equality/bound direction probes must meet existing
scaled1e-7 AND relative/floor1e-10 threshold1e-5 gates. Record dense errors,
including coordinate/block outliers, before substitutions. Preserve and
diagnose reference-only errors; never blame substitutions automatically.
An independent consumer and a large deliberate composed-gradient negative
control must demonstrate that downstream probes exercise changed operators.
Dummy reconstruction lambda is NOT an original dual-residual diagnostic.

## Matched merge inputs and distribution

For each pair reconstruct matching intervals under balanced production scan.
Initial bounded change sampling: merges(depth1,left13,right14),
(depth2,left5,right7),(depth3,left9,right13). One-based record indices,
terminal25. Validate interval identities and mutation/read dependencies.
Use H=left.Mbb+right.Maa, B_L=left.Mab,B_R=right.Mab,
h=left.mb+right.ma (PintElem/PintElem locations). Keep raw asymmetry;
production uses ORIGINAL-coordinate UPPER completion of H.
Report ALL3 matrices and vector, plus complete value blocks/vectors where
practical. Absolute/Frobenius/max-coordinate error, strict-relative with
explicit zero flags and unit/floor scaling; row/column and battery/feeder
block quantiles, top outliers, diagonal/offdiagonal energy/concentration.
Record singular spectra/numerical ranks of DELTAS separately from originals,
actual-upper conditioning/inertia and applied regularization metadata.
Small magnitude, temporal stability and downstream importance stay distinct.

## Initially24 interventions, all offline

Use the interior merge(depth2,left5,right7,interval5:8) for influence first:
its record is subsequently read by prefix merges, avoiding a trivial
terminal value that no stage consumer uses. Four pairs x six substitutions:
retain previous full H, B_L, B_R or h separately (16 total);
then omit ONLY diagonal H change or ONLY offdiagonal B_L change (8 total).
Keep every other current merge input/local output block/vector, current
stage S/network factors/barrier/nominal and true terminal/future current.
These mixed-input substitutions may not be internally consistent Newton
systems; label them SENSITIVITY DIAGNOSTICS, not deployable reuse schemes.

Apply at the selected operation, then execute EVERY remaining production
scan round with dependency tracking. Trace exact affected suffixes and
recover the battery trajectory and full controls/equalities/bound multipliers
over ALL24 stages, including forward propagation beyond changed suffixes.
Do not regenerate an artificial terminal at a truncated window.
Compare intermediate output fields, full suffix gradients/Hessians,
their actions on saved actual SOC and additional normalized probes, then
current-stage reconstruction. Record first amplification location and
max coordinate/feeder outliers, not just global Frobenius norms.
If an intervention fails factorization, record it and preserve dense fallback
as an analysis outcome; do not silently replace it then call reuse harmless.

Conditioning/SVD/change-detection/construction/substitution/reconstruction
costs must be recorded as diagnostic costs, not clean Phase B benchmarks.
Assess which measured expensive operations a genuinely consistent update
might avoid, charging how its change would be detected/represented.
Stable-but-cheap, expensive-but-sensitive, promising candidate, or no useful
structure are all acceptable findings. No automatic implementation follows.

## Finish

Resource limits: sequential64CPU300GiB Gold6438Y+, bounded capture40minprocess/
90minbatch and bounded frozen replay<=2h initially; BLAS/Julia/OMP1,
private first depot plus readonly old cache, occupancyguard/no-requeue.
Implement/unit-check selectors and dependency propagation on tiny synthetic
arrays before cluster replay. Strictly verify capture physical/objective/
OpenDSS status, eight-state/all-stage coverage and hashes.
Archive diagnostic code, raw spectra/change/outlier/direction/condition/timing
CSVs, provenance/input/source/capture hashes, hardware/job accounting,
negative controls and failures. Binary captures stay immutable on cluster.
Commit/push concise final assessment and recommend only the specific next
experiment supported by evidence. Keep the queue heartbeat through final
collection/push; no handoff merely because a related harness is needed.

