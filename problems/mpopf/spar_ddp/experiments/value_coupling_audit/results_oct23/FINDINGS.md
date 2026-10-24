# Bounded block-scaled follow-up: validated October23 assessment

The controlled follow-up is complete. Improved cross-block error control
substantially improves the original1020 late-window directions, with useful
coefficient compression. It does NOT yet supply an accurate general-purpose
compressed solver or a faster implementation. No optimization comparison was
launched: medium late missed the predeclared accuracy gate, and warmed total
cost was higher than dense at both qualifying points. This is a specific
limitation of the tested representation/control/implementation, not evidence
that value-function compression is inherently unsuccessful.

Production defaults, physical batteries/constraints, existing benchmark rows
and the paper are unchanged. No battery-group dispatch or larger campaign.

## Scope, completion and reference checks

Job31368755, immutable overlaye2a86a50 in value_block_followup_oct23_v4,
COMPLETED0:0 in3m32 on Gold6438Y+,64CPU300GiB,Julia1.12.4/BLAS1/OMP1.
Strict collection validates all three requested frozen points: medium late and
early fullT24 suffixes, and original1020 large late first4period window with
the TRUE production-dense future. The latter is NOT a fully compressed T24
trajectory. All saved dense full-control/equality probes reproduce exactly0;
ordinary versus production value errors are below1e-7. Overlay, inherited
source, legacy source, input and capture file SHA checks pass; stderr is empty.
All eight accuracy trials retain the recorded upper-completed elimination
inertia classification relative to dense. Numerical eigenvalue thresholds
are diagnostic, not interval certificates or convergence proofs.

Failed setup/measurement attempts and compiler-cap partial evidence remain
preserved (RUN.md). The fixed recursive payload container avoids specializing
on entire nested tree/rank-choice types. Synthetic four-budget and real
medium-late leaf0/1e-6 comparisons against pinnedf8f60d39 pass bit-identical
matrices, ranks, error/fallback and coefficient counts. No numerical threshold,
SVD arithmetic or physical input changed. The earlier22minute compiler-cap
attempt contributes NO runtime benchmark observations.

## Accuracy and storage

Raw Maa/Mab/Mbb are compressed independently using each directed hierarchy
cross block's own Frobenius scale. Local blocks and both vectors stay exact;
raw skew is not projected away. Conditioned tolerance is base/sqrt(kappa),
where kappa is the Inf-condition of the actual upper-completed battery-state
elimination matrix; initial leaves use Maa+Mbb as an explicit heuristic.

Direction gates are norm(delta)/max(1,norm(reference))<1e-7 AND
norm(delta)/max(1e-10,norm(reference))<1e-5 over battery paths and selected
controls/equality/bound-multiplier zero/sine probes. Candidate selection uses
dense direction comparisons OFFLINE, not a deployable free oracle.

Original1020 large late: block-scaled1e-6 already passes with maximum scaled
direction error1.622e-9. Conditioned1e-6 passes at6.133e-10 (relative-floor
error7.182e-10), maximum cross rank30 and coefficient ratio0.0344112, about
29.06x smaller. This preserves nonzero cross-group interactions rather than
dropping them. The prior joint-global1e-8 prototype had0.001716 scaled error
on this same late window; the new result is encouraging, but combines changed
budget scaling/raw representation and is not an isolated attribution study.

Medium early: block-scaled1e-6 misses the scaled gate at4.966e-7. Conditioned
1e-6 passes at3.996e-8 (relative-floor1.119e-7), maximum rank9 and coefficient
ratio0.244300, about4.09x smaller.

Medium late: none of the four trials passes. Maximum scaled errors for plain
block1e-6, conditioned1e-6/1e-9/1e-12 are respectively2.382e-7,2.237e-7,
9.559e-7,3.793e-6. Payload ratios grow0.19216/0.20918/0.27431/0.34188 and
maximum ranks11/11/17/109. Tighter truncation does NOT monotonically improve
the recovered direction. The plain trial's worst field is stage1 sine-probe
equality direction; conditioned1e-6's worst is upper-bound multiplier direction
(absolute2.237e-7, reference norm0.2581). The tighter trials' worst fields are
again equality directions. This is not a demonstrated optimization failure;
it fails the explicitly chosen conservative frozen-direction threshold.
Raw dense remains the fallback, not a compression success.

The current battery-state condition heuristic does not capture all amplification
through full-stage policy reconstruction. Roundoff/conditioning sensitivity
is consistent with the nonmonotonic errors, but its mechanism has not been
isolated. Do not assert that ever-smaller SVD budgets alone solve it. No whole
element raw reconstruction fallback fired in these eight trials; unprofitable
individual cross blocks still use dense storage (rank_stats.csv records it).

Coefficient ratios are arithmetic means over represented PintElem coefficient
counts (leaf/intermediate elements), not the old joint prototype's median
convention. They exclude tree/index/object overhead and dense closed PintVal
suffixes. They are NOT resident-RAM ratios or whole-PhaseB storage ratios.

## Warmed costs: separate from cold trials and Phase B

Only the two passing conditioned1e-6 candidates were timed after full-path
warmup. Each has three identical-input dense/candidate pairs IN ONE PROCESS,
alternating order dense-first/candidate-first/dense-first. Direction checks,
file I/O and explicit boundary GC are outside timers. Condition estimation,
representation construction, materialization, dense combination, allocation,
GC, workspaces/deepcopy and in-memory timing-record overhead are charged.
The offline candidate-selection oracle is not included in warmed execution.

Large late: dense median0.891338s; candidate median5.237608s (range5.231897--
5.590300), or5.88x the dense time. Component medians are approximately
condition0.7443s, construction3.3566s, materialization0.1928s and combination
0.9019s. Candidate median cumulative allocation2.709GB versus dense0.1333GB;
candidate total measured GC varies0.0484--0.4060s versus dense0. Three repeated
direction checks remain at6.133e-10 scaled error and the same ranks/payload.

Medium early: dense median0.236605s; candidate median2.079740s (range1.727441--
2.080953), or8.79x the dense time. Component medians are condition0.1844s,
construction1.2087s, materialization0.1281s and combination0.2355s. Candidate
median cumulative allocation1.357GB versus dense0.03783GB; GC varies0.0335--
0.3856s versus dense0. The repeat spread is retained, not hidden as compiler
cost or replaced by the fastest observation. Three repetitions establish this
bounded component comparison, not universal timing/noise characteristics.
Component medians need not sum to the median whole-call time.

Medium late has NO accepted candidate and NO warmed comparison; do not use its
cold accuracy-trial times. Full batch Slurm sampled MaxRSS2548712KiB covers
setup/reference/audits and both representations, not candidate-only peak RAM.
Cumulative allocations above are not resident memory. The prototype holds
compact payload AND dense expansion, and still performs dense production
combinations. This is NOT genuinely compact computation, parallel PhaseB wall
timing, communication-only timing, or an end-to-end solve speedup measurement.

## Decision and smallest next discriminator

Do not run a full optimizer with this prototype: cost evidence fails even at
the accurate points, while medium late needs a more reliable direction-error
strategy. Feasibility, objective agreement, iteration counts, feasible
near-optimality and full-convergence times for a compressed optimizer remain
UNTESTED, not failed or inferred from the dense capture's OpenDSS checks.

The concrete next engineering target is a factor-aware compact combination
that avoids repeated SVD discovery plus dense expansion, with reuse of existing
separator factors and a cheaper condition estimate/reuse of factorizations.
Retaining every physical battery and strong local block remains mandatory;
fixed low rank or horizon-independent separator rank is not established.
At the same medium-late frozen point, error control must additionally account
for full-stage/Schur-direction sensitivity and explicitly preserve raw dense
fallback where accuracy cannot be certified. Test these kernels on the same
saved points before any further solver campaign. Only an accuracy AND total-
cost pass should trigger bounded convergence comparison, then separate clean
near-optimal/full timings including extra iterations and representation cost.
No such follow-on kernel campaign is automatically launched by this report.

## Evidence

v4_evidence.tgz preserves every raw CSV/log/source overlay/hash/hardware and
accounting record, SHA256
e84955b1d587e6d7dff3cee40cbce7d7da42dc8dd93c8c592e946d62c8e26d01.
Reproduce with analyze_followup.py EXTRACTED_ROOT OUTPUT (strict, no --partial).
analysis/validation.json reports3complete points and2warmed points, not three
compressed optimizer successes. accuracy.csv,field_errors.csv,rank_stats.csv,
warmed.csv,components.csv,medians.csv and evidence.sha256 retain all outcomes.
