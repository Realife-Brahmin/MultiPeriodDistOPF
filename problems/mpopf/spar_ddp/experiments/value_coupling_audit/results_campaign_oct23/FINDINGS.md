# Decision: revise the consumption architecture; do not integrate this candidate

The coverage gap is closed and direct generation/composition is accurate on
the complete four-stage workload. It does NOT provide a worthwhile fresh-input
cost benefit. Stop this bounded implementation campaign here; do not promote it
to a compressed optimizer or automatically extend cases/horizons. This is a
specific implementation/consumption bottleneck, not a disproof of compression.

## Verified inputs and accuracy

New capture31369462 reproduced the original1020-battery large10k late Newton
state BIT-EXACTLY: all24 leaves, original overlapping nominal primal/dual state,
k14/j8/mu, and all saved full control/equality probes. New full elimination
records1:4 close the old missing2:4 gap. Original-model physical residual
6.82e-13 and independent OpenDSS discrepancy8.18e-9pu/no violations validate
the dense capture, NOT an unrun compressed optimizer. Both failed capture
attempts and their exact setup/checker causes remain archived separately.

Job31369517/overlay3f1e2e06 COMPLETED0:0 in3m10 on Gold6438Y+,64CPU300GiB,
Julia1.12.4/BLAS1/OMP1, no competing user jobs. Strict collector validates
source/runtime/capture/input hashes, synthetic units, full stdout/stderr,
accounting, actual four-stage coverage, all repeats and zero compilation AND
recompilation in EVERY measured phase.

All FOUR initial factors are generated from the matching elimination records,
not dense boundary SVD. Poisoning saved Sh0 leaves generation unchanged. All
three raw blocks/two vectors are checked; raw skew remains explicit. Original
coordinate upper completion is retained. All four closed value operators are
compared on the full coordinate basis outside execution timers; timed actions
use three probes. All four stages consume zero/sine full control, equality and
bound-multiplier directions, plus the battery trajectory. No physical battery,
constraint, location, strong local block or full vector was removed.

Maximum scaled value error1.82014e-8; full-direction error8.29262e-8 and
relative/floor1e-10 error1.01104e-6 pass unchanged1e-7/1e-5 gates. All four
original-upper inertia checks match. Corrupting the composed gradient yields
1.75761e7 scaled error and is rejected. Every timed repetition passes again.
The independent consumer with DENSE future already has8.292617e-8 error:
the limited margin is substantially a consumer/reference sensitivity, not
evidence of compression alone. In particular stage2 sine equality recovery
has8.22451e-8 scaled error. Do not infer robust late-iterate accuracy from
passing this single point, or claim tighter SVD alone would fix sensitivity.

## Complete fresh-input warmed costs

THREE alternating identical-input pairs in ONE process. No numerical factors
are reused across Newton iterates or repetitions. Both methods reconstruct the
true incoming V5 from20 frozen tail leaves, preserving GLOBAL production scan
association. Candidate additionally pays its boundary factor discovery. The
four-period window is NOT a fully compact24-period horizon.

Median selected-window fresh total, INCLUDING generation, discovery,
combinations, condition probes, all actions/consumers and cleanup:
dense2.66357s versus compact2.70074s (compact1.40% slower).
Whole charged median including true-future construction:
dense8.91349s versus compact8.91469s (effectively tied, not a speedup).
Paired whole totals: dense/compact8.91349/8.83259,
8.97365/8.91469,8.83041/8.95223s. Thus not all pairs improve; the predeclared
>=10% window saving and consistently lower whole-total gate FAILS.

Phase medians (dense vs compact where applicable): initial generation
0.08806/0.16003s; candidate incoming discovery0.13553s;
value construction0.89438/0.26694s; candidate condition probes0.01145s;
all value actions0.00328/0.09056s; all-stage consumers1.17704/1.49504s;
cleanup0.50753/0.51954s. Common incoming construction6.24992/6.20918s.
Per-phase medians need not sum to median totals. Setup INCLUDING incoming
construction/generation/discovery6.33799/6.50951s; remaining charged recurring
work2.57551/2.39911s. Neither category is omitted from fresh totals.
Without explicit cleanup, whole medians8.41269/8.36858s still do not establish
a worthwhile selected-window benefit or satisfy the predeclared gate.

The important mechanism: direct generation removes the previous~1.42s dense
factor discovery bottleneck, and compact value construction is~3.35x faster.
But this saves only~0.627s here; slower full-stage consumption (~0.318s), value
actions (~0.087s), generation/discovery/conditioning and cleanup absorb it.
The earlier stage1-only consumer result was insufficient coverage, not a
whole-workload benefit that can be extrapolated to all stages.

Cumulative allocation1.756GB dense versus2.565GB compact (~46% higher),
including common tail. Median recorded in-phase plus explicit GC0.559/0.679s.
All allocations, recompressions and dense CROSS-block fallbacks are charged.
There is no algorithm-required dense boundary export; dense validation-only
reconstructions do not feed execution and are outside timers.

## Structure, memory and scope

Generated stage separator width6, retained initial rank3, local elimination
dimension20; each stage performs1884 small-core SVDs (max core8) and966--972
dense cross fallbacks. Updates record6344 core SVDs, max rank18/core28,
Woodbury core36, max local dense60 and4860 cross fallbacks. These are observed
structure, NOT fixed-rank/horizon-independent bounds or closed separator algebra.

Retained suffix dimensions1020/2040/3060; coefficient bytes1,157,600/
3,240,960/6,072,000. Initial coefficients3,370,384bytes; whole object graph
17,318,712bytes. These are NOT resident-memory savings. Sampled batch MaxRSS
3,776,588KiB includes BOTH methods, captures, references and validation.
The stage consumers still use saved phaseA/network factors with dense Sh0;
generator independence does not remove that matrix from the whole solver.
Incoming future remains constructed densely, not supplied free. Structural
topology setup0.869s is separately recorded COLD (includes possible JIT), not
a warmed benchmark; adding first-use setup cannot create a missing fresh win.

## Recommendation and stopping rationale

REVISE before any optimizer integration: the strongest next design question
is a reusable factor-backed full-stage consumer that avoids repeated iterative
Schur actions and recovers sensitive equality/bound directions economically.
That is a different integration/consumer design, not another narrow tolerance
or discovery tweak. Known-factor generation is now measured, not missing.
No single remaining discovery cost dominates as before; a blind retry or
another small timing sweep would not establish the required benefit.

The current accurate path passes only near the direction threshold, allocates
more and misses the declared total-cost gate despite favorable low-rank
large10k structure. This supports stopping this bounded implementation,
preserving it as an accurate compact prototype. Medium early/late, longer
composition and opt-in optimization did NOT pass their prerequisite and were
not run. No compressed feasibility/objective/iterations, near-optimality,
full-convergence, parallel PhaseB wall, communication-only or optimizer-speedup
result is claimed. Published rows, production defaults and paper unchanged.

## Reproduction and provenance

Run analyze_campaign.py EXTRACTED_REPLAY_ROOT NEW_OUTPUT_DIRECTORY.
Strict outputs in analysis/ include validation.json, totals.csv, phases.csv
and hashes. replay_evidence.tgz SHA256
491720c3171d8f8650069e34bb4f031df70a66bec91772ba6e9144ebff9ec40f
contains complete replay source, logs, raw results, manifests and accounting,
excluding only the private dependency depot. Capture binary hashes remain
pinned to immutable cluster paths; capture_verified_evidence.tgz retains
textual/metadata/runtime/physical/OpenDSS evidence. Failed snapshots are not
overwritten or relabeled. No expensive jobs remain active at collection.
