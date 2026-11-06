# Exact scan barrier cache: correct extension, negative charged-cost result

The two-sided scan obeys the established barrier-only curvature invariance
throughout the captured fullT24 horizon. An exact process-parallel cache works,
including current gradients and complete recovery, but this conservative pilot
should NOT be integrated into the optimizer: final five-pair charged median
5.204454s versus fresh3.636563s,43.1% slower. The bottleneck is exact content-
fingerprint detection/allocation and cache management, not matrix reuse accuracy.

## Existing optimization and source conditions

This extends the regular-DDP optimization in SPARDDP.jl/src/reuse_pass.jl;
it is not a new mathematical claim or a repeat of that sequential study.
Its solver/k/Hessian bookkeeping and barrier-only control-flow reasoning apply;
its sequential factors are not interchangeable with new scan merge factors.
procs_pass.jl explicitly omitted that sequential reuse path.

SOURCE_ARGUMENT.md traces the actual dependencies. With fixed x/u/phi/zl/zu,
slacks/bounds, model, dynamics/derivatives, terminal curvature and execution/
regularization path, Sigma=z/slack is independent of mu. Changing mu modifies
mubar, gradients/RHS and complementarity recovery, not those matrices.
The pinned _tp_pre/tree path has no added regularization retry/shift; failed
Cholesky uses unshifted LU. This source fact establishes the supported applied
unshifted path; the subsequently reset reg_last field is NOT evidence for it.
Unknown or changed regularization disables reuse. Accepted steps, hybrid/
Hessian changes, matrix/factor/layout/model/source/backend or terminal-curvature
changes invalidate. Soft terminal-SOC curvature is in the last leaf; the true
scan terminal25 has zero V/v, not an artificial truncated-window boundary.

## Complete global evidence and validation

Original1020-battery large10k, fixed original balanced/in-place/original-UPPER
scan, same immutable k4/j0->j1 capture (mu1->.2). Both fresh leaf reassemblies
match saved coefficients and full saved live probes. ALL24leaf matrices,
terminal curvature and EVERY one of68 merges' matrix inputs AND outputs are
bit-identical:520 data-driven matrix comparisons plus the explicit zero
terminal-curvature invariant (521 rows overall), including raw asymmetric inputs
and external diagonal blocks. Exact stage/factor signatures also match.
All24suffix gradients/Hessians, value actions, complete battery forward path,
full controls, equalities and lower/upper-bound directions pass unchanged
scaled1e-7 AND relative/floor1e-10 threshold1e-5 gates. Correct population,
hit, warmup and every measured output have ZERO numerical differences.

Seven controls: an actual ordinary accepted k0->k1 pair, changed Sh0 input,
missing merge entry, unknown applied regularization, model/execution context,
and changed terminal curvature all invalidate and reproduce fresh outputs.
Deliberately stale gradients are rejected: max scaled trajectory error1.388.
A source-kernel unit additionally forces failed Cholesky and verifies unchanged
fresh LU fallback with no stored entry. No nearly-unchanged eligibility rule.
Raw absolute/scaled/relative errors and coordinate outliers are retained.
Current reconstruction uses current S/barrier records and rebuilt policies;
dummy lambda's reported scalar is not claimed as an original dual residual.

## Minimal dense-BLAS layout

Each of44 two-boundary entries stores U,Y1,Y2; each of24 closing entries
stores U,Y1,V. Exactly204dense1020x1020 matrices:1,697,932,800numeric bytes
(1.58GiB); worker cache object graphs1,697,959,872bytes, plus a small coordinator
stamp and allocator overhead. No leaf/internal curvature blocks, matrix inputs,
gradients or old feedforward policies are retained. U is kept in dense layout
for the original triangular kernel, not asserted to be a universal packed lower
bound. This is half the initial prototype's numeric storage.

An eligible complete hit needs no intermediate matrix materialization: current
ma/mb recurrences use cached U/Y1/Y2, and every final closed value restores V.
It skips68 Cholesky factorizations,112matrix triangular solves and156matrix
products;68current vector triangular solves and112gradient matvecs remain.
Counts follow the verified DAG/entry types, not hardware performance counters.
All current leaf gradients, h, both gradient updates, suffix gradients,
terminal gradient, policies and full recovery are refreshed.

Whole-cache completeness is checked before reuse. A missing/failed entry
makes the whole pass fresh. Unexpected mid-pass entry loss abandons speculative
work and restarts from current leaves with unchanged fresh kernels, never
fresh-combining unmaterialized scratch curvature. Current-source fingerprints
are exact and conservative, not free trusted version counters.

## Fully charged timing and memory

Final job31609146 completed0:0 in7m11 on cn168, Gold6438Y+,64physical CPUs,
300GiB allocation;24processes x BLAS1 plus coordinator on disjoint physical
cores. Launch occupancy guard and all-users node accounting show only the two
sequential pilot jobs in the study window. Same balanced merge assignment,
shared storage and actual process round synchronization as the captured scan.

Both full workloads warmed, then FIVE alternating baseline/candidate pairs.
Every clean coordinator AND worker call has zero compilation/recompilation.
Fresh first+fresh second versus fresh INCLUDING population+eligible hit.
Median charged totals3.636563/5.204454s; ranges3.307887-3.711062 and
4.845163-5.774132. Every paired candidate is slower; no confidence interval
or production speedup is claimed. First-pass medians1.871550/2.719881s;
second-pass1.735257/1.981728s; retirement/GC.073443/.500699s. Component
medians need not sum to the median total.

On the HIT pass, scan wall falls.759100->.016196s (~46.9x scan-only), but
exact fingerprint detection costs.891971s. Population detection costs.905952s;
population scan including allocation/writes.800298s versus fresh.755416s.
Current consumer medians on the second pass.697852/.827043s and forward
recovery.266034/.265592s remain charged. No free cache preparation, lookup,
copies or synchronization. Complete worker allocation deltas9.438/22.421GB;
summed worker GC4.632/13.666s overlaps across processes and is NOT additive
wall time. Allocation counters include small counter-RPC overhead. Guard
serialization materializes byte buffers; per-operation allocation attribution
was not independently profiled. Batch MaxRSS55,180,452KiB includes processes,
shared mappings, captures, references/JIT/diagnostics and is not a private cache
increment. Retained cache bytes are the direct storage measure.

Common frozen inputs are current pre-phase-C stage records. Their loading and
original phaseA derivative/network-factor construction are outside BOTH clocks.
Timed work includes leaf assembly/gradient refresh, scan, current policy
construction, full forward directions, all cache detection/population and
retirement GC. Validation/probes/outliers are outside clean clocks. This is a
complete frozen two-pass workload, NOT a new whole optimizer solve. Cold
package/runtime/process setup is separately logged, not a warmed optimizer time.

Initial job31608936 is preserved in results_v1: exact but over-retained leaf/
intermediate outputs (3.396GB), charged5.436s versus3.251s. V2 removed those
unneeded copies and added direct population/full-output and missing-entry
checks. Versions are not pooled or compared as a cross-batch timing win; final
conclusion uses V2's fresh controls. No failed scientific gates were loosened.

## Frequency, solve relevance and recommendation

Two existing matched original24x1/BLAS1/balancedT24 full solves have23passes,
8 source-eligible transitions and15fresh passes. Eligibility requires prior
outcome=barrier_update, consecutive passes, same k, next j, decreasing mu,
no accepted step between, fixed source/backend/Hessian/model and unshifted
route, and no configured hybrid switch. This is more than counting mu drops.
Only the representative k4 transition has complete captured numerical keys;
we do not pretend logs fingerprint the other seven states. Those remain
source-conditional opportunities, with successful complete cache factors also
required. The projection baseline is this original kernel, NOT latest grouped/
right/BLAS2 results or a different horizon.

Eligible fresh scan work is6.007599/6.016921s of22.595/22.746s optimizer time:
even making ALL eligible scan work free has only~26.5% whole-solve reduction
ceiling. Under the explicit optimistic assumptions that all8transitions pass
exact keys/Cholesky and representative costs apply throughout, the measured
hit scan ratio gives~5.884s gross scan saving. Charging population scan writes
on15fresh passes adds.673s; exact detection on15fresh+8hit passes adds20.725s.
Net scan/detection balance is about15.51s WORSE, before additional GC/current-
consumer effects. projections.csv retains both source-run scenarios. This is
an accounting estimate, not a measured optimized-solve time or forecast.

REVISE, do not integrate this implementation. The supported next experiment is
exact source-controlled curvature epochs/narrow state-version checks, using
the established solver bookkeeping, instead of serializing whole stage/factor
objects for eligibility. Mutation/invalidation must remain conservative and
negative-tested; do not substitute tolerance-based equality or grant detection
free. Re-run the same bounded charged pilot if that design is pursued before
optimizer integration. Matrix invariance and vector refresh are now verified;
another compression campaign or ordinary-iterate reuse is not supported here.

Code, source/runtime/input SHA,521global comparisons, all current errors/
outliers, seven controls, five pairs, storage, compiler/allocation counters,
affinity, hardware, accounting and both implementation outcomes are archived.
Run collect.py EXTRACTED_ROOT NEW_OUT, frequency.py ROOT/frequency NEW_OUT,
then project_benefit.py ANALYSIS FREQUENCY_ANALYSIS NEW_OUT. The raw archives
include pinned lib helpers; prepare.sh assembles a NEW immutable snapshot from
existing captures/harnesses. Binary captures remain on Kamiak. No new capture,
horizon/network/thread sweep, approximate optimizer, production default,
paper benchmark, MSOPF change or full optimizer integration occurred.
