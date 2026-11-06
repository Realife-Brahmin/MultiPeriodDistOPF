# T96 24x4 versus 48x2: validated complete solves

48x2 is faster overall in this paired experiment. Two-repeat warmed optimizer
medians are 26.7625 s (48x2; range26.715-26.810) and29.855 s
(24x4; range29.715-29.995): 24x4 takes11.56% longer. The ranges do not
overlap, so the predeclared overlap-triggered third pair is unnecessary.
Two repeats establish this observed comparison, not a confidence interval.

Phase medians in seconds (48x2 ->24x4): A4.036786 ->7.510736,
B16.557723 ->13.050924, C3.392792 ->4.777595.
Phase B decreases21.18% (3.506799 s), while A increases3.473950 s and
C increases1.384804 s. A/B/C do not exhaust the complete optimizer timer.
Reducing worker count while increasing BLAS threads benefits the dense merge
phase but costs stage parallelism; this experiment does not isolate those
mechanisms further. It supports retaining48x2 for this T96 configuration,
not a general worker/thread optimum or a centralized speedup claim.

## Matching and scientific validation

All four solves COMPLETE, T96.finished, Slurm31580246 COMPLETED0:0,
34m38s. All ran on cn179, Xeon6760P,128physical CPUs,600GiB allocation.
The user approved any matching node; the initial cn178 restriction and its
removal are preserved in the archive's scheduling records. Both settings
use97physical cores including the single-thread coordinator. Actual OpenBLAS
worker thread counts2/4 and OS masks were audited after setup/solve; each
team stays within one socket and groups use disjoint physical cores.
Every-merge and phase thread assertions passed. No expensive job overlapped.

Same original1020 individual batteries, frozen T96 input/source/start hashes,
Cb2.4e-5/gamma420, grouped b2/right kernel, profiles, screening, pass GC and
unchanged tolerances. Source/input/harness before/after checks pass.
Both arms/repeats have16accepted steps/solver iterations,25backward passes,
8mu reductions and0backtracks. All25 objective/physical-feasibility trace
keys match; objective differences are zero and maximum original-primal
trace difference4.55e-13. Same-repeat full solutions satisfy unchanged
same-iteration/scaled state/control1e-10 checks: state2.89e-15,
control3.52e-16. The48x2 repeats are bit-identical. Original numerical
primal/dual/complementarity/bound gates1e-5 pass; independent physical
residual<=1.14e-12, objective2704807.0669342116 in all four runs.
Independent OpenDSS reports zero lower/upper voltage violations.
Strict collector validation=true, no missing evidence or equivalence failures.

## Timer and memory scope

Coordinator model builds209.294-209.954 s, complete warmups107.925-119.200 s,
and whole per-run process wall504-551 s are separate from the26.7-30.0 s
warmed optimizer. Whole-process wall includes setup, warmup and validation;
it is not a warmed solve time. Full raw passes/trials/directions/quality/
iteration timing and worker audits are retained. analysis/gc.csv reports
summed worker GC and sums of per-pass worker maxima; neither is an additive
wall-time attribution and no GC causality is inferred.

48x2 coordinator lifetime peak RSS2964-3142MiB, individual worker max
2294-2332MiB, sum of worker lifetime maxima103249-103553MiB.
24x4 coordinator3191-3195MiB, individual worker3084-3091MiB,
summed maxima71432-71475MiB. These RSS counters include setup/warmup and
shared mappings; sums are neither simultaneous nor private RAM measurements.
Slurm accounting is also preserved with its reporting limitations.

## Reproduction and evidence

Run collect.py EXTRACTED_ROOT OUT without --partial. results/analysis retains
per-run/median/solution/trajectory/validation outputs and GC aggregates.
The raw archive contains logs/CSV/audits/masks/hardware/manifests/solution
hashes/accounting/scheduling; binaries stay in the immutable cluster snapshot
/home/aryan.jha/kamiak_jobs/dense_merge_teams96_nov05.
No fresh centralized solve, production-default, paper or main benchmark change.
The numerical influence audit continues separately; this result does not
complete that audit or authorize a compression implementation.
