# Full-solve value grouping diagnostic, October22

Job31356155 COMPLETED34m12s/exit0, cn164 Gold6438Y+,64CPUs,300GiB,
Julia1.12.4/OMP1/BLAS1. Sourcebae4f9d7. All3 IEEE123T3 smokes and all9
T24 full cases passed convergence/original residual/bounds, independent
physical objective/residual re-evaluation, and OpenDSS. The24 existing raw
numerical/OpenDSS validators were rerun locally and passed. No paper changes.

## What was actually tested

All physical batteries retained:51/249/1020 for small/medium/large. No changed
loads/PV/network, no battery aggregation or altered ratings. Dense is the
unchanged scan; root operates independently on branches immediately below
the root; chunks16 partitions sorted battery bus indices into groups of at
most16, deliberately cutting possible strong interactions. Approximation is
ONLY in the two-sided scan quadratic, not the physical model. Original stage
systems still use exact curvature and the approximate future in Phase C.

All modes use same periodic exports,physical starting values,screening,cost
coefficients and soft terminal penalty;24 workers atT24/two atT3. Every arm
fully warmed. One exploratory solve per mode, not two-repeat benchmarks.
Source/input manifests and hardware/affinity are in the raw archive.

## Outcome: topology matters, and approximation can cost iterations

Small: root partition has ONE group of51 batteries; it therefore makes no
approximation. Dense/root both22 accepted iterations and32 backward passes.
The chunks16 stress partition has4 groups; it converged in31 iterations and
41 passes (+41% accepted iterations), objective relative difference1.05e-10
versus dense. Original physical residual2.60e-16.

Medium: root again has ONE group of249, not a nontrivial approximation.
Dense/root both52 iterations and61 passes. The16 stress groups converged in
76 iterations and85 passes (+46%), objective relative difference4.52e-7
(0.0000453%), original physical residual7.65e-16.

Original large10k: root produces102 groups of10. Dense/root both14 iterations
and23 passes; final objective difference1.74e-16 relative. Their printed
iteration/objective/primal/dual/complementarity trajectories agree to2e-16
relative under the collector's max(1,abs(reference)) scaling. This supports
the network-specific separability concern; it is NOT a general exactness proof.
Chunks16 creates64 groups and converged in21 iterations/30 passes (+50%),
objective relative difference1.12e-10, original physical residual6.82e-13.

All objective comparisons are against OWN same-system dense DDP reference,
not a newly run centralized solution or a claim of global optimality. All
grouped arms meet the0.1% agreement threshold, with no changed tolerances.
The small/medium root arms are controls for grouping overhead, not evidence
that a genuinely compressed representation works on those networks.

## Value-function errors: not universally negligible

The dense same-iterate ordinary-recursion audit passed<1e-7; largest printed
error over all dense cases about1.45e-14. Stress approximation first-pass
Vxx/Vx errors (relative infinity norms with denominator max(1,norm(reference)))
are about0.00191/0.00540 small,0.19465/0.38083 medium,0.00290/0.00107 large.
Medium atiteration5 had0.05052/0.02441 errors. Thus favorable final objectives
do NOT mean the value functions are identical. In particular the medium case
has significant first-pass errors and extra iterations, unlike the tiny
cross-root interactions in the original large10k first-pass replay.

IMPORTANT AUDIT COVERAGE LIMITATION: the harness selected k0/5/10/15 only
when barrier_iteration j==0. That condition excludes later barrier stages.
The actual saved/printed coverage is k0 only on small/large, and k0/5 on
medium. No late-barrier or final-point value-function error is established.
The full solves genuinely completed, but the intended early/mid/late audit
was incomplete. Before claiming structure persists, a separate immutable
diagnostic must audit selected k across all relevant j, or re-evaluate saved
states independently. Do not silently edit this completed snapshot or infer
later accuracy from endpoint convergence. No additional campaign auto-submitted.

## Timing and memory: diagnostic only

Audit/capture work occurs inside the printed optimizer clock, kernels use
dense storage, and this is one repetition. Raw logged clocks are0.279/0.281/
0.356s small,4.436/4.811/4.631s medium,25.028/8.019/10.095s large for
dense/root/chunks16. These are NOT clean benchmark speedups. Preserve them
as diagnostic observations only; no main-paper replacement rows are justified.
Scan totals also include reference-build overhead. Later audits were skipped
as above, so these numbers cannot be extrapolated to a fully audited run.

Shared arrays remain full-sized and group subblocks are copied into temporary
arrays: no compact-storage backend or memory-saving claim. Slurm batch MaxRSS
39417912KiB (~37.6GiB) covers the campaign allocation accounting, not a
backend-specific or optimizer-only simultaneous peak. Per-process raw memory
and pass/GC clocks are retained in logs.

## Decision / next smallest step

Natural-root splitting is promising specifically where graph/model separators
actually decouple batteries. First establish the exact conditions (fixed root
voltage, absence of aggregate nonlinear costs/shared limits, separator size),
rather than labeling every root split exact. Stress partition shows a genuinely
inexact value representation can still converge here, but may need~1.4--1.5x
as many accepted iterations. Neither result establishes broad reliability.

Next: fix the audit sampling condition in a NEW diagnostic snapshot and use
early/mid/late original-state captures to test structure/error/inertia and
direction changes. Only then design a compact exact separator backend or
adaptive error-controlled compression and run clean repeated timing comparisons.
No new physical test network, horizon escalation or centralized retuning yet.

## Evidence

Cluster immutable snapshot:
/home/aryan.jha/kamiak_jobs/value_group_solves_oct22.
All57 saved dispatch/capture arrays remain there (~285MB results). Their SHA256
manifest is saved_arrays.sha256. Binary arrays are not included in Git; this
is not an additional off-cluster backup of them. evidence.tgz includes raw
text logs, group outcomes, objective comparison, source/input manifests,
hardware and gate markers; SHA256:
41d64861f5ca38999685656ccd7f74a2fcbeeebd83dc584589fa4b760d38fe23.
Slurm completion/accounting is in slurm.psv. Extract to an empty directory;
run `python ../analyze.py EXTRACTED REPO_ROOT OUTPUT` to reproduce per_run.csv,
reference_checks.csv,passes.csv,validation.json and rerun24 raw validators.
campaign.finished alone was not used as evidence that all approximations worked.
