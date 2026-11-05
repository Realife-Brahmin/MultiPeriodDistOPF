# Matched 128-core and T24 follow-ups

Authorized November5: complete original large10k1020-battery solves.
Production defaults, published benchmark rows, paper and MSOPF stay unchanged.

Run sequentially after31575217 in two new immutable snapshots:
* dense_merge_T24_nov05: Gold6438Y+,64cores300GiB4h. T24 original
  balanced/current24x1 versus groupedb2/right24x2, plus MA57 and MA97.
  Two repeats each, reversed order. min(T,workers) means24x2, not31x2.
* dense_merge_core128_nov05: Xeon6760P,128physicalcores600GiB6h.
  T96 groupedb2/right63x1,96x1,48x2, plus MA57 and MA97.
  Two repeats each, reversed order within the SAME node/allocation.
  Centralized BLAS/OpenMP stays1 per inherited matched protocol;
  this is not a centralized threading sweep.

Reserve coordinator core; pin disjoint physical-core groups with each BLAS
team within one socket. Record actual OS-thread masks, topology, backend and
phase/merge thread assertions. A cross-generation comparison cannot isolate
core count. One worker per timestep is a hypothesis, not an established optimum.

Reuse exact validated numerical source and original frozen T24/T96 inputs.
Only centralized runtime audit logging and the allocator core budget are new.
Full warmup, physical start, profiles/screening/Cb2.4e-5/gamma420/passGC and
stopping criteria stay inherited. All arms require original numerical,
physical/objective/bound and independent OpenDSS checks. Preserve failures.

Score complete warmed optimizer time separately from setup/warmup/processwall.
Also score objective gap<=0.1% AND original feasibility<=1e-5 against mutually
agreeing fresh central references; use faster MA57/MA97 two-repeat MEDIAN.
The historical0.86 ratio needs >14% DDP time reduction if central time stays
fixed. Do not recolor old rows or conflate near-optimality/full convergence.

Before final reporting: all COMPLETE markers, clean stderr, successful Slurm,
source/input/harness hashes, finite full traces, physical residual<=1e-5,
central objectives agree<=1e-5relative, every final objective<=0.1%gap,
matched physical initial objectives, identical DDP start-array hashes.
Use repaired compare_solutions.jl for DDP candidates versus same-repeat DDP
reference: same iterations and scaled state/control<=1e-10. Preserve any false
gate; do not automatically promote. Verify actual thread masks against pinned
physical-core maps, not only requested settings. Add third paired repeat only
if overlapping ranges require it. No automatic further tuning.

Keep raw logs/audits/solution hashes/manifests/accounting, binaries on cluster.
Safe repairs only new immutable attempts, never edits to queued/running source.
afterok dependencies serialize expensive jobs; inspect predecessor failures.
Commit/push final evidence and assessment before removing completion heartbeat.

