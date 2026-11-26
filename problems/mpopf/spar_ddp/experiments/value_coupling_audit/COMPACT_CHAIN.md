# Bounded compact-to-compact continuation, October 23

Base: ef43e34e. Only the immutable original1020 large10k late T24 capture,
first four periods, true dense future at period5. No optimizer, medium case,
paper, benchmark-row or production-default changes.

## Representation and checks

Reuse separate raw Maa/Mab/Mbb topology factors (input budget1e-9) and full
vectors/local blocks. First merge3/4 returns CKSchur. Prepending period2 takes
that object directly, not its dense export. It retains the two eliminated
time slices in a spatially hierarchical block-tridiagonal operator, refactors
that enlarged operator and provides all three boundary actions. Closing with
the true future adds a third interior. This is retained-interior composition,
NOT proof of fixed-size/rank boundary closure. Dimensions grow1020/2040/3060.
The previous interface representation is reused; its LU is NOT claimed reused
inside the enlarged factorization. Initial factor generation still uses SVD
on the frozen dense stage matrices and dense true-future boundary.

Original-coordinate upper completion occurs explicitly on the ORIGINAL raw
diagonal sums, never by averaging their skew. Lifted order is time-major in
original coordinates and spatially permuted only for storage/factorization.
Raw cross matrices and transposes remain distinct. A Schur action is evaluated
algebraically; floating-point differences from production's copytri!/Cholesky
association are measured, not presumed zero. Validation compares the three
quadratic blocks/two vectors at both merges and the final value with both a
matching dense chain and the unchanged production suffix. No intermediate
boundary export feeds execution.

QR/core-SVD recompression budget1e-12, growing ranks and full dense cross-block
fallback are charged. Local factor cap64 normally; rank-growth fallback may
factor a bounded <=192 local submatrix, explicitly logged/charged. A larger
rank-growth rejection is preserved, not silently densified globally.

Full-stage consumer uses saved network elimination factors and the raw dense
stage Sh0 (existing input, not a reconstructed future). Its grouped-coordinate
upper triangle matches production. A bounded matrix-free diagonally
preconditioned CG solves the energy Schur system using future operator actions.
Tolerance2e-13, true residual recheck, max300 iterations; non-SPD/breakdown/cap
rejects the consumer. This tolerance is NOT a replacement for the unchanged
direction gates1e-7 unit-floor AND1e-5 relative-floor1e-10 on battery/control/
equality/bound directions. The consumer must first pass using the actual
unchanged dense production future. Corrupting the composed gradient must fail
the same direction gates. No network-direction claim from battery-only checks.

## Cost experiment and limits

After accuracy passes, warm exact paths and run three alternating dense/
compact identical-input pairs in one process. Charge fresh discovery (all
three raw elements plus terminal factors), both compact merges, enlarged
factorizations/terminal closure, condition probes, downstream directions,
allocations and GC. Dense comparison uses the same two merges/closure and
full policy plus needed battery suffix consumers. No unused suffix1 work.
Report fresh totals and conditional within-iterate reuse totals separately;
the latter subtracts only discovery, not updates or downstream consumption.
Topology setup is separately recorded cold, not benchmarked. Explicit boundary
GC is separately charged to both paths. Validation-only dense exports and
eigenvalue checks stay outside execution timers and never feed execution.
Necessary algorithm conversions would be charged; the intended path has none.

Coefficient counts, object-graph summarysize and whole-process sampled RSS are
different measures. Both dense/compact/reference paths coexist. Component
totals are NOT parallel PhaseB wall or solve timings; no memory/solver win is
inferred. Direct initial factors from network elimination remain unresolved.
No reuse across changing Newton iterates is assumed.

## Run policy

New immutable small overlay snapshot; read-only existing numerical/capture
roots, private first depot/read-only prior caches. Gold6438Y+64CPU300GiB,
Julia/BLAS/OMP1, occupancy guard, no-requeue; batch35min/process25min.
No overlapping expensive jobs. Pin source/inherited/input/leaf/stage hashes,
run synthetic raw-skew/interleaving/chain/terminal unit gates on node first.
Preserve every rejection/log, even if a setup-only replacement is necessary.
Assess bounded accuracy/cost outcome, commit/push raw evidence and findings;
do not automatically extend to medium or a full optimizer.
