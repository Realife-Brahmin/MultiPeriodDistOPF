# Bounded genuinely compact component test (October23)

Continues68b2b613. No production defaults, benchmark rows, paper, batteries,
constraints or dispatch model changes. Start ONLY original1020 large late,
immutable capture leaves/state/factors. No automatic medium/optimizer campaign.

## Exact tested scope

One leaf1/leaf2 Schur merge, corresponding to the first scan combination in
the existing four-period frozen window with true dense future. It returns
three factored Schur operators and two explicit vectors, not dense matrices.
HODLR matrix products and recursive Woodbury solves operate on stored factors.
The eliminated matrix and output diagonal bases follow production's ORIGINAL
coordinate upper completion. Raw input asymmetry remains stored independently;
Mab is never symmetrized. Interleaved topology clusters use masked factors to
implement original-coordinate upper completion, with possible rank growth.

Only local blocks of size<=64 may be assembled/factored densely. Larger node
Woodbury-core dimension>=node dimension explicitly rejects the compact path
and records a dense fallback, instead of silently reconstructing a full nB
matrix. Local strong interactions, full vectors and all physical variables
remain. Initial off-diagonal compression uses own-block relative tolerance1e-9;
factor-sum/upper-mask recompression uses QR and a small core SVD at1e-12.
These are two fixed controls, not a broad rank/tolerance tuning campaign.
No fixed rank or horizon-independent rank is assumed.

The algebra is recursive block-diagonal plus low-rank Woodbury updates. The
HODLR solver family is established; see the authors' implementation overview:
https://hodlrlib.readthedocs.io/en/latest/ . This is an independently derived
small Julia diagnostic, not a claimed new algorithm or copied implementation.
No general low-complexity guarantee is inferred from this system's ranks.

## Reuse and deferred work: no free lazy-output speedup

Topology can be reused for a fixed network/battery set; its setup time is
reported. Initial six raw blocks still need factor discovery (dense cross-block
SVD), and that full cost is measured. Within the merge, factors are reused and
updated through QR/small-core SVD; no dense output is constructed in the kernel.
The eliminated factorization is reused for gradients, condition probes and
Schur actions. These caches are valid ONLY while the iterate's matrices are
unchanged. No free reuse across optimization iterations is assumed.

The prototype does NOT implement closed compact-to-compact scan composition.
Schur output actions remain factored operator chains. Today's next scan merge
and stage policy code need dense blocks. Therefore a separate EXPLICIT legacy
export applies these operators to batches of16 columns to build dense matrices;
its entire cost and allocation are charged. A hybrid four-period replay then
uses this replacement for only the first merge; all remaining operations are
unchanged dense production operations. This is a component test, not a fully
compact scan. Cheap operator-object creation cannot justify medium/optimizer
extension without the input discovery, deferred-consumer and total cost checks.

Direct reuse of ORIGINAL stage tree-KKT separator factors as a generator of
temporal boundary factors is not implemented by this prototype. The existing
separator topology and discovered boundary factors are reused. Their interface
width need not stay fixed through time (SEPARATOR_ASSESSMENT.md). Report any
remaining initial discovery cost as a specific unresolved obstacle, not as
evidence that separator-factor reuse was measured or disproved.

## Gates and timing boundaries

First require saved dense full-control/equality probes and ordinary/production
value gate<1e-7. Verify all three merged blocks and two vectors against the
actual production in-place merge, all resulting hybrid suffix values, operator
actions and compact solve residual<1e-7 scaled. Recover battery paths and selected
saved full control/equality/bound directions using the unchanged1e-7 scaled AND
1e-5 relative-with1e-10-floor gates. Full matrix reconstruction is permissible
ONLY as explicit validation/legacy-export work outside the compact kernel.

If accuracy passes, fully warm paths, then three alternating dense/compact
pairs in the SAME process on identical data. Record fresh input factor discovery,
factor updates, approximate condition probing, operator actions, dense export,
legacy downstream replay and explicit cleanup GC separately, including all
allocations/within-phase GC. Dense uses matching first production merge and
the same downstream continuation. Its scratch allocation is recorded. File
I/O and accuracy checks stay outside timers; initial cold clocks are not warmed
results. Condition probing is NOT a certified norm/condition bound, nor a
substitute for direction accuracy. No expensive profiler inside timers.

Assess both recurring cached-input kernel cost and full fresh-input/legacy
pipeline cost; never omit construction/export from the overall conclusion.
Coefficient bytes, object-graph size and whole-job MaxRSS are different; raw
reference and compact objects coexist in validation, so there is no isolated
resident-RAM comparison. Component time is not parallel PhaseB wall time, and
hybrid replay time is not full-solve performance. Full optimizer quality,
feasibility/OpenDSS, iteration, near-optimality/full-convergence timing remain
untested until a general accuracy-and-total-cost pass justifies that next gate.

One job only: Gold6438Y+,64CPU300GiB,Julia/BLAS/OMP1,occupancy guard,35minute
batch/25minute process cap,no requeue,private first depot and immutable reused
source/captures. Preserve ALL failed attempts; repairs in new snapshots only.
