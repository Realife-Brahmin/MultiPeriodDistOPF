# Phase B representation: initial assessment, 2026-10-21

Exploratory only. No optimizer, queued job, frozen benchmark source, solver
default, or paper was modified. Battery-count PR173 remains the benchmark
workstream; this note is a separate follow-up, not a new benchmark result.

## Demonstrated implementation

`examples/power_system/pint_scan.jl` defines `PintElem` with dense Float64
`Maa`, `Mab`, `Mbb` (each nB by nB), and `ma`, `mb` (each nB).
They represent the quadratic in entering and exiting battery states. Constant
cost is omitted by this derivative representation. Maa and Mbb are symmetric;
Mab is a cross-boundary block and need not be symmetric.

`paired_pass.jl::_pp_combine!` forms K = left.Mbb + right.Maa, factors it,
performs two triangular matrix solves, two symmetric rank updates and one
dense matrix product for a full two-sided combination. A closing combination
with PintVal does less work. Cholesky failure uses the existing LU fallback.
This is O(nB^3) dense work per full combination, not just message passing.
Persistent PairWork has another three dense matrices and one vector per worker.

`procs_pass.jl` maps three matrices and two vectors per period in shared
memory on one node. Reading shared data is not network communication.
The balanced scheduler can assign a combination away from its stage owner;
round completion is a barrier. Explicit-fetch experiments additionally incur
serialization, transfer and source-worker availability. Neither rounds_s nor
summed task/fetch clocks identify communication-only time.

Payload per period is 8(3 nB^2 + 2 nB) bytes: at nB=1020, 24,985,920 bytes;
T96 requires 2,398,648,320 bytes for these shared slots alone. This excludes
network models, factors, private scratch and all other process memory.
Storing only one triangle of Maa/Mbb would reduce payload to
8(2 nB^2 + 3 nB) = 16,670,880 bytes, about 33.3% less, exactly. That alone
does not remove cubic work; packing/unpacking or packed BLAS may negate gains.

The scan computes all remaining-future descriptions, not a single reduction:
the current schedule makes 368 combinations at T96 and 1426 at T288.

## Existing original-1020 capture inspected locally

Read-only copies of stage_1 through stage_4 from
`/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/capture` were inspected using
Julia 1.12.7 on ETRL309, BLAS1. All four SHA256 hashes match the saved capture
manifest. No new cluster computation/capture was launched. The original
capture input's Bset and parent map match the unchanged original-battery T3
metadata used in the present campaign. The capture is historical first-pass
data, not evidence about the final campaign's later iterations.

Battery rows partition into 102 topology groups of ten. A diagnostic split
places 51 whole groups on either side (510 rows each). Production
`pint_combine` was used to construct intervals 0:2 and 0:4 from the four
captured leaves; this was untimed analysis, not a solver performance test.
Generated data and reproducible inspection script are in
`experiments/value_representation/`.

Findings in this limited sample:

- All three blocks are 1020 by 1020. Leaf 0:1 has 350570 exact nonzeros per
  block, leaf 1:2 has 109180; compositions have 1040400 (fully populated).
  Exact nonzero counts are not meaningful numerical-rank or coupling tests.
- Maa/Mbb relative symmetry error is below 3e-18 in these leaves and zero
  after the production symmetric updates. Mab's symmetry error reaches
  1.1e-7 in 0:4: do not symmetrize Mab.
- The selected cross-group block has Frobenius norm only about 1.6e-21 to
  1.5e-20 of its parent matrix. Its leaf spectrum is essentially rank one.
  At relative-to-that-small-block threshold 1e-8, apparent ranks become
  4/2/4 for Maa/Mab/Mbb at 0:2 and 8/7/8 at 0:4.
- The sampled ten-battery within-group block remains rank ten at every
  tested threshold (1e-12 through 1e-4 relative to its largest singular value).

These exceptionally small off-group terms make relative spectral ranks
sensitive to rounding/cancellation. They are NOT evidence of physically
significant rank growth or a universal compression tolerance. Nor is the
global matrix itself low rank: strong local/full-rank terms must be retained.
Only one partition, four leaves, scan depths zero through two, and the first
optimization pass were inspected. No across-iteration claim is supported.
Dropping tiny Hessian entries without assessing scale/conditioning is unsafe.

## Exact candidates and their assumptions

1. Symmetric-storage reduction above is unconditional for Maa/Mbb, not Mab.
   It preserves all physical variables, but needs a compatible efficient kernel.
2. Block-separable problems: if all three blocks share an exactly block-diagonal
   partition and K is invertible blockwise, composition stays block diagonal.
   This is a closure statement, not a claim that our model is exactly separable.
3. Small network separators: `tree_kkt.jl::_tree_kkt_prepare` and
   `_tree_kkt_block0` already form local group blocks plus a root-mediated
   correction before materializing dense Sh0. This motivates preserving that
   representation at leaf construction. The older structured solver comments
   explicitly assume diagonal Hessians in one path; they are not proof of
   exact closure for today's exact-curvature scan.

The general principle is a graph separator, not a hardcoded count of 102
feeders. Larger separators, meshed networks and global constraints can increase
or destroy the compactness. Grouping is a matrix partition; all individual
battery states, locations, dynamics and limits remain.

For block diagonal plus low rank, addition concatenates bases, inverses can
use Woodbury if the requisite factors are nonsingular, and Schur products
introduce additional bases. Closure of the representation class does NOT
mean closure at a fixed small rank: rank can grow with interval length and
eventually saturate at nB. Retaining separator variables for every period
also grows the interface with interval length. A fixed rank-three temporal
representation is not established by a small per-period root separator.

Storage for one symmetric structured block is O(sum b_g^2 + nB*r + r^2)
versus O(nB^2); a generic nonsymmetric low-rank correction needs both left
and right bases. Factor/solve cost is roughly O(sum b_g^3 + nB*r^2 + r^3),
plus applying bases and any multiple-right-hand-side work. These are candidate
counts, not measured gains; rank growth and conversion costs must be included.

## Approximation candidates and safeguards

Preserve dense local blocks and compress only off-group/separator interactions,
possibly hierarchically using automatically detected graph separators. Compare
spectral and Frobenius tail errors in scaled coordinates, absolute magnitude
relative to the full matrix, and action on representative vectors. Select ranks
from an error budget, not a system-specific battery count or arbitrary cutoff.

Schur inversion can amplify error: delta(K^-1) is approximately
-K^-1 deltaK K^-1. Small local compression error need not imply small value
gradient, Newton direction, objective or feasibility errors. Maintain symmetry
of the full quadratic and check conditioning/inertia; do not assume a truncated
cross block preserves positive definiteness. A changed representation is an
inexact Newton step, even though physical constraints and states remain.
Full-solve feasibility, line-search/filter behavior, iterations and OpenDSS
must subsequently be checked, with dense fallback and an explicit error policy.

## Smallest next controlled experiment (proposed, not submitted)

First use existing first-pass captures only: inspect four to eight adjacent
periods and the same dense compositions at successive depths. Compare an
exact block-plus-separator implementation with the dense result, including
all three matrices and both vectors. Establish algebraic closure and rank/
storage growth before proposing a production backend. Then test off-group
truncation separately, reporting errors against dense, never calling it exact.

For performance isolation, use identical frozen operands and worker placement:
(a) one-worker local kernel, (b) balanced shared-memory execution, (c) preloaded
private operands with matched task/barrier schedule. Include matched no-op
dispatch/barrier controls. Time K construction, factor, triangular solves,
rank updates/product, copies, allocation/GC, dispatch and barrier wall time.
Memory/cache effects cannot be obtained by simply subtracting overlapping
task times; reset policy and instrumentation overhead need controls. Report
both per-pass and per-round distributions, not just accumulated Phase B.

Only if the historical first-pass evidence survives those checks, obtain
separate untimed early/mid/late captures from the current original-1020
configuration. Do not inject profiling into the running battery-count campaign.
Choose optimization passes spanning barrier changes, not merely consecutive
iterations. Record accepted iterations and total backward-pass count alongside
per-pass costs, because retries/barrier passes can change accumulated time.

Current priority: exact separator-aware leaf representation and explicit
rank-growth diagnosis before whole-matrix SVD compression or physical battery
aggregation. No new simulations, capture jobs, production implementation or
paper changes were made for this assessment.
