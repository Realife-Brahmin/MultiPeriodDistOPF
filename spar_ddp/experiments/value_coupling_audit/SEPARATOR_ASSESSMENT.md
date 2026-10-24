# Exact separator possibility versus implemented approximation

This is an algebraic assessment, not a claim that an exact compact temporal
kernel has been implemented. It applies to general spatial partitions, not
just the original large10k root feeders. All battery variables are retained.

## What would make an exact representation possible?

Suppose a partition of the FULL time-expanded quadratic/KKT system has left
variables L, right variables R and separator variables S. If, conditional on
S, there is no L--R coupling, eliminating S gives the cross term

    H_LR(reduced) = -H_LS * inverse(H_SS) * H_SR.

This factorization is exact when the chosen elimination is well defined, and
its rank is at most the number of separator variables (also bounded by the
two block dimensions). For an indefinite KKT system the same block-elimination
identity applies, but invertibility/pivoting must actually be checked. It does
not grant a positive-definiteness assumption. If direct L--R edges remain,
their contribution must also be retained; the displayed factor alone is then
not the complete interaction.

The existing tree network solver provides a useful starting point: it validates
the radial stage KKT pattern and stores local own/interface indices, factors
and their coupling positions. A general hierarchy can preserve these spatial
separator factors instead of discovering all coupling again with dense SVD.
However, a battery-weighted graph partition is not by itself a proof that a
particular three-entry parent interface separates the entire temporal system.
Separator batteries, dynamics, equality multipliers and any other links must
be included in the time-expanded cut. Our replay records the actual stage
interfaces; it does NOT yet construct and verify that full temporal cut.

## Why the single-period width is not a fixed scan rank

Across an interval of L periods, separator responses from different periods
can be independent. A valid exact interface can therefore grow with L rather
than stay at the one-period size. One generic conservative bound is the size
of the complete time-expanded separator, capped by the boundary block sizes;
no horizon-independent bound has been established here. Composition requires
updates to separator factors, local solves and their Schur complements, not
just forwarding an unchanged fixed-rank object. Fill/rank growth and stable
indefinite elimination are essential parts of any exact-kernel design.

Observed check, from completed medium T24 frozen points (not a benchmark):
at the TOP hierarchy split, the early-point numerical rank using sigma >
1e-12*sigma_max grows across scan depths0/1/2/3/4 as3/6/12/20/26 (maxima over
the intervals at each depth). Barrier-before and after give3/6/12/20/23.
This is descriptive floating-point rank, not proof of exact algebraic rank,
but directly contradicts treating the observed single-period rank3 as fixed.
The raw spectra and interval/depth indices are retained by the replay.

## What the implementation does today

Zero-truncation mode reproduces the full quadratic to rounding, falling back
to dense cross blocks when exact numerical factors would cost more. That is
an exactness CONTROL, not a compact exact algorithm: roundoff-scale singular
values frequently prevent useful zero-truncation rank reduction. Even its
storage ratio includes a symmetry-storage convention and is not evidence of
network-separator compression alone.

Adaptive modes approximate hierarchy OFF-DIAGONAL blocks of the JOINT two-sided
Hessian, preserve full local blocks and both vectors, and recompress after
element composition. They retain cross-group interactions rather than delete
them. Budget1e-8 refers to each element's joint Frobenius reconstruction, not
to every individual block, final suffix, direction, or optimization residual.
These downstream errors and KKT inertia are independently measured. The
current prototype materializes dense matrices and still uses dense production
combination, so neither an exact compact kernel nor a compressed fast solver
has been delivered. No production or paper result is changed.

## Next smallest computational discriminator

After all frozen accuracy checks, first isolate warmed construction and
combination on the SAME saved medium early/late inputs, including allocation,
GC and materialization. Current cold serial clocks mix compilation with SVD/
object construction and cannot establish the cost of a warmed compact kernel.
If storage and accuracy justify continuing, implement a separator/factor-aware
or hierarchical low-rank combination with dense fallback, then compare its
warmed total construction+composition cost against dense. Only a candidate
with justified error/inertia AND cost checks should enter a convergence test.
An adaptive per-element budget may need tightening near convergence; it must
be judged using original physical feasibility and the total iteration count,
not the Hessian reconstruction tolerance alone.
