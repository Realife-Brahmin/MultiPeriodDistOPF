# October 23: bounded block-scaled / conditioning-aware follow-up

Authorized follow-up to 25dc3e46, not a production or paper change. No dispatch
aggregation. Reuse immutable medium late, original1020 large late, medium early
control. Medium replays all24 suffixes; large uses first4 periods with the true
dense future. No new captures, optimizer trajectories or benchmark-row edits.

Raw Maa/Mab/Mbb are represented independently. Hierarchy local blocks and both
vectors remain exact, including raw floating-point skew. Each directed hierarchy
cross block gets a relative Frobenius budget based on ITS OWN norm, rather than
the large joint-Hessian norm. Adaptive SVD retains coupling; dense blocks are
kept whenever factor storage is not smaller. This intentionally stores upper
and lower factors separately instead of assuming bitwise raw symmetry.

Trials: block-scaled1e-6 control, then conditioning-aware base1e-6,1e-9,1e-12,
stopping at the first accurate compressed conditioning-aware candidate. Actual
upper-completed elimination K condition is computed explicitly in Inf norm;
effective tolerance is base/sqrt(max(1,kappa)), or zero/raw on nonfinite kappa.
Initial leaves use Maa+Mbb as a conditioning heuristic, NOT a proof about future
eliminations. Condition-control costs are charged and reported separately.
Failed reconstruction uses raw element fallback at the unchanged bound.

First reproduce saved dense full-control/equality probes<1e-7, and ordinary
versus production values<1e-7. Candidate gates require maximum recovered
direction error norm(delta)/max(1,norm(ref))<1e-7 AND
norm(delta)/max(1e-10,norm(ref))<1e-5 across battery paths and selected full
controls/equality/bound-multiplier probes. This stricter relative floor is
explicit, not an original physical residual. Every rejected trial is logged.
Direction checks against dense are an OFFLINE selection oracle; they are not
free or available inside an eventual compressed optimizer. No feasibility or
convergence claim follows from a frozen-direction pass.

Only after a candidate passes: warm chosen paths, then alternate dense/candidate
order across3 same-input repetitions IN ONE PROCESS per point. Construction,
materialization, condition control, combination, allocation/GC are distinct;
whole-call time additionally includes deepcopy, workspaces and in-memory record
overhead. File I/O, direction checks and explicit boundary GC are outside warm
timers. Warmup and exploratory accuracy trials are NOT warmed results. No claim
that3repeats eliminate all system noise. Dense baseline uses production_values.

Prototype still holds payload AND dense expansion, closed values stay dense;
coefficient ratios are not resident RAM reductions, component times are not
parallel PhaseB wall times. This is not a compact computation kernel. Measure
current implementation fairly, not infer the cost of a future compact kernel.
Only accuracy AND total-cost evidence can justify a bounded optimizer comparison;
otherwise report the precise limitation without declaring compression itself
unsuccessful or expanding the campaign. Any optimization then requires original
feasibility/objective/OpenDSS plus near-optimal and full times separately.

One sequential job: Gold6438Y+,64CPU300GiB,occupancy guard, Julia/BLAS/OMP1,
75minute batch cap /22minute independent point caps. Private first depot, old
numerical snapshot/caches read-only. Check quota first; no full source copies
or deletion of old evidence. New immutable overlays only for setup repairs.

The initial implementation used recursively nested NamedTuple payload types,
which specialize on tree shape and dense/rank choices. A fixed recursive node
implementation preserves operation order, arrays, budgets/ranks and coefficient
counts while bounding compiler specialization. followup_equivalence.jl compares
both implementations (synthetic all four budgets, real medium-late leaf at zero
and1e-6), requiring byte-identical arrays, ranks, errors and coefficient counts.
The old code is read from the immutable f8f60d39 snapshot on Kamiak. This is an
implementation repair, not a scientific tolerance or algorithm change.
