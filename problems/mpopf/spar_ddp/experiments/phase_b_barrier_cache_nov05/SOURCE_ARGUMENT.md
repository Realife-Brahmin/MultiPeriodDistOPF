# Source-level argument and conservative eligibility

Existing optimization: SPARDDP.jl/src/reuse_pass.jl _reuse_candidate and
_backward_pass_reuse!. ReusePass is tied to solver identity/k/Hessian model
and surviving factor policies; backward_pass.jl populates it only on supported
successful unregularized battery-block passes. It re-solves new mu-dependent
RHSs through saved sequential stage factors. Those are NOT scan merge factors.
procs_pass.jl explicitly omits this path and its pass-GC clears policies.

solve.jl barrier branch checks inner tolerance, reduces mu/increments j,
resets filter and CONTINUES without forward_pass!/updating nominal x/u/phi/z.
A hybrid switch can accompany this branch and MUST invalidate. Accepted
steps advance k and nominal primal/dual, hence invalidate even if some arrays
happen to look close. Changing mu by itself is not sufficient eligibility.

In _tp_pre: slack u-l / upper-u and Sigma_L=zl/slackL,
Sigma_U=zu/slackU determine curvature. These depend on current primal/dual/
bounds, NOT explicitly on mu. The barrier RHS mubar=mu*(inv_uu-inv_ul)
and complementarity updates DO depend on mu. Hessian derivatives depend on
x/u and constraint multipliers phi; fx/fu and Jacobians depend on model/x/u.
Supported linear dynamics have no multiplier-dependent second derivatives.
K and network factors derive from those fixed matrices. _tp_pre contains no
reg input or regularization retry; network tree factor/finish0/finishV use
unshifted K, Cholesky then LU fallback, not a hidden diagonal shift. Any fixed
Hessian floor is part of its model/path and numerical K and must be unchanged.
The generic backward_pass regularization loop is a DIFFERENT execution path.
Pinned capture _tp_pre/procs hook success therefore establishes the applied
unshifted path from source, not from d.reg_last (which procs resets afterward).
Unknown/different regularization or a fallback into another driver invalidates.

_pp_elem and _pp_elem!: leaf Mbb is diag(1/f)*Sh0*diag(1/f) in battery order;
Mab uses Rhat and that block; Maa also adds Rhat transpose and D. Rhat/D,
Sh0, f and static elimination factors depend only on the matrices above.
Leaf ma/mb use rhat,dvec, current S.b/network RHS/gradient and thus change
with mu. Terminal soft-SOC gamma is fixed, its curvature is in the FINAL
leaf's objective Hessian; its gradient depends on terminal nominal state.
Actual scan terminal25 is zero V/v, not a zero artificial window boundary.
Model/gamma/terminal coefficients/state changes invalidate as appropriate.

Induction over EXACT same scan DAG: unchanged leaf curvature and terminalV
imply unchanged H=x.Mbb+y.Maa (or y.V), BL=x.Mab,BR=y.Mab and each merged
curvature. Cholesky uses Symmetric(H,:U) in ORIGINAL coordinates, preserving
raw input arrays rather than averaging asymmetry. U' solves produce Y1/Y2;
curvature uses SYRK/GEMM. Given exact invariance, retaining U,Y1,Y2/output
blocks avoids those same matrix operations. Current h=x.mb+y.ma(or y.v)
still requires U'\h and both gradient matvec updates. No gradient, h or
feedforward/bound update is retained. Every node is included, not sampled.

Eligibility signatures include x/u/phi/zl/zu, bounds/masks, derivative and
KKT coefficient arrays, factor/layout/elimination data, dynamics, terminalV,
solver/model/source/backend/execution/Hessian/regularization context; k and
j/mu chronology are separately checked. Exact bit-sensitive signatures (no
near-equality tolerance) conservatively invalidate the WHOLE cache when any
relevant item differs. A changed factor with identical final leaf curvature
still invalidates. Unsupported applied-reg proof disables reuse. Signatures
are detection work and fully charged, not free trusted version counters.

This argument has supported-config assumptions. Global captured comparisons
and negative controls test its implementation; they do not generalize it to
arbitrary mu-dependent objectives, alternate dynamics or changing models.

## Minimal V2 storage and safe miss behavior

V1 exact prototype retained every intermediate output and leaf matrix so a
local miss could run fresh at any node. Fully charged cost was negative and
its3.396GB payload was unnecessary for a complete hit pass; it is preserved.
V2 uses an all-or-nothing cache, not approximate reuse. Each of68 entries has
U,Y1,and either Y2 (44two-boundary merges) OR its FINALclosedV (24value merges).
Thus204dense1020x1020 matrices,1,697,932,800 numeric bytes, not408matrices.
No gradients, leaf matrix cache, matrix inputs or internal Maa/Mab/Mbb outputs
are retained. This is the minimal dense-BLAS layout for these selected kernels,
not a universal packed-storage lower bound: unused lower U could be packed
only with a different storage/solve interface and additional conversion cost.

On a verified complete hit, leaf and intermediate curvature need not be
materialized: vector recurrences only read current ma/mb, U,Y1,Y2. Every
closed value restores its cached V, so ALL24suffix gradients/Hessians and
current consumers are fully available. All68 vector solves/112gradient
matvecs still execute, followed by unchanged current-stage reconstruction.
Before reuse require68entries and exact global eligibility. A missing entry
makes the ENTIREpass fresh, preserving all original curvature operands. An
unexpected mid-round loss also abandons speculative vector work and restarts
from current leaves with unchanged fresh kernels; never fresh-combine an
unmaterialized intermediate. Failed first-pass factorizations leave missing
entries and disable future whole-pass reuse. Both population and hit outputs
are explicitly checked against full fresh replay before measured pairs.
