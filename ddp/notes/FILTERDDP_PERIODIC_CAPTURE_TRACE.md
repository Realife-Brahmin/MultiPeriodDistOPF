# Periodic backward-pass array capture (ieee123C_1ph, T=3)

Diagnostic instrumentation for inspecting how every matrix/vector FilterDDP
carries per stage evolves over outer DDP iterations, and specifically whether
the multi-right-hand-side (multi-RHS) KKT solve — the dominant per-stage cost
— has exploitable structure. This is exploratory, not part of the paper's
optimization stack, and not a timing benchmark: instrumentation overhead is
not measured or reported anywhere else.

## Setup

`ddp/patches/periodic_capture_instrumentation.patch` adds two env-var-gated
hooks to `backward_pass!` (apply after the full patch stack, same as
`factor_bound_sensitivities.patch`):

- `FILTERDDP_PERIODIC_CAPTURE_DIR`: if set, every `FILTERDDP_PERIODIC_
  CAPTURE_STRIDE`-th outer iteration (default stride 5, gated on `data.k %
  stride == 0`), every stage serializes one `.jls` file containing: `K` (the
  sparse KKT matrix), `rhs_multi` (the true pre-solve multi-RHS block, saved
  *before* `ldiv!` overwrites it in place — this bit matters, see below),
  `kkt_solution` (the same block after solving), `alpha`/`beta`/`psi`/`omega`,
  the factored bound-dual sensitivity scale vectors `Sigma_L`/`Sigma_U`, and
  the incoming value-function curvature/gradient `Vxx_incoming`/`Vx_incoming`.
- `FILTERDDP_PERIODIC_CAPTURE_STRIDE`: the stride (default `5`).

Run reproduced with:

```bash
FILTERDDP_PERIODIC_CAPTURE_DIR="ddp/results/network_filterddp/periodic_capture_ieee123C_1ph_T3" \
FILTERDDP_PERIODIC_CAPTURE_STRIDE=5 \
julia --startup-file=no --project=envs/ddp2026 \
  ddp/examples/power_system/ieee123c_filterddp.jl ieee123C_1ph 3 solve quiet
```

Solved cleanly in 48 iterations (identical to the unmodified stack: same
objective and residuals to all printed digits), confirming the instrumentation
does not perturb the algorithm. 30 captures (10 sampled iterations x 3
stages), analyzed with `ddp/examples/power_system/analyze_periodic_capture.jl`,
which writes `periodic_capture_summary.csv` alongside the raw captures.

**The raw `.jls` captures (54 MiB for this one T=3 run) are not committed** —
they are exactly reproducible from the patch and script above, and committing
periodic captures at larger horizons/systems would scale to hundreds of MiB
or more. Only the summary CSV, the patch, and this note are tracked. Ask if
you want a specific run's raw captures committed anyway.

## A capture bug caught and fixed mid-investigation

The first version of the instrumentation placed the whole capture block after
`kkt_solution = rhs` — but that line is `ldiv!(F, rhs)` executed *in place*,
so by the time the field `rhs_multi=copy(rhs)` ran, `rhs` no longer held the
right-hand side; it held the solution. This produced a nonsensical result
(an apparent multi-RHS "rank 1-2" collapse driven entirely by the feedforward
column's enormous relative magnitude once mixed with the mislabeled block,
while beta/omega showed near-full rank — impossible for a genuine RHS, since
a nonsingular solve preserves rank exactly, so beta/omega's rank can only be
less than or equal to the true input RHS's rank, never higher). Fixed by
snapshotting `rhs` immediately after it is assembled, before the solve.

## Findings

Shapes are stage-invariant for this network model (`nx, nu, nc` are backward-
pass type parameters, fixed across the horizon): `K` is `1353x1353` sparse
(`nu+nc = 791+562`), density ~0.30-0.44%; the multi-RHS state-sensitivity
block is `1353x51` (`nx=51` columns); `beta` is `791x51`; `omega` is `562x51`;
`Sigma_L`/`Sigma_U` are length-791 vectors; `Vxx_incoming` is `51x51`;
`Vx_incoming` is length-51.

**The terminal stage (t=3, t=N) is exactly, structurally different, not just
numerically.** Its multi-RHS state-sensitivity block has *every* singular
value equal to exactly `1.0`, at every sampled iteration. This is not
approximate compressibility — it is because at `t=N` the incoming value
curvature `Vxx` is zero on the very first pass through this stage each outer
iteration (there is no stage beyond the horizon), so the `B_active` term the
RHS would otherwise carry vanishes, leaving the RHS built entirely from `-cx`
— and the only state-dependence in the last stage's constraints is the
per-battery energy-slack equality `x_b - dt*u_pb - emin - slack_b = 0`, whose
Jacobian w.r.t. `x` is (up to zero rows) exactly `-I`. **This terminal-stage
KKT solve should be exploiting the RHS's exact identity structure directly
(row selection, not a dense multi-RHS solve at all) rather than treating it
as 51 generic dense columns** — this is a free, exact simplification, not an
approximation.

**Stages 1-2 (the ones actually carrying a nonzero incoming value message) are
not meaningfully rank-deficient early on**, and mildly compressible later:

- At the cold start (`iter=0`) and through early iterations (`iter<=15`), the
  state-sensitivity RHS block needs essentially its full 51-dimensional basis
  even at a loose 10% relative-Frobenius tolerance (`rhs_rank_10pct` = 48-51).
  No exploitable structure here.
- From roughly `iter=20` onward — as the barrier parameter shrinks and bound
  constraints become active/inactive — the RHS starts to compress: by
  `iter=35-45`, `rhs_rank_10pct` drops to 10, and `rhs_rank_1pct` to ~44-49
  (condition number climbing from ~4x at the start to 150-280x by the end).
  This tracks convergence, not a fixed property of the problem: a compression
  scheme would need to adapt its rank budget per iteration, not use one fixed
  rank for the whole solve.
- `omega` (the equality-multiplier sensitivity) compresses noticeably more,
  and earlier, than `beta` (the control sensitivity): by `iter=25`,
  `omega_rank_5pct` is down to 2-4 (out of 51) at stage 1, while
  `beta_rank_5pct` is still 46-48. `beta` stays close to full rank throughout
  this run. If a compressed/reduced-column KKT solve is worth pursuing, it
  looks like a much better candidate for `omega` than for `beta` — a uniform
  low-rank treatment of the whole `[beta; omega]` stack would be limited by
  `beta`'s near-full rank.

## Caveats

- Single system (`ieee123C_1ph`), single horizon (`T=3`), one stride (5).
  Whether the mid-to-late-iteration RHS compression and the beta/omega
  asymmetry hold at larger `T`, or at `ieee2522C_1ph`/`large10kC_1ph` where
  the multi-RHS solve is the documented bottleneck, is untested.
- "Rank" here is the smallest truncation whose relative Frobenius-norm error
  is within tolerance — a data-analysis proxy for compressibility, not a
  demonstrated algorithmic speedup. Actually exploiting either finding (exact
  terminal-stage structure, or iteration-adaptive omega compression) would
  require its own implementation and correctness check against the existing
  reproducibility evidence in `ddp/README_FILTERDDP_EXPERIMENT.md`.
