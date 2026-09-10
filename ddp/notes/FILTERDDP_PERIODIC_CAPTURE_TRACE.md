# Periodic backward-pass array capture (T=3, all three systems)

Diagnostic instrumentation for inspecting how every matrix/vector FilterDDP
carries per stage evolves over outer DDP iterations, and specifically whether
the multi-right-hand-side (multi-RHS) KKT solve — the dominant per-stage cost
— has exploitable structure. This is exploratory, not part of the paper's
optimization stack, and not a timing benchmark: instrumentation overhead is
not measured or reported anywhere else.

## Setup

`ddp/patches/periodic_capture_instrumentation.patch` adds env-var-gated hooks
to `backward_pass!` (apply after the full patch stack, same as
`factor_bound_sensitivities.patch`):

- `FILTERDDP_PERIODIC_CAPTURE_DIR`: if set, every `FILTERDDP_PERIODIC_
  CAPTURE_STRIDE`-th outer iteration (default stride 5, gated on `data.k %
  stride == 0`), every stage captures. Two modes:
  - **Raw dump (default)**: one `.jls` file per (iteration, stage) containing
    `K` (the sparse KKT matrix), `rhs_multi` (the true pre-solve multi-RHS
    block, saved *before* `ldiv!` overwrites it in place — this bit matters,
    see below), `kkt_solution` (the same block after solving),
    `alpha`/`beta`/`psi`/`omega`, the factored bound-dual sensitivity scale
    vectors `Sigma_L`/`Sigma_U`, and the incoming value-function
    curvature/gradient `Vxx_incoming`/`Vx_incoming`. Used for `ieee123C_1ph`
    and `ieee2522C_1ph`.
  - **Inline stats** (`FILTERDDP_PERIODIC_CAPTURE_INLINE_STATS=1`): computes
    the same shape/SVD-rank statistics `analyze_periodic_capture.jl` would
    compute offline, right where the arrays already live in memory, and
    appends one small CSV row directly to `periodic_capture_summary.csv` —
    never serializing the arrays themselves. **Required at `large10kC_1ph`
    scale**, see incident below. Verified byte-identical to the raw-dump
    pipeline's output on `ieee123C_1ph` before trusting it at scale (same
    run, same iteration/stage, same 30 rows, only differing in row order).
- `FILTERDDP_PERIODIC_CAPTURE_STRIDE`: the stride (default `5`).

Reproduce e.g. ieee123 with:

```bash
FILTERDDP_PERIODIC_CAPTURE_DIR="ddp/results/network_filterddp/periodic_capture_ieee123C_1ph_T3" \
FILTERDDP_PERIODIC_CAPTURE_STRIDE=5 \
julia --startup-file=no --project=envs/ddp2026 \
  ddp/examples/power_system/ieee123c_filterddp.jl ieee123C_1ph 3 solve quiet
```

or large10k with `FILTERDDP_PERIODIC_CAPTURE_INLINE_STATS=1` added.

All three T=3 runs solved cleanly, matching each system's known iteration
count from the paper's own table (ieee123: 48, ieee2522: 56, large10k: 115)
and objective (all three agree with the centralized reference to 1e-3 to
1e-4 relative), confirming the instrumentation does not perturb the
algorithm.

**Raw `.jls` captures are not committed** — they are exactly reproducible
from the patch and script above. Only the summary CSVs, the patch, and this
note are tracked.

## Disk-space incident (large10k, raw-dump mode)

The raw-dump mode's per-snapshot cost is `O(nu*nx + nc*nx)` dense doubles
(`beta`, `omega`, `rhs_multi`, `kkt_solution` are all dense at this size).
That's KB-MB at ieee123/ieee2522 scale but ~2.4 GiB *per stage per sampled
iteration* at large10k scale (`nu~54665, nc~42303, nx=1020`) — not
anticipated before running it. The first large10k attempt (3 stages x ~24
samples) filled the disk to 99% (5.8 GiB free out of 476 GiB) before the
solve finished, truncating its last capture file (caught as an `EOFError`
on read). All raw captures from that attempt, and the ieee123/ieee2522 raw
captures (54 MiB and 4.9 GiB respectively — fine individually, but all
already fully analyzed and no longer needed), were deleted; only the CSV
summaries survive. The inline-stats mode above was added specifically to
avoid repeating this, and was cross-checked against the raw-dump mode's own
output on ieee123 before being trusted for the large10k rerun.

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

Shapes are stage-invariant within each system (`nx, nu, nc` are
backward-pass type parameters, fixed across the horizon):

| System | nx | nu | nc | K shape | outer iterations (T=3) |
|---|---|---|---|---|---|
| ieee123C_1ph | 51 | 791 | 562 | 1353x1353, ~0.3-0.4% dense | 48 |
| ieee2522C_1ph | 250 | 13358 | 10337 | 23695x23695, ~0.02% dense | 56 |
| large10kC_1ph | 1020 | 54665 | 42303 | 96968x96968, ~0.015% dense | 115 |

**The terminal stage (t=N) is exactly, structurally different, not just
numerically — confirmed on all three systems, including large10k.** Its
multi-RHS state-sensitivity block has *every* singular value equal to
exactly `1.0`, at every sampled iteration, at all three scales. This is not
approximate compressibility — it is because at `t=N` the incoming value
curvature `Vxx` is zero on the very first pass through this stage each outer
iteration (there is no stage beyond the horizon), so the `B_active` term the
RHS would otherwise carry vanishes, leaving the RHS built entirely from `-cx`
— and the only state-dependence in the last stage's constraints is the
per-battery energy-slack equality `x_b - dt*u_pb - emin - slack_b = 0`,
whose Jacobian w.r.t. `x` is (up to zero rows) exactly `-I`. **This
terminal-stage KKT solve should be exploiting the RHS's exact identity
structure directly (row selection, not a dense multi-RHS solve at all)
rather than treating it as `nx` generic dense columns** — this is a free,
exact simplification, not an approximation, and it holds at every scale
tested so far, independent of everything below.

**Stages 1-2 (the ones actually carrying a nonzero incoming value message):
the interior-stage picture does *not* scale monotonically with network
size.** ieee2522 looked like the start of a "compresses more at larger
scale" trend; large10k breaks it.

- **RHS-block compression is pronounced on ieee2522 and nearly absent on
  ieee123 and large10k.** ieee123: `rhs_rank_10pct` drops from ~50 (out of
  51) to 10 by `iter=35-45`. ieee2522: the sharpest case, dropping to 7-9
  (out of 250, under 4%) by `iter=30-55`, condition number climbing to
  1200-2400x. large10k: **the RHS barely compresses at all** —
  `rhs_rank_5pct` stays at 920-1018 (out of 1020, 90-99%) for the *entire*
  115-iteration run, and `rhs_rank_1pct` never drops below 1010. Whatever
  makes ieee2522's RHS collapse late in the solve is not a generic
  large-network effect; it is not present at 4x ieee2522's state dimension.
- **`omega` compresses more than `beta` on all three systems, but the *depth*
  of that gap is not monotonic in scale either.** ieee123: `omega_rank_5pct`
  bottoms out around 2-4 (out of 51, near-total) then partially recovers to
  14-23; `beta_rank_5pct` stays 28-50. ieee2522: the most dramatic case —
  `omega_rank_5pct` down to 13-17 (out of 250, ~6%) against `beta_rank_5pct`
  of 48-150 (20-60%). large10k: `omega_rank_5pct` only reaches 165-365 (out
  of 1020, 16-36%) — a real gap below `beta`, but nowhere near ieee2522's
  depth. Across all three, `omega` never compresses less than `beta` at the
  same iteration; the ordering is robust even though the ieee2522 magnitude
  does not generalize.
- **large10k shows a late-iteration recovery the shorter runs never reach.**
  `beta_rank_5pct` falls monotonically from 982 (96%) at `iter=0` to a
  minimum of 562 (55%) at `iter=80`, then *rises back* to 936 (92%) by
  `iter=115` — a U-shape, not a one-way decay. This coincides with the
  barrier parameter `mu` collapsing from `1e0` toward its floor (`1e-8` by
  the final iterations, per the CSV's `barrier_mu` column): the late-solve
  regime the shorter ieee123/ieee2522 runs (48 and 56 iterations) never
  reach at all. Any compression scheme would need to know it is in this
  terminal-convergence phase and back off, not just use a fixed rank budget
  or a monotonically shrinking one.

## Caveats

- All three systems covered at `T=3` only; whether any of the interior-stage
  patterns above hold at larger `T` is untested. Given they already don't
  generalize cleanly across *systems* at fixed `T=3`, extrapolating to
  larger `T` is speculative, not a natural next assumption.
- The one pattern that *is* robust across all three systems and worth
  building on is the terminal-stage exact-identity structure — everything
  else in the interior-stage findings should be read as system-specific
  data points, not a scaling law.
- "Rank" here is the smallest truncation whose relative Frobenius-norm error
  is within tolerance — a data-analysis proxy for compressibility, not a
  demonstrated algorithmic speedup. Actually exploiting either finding (exact
  terminal-stage structure, or iteration-adaptive omega compression) would
  require its own implementation and correctness check against the existing
  reproducibility evidence in `ddp/README_FILTERDDP_EXPERIMENT.md`.
