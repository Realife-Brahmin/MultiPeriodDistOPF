# Parallel in time, and what profiling it turned up (2026-09-30)

Agenda item for the meeting of 2026-10-07: Sundermann et al. and tADMM both
decouple the time steps and treat the coupling separately. Can FilterDDP?
Branch `ddp-parallel-in-time-sep30`. All runs are in the Table II
configuration (diagonal Hessian, blocked solve, factor-backed policy, exact
assembly rewrites, per-system `C_B`, 8 Julia threads except where stated).

## 1. What in a stage waits for stage t+1

Per-category time of the Table II run, large10k `T=24` (3986 s):

| part | time | waits for `t+1`? |
|---|---|---|
| derivatives | 307 s | no |
| KKT assembly | 458 s | no, except the `n_B` battery-power diagonal entries from `V_xx` |
| factorization | 528 s | only through those `n_B` entries |
| multi-column solve (8 threads) | 855 s | yes: right-hand side carries `V_x`, `V_xx` |
| algebra + value update | 587 s | yes |
| forward pass | 1052 s | sequential rollout |

Two tests follow. **A (exact):** evaluate everything that does not wait for
`t+1` for all stages in parallel before the sweep. **B (approximate, the
Sundermann analog):** lag the only coupling entries by one sweep, so that all
`T` factorizations are independent too.

## 2. Profiling test A found a plain inefficiency instead

Timing each callback of one large10k stage (`nu = 54665`, `nc = 42303`):
the residual callback `c(x, u)` takes **104 ms**; every other callback
takes under 2 ms. It pushed every residual into a `Vector{Any}` through Dict
lookups, and FilterDDP calls it at every stage of every backward sweep and
every line-search rollout -- roughly 900 of the 3986 s.

`FILTERDDP_TYPED_EQUATIONS=1` (`typed_equations` in `ieee123c_filterddp.jl`)
resolves indices once and repeats the same floating-point operations in the
same order. `check_typed_equations.jl`: **0 bitwise mismatches** at every
stage on all three systems; callback 142x (ieee123), 202x (med2522), 510x
(large10k, 99 -> 0.19 ms) faster. Full runs follow the identical iteration
trace:

| case | iterations | before | typed | change |
|---|---|---|---|---|
| ieee123 `T=6` (1 thread) | 67 = 67 | 23.2 s | 19.0 s | -18% |
| med2522 `T=24` | 77 = 77 | 374.4 s | 284.1 s | -24% |
| large10k `T=24` | 98 = 98 | 3986.1 s | clean rerun PENDING | |

The first large10k run (r1) followed the identical trace in 3478.5 s, but
overlapped the test-B runs from 12:48, so its time is not reported; the
clean repeat is r2.

med2522: derivatives 45.8 -> 9.9 s, forward pass 41.4 -> 11.1 s.

**Consequence for test A:** its budget was mostly this callback. With it at
0.2 ms, what remains parallelizable exactly is a few ms per stage of
callbacks and barrier terms, so test A is not built.

## 3. Test B: lagged battery curvature

`FILTERDDP_LAGGED_CURVATURE=1` (`_lagged_prefactor` in `backward_pass.jl`):
the battery-power diagonal `diag(f_u' V_xx f_u)` of each stage matrix comes
from the previous sweep; all `T` stage matrices are assembled and factorized
in parallel before the sweep (BLAS pinned to one thread meanwhile).
Right-hand sides, gains and the value recursion stay exact; only the matrix
lags. The first sweep of a solve has nothing to lag and runs sequentially.

| case | with the lag | without |
|---|---|---|
| ieee123 `T=6` | near-optimal at iteration 72 | 67 |
| med2522 `T=24` | **fails**: line search, iteration 55, `mu = 3.2e-4` | near-optimal at 77 |
| large10k `T=24` | **fails**: line search, iteration 83, `mu = 8e-3`, objective 0.12% above Ipopt, residual `6.7e-3` | near-optimal at 98 |

Both failures are the filter line search shrinking the step to `1e-16`: the
direction computed with the lagged matrix is not acceptable to the filter.
The lag is small in norm (that curvature term moves 0.3-1.6% per iteration,
KKT_EVOLUTION.md), but the gains, the value recursion and the forward-pass
policy then come from slightly different quadratic models, and the filter
does not tolerate it. **Negative result: test B is closed.** Timings of these
runs are not reported: by a queuing error they ran concurrently with the
large10k baseline. Even ideally the gain was bounded: factorization is 29% of
a med2522 backward sweep, and the parallel prefactor took 0.63 s per sweep
against 1.02 s sequential, under that contention.

Why Sundermann et al. and tADMM can decouple time and FilterDDP cannot in
this way: Sundermann et al. only precondition with the decoupled matrix and
recover the exact Newton step with GMRES (affordable at one right-hand side
per iteration); tADMM never forms a Newton step across time at all, and pays
in consensus iterations. Test B does neither: it takes an inexact step and
hands it to a line search built for exact ones.

Not tested: lagging `V_xx` *consistently* at stage `t` (in the gains and the
value recursion too, not only the matrix). Every stage's 1021-column solve
would then be independent and could run in parallel, with only the one
`alpha` column sequential. The curvature would then propagate one stage per
sweep, which is the staleness that stalled the user's first-order DDP
(CLAUDE.md), so it is not expected to help.

## 4. The lead this exposed: the value update needs only the battery block of `K^{-1}`

With the factor-backed policy the stage backward pass solves
`K [α β; ψ ω] = -[Q̂u B̃; c c_x]` for `1 + n_x` columns (1021 at large10k),
but then uses:

- `α`, `ψ` in full (feedforward, dual updates);
- `β` only in the `n_B` battery-power rows: `V_xx = C + β_B' B_active + ω' c_x`;
- `ω` only in the `n_B` energy-constraint rows, the only nonzero rows of `c_x`;
- `β' Q̂u + ω' c` in `V_x`, which by symmetry of `K` equals `B̃' α + c_x' ψ`
  and so needs only `α`, `ψ`;
- `β δx` in the forward pass, recomputed per rollout from the stored factor.

So of the `(n_u + n_c) x n_x` solution (96,968 x 1,020 at large10k) the
backward pass consumes `2 n_B x n_x` (2,040 x 1,020), about 2%. The
right-hand side is nonzero in those same rows, so what is needed is the
`2 n_B x 2 n_B` block of `K^{-1}` on `E` = {battery powers, energy rows}. That
block is the inverse of the Schur complement of `K` onto `E`: the network KKT
Kron-reduced onto the battery buses, i.e. how each battery's marginal price
responds to every other battery's injection.

In DDP terms: `Q_ux = f_u' V_xx f_x` is nonzero only in the battery-power
rows, so `V_xx = Q_xx - Q_ux' Q_uu^{-1} Q_ux` (constrained) needs only the
battery x battery block of the constrained `Q_uu^{-1}`. This is where the
MPOPF-specific derivation should land.

How it could be computed, untested: a factorization with `E` ordered last
leaves the Schur complement as its final dense block (MUMPS `ICNTL(19)`,
PARDISO; UMFPACK does not expose it), at roughly the cost of a factorization
plus a dense 2,040 x 2,040 block, against today's factorization plus a
1,021-column solve.
