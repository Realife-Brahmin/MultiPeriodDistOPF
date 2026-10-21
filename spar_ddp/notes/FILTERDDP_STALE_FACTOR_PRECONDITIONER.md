# Can a stale KKT factorisation be reused? No -- and now for a reason, not just empirically

Asked by R. Gupta for the 2026-09-18 meeting: *"Biggest bottleneck is solving
AX=B for every time-step at every iteration. Any scope of being able to solve for
X = Ainv*B? Assuming that Ainv doesn't change much?"*

**The premise fails. `A` does not change a little between interior-point
iterations -- it changes by 100% to 240% in Frobenius norm after a SINGLE
iteration, and by a factor of 1e9 over twenty.** Every scheme that assumes a
slowly-varying `A` is therefore ruled out at once, including the ones that were
still open.

## Why this needed a second experiment

`FILTERDDP_FREEZE_KKT=N` already showed that reusing `lu(K_j)` *as the solve*
diverges (commit `6a32f3ed`: N=1 converges in 50 iterations, N=2 fails at 300
with `du_inf` 3.6e+02). But that test throws the true `K` away, so it could not
distinguish "the factor is stale" from "the approximation is unusable". The sound
version keeps the true `K` in the residual and uses the stale factor only to
accelerate:

    x <- x + M^{-1} (b - K_k x),     M = lu(K_j),  j < k

This converges to the **exact** solution whenever `rho(I - M^{-1} K_k) < 1`, and
costs one stale triangular solve per step instead of a fresh factorisation. It is
the version of the user's question that could actually have worked.

**Break-even.** `k` refinement steps beat a fresh factorisation when
`k * solve < factor + solve`, i.e. `k < 1 + F/S`. Measured per-stage `F/S`:
ieee123 T=3 `13.06/4.82 = 2.71`, large10k T=6 `12.88/4.58 = 2.81`. So refinement
must converge in **3 steps or fewer** to pay for itself.

Refinement, not GMRES, is the relevant scheme: the backward pass solves against
`nx+1` right-hand sides simultaneously (52 columns on ieee123, 1021 on large10k),
and refinement handles the whole block in one application of `M^{-1}` where
block-GMRES would not. Per-column GMRES is reported anyway, as the best case any
Krylov method could reach with this preconditioner.

## Measured: ieee123 T=3 periodic, `C_B = 1e-3`, stage 1, all 51 iterations captured

Script `ddp/examples/power_system/stale_factor_preconditioner.jl`; capture via
`FILTERDDP_PERIODIC_CAPTURE_DIR` with `STRIDE=1`.

Short lags, every anchor (lag = iterations between the factorisation and the solve):

| anchor | lag | `‖ΔK‖/‖K‖` | `‖ΔΣ‖/‖Σ‖` | `rho(I-M⁻¹K)` | refine steps | GMRES steps |
|---|---|---|---|---|---|---|
| 0 | 1 | 0.166 | 1.245 | 0.773 | >30 | 17 |
| 10 | 1 | 1.336 | 1.394 | 0.757 | >30 | 17 |
| 20 | 1 | 2.363 | 2.363 | 0.159 | 12 | 6 |
| 30 | 1 | 2.293 | 2.293 | 0.273 | 12 | 5 |
| 40 | 1 | 1.132 | 1.132 | **1.426** | diverges | 6 |
| 20 | 2 | 6.722 | 6.722 | **1.315** | diverges | 14 |
| 30 | 2 | 10.760 | 10.760 | 0.457 | 19 | 7 |
| 30 | 3 | 26.659 | 26.659 | **1.982** | diverges | 14 |

Long lags, anchor 20: `‖ΔK‖/‖K‖` reaches 4.6e+04 at lag 10 and 1.3e+09 at lag 20,
with `rho` at 36 and 1486. There is no lag at which the stale factor is close.

**The best result anywhere in the sweep is 12 refinement steps at lag 1** -- four
times past the 3-step break-even. Even GMRES, which cannot be used on the
multi-column block, needs 5-6 steps at its best and so would not pay either.

## The mechanism

`‖ΔΣ‖/‖Σ‖` equals `‖ΔK‖/‖K‖` to four significant figures at anchors 20, 30 and
40. The interior-point barrier terms `Sigma = z/s` are not a perturbation of `K`;
they **are** `K`'s norm, and they grow without bound as iterates approach their
bounds. This is exactly why the fast-decoupled-power-flow analogy breaks: FDPF
can freeze `B'` because susceptance is constant by construction, whereas `Sigma`
is a function of the iterate.

It also explains the independent large10k T=6 observation that per-sweep
factorisation time rises 38x (3.23 s -> 126.01 s) **while `mu` is still pinned at
1.0** for roughly the first 100 of 139 sweeps
(`iteration_timing_large10kC_1ph_T6.csv`). What moves `K` is the slacks
approaching bounds, not the barrier parameter descending.

## What is now closed

- Freezing the factorisation (already known: diverges).
- Stale factor as a **preconditioner** for iterative refinement -- 12 steps at
  best against a 3-step break-even.
- Stale factor as a **Krylov preconditioner** -- 5 steps at best, and unusable on
  the multi-RHS block regardless.
- **Refactor on `mu` change.** `mu` is constant for ~100 of 139 sweeps at
  large10k T=6 while `K` is changing fastest, so this trigger would freeze the
  factor precisely when it must not be frozen.
- Forming `A^{-1}` explicitly. At ieee123 it costs ~1353 solves (~0.044 s)
  against a 0.0049 s factorisation, and the dense `A^{-1}B` GEMM is slower than
  the sparse triangular solve; at large10k the inverse is ~75 GB.

Previously closed elsewhere, same question: RHS blocking (all-at-once 3.127 s is
fastest, widths 1-512 span 3.192-3.384 s, `FILTERDDP_KKT_RHS_BLOCK_BENCHMARK.md`);
parallel UMFPACK over RHS columns (7.59% at best,
`FILTERDDP_PARALLEL_UMFPACK_RHS.md`); reusing the sensitivity map `beta` instead
of the factor (`beta` needs rank 51 of 51 for 1% error,
`IPOPT_SENSITIVITY_REUSE_AND_MAP_COMPRESSION.md`); symbolic-factorisation caching
(slower: 295.3 s vs 233.1 s on ieee2522 T=3, `../results/network_filterddp/README.md:177`).

## What is still open

- A **slack-triggered** refactorisation (refactor when `‖Sigma‖` moves more than
  x%) is not tested, but the numbers above make it unpromising: `‖ΔΣ‖/‖Σ‖` exceeds
  1 within a single iteration for most anchors, so the trigger would fire every
  iteration and reduce to the default.
- Only ieee123 stage 1 is measured. The mechanism is instance-independent
  (`Sigma` behaves this way in any interior-point method) but the specific step
  counts are not.
- The route this leaves open is a **cheaper factorisation**, not a reused one:
  `FILTERDDP_DIAGONAL_HESSIAN.md` cuts per-stage factorisation 63% on ieee2522 by
  changing what is factorised rather than how often.
