# Diagonalising the full-space stage Hessian

Asked by R. Gupta for the 2026-09-18 meeting: *"Hessian -- can we not just use a
diagonalized matrix? Does that converge? Time benefit?"*

**Answer: yes it converges, to the same optimum, and it is a real win that grows
with system size, confirmed up to large10k (1020 batteries, `nu = 54665`).** On
large10k it cuts per-stage factorisation by 90% and total wall by **46-57%**
(two independent runs, see below), paying 16.5% more iterations for it -- in line with the 16-20% iteration penalty
already seen at ieee123 and ieee2522, so the penalty does not worsen with scale.

This is the **full-space** question. Earlier cheap-curvature results in this repo
(`battery_only`, Nystrom low-rank, floored Nystrom) are all **reduced-space** and
all stalled; they do not transfer, and conflating the two would give the wrong
answer here. See `REDUCED_SPACE_INNER_OPF_FEASIBILITY.md` for those.

## What is replaced

The backward pass factorises, per stage per iteration,

    K = [ H   cu' ]        H = luu + diag(Sigma_L + Sigma_U) + fu' Vxx fu + fuu + cuu + reg*I
        [ cu   0  ]

(`ddp/DDP4OPF.jl/src/backward_pass.jl`, H assembled at 178-209, K at 231).
`FILTERDDP_DIAG_HESSIAN=1` replaces H by `diag(H)`, leaving the constraint
Jacobian `cu` untouched. `FILTERDDP_DIAG_HESSIAN_FLOOR` (default `1e-8`) floors
each diagonal entry from below, which is **required**, not cosmetic: at the
terminal stage `nnz(H) = 536` against `nu = 791`, so line flows, voltages and
currents have no stored diagonal entry at all and `diag(H)` is singular without
it.

## The structural prediction, and why it was wrong

On ieee123 T=3 roughly 83% of `nnz(H)` is the dense `nB x nB` block `fu' Vxx fu`
(`IEEE123_FILTERDDP_IPOPT_MATRIX_COMPARISON.md:34-47`). The network part of H is
already nearly diagonal. So diagonalising H does essentially one thing: it
**deletes `Vxx`** -- the intertemporal curvature that makes this a second-order
method. The expectation was that this would cost convergence outright, as it does
in the reduced space and as removing `Vxx` does in the user's own first-order DDP
(period-two limit cycle, `stage12_user_ddp.log`).

It does not. The iterate count rises by 20% (ieee123) and 16% (ieee2522) and the
objective is unchanged to 9-10 significant figures. The reading is that in an
interior-point method the **diagonal barrier term `Sigma = z/s` dominates the
stage curvature**, so discarding the off-diagonal part perturbs the Newton
direction much less than its share of `nnz` suggests.

## Measured, `T = 3` periodic exports, `C_B = 1e-3`, `FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8`

Timing and nnz diagnostics on for both arms, so wall times are directly
comparable to each other but not to uninstrumented runs.

| | ieee123 exact | ieee123 diag | ieee2522 exact | ieee2522 diag |
|---|---|---|---|---|
| wall | 16.970 s | 15.312 s (-9.8%) | 112.052 s | **100.315 s (-10.5%)** |
| iterations | 50 | 60 (+20%) | 56 | 65 (+16%) |
| objective | 2688.984189773 | 2688.984194219 | 8238.431077492 | 8238.431089973 |
| per-stage factorisation | 13.05 ms | **7.71 ms (-41%)** | 157.18 ms | **57.57 ms (-63%)** |
| total factorisation | 2.39 s (43.7%) | 1.64 s (37.1%) | 31.59 s (34.3%) | 13.13 s (18.7%) |
| total RHS solve | 0.48 s (8.8%) | 0.36 s (8.1%) | 43.02 s (46.6%) | 30.35 s (43.2%) |
| mean `nnz(K)` | 7226 | 5269 (-27%) | 142339 | 96025 (-33%) |
| mean `nnz(LU)` | 18545 | 10981 (-41%) | 369777 | 213259 (-42%) |
| fill ratio | 2.567 | 2.084 | 2.598 | 2.221 |

Two things scale the right way: the per-stage factorisation saving grows with
size (41% -> 63%), and the RHS solve gets cheaper too (-29% on ieee2522), which
the structural argument did not predict -- fewer `LU` nonzeros means cheaper
triangular solves against all `nx+1` columns, not just a cheaper factorisation.

## The floor matters, and not monotonically

ieee123 T=3, same instance:

| floor | iterations | wall | outcome |
|---|---|---|---|
| `1e-8` | 60 | 14.34 s | converged, objective matches to 9 digits |
| `1e-4` | 58 | 14.33 s | converged, objective off by 1.9e-04 relative |
| `1e-1` | 200 | 22.54 s | **failed**, `du_inf` stalls at 3.3e-07 |
| `1e1`  | 200 | 17.85 s | **failed**, `du_inf` 5.2e-07, alpha collapses to 4.9e-04 |

A floor large enough to act as damping destroys convergence rather than
stabilising it. Keep it at the smallest value that removes the structural zeros.

## large10k, measured 2026-09-18

`large10kC_1ph`, `T = 3`, periodic profile, `C_B = 1e-3`, same
`FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8` as the other two systems:

| | large10k exact | large10k diag |
|---|---|---|
| status | converged (0) | converged (0) |
| iterations | 103 | 120 (**+16.5%**) |
| objective | 3124025.510478886 | 3124025.509532730 |
| obj rel diff | -- | **3.03e-10** |
| max equality residual | 2.083e-09 | 9.391e-10 |
| per-stage factorisation | 1669.76 ms | 164.53 ms (-90.1%) |
| wall | 3434.995 s | 1467.627 s (-57.3%) |

**Converges cleanly, to the same optimum, at the largest tested scale.** The
iteration penalty (+16.5%) lands squarely inside the range already seen at
ieee123 (+20%) and ieee2522 (+16%) -- it does not get worse with size, which
was an open question. The earlier extrapolation in this note guessed roughly
35% off wall from the ieee2522 arm; the measured reduction is 57.3%, well
past that guess, because per-stage factorisation collapsed by 90% rather than
scaling with the ieee2522 arm's 63%. `mean_nnz_LU` fell from 1.78M to 868K
(-51.2%).

**A second, independent run** (`run_agenda_pipeline.sh full`, 2026-09-17
18:20-19:28, idle machine) reproduced the same trajectory exactly -- 103 and 120
iterations, objectives identical to every printed digit, same `nnz` -- but with
different wall times:

| large10k | exact wall | diag wall | wall saving | per-stage factor |
|---|---|---|---|---|
| run above (2026-09-18) | 3434.995 s | 1467.627 s | -57.3% | 1669.76 -> 164.53 ms (-90.1%) |
| pipeline (2026-09-17) | 2623.629 s | 1409.841 s | -46.3% | 2168.62 -> 208.34 ms (-90.4%) |

The diag arm is stable across runs (1410 vs 1468 s, 4%); the exact arm is not
(2624 vs 3435 s, 31%). Since both trajectories are bit-identical, that spread is
machine timing noise concentrated in the arm that spends most of its time
factorising. **Quote the saving as 46-57%, not as a single number.** The
per-stage factorisation reduction (-90%) is robust across both.

**Why the benefit scales so steeply:** the deleted block `fu' Vxx fu` is dense
`nB x nB`, so it grows as the SQUARE of the battery count while the network part
of `K` grows linearly. At large10k it is ~65% of `nnz(K)` by itself (1.11M ->
0.39M mean); on ieee123 diagonalising removes 27%. After diagonalisation the
factorisation is no longer the bottleneck at all -- the `nx+1`-column RHS solve
is 54% of the diag arm's wall (764 s of 1410 s in the pipeline run).

## What this does not yet establish

- Only `T = 3`. The iteration penalty comes from discarding `Vxx`, whose
  influence should grow with the horizon, so the 16-20% figure is a lower bound
  on what longer horizons would pay.
- Only `C_B = 1e-3`. `C_B` enters `luu` as a perfectly-conditioned diagonal, so
  it lands entirely inside the part that survives diagonalisation -- a smaller
  `C_B` should make the approximation worse, in the same way it breaks every
  cheap reduced-space curvature model.
- No block-diagonal variant was tried (keep the `nB x nB` `Vxx` block, diagonalise
  the rest). Since the network part of H is already nearly diagonal that variant
  would be close to exact and save almost nothing -- but it would isolate how much
  of the 16-20% iteration penalty is specifically `Vxx`.
