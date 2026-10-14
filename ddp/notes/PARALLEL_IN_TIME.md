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
| large10k `T=24` | 98 = 98 | 3986.1 s | 2675.0 s | -33% |

All nine Table II cells, clean evening session (2026-09-30,
`run_table2_typed.sh`), time to near-optimality, same iterations and
objectives as Table II: ieee123 `T=6/24/96` 18.0 / 26.8 / 74.6 s (-13 / -16 /
-22%), med2522 85.1 / 293.3 / 1355.9 s (-17 / -21 / -23%), large10k `T=6/24/48`
842.9 / 2672.6 / 4491.6 s (-31 / -33 / -29%). These replaced Tables II and
III in the TPEC paper. Runs repeated after the user left were 4-6% slower
than the same runs while they were at the machine (identical traces),
probably the Balanced power plan; all table values come from the one
session.

The first large10k run (r1) followed the identical trace in 3478.5 s, but
overlapped the test-B runs from 12:48, so its time is not reported; the
clean repeat is r2.

med2522: derivatives 45.8 -> 9.9 s, forward pass 41.4 -> 11.1 s. large10k
(r2, background load 0.33 cores): derivatives 307 -> 48 s, forward pass
1052 -> 111 s, every other category within 4%; near-optimality at 2672.6 s
against Table II's 3983.5 s.

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

### Tested (2026-09-30): exact, and break-even with MUMPS off the shelf

`battery_block_check.jl` on 20 captured stage systems: `E` has exactly
`2 n_x` rows (102 of 1,353; 498 of 23,691; 2,040 of 96,968), the Schur-
complement rows match the full solve to <= 3e-11, and the `V_x` identity
holds to <= 5e-10.

`FILTERDDP_BATTERY_SCHUR=1` computes the value increments from MUMPS's Schur
complement (`battery_schur_hook.jl`; UMFPACK has no Schur interface). It skips
the 1,021-column solve at every stage except the terminal one, where the soft
terminal-SOC penalty makes `l_ux` nonzero. UMFPACK still factors each stage
for the feedforward column and the forward-pass policy. Full runs reach
near-optimality **at the same iteration with the same objective** as the full
solve: ieee123 `T=6` (67), med2522 `T=24` (77, every digit), large10k `T=6`
(101, every digit, also with the analysis reused).

The stage KKT pattern is constant (393,075 nonzeros in all 2,520 large10k
stage evaluations), so MUMPS's analysis runs once and later calls only
refactor. Per large10k stage, informal (machine in use):

| per stage | full solve (Table II) | Schur, fresh analysis | Schur, analysis reused |
|---|---|---|---|
| 1021-col solve, or Schur + dense solve | 0.37 s | 0.98 s | 0.60 s |
| right-hand-side assembly | 0.22 s | 0.07 s | 0.07 s |
| value update | 0.15 s | 0.10 s | 0.11 s |

Clean evening run, large10k `T=24` (background load 0.30 cores): near-optimal
at iteration 98 with the same objective as the full solve, in 2,994.5 s
against 2,672.6 s, i.e. **12% slower**. Per stage: Schur + dense solve
0.74 s against the blocked solve's 0.33 s; assembly 0.17 -> 0.02 s, value
update 0.09 -> 0.04 s. At `T=6` (informal) it was about break-even.

So off-the-shelf MUMPS makes the battery-block route break-even at best.

**Status (2026-09-30, closed with the user):** exact but not faster. UMFPACK
cannot help: its permutations are known, but it pivots the battery rows into
the middle of its elimination order, so reading off only those rows still
needs essentially full triangular sweeps; only an ordering that puts them
last (MUMPS's Schur option) makes them cheap, at the price of MUMPS's slower
factorization and a dense 2,040 x 2,040 block. A single MUMPS factorization
serving the step and the forward pass too would remove the second
factorization but is still estimated at parity (0.33 + 0.23 s against
0.20 + 0.33 s per large10k stage). The hand-written leaves-to-substation
(Kron) sweep is parked as a discussion item for the 2026-10-07 meeting. Of its
0.60 s, 0.19 s is the dense 2,040 x 2,040 LU, and the separate UMFPACK
factorization (0.2 s) is still paid. Next: one factorization for both (MUMPS
reduced right-hand sides, `ICNTL(26)`), or a tree elimination that
Kron-reduces the radial stage network onto the battery rows directly.

## 5. The radial tree solver (2026-09-30/10-01): the battery block done cheaply

The battery-block study (Section 4) was closed as "exact but not faster" with
general solvers. Exploiting the feeder tree by hand changes that.

**Structure.** Group the stage KKT by bus: a non-root bus owns its incoming
line's `P, Q, ell`, SOC slack, its voltage, its DER control, and its
P-balance, Q-balance, voltage-drop and SOC rows (at most 10 unknowns). Checked
against the sparsity pattern on all three systems, each group couples only to
its parent, through 3 unknowns (the parent's two balance rows and voltage).
Eliminating groups from the leaves to the substation is block Gaussian
elimination on the tree: a small dense LU per bus and a 3x3 update to the
parent, no fill.

**The battery block is feeder blocks plus a low-rank term.** The substation is
the only place feeders meet, through 3 root unknowns. With the root last,
`S = D - U inv(M) U'`: `D` block diagonal by feeder (large10k: 102 blocks of
about 20), `M` the 6x6 root block, `U` of rank at most 6. `S \ R` is then
Woodbury with one small solve per feeder; no 2,040 x 2,040 matrix is formed.
The same factorization solves full systems `K \ b` with two tree sweeps, so
it replaces the stage's sparse LU for the feedforward column, the battery
rows and the forward-pass policy (`tree_kkt.jl`, `FILTERDDP_TREE_KKT=1`).

**Exact.** On production-config captures: `S` against MUMPS 1e-16 (ieee123),
3e-11 (med2522), 3e-13 (large10k); battery rows against a full UMFPACK solve
5e-16, 3e-12, 8e-16; one-column solves 6e-14, 7e-12, 4e-15, with residuals no
larger than UMFPACK's.

**Two exact rewrites it exposed** (`FILTERDDP_STRUCTURED_DYNAMICS=1`): `f_x` is
the identity and `f_u` is `-dt` times a selector, so the four dense `n_x^3`
products per stage with them are replaced by sparse ones; and the soft
terminal SOC penalty puts `l_ux` only in the battery-power rows, so the
terminal stage is structured too (it had been building a dense `n_u x n_x`
block, an UMFPACK factorization and the full 1021-column solve: 44% of a
`T=6` backward pass once everything else was fast).

**Full runs, large10k `T=6`** (same iterate path: near-optimal at iteration
101, identical objective):

| configuration | time | per stage (backward) |
|---|---|---|
| Table II on 2026-09-30 morning | 1222.9 s | |
| typed residuals | 842.9 s | 1.19 s |
| + tree solver | 679.3 s | |
| + sparse dynamics products | 491.8 s | 0.66 s |
| + structured terminal stage | **290.5 s** | 0.36 s |

That is 10.5x Ipopt with MA57 (27.6 s), against 44x at the start of the day.
Per stage now: assembly 0.008, tree build 0.071, solve 0.166, derivatives
0.028, algebra 0.067, update 0.016 s; forward pass 54 s of the 290.

**Limits.** med2522 is a single feeder, so `D` is one 498 x 498 block and the
tree solver is no faster than UMFPACK there (110.0 against 85.1 s at `T=6`);
it needs the same low-rank idea applied at branch points inside a feeder.
With every solve on the tree solver the stopping iteration can move by one
where a criterion is borderline (ieee123 `T=6`: 68 against 67, primal
residual 9.3e-7 against a 1e-6 threshold); the solves agree to rounding.

### Clean overnight runs, all nine Table II cells (2026-10-01)

Time to near-optimality; same iteration counts as Table II except ieee123
`T=6` with the tree solver (68 against 67). `sd` = structured dynamics with
UMFPACK; `tree` = structured dynamics + tree solver, after two fixes found by
profiling med2522: the downward sweep computes only the rows a bus's
children and battery read (it had been allocating a 10 x 252 matrix per
bus), and the tree solver's dense steps run on one BLAS thread (a 498 x 498
LU takes 22 ms on OpenBLAS's ten threads against 3 ms on one).

| case | typed | sd | tree | tree vs typed | x Ipopt MUMPS / HSL |
|---|---|---|---|---|---|
| ieee123 `T=6` | 18.0 | **17.9** | 23.2 | +29% | 55 / 95 (sd) |
| ieee123 `T=24` | 26.8 | **25.7** | 35.7 | +33% | 21 / 35 (sd) |
| ieee123 `T=96` | 74.6 | **70.4** | 107.9 | +45% | 3.1 / 18 (sd) |
| med2522 `T=6` | 85.1 | 76.9 | **57.5** | -32% | 8.4 / 12.9 |
| med2522 `T=24` | 293.3 | 267.1 | **177.8** | -39% | 4.3 / 7.3 |
| med2522 `T=96` | 1355.9 | 1263.3 | **811.7** | -40% | 3.7 / 6.4 |
| large10k `T=6` | 842.9 | 645.0 | **259.7** | -69% | 5.0 / 9.4 |
| large10k `T=24` | 2672.6 | 2363.4 | **924.2** | -65% | 3.4 / 7.0 |
| large10k `T=48` | 4491.6 | 4061.8 | **1537.0** | -66% | 2.5 / 5.7 |

So the tree solver is not only a many-feeder effect: med2522 is one feeder
with 249 batteries and gains 32-40%. At that size the feeder's dense block
(498 x 498) is cheap; what had made the first version slower there was
overhead. Only ieee123 (128 buses), where UMFPACK needs 4 ms per stage,
loses. Per stage on captures: med2522 0.028 s against ~0.10 s for UMFPACK's
factorization and blocked solve, large10k 0.076 s against ~0.53 s.

### Robustness checks (2026-10-01, machine in casual use: timings indicative)

Run to FilterDDP's strict tolerance instead of stopping at near-optimality,
structured dynamics with UMFPACK (`sd`) against the tree solver:

| case | sd | tree |
|---|---|---|
| med2522 `T=6` | 81 iterations, 91.5 s | 81 iterations, 66.9 s, identical objective and residuals |
| med2522 `T=24` | 93 iterations, 327.2 s | 94 iterations, 214.3 s, objective equal to 1e-8 |
| ieee123 `T=24` (8 threads) | 93 iterations, 41.8 s | 93 iterations, 38.6 s, identical objective |

So the tree solver stays accurate through the whole barrier sequence.

`C_B = 0` with the tree solver fails exactly as the UMFPACK runs did:
large10k `T=6` line-search failure at iteration 2, med2522 `T=96` stalls at
primal `1.2e-4` (iteration 93). Those failures belong to the diagonal
Hessian, not to the linear solver.

### Long horizons, one thread, repeats and memory (clean queue, 2026-10-01/02)

`run_tree_long_horizon.sh`, tree solver + structured dynamics, quiet machine.

**Repeats** (med2522, time to near-optimality): 58.3 / 178.3 / 813.9 s against
57.5 / 177.8 / 811.7 s in the first run: reproducible to 1.4%.

**One Julia thread:** large10k `T=6` 304.8 s (259.7 on eight), med2522 `T=24`
215.2 s (177.8): 17-21% slower, same iterations.

**Longer horizons**, against the HSL runs on the same instances (peak
resident memory in GiB):

| case | DDP iters | DDP time | DDP peak | Ipopt HSL time | Ipopt peak | ratio |
|---|---|---|---|---|---|---|
| med2522 `T=96` | 96 | 813.9 s | 2.5 | 127.3 s (MA57) | | 6.4 |
| med2522 `T=192` | 108 | 1904.7 s | 4.0 | 276.7 s (MA57) | 4.6 | 6.9 |
| med2522 `T=384` | 114 | 4744.7 s | 7.2 | 589.8 s (MA57) | 9.0 | 8.0 |
| med2522 `T=1152` | 130 | 33364.4 s | 19.7 | 2045.5 s (MA57) | 19.9 | 16.3 |
| large10k `T=96` | 75 | 2732.9 s | 7.3 | 926.5 s (MA97) | 9.4 | 2.9 |
| large10k `T=192` | 120 | 9165.6 s | 13.9 | 1432.9 s (MA97) | 19.0 | 6.4 |

Against MA57, large10k is 1.46x (`T=96`) and 1.34x (`T=192`).

**med2522 `T=1536`** (MA57 and MA97 both ran out of memory there): FilterDDP
**failed**, a line-search failure at iteration 33 (`mu = 4e-2`, steps of
`6e-5` with the regularization active from about iteration 30), at 22.5 GiB
peak. So a solution beyond the centralized memory limit is still not
demonstrated.

**What this shows.**
- No memory advantage yet. Retained storage grows linearly, about 17.5 MiB
  per med2522 stage and 74 MiB per large10k stage, so at med2522 `T=1152`
  FilterDDP needs as much memory as Ipopt (19.7 against 19.9 GiB); at
  large10k it needs about 25% less. Most of it is avoidable: every stage
  keeps a tree solver made of thousands of small matrices (object overhead
  ~2.5 kB per bus), its own copy of `K`'s values and its own work buffers.
- The gap to HSL does not keep narrowing. On med2522 it widens with horizon
  (6.4, 6.9, 8.0, 16.3): iterations grow (96 to 130) and the per-stage time
  rises from 0.068 to 0.157 s at `T=1152`, with the forward pass going from
  15% to 25% of the run, consistent with the 20 GiB working set. On large10k
  the ratio moves with the two solvers' iteration counts (DDP 75 then 120,
  Ipopt 108 then 65).

## 6. Agenda of 2026-10-07: singular values of the stage objects

`sensitivity_spectra.jl` on production-config captures (stage 2 of `T=6`,
one early and one late iteration where available). Rank needed to capture
each matrix to 10% / 1% in Frobenius norm:

| object | ieee123 (51) | med2522 (249) | large10k (1020) |
|---|---|---|---|
| `V_xx` incoming | 51 / 51 | 244-247 / 249 | 1005 / 1020 |
| `V_xx` minus its diagonal | 23-24 / 46-47 | 15-69 / 160-209 | 703 / 952 |
| `beta`, all controls | 2-50 / 47-51 | 6-74 / 190-247 | 995 / 1020 |
| `beta_B`, battery rows | 51 / 51 | 245 / 249 | 1005 / 1020 |

`V_xx` has a flat spectrum (largest / smallest singular value 1.1-7) and is
nearly diagonal: its off-diagonal part is 0.7-1.0% of its Frobenius norm on
ieee123, 1.3% (late) to 29.5% (early, `mu = 0.2`) on med2522, 11% on
large10k. The off-diagonal part is not low-rank. The full sensitivity map is
dominated by a few directions only early in the solve and only to 10%.

So nothing here is compressible by rank. The structure is diagonal dominance
of the battery block, which is why the diagonal Hessian works, and what it
misses is a small full-rank remainder: carrying the full `V_xx` in the
battery block (BATTERY_BLOCK_REDUCTION_EXPLAINED.md, Section 6) would capture
it; a diagonal-plus-low-rank correction would not.

**Wider sample (2026-10-04):** a middle stage of `T=24` at an early, middle
and late iteration (ieee123, med2522) and large10k `T=6` stage 3 late, all in
the current configuration (per-system `C_B`). Rank for 10% / 1%:

| object | ieee123 (51) | med2522 (249) | large10k late (1020) |
|---|---|---|---|
| `V_xx` incoming | 51 / 51 | 243-247 / 249 | 1010 / 1020 |
| `V_xx` minus its diagonal | 14-22 / 45-47 | 37-56 / 168-195 | 659 / 933 |
| `beta`, all controls | 2-42 / 46-51 | 11-16 / 227-231 | 386 / 929 |
| `beta_B`, battery rows | 51 / 51 | 244-247 / 249 | 1010 / 1020 |

Off-diagonal share of `V_xx`: ieee123 1.1-1.2% throughout; med2522 26.7%
early (`mu = 4e-2`), 4.1% mid, 3.4% late; large10k 1.8% late (11% early).
Same conclusion: nothing is low-rank, `V_xx` is nearly diagonal, and it is
least diagonal early in the barrier phase on med2522.

## 7. Does the diagonal Hessian degrade with horizon? (2026-10-04)

`run_diag_vs_exact_horizon.sh`: ieee123, UMFPACK, structured dynamics, run to
strict tolerance; both arms converge everywhere.

| `T` | diagonal | exact | extra | time diag / exact |
|---|---|---|---|---|
| 96 | 129 | 109 | +18% | 76 / 98 s |
| 384 | 155 | 153 | +1% | 300 / 437 s |
| 1536 | 309 | 169 | **+83%** | 2473 / 2846 s |

Not monotone, but at `T=1536` the diagonal Hessian needs nearly twice the
iterations, 50 of them at the barrier floor against 20; the exact count grows
only mildly with horizon. Objectives agree to 1e-8. So the diagonal
approximation is the weak point at long horizons, on this system. (It does
not by itself explain the med2522 `T=1536` failure, which stopped at
iteration 33, far earlier.) This is the case for routing the full `V_xx`
through the battery block of the tree solver.

## 8. Exact stage Hessian through the battery block (2026-10-05)

`f_uᵀ V_xx f_u` is dense but lies entirely in the `P_B` x `P_B` block the tree
solver keeps, so with `FILTERDDP_TREE_KKT=1` and no `FILTERDDP_DIAG_HESSIAN`
it is not assembled into `K`: the network elimination is the diagonal case's,
and the dense block is added to the battery Schur complement, which is then
factored densely (2 `n_B` x 2 `n_B`). No sparse fill. Reproduces the original
exact-Hessian path on ieee123 `T=96` (109 = 109 iterations, identical
objective).

Clean runs, time to near-optimality, tree solver in both columns:

| case | diagonal: iters, time | exact: iters, time | change |
|---|---|---|---|
| ieee123 `T=6` | 68, 23.2 s | 56, 22.7 s | -2% |
| ieee123 `T=24` | 85, 35.7 s | 85, 34.9 s | -2% |
| ieee123 `T=96` | 120, 107.9 s | 100, 86.3 s | -20% |
| med2522 `T=6` | 68, 57.5 s | 45, 46.0 s | -20% |
| med2522 `T=24` | 77, 177.8 s | 67, 154.6 s | -13% |
| med2522 `T=96` | 96, 811.7 s | 83, 687.9 s | -15% |
| large10k `T=6` | 101, 259.7 s | 97, 362.3 s | +39% |
| large10k `T=24` | 98, 924.2 s | 90, 1270.9 s | +38% |

Per stage: med2522 0.064 s against 0.067 s (the 498 x 498 block is free);
large10k 0.47 s against 0.31 s (the 2,040 x 2,040 block costs 0.17 s) for 4-8
fewer iterations. So exact curvature pays on med2522 (5.4-10.3x HSL) and not,
as implemented, on large10k. Eliminating the energy rows analytically would
halve that block's dimension (about 8x less factorization work).

### The same comparisons at near-optimality, and longer horizons (2026-10-05)

The ieee123 horizon test of Section 7 ran to strict tolerance because no
Ipopt reference existed at `T=384`, 1536. With MA57 references now run
(15.3 s and 60.3 s) and `near_opt_posthoc.py` reading the per-iteration logs
(checked against a run that stopped at near-optimality itself: iteration 77
and 180.6 s against 77 and 178.3 s), iteration and time to near-optimality:

| ieee123 `T` | diagonal (UMFPACK) | exact (UMFPACK) | exact (tree) | diagonal's extra iterations |
|---|---|---|---|---|
| 96 | 120, 72 s | 100, 91 s | | +20% |
| 384 | 142, 277 s | 140, 401 s | | +1% |
| 1536 | 247, 1970 s | 148, 2511 s | 148, 1633 s | +67% |

At `T=1536` the exact Hessian through the tree solver is the fastest of the
three (27x MA57 against 33x for the diagonal), with the old exact path's
iterations and objective. These timings are single runs on one thread.

med2522 at longer horizons, clean runs, tree solver, near-optimality:

| `T` | diagonal: iters, time, x MA57 | exact: iters, time, x MA57 | change |
|---|---|---|---|
| 96 | 96, 811.7 s, 6.4 | 83, 687.9 s, 5.4 | -15% |
| 192 | 108, 1904.7 s, 6.9 | 92, 1595.4 s, 5.8 | -16% |
| 384 | 114, 4744.7 s, 8.0 | 100, 4044.1 s, 6.9 | -15% |

The exact Hessian saves 13-16 iterations at each horizon at no extra cost per
stage, and the same peak memory (4.1 and 7.3 GiB). The ratio to MA57 still
rises with horizon (5.4, 5.8, 6.9): iterations keep growing in both arms.

## 9. Parallel in time, exactly: what in a stage waits for `t+1` (2026-10-06)

Test B (Section 3) decoupled the periods by lagging the battery curvature and
failed. This section does it without approximation, in the configuration of
`CONSTRAINT_SCREENING.md` (screened model, power-flow start, exact step,
diagonal Hessian, tree solver).

### The split

With the tree solver the value function of stage `t+1` reaches stage `t` in
two places only, both inside the block `E` the solver keeps (battery powers
and energy rows): the battery-power diagonal of the stage matrix, and the
battery-power rows of the right-hand side. So a stage's backward work falls
into three parts:

| part | what | needs |
|---|---|---|
| pre | derivatives, barrier terms, KKT assembly, network factorization and Schur preparation (`_tree_kkt_prepare`), first network solve of the feedforward column | nothing from `t+1` |
| seq | products with `V_xx`, the battery block (assembly, factorization, feedback rows), the value update | the value function of `t+1` |
| post | second network solve (back-substitution), multiplier updates, stored policy | this stage's battery solution |

The forward pass splits the same way. The state depends on the controls only
through the battery powers, and their feedback rows come out of the backward
sweep, so the states of all periods follow from a recursion of `n_B x n_x`
products; each period's full policy solve, its directions, the ratio test and
the trial evaluations are then independent.

With one worker per period a pass costs `max_t(pre) + sum_t(seq) + max_t(post)`
instead of the sum of all three.

### Checks

- `tree_kkt_phase_check.jl` overwrites every entry of `K_EE`, and the `E` and
  energy-slack entries of the right-hand side, with NaN: the network factor,
  the preparation and the first network solve come out bit-identical on
  captured stage matrices of all three systems (153, 747 and 3,060 entries
  destroyed), and the finish step does read them.
- The instrumented run (`FILTERDDP_PARSIM=1`) is the ordinary sequential run
  with timers: it reproduces the committed ieee123 `T=6` run to the last digit
  of the objective and of every residual line.
- The state recursion is run beside the actual rollout: largest difference
  `1e-14`.

### The accounting, and what it assumes

`parsim_from_log.py` replaces sums over periods by maxima for the pre and post
parts and the per-period forward work, and keeps everything else as measured.
This is the accounting of the tADMM paper: per-period work is run in sequence
and the slowest period is counted. It assumes one core per period for the
independent parts (they are single-threaded here) and the machine's threads
for the sequential sweep, and it ignores communication and memory contention.
The maxima include whatever garbage-collection pause lands in one period (on
med2522 the largest pre is about 0.2 s against a mean of 0.04 s), so they are
pessimistic. Ipopt, in every comparison below, is a single process.

### Exact rewrites found on the way

The first count showed that the sequential part was most of a stage, and that
almost none of it was the battery-block factorization (1.3 ms of 197 ms on
large10k). Three rewrites followed. The first two change no floating-point
operation (identical residual lines on every pass); the third changes the
elimination order (iterates agree to rounding: `compare_runs.py`, residuals
within `1e-8` relative, objectives identical to 16 digits).

| flag | what | large10k stage |
|---|---|---|
| `FILTERDDP_LEAN_VALUE=1` | `C = l_xx + V_xx + f_xx` without products by `f_x = I`; `B` as scaled rows of `V_xx`; the feedback right-hand side written in the solver's grouped order into one reused workspace and solved in place, feeder blocks on the Julia threads | seq 197 -> 81 ms |
| (same flag) | the diagonal stage Hessian added up as a dense vector: the code summed the sparse matrices, took `diag()` (a SparseVector), floored it entry by entry and rebuilt a sparse diagonal | pre 165 -> 105 ms |
| `FILTERDDP_HALF_BLOCK=1` | the energy rows of `E` have a diagonal block once the energy slacks are eliminated, so they are eliminated too: the battery block is `n_B` instead of `2 n_B`, and the value increment is `-X_P' Rhat - cxE' inv(S_NN) cxE` with one `n_B x n_x` right-hand side (agenda 4) | seq 81 -> 57 ms |

On med2522 the same three give pre 47 -> 40 ms and seq 15 -> 9 ms.

### Results (2026-10-06, idle machine, compilation excluded)

Per stage, in ms, with all three rewrites:

| | pre | seq | post |
|---|---|---|---|
| large10k | 108 | 57 | 16 |
| med2522 | 40 | 9 | 3.5 |
| ieee123 | 3 | 2 | 0.2 |

Seconds to near-optimality. "Sequential, 10-05" is the configuration committed
on 2026-10-05 (the table in the paper); Ipopt is the faster of MA57 and MA97
on the same screened model and start.

| case | sequential, 10-05 | sequential, now | one worker per period | Ipopt HSL |
|---|---|---|---|---|
| large10k `T=6` | 28.8 | 16.1 | 7.8 | 5.0 |
| large10k `T=24` | 128.9 | 73.3 | 26.8 | 30.2 |
| large10k `T=48` | 272.1 | 154.9 | 54.3 | 65.3 |
| med2522 `T=6` | 18.0 | 13.2 | 5.6 | 3.1 |
| med2522 `T=24` | 82.5 | 62.5 | 17.6 | 18.7 |
| med2522 `T=96` | 418.8 | 319.2 | 67.2 | 90.0 |
| ieee123 `T=6` | 1.0 | 1.4 | 1.0 | 0.13 |
| ieee123 `T=24` | 2.9 | 3.6 | 1.8 | 0.45 |
| ieee123 `T=96` | 32.4 | 37.1 | 19.0 | 2.04 |

(The 10-05 column for ieee123 used UMFPACK, which is faster than the tree
solver on a feeder that small; the other columns use the tree solver.)

Sequential FilterDDP is 43-44% faster than on 2026-10-05 on large10k and
24-27% on med2522, with the same iterates: 2.4 times Ipopt-HSL at large10k
`T=24`, 48 (4.2-4.3 before) and 3.3-3.5 times at med2522 `T=24`, 96 (4.4-4.7).
Counted with one worker per period it is ahead of Ipopt-HSL from `T=24` up on
both larger feeders (0.75-0.94 of its time) and behind at `T=6`.

Before the halved block the same count was 35.3 / 70.7 s (large10k `T=24`, 48)
and 23.4 / 95.3 s (med2522 `T=24`, 96): logs in `nine_cells_lean/`.

### The exact Hessian in the same configuration

Four cells with the exact battery-block Hessian, lean path and halved block
(logs `*_exact_half.log`). These runs stop at a gap of `1e-5` to `1e-7`,
because with exact curvature feasibility arrives last, so the comparison
column is the diagonal Hessian carried to a gap of `1e-5` (its iteration count
there from the strict runs of 2026-10-05, its time scaled from the run above).

| case | exact: iterations, sequential, one worker per period | diagonal to a gap of 1e-5 (estimate) |
|---|---|---|
| large10k `T=6` | 12, 47.7 s, 34.4 s | 17 it, about 27 s |
| large10k `T=24` | 12, 187.7 s, 126.4 s | 23 it, about 153 s |
| med2522 `T=24` | 52, 93.2 s, 30.0 s | 62 it, about 99 s |
| med2522 `T=96` | 66, 465.2 s, 121.9 s | 77 it, about 500 s |

On med2522 the two are level at equal accuracy, as before. On large10k the
diagonal Hessian is now ahead even at equal accuracy: its sequential part is
57 ms per stage against 246 ms for the exact one (dense `n_B` block: 123 ms to
assemble and factor, 62 ms for the feedback rows). The rewrites of this
section went into the diagonal path first; the exact path still assembles the
dense block at full size before reducing it. Against 2026-10-05 the exact
runs themselves are 19% (large10k) and 5-24% (med2522) faster.

### What limits the count

The sequential sweep: `T` times 57 ms (large10k) or 9 ms (med2522) per pass,
about 80% of the count at the longer horizons. What is left in it is dense
`n_B`-sized algebra on the value Hessian: on large10k mostly one
`n_B x n_B x n_x` product and the memory traffic around it.

### The sweep in parallel: what it would take (assessment, nothing built)

After the network is eliminated, what the sweep solves is a linear-quadratic
problem in the battery energies alone. Per stage, with `x = B^{t-1}`,
`p = P_B^t` and `x+ = x - dt p`:

- the network leaves a quadratic in the battery powers, `1/2 dp' N_t dp + g' dp`,
  where `N_t` is the stage matrix reduced onto the battery powers (the halved
  block without the value term): dense `n_B x n_B`;
- the energy window leaves a diagonal quadratic in the next state,
  `1/2 dx+' Sigma_t dx+` (this is what eliminating the energy rows produces);
- the dynamics are `dx+ = dx - dt dp`.

With the exact battery-block Hessian the value Hessian then obeys

    P_{t-1} = W - dt^2 W (N_t + dt^2 W)^{-1} W,      W = Sigma_t + P_t,

that is, `inv(P_{t-1}) = inv(Sigma_t + P_t) + dt^2 inv(N_t)`: a Riccati
recursion with identity dynamics. Compositions of such linear-fractional maps
are again of that form, so the `T` steps can be combined pairwise in about
`2 log2 T` rounds (the conditional-value-function scan of Sarkka and
Garcia-Fernandez for LQ control is the numerically stable way; multiplying
chain matrices is not). The diagonal Hessian does not fit: it puts only
`diag(P_t)` into the stage matrix, which is not a linear-fractional map of
`P_t`.

Rough cost, from the measured sizes. One sweep step is now about 3 dense
`n_B^3` operations (57 ms on large10k with ten BLAS threads, 9 ms on med2522);
one pairwise combination is about 20 (two solves with `n_B` right-hand sides
and several products).

| | sweep per pass | scan per pass, one core per worker | scan, ten cores per worker |
|---|---|---|---|
| med2522 `T=24` | 0.22 s | 0.07 s | |
| med2522 `T=96` | 0.86 s | 0.10 s | |
| large10k `T=24` | 1.4 s | about 4.5 s | about 0.6 s |
| large10k `T=48` | 2.7 s | about 5 s | about 0.7 s |
| large10k `T=96` | 5.5 s | about 6 s | about 0.8 s |

So it would pay clearly on med2522 (the count at `T=96` would go from 67 s to
roughly 25 s, against 90 s for Ipopt), and on large10k only if each worker has
several cores or the horizon is long: with 1,020 batteries one combination
costs as much as a third of a 24-period sweep. It also changes the algorithm
to the exact battery-block Hessian, whose iteration counts differ (fewer at
tight accuracy, more at the 0.5% stop on med2522). These are estimates from
operation counts, not measurements. Days of work; a decision for the user.

### An observation about large10k, not used

`FILTERDDP_VXX_SEPARABILITY=1` measures how much of the value-Hessian update
couples batteries on different feeders: at most `7e-10` in any entry against
828 within a feeder (Frobenius ratio `2e-7` at worst over the 78 stage solves
of a `T=6` run). The 102 feeders meet only at the substation, whose voltage is
fixed and whose import has a linear cost and, after screening, no bound; what
couples them is the `1e-8` Hessian floor on the substation variables. So
large10k is, to that accuracy, 102 independent feeder problems. A value update
by feeder would make its sequential part nearly free, but it would do nothing
for a single feeder such as med2522 and is a property of this test system, so
it was not built (user, 2026-10-06). It is worth remembering when large10k
results are read: Ipopt could equally be run one feeder at a time.

### Measurement lessons

- **`NO_COLOR` defeats the warm-up.** With it in the environment Julia wraps
  stdout in an IOContext; prints compiled during the warm-up (a plain
  IOStream) are compiled again in the first pass of the timed solve: 0.7 s in
  an ordinary run, 1.5-2.3 s with the accounting lines. Job queues started
  from PowerShell carry it. The run scripts now unset it and the driver prints
  `WARMUP_STREAM_MISMATCH` if it happens. The overnight queue of 2026-10-05/06
  carries about 0.7 s per run, which matters only at `T=6`.
- A `@printf` with more than 32 arguments is not specialised.
- A queued batch loads the solver source afresh for every run: no source edits
  while one is running.

### Files

| | |
|---|---|
| `DDP4OPF.jl/src/backward_pass.jl`, `affine_linesearch.jl` | the timers (`FILTERDDP_PARSIM`), the lean value path, the separability diagnostic |
| `tree_kkt.jl` | `_tree_kkt_prepare` / `_tree_kkt_finish`, grouped feedback rows, the halved block |
| `tree_kkt_phase_check.jl` | the NaN test |
| `run_parsim.sh`, `parsim_from_log.py`, `compare_runs.py` | one run, the count, the iterate comparison |
| `results/parallel_in_time/nine_cells_half/`, `nine_cells_lean/` | logs and summaries of the nine cells |

## 10. The sweep by pairwise combination, and a lean exact step (2026-10-10/11)

Asked by the user after reading Sarkka and Garcia-Fernandez, "Temporal
Parallelization of Dynamic Programming and Linear Quadratic Control" (IEEE TAC
68(2), 2023; arXiv 2104.03186): implement it for the three systems at
`T = 6, 24, 96`, with the exact battery-block Hessian. The assessment in
Section 9 ("what it would take") is superseded by the measurements here.

### What was built

`pint_scan.jl` (`FILTERDDP_PINT_SCAN=1`, `run_pint_scan.sh`). A diagnostic
inside an ordinary exact-Hessian solve: on the selected backward passes every
stage hands over its battery-level data, and the value recursion is redone
from that data in two ways, the ordinary sweep and the pairwise combination.
The solve itself is unchanged.

With `a` the state before a stage, `dp` its battery powers and `b = a + f .* dp`
the state after it (`f = -dt`), the stage leaves

    q(a, dp) = 1/2 dp' S0 dp - (rhat + Rhat a)' dp + 1/2 a' D a + d' a,

where `S0` is the halved battery block without the curvature term and none of
the five objects depends on the next stage. The sweep is
`V_new(a) = min_dp q(a, dp) + V(a + f .* dp)`. The pairwise scheme writes each
stage as a function of both end states, `G(a, b) = q(a, (b - a) ./ f)`, combines
neighbours by minimising over the state between them (one Cholesky of size
`n_B`, two triangular solves with `n_B` right-hand sides, three products), and
does the combinations level by level: in every block of `2^(r+1)` positions
each entry of the left half is combined with the first entry of the right
half. `ceil(log2(T+1))` levels: 3, 5, 7 for `T = 6, 24, 96`. A combination whose
right operand already reaches the end of the horizon costs as much as a sweep
step.

Because `f` is invertible here (one battery power per battery energy), the
two-ended function is an ordinary quadratic in `(a, b)`; the dual
parameterisation of the paper, needed when there are fewer controls than
states, is not.

### It is exact

57 backward passes over the nine cases (`pint_*_scan2.log`): the value
functions of the pairwise scheme agree with the solver's to `7e-10` relative
or better on 56 of them, `2.5e-7` on one late med2522 `T=96` pass (`mu` near
its floor). The same holds for the sweep redone from the stage data, so the
difference is conditioning, not the scheme. No Cholesky failure in any
combination or step.

### Step B only: the count

One worker per period, one BLAS thread each. *sweep*: the value recursion in
sequence from the stage data. *paired*: slowest combination of each level,
plus the slowest stage function and the slowest final battery solve.
Milliseconds per backward pass, medians over the scanned passes; the garbage
collector is kept out of the timed operations (it had inflated the slowest
combination by 4 to 20 times in the first batch).

| | `T` | sweep | paired | times faster | one combination / one step |
|---|---|---|---|---|---|
| ieee123 (51) | 6 | 0.37 | 0.45 | 0.8 | 0.09 / 0.05 |
| | 24 | 1.33 | 0.77 | 1.7 | |
| | 96 | 5.1 | 1.1 | 4.6 | |
| med2522 (249) | 6 | 11.3 | 13.1 | 0.8 | 3.8 / 1.8 |
| | 24 | 44.2 | 24.5 | 1.8 | |
| | 96 | 192 | 38.2 | 5.1 | |
| large10k (1,020) | 6 | 348 | 413 | 0.8 | 129 / 58 |
| | 24 | 1,364 | 705 | 2.0 | |
| | 96 | 5,615 | 1,188 | 4.7 | |

With the typical (median) combination instead of the slowest: 0.9-1.2, 2.0-2.6
and 5.9-7.1. A combination costs a little over two steps at every size, so
the ratio is set by the horizon, not by the number of batteries. The paper
reports losing its advantage at state dimension 256-512; there the sequential
baseline runs on a GPU that already parallelises each step.

Not measured: a solve that uses the scheme. It is a count on step B from a
serial run, under the High performance power plan.

### A lean exact step B (in the solver)

The sweep redone from the stage data turned out 2.5 to 4.5 times cheaper than
the sequential part of the solver under the exact Hessian (14 times on
ieee123). The solver built the battery block at `2 n_B` with the curvature
inside, reduced it and factored it, all in the sweep. But

    Shat = Shat0 + Vc,   Shat0 = blockdiag_g(D_PP - D_PN inv(D_NN) D_NP) - U_P inv(M) U_P',

and `Shat0` reads nothing from the next stage. `FILTERDDP_LEAN_EXACT=1` (with
`FILTERDDP_HALF_BLOCK=1`): `_tree_kkt_finish0` builds `Shat0` before the sweep
at its final size; the sweep adds `Vc`, factors by Cholesky (LU if that
fails) and forms the value increment from one triangular solve and one
symmetric product. The battery rows of the feedback are not needed by the
sweep; under `FILTERDDP_PARSIM` they are formed and timed with the
back-substitution of the stage. The battery curvature is formed entry by
entry (`_lean_curvature`) instead of by a sparse-dense-sparse product.

Same iterates as the plain exact path (`compare_runs.py`): ieee123 `T=6`,
med2522 `T=24`, large10k `T=6`, residuals within `1.5e-9` relative, objectives
identical to 16 digits.

| exact Hessian | step B per stage | sequential | one worker per period |
|---|---|---|---|
| large10k `T=6`, plain | 165 ms | 34.2 s | 23.7 s |
| large10k `T=6`, lean | 58 ms | 24.4 s | 12.3 s |
| med2522 `T=24`, plain (2026-10-06, Balanced plan) | 13.6 ms | 93.2 s | 30.0 s |
| med2522 `T=24`, lean | 3.7 ms | 64.4 s | 14.3 s |

(Check runs, near-optimality stop, High performance plan; the clean figures
are those of the batch below.) With it step B of the exact path costs about
what that of the diagonal path does (57 ms on large10k), so what remains in
favour of the diagonal Hessian is its iteration count on med2522 at the 0.5%
stop.

### Power plan

Same configuration, same instance, same iterates, backward pass per iteration:
med2522 `T=24` 1.68 s (2026-10-05) and 1.62 s (2026-10-06) under Balanced,
1.28 s under High performance (2026-10-09); large10k `T=6` 1.98 against
1.70 s. One run each on different days, load not logged for the last, so
indicative: 14-24%. Every time in the TPEC table was measured under Balanced,
for both solvers.

### Stopping rules on one footing

- `IPOPT_TRACE=1` (`centralized_ipopt_matched.jl`): one `IPOPT_ITER` line per
  iteration through the intermediate callback of Ipopt (clock, objective,
  primal and dual infeasibility, barrier parameter). `ipopt_trace_posthoc.py`
  reads off when a rule was first met, as `near_opt_posthoc.py` does for
  FilterDDP. Same iterations and objective with and without it.
- **Two primal infeasibilities, and they differ a great deal.** The iteration
  table of Ipopt prints the violation of the constraints as posed: zero for
  the cone whenever `P^2 + Q^2 <= v ell` holds. Its callback gives the
  internal one, which counts an inequality through its slack. The
  `primal_inf` of FilterDDP is of the internal kind (the cone row carries an
  explicit slack). By the internal measure Ipopt reaches the 0.5% point at
  83-84% of its solve time on ieee123 `T=24` and med2522 `T=6`; by its printed
  one, at 64-65%. FilterDDP now logs both (`FILTERDDP_ORIGINAL_PRIMAL`, the
  as-posed one through `install_original_primal_hook`); `near_opt_posthoc.py
  --primal=original` uses it. Which one the paper applies is the call of the
  user; it must be the same for both solvers.
- User, 2026-10-10: near-optimality always means the objective criterion and
  the infeasibility criterion together; FilterDDP runs go to a tolerance of
  `1e-5` and centralized runs to `1e-6`, with the objective and infeasibility
  of every iteration logged, and the reported criterion is chosen afterwards.

### The OpenDSS check belongs to the run

`run_parsim.sh` now keeps the dispatch of each run and hands it to OpenDSS
before anything else starts (`OPENDSS_CHECK=1`, default). It had skipped the
solution write; there is no timing conflict in checking after the timer has
stopped. User, 2026-10-11: every dispatch goes through OpenDSS, the
centralized ones included.

### Clean batch (2026-10-11, in progress when this was written)

`run_clean_hp_batch.sh`: Ipopt (MA57, MA97, MUMPS; `tol = 1e-6`; trace) and
FilterDDP (tolerance `1e-5`; diagonal Hessian, and exact Hessian with the lean
step B) on ieee123 and med2522 at `T = 6, 24, 96` and large10k at
`T = 6, 24, 48, 96`, High performance plan, quiet wait and background load
logged per case, full warm-up. Logs: `results/ipopt_hsl/logs/*_hp_*`,
`results/parallel_in_time/logs/*_hp_tol1e5.log`, `results/opendss_check/parsim_*`.

### Files

| | |
|---|---|
| `pint_scan.jl`, `run_pint_scan.sh`, `pint_scan_summary.py` | the pairwise scheme, its launcher, the summary |
| `tree_kkt.jl` | `tree_kkt_block0`, `tree_kkt_reduced_rhs` (stage data for the diagnostic); `_tree_kkt_finish0` / `_tree_kkt_finishV` (lean exact step) |
| `DDP4OPF.jl/src/backward_pass.jl` | `PINT_SCAN_HOOK`, `PRIMAL_ORIGINAL_HOOK`, `_lean_curvature` |
| `centralized_ipopt_matched.jl`, `ipopt_trace_posthoc.py` | the Ipopt trace |
| `run_clean_hp_batch.sh`, `run_clean_hp_opendss.sh` | the clean batch; dispatches at a looser stop |
| `results/pint_scan/logs/` | `*_scan.log` (first batch), `*_scan2.log` (collector kept out of the timings) |

## 11. Clean batch results, reuse after a barrier update, barrier schedules, the passes on threads (2026-10-11)

Branch `ddp-oct10-new-tests`.

### The clean batch

All 20 FilterDDP runs converged and all 20 dispatches passed OpenDSS (no
voltage outside its limits; voltages within `1e-5` pu of the optimizer's).
Time to the 0.5% point in seconds, primal infeasibility below `1e-6` (ieee123)
or `1e-5` (med2522, large10k), each solver's own (internal) infeasibility.
"Count" is one worker per period: one core each, no communication; Ipopt is
one process. `clean_hp_table.py` prints this and the other rules.

| case | Ipopt MA57 | DDP diagonal, sequential | DDP diagonal, count | DDP exact, count |
|---|---|---|---|---|
| ieee123 `T=6` | 0.10 | 0.33 | 0.19 | 0.16 |
| ieee123 `T=24` | 0.33 | 1.31 | 0.63 | not met |
| ieee123 `T=96` | 1.68 | 16.5 | 9.38 | 4.34 |
| med2522 `T=6` | 2.62 | 10.3 | 4.56 | 5.79 |
| med2522 `T=24` | 15.2 | 47.9 | 14.3 | 14.8 |
| med2522 `T=96` | 76.5 | 245 | 55.2 | 49.7 |
| large10k `T=6` | 5.16 | 16.3 | 7.56 | 13.5 |
| large10k `T=24` | 22.0 | 57.5 | 19.3 | 39.3 |
| large10k `T=48` | 47.2 | 125 | 37.0 | 70.9 |
| large10k `T=96` | 101 | 263 | 67.1 | 150 |

- **large10k's threshold is `1e-5` since this batch** (user, 2026-10-11; it was
  `1e-4`). At `1e-4` the Ipopt column reads 4.83 / 20.5 / 44.3 / 84.1 and the
  diagonal count 5.98 at `T=6`, unchanged otherwise: Ipopt's own infeasibility
  falls below `1e-5` only near the end of its solve, the FilterDDP iterate at
  the 0.5% point is already below it from `T=24` up.
- **The table depends on the infeasibility measure.** With the violation of
  the constraints as posed (`--primal=original`) Ipopt qualifies much earlier
  (large10k `T=96`: 32.5 s, iteration 3) and FilterDDP's count is ahead only
  on med2522 with the exact Hessian (`T=24`: 9.0 against 10.3 s; `T=96`: 31.8
  against 45.2 s). On large10k the power-flow start, batteries idle, is
  feasible as posed to `1.2e-6` and 3 to 4% from the optimum; a few
  iterations bring the objective within 0.5%, and by that measure an iterate
  whose inequalities hold counts as feasible whatever its slack residuals.
- **Run to tolerance, the exact Hessian is faster than the diagonal one in all
  ten cases** (large10k `T=96`: 505 s and 16 iterations against 1085 s and
  58), and sequential FilterDDP is 3 to 5 times slower than MA57 (0.10 / 0.36
  / 1.90, 2.93 / 17.2 / 83.5, 5.17 / 25.1 / 53.0 / 113 s).
- **Why the exact Hessian qualifies late by the internal measure.** It leaves
  each barrier level after three or four iterations, as soon as the level's
  own tolerance (ten times the barrier parameter) is met, and the first steps
  at the next level are long again; the residual of the cone rows, which is
  second order in the step, has no time to settle. The diagonal Hessian stays
  longer at each level and ends it with a residual of `1e-7` or less. On
  ieee123 `T=24` the exact run ends at `3.7e-6` against the threshold of
  `1e-6` ("not met"): a FilterDDP tolerance of `1e-5` is not tight enough for
  that threshold with the exact Hessian.
- **FilterDDP at tolerance `1e-5` stops with the barrier parameter at `1e-6`**
  (a tenth of the tolerance), which leaves its objective above Ipopt's by up
  to `2.7e-5` (med2522 `T=96`): a `1e-5` objective rule cannot be met there.

### No new factorization after a barrier update (`FILTERDDP_REUSE_BARRIER=1`)

When a barrier problem is solved, `solve!` lowers the barrier parameter and
repeats the backward pass at the same trajectory: 8 times in a run to `1e-5`
with the default schedule, 10 to 27% of the passes. No stage matrix depends
on the barrier parameter (the barrier terms of the Hessian are `z ./ slack`,
and the value Hessian is built from the stage matrices), only the right-hand
side of the step does. `reuse_pass.jl` therefore keeps every stage's factor
for the repeated pass: one solve per stage and vector work, no derivatives,
no assembly, no factorization. A reused pass takes 6 ms against 43 ms on
ieee123 `T=24` and about 0.2 s against 1.0 s on large10k `T=6`.

The iterates are those of the full pass exactly: on ieee123 `T=6/24/96`,
med2522 `T=6/24` and large10k `T=6/24`, both Hessians, the runs follow the
clean batch step for step with objectives equal to all 16 digits
(`compare_runs.py`, which now matches residual lines by iteration).

### Barrier schedule and hybrid Hessian (132 untimed runs)

`run_barrier_sweep.sh`, `barrier_sweep_summary.py`. The solver's schedule is
`mu <- max(tol/10, min(kappa_mu mu, mu^theta_mu))` from `mu_init`, a level
being left when its error is below `kappa_eps mu`: defaults (1, 0.2, 1.2, 10),
9 levels down to `1e-6`; Ipopt's defaults (0.1, 0.2, 1.5, 10), 6 levels.
Iterations to tolerance summed over ieee123 `T=24/96`, med2522 `T=6/24`,
large10k `T=6/24`:

| schedule | diagonal | exact |
|---|---|---|
| base (1, 0.2, 1.2, 10) | 252 | 165 |
| Ipopt's (0.1, 0.2, 1.5, 10) | 187 | 157 |
| `mu_init` 0.1 alone | 212 | 148 |
| `mu_init` 0.01, `theta_mu` 1.5 | 188 | 158 |
| `theta_mu` 1.5 alone | 271 | 168 |
| `kappa_eps` 30 | 225 | 163 |

The largest single effect is ieee123 `T=96` with the diagonal Hessian: 56
iterations at base, 22 with Ipopt's schedule (at base it spends 32 iterations
at two late levels). To the 0.5% rule with the internal infeasibility the
picture is mixed for the diagonal Hessian: Ipopt's schedule gives
20 / 21 / 38 / 39 / 15 / 14 iterations against 19 / 55 / 33 / 39 / 13 / 11 at
base. It repairs ieee123 `T=96` and costs two or three iterations on
large10k, where the base schedule is the best one tried. No schedule is best
everywhere.

The hybrid (diagonal Hessian until the barrier parameter falls to a value,
exact from then on; `FILTERDDP_HYBRID_MU`) ends with the exact Hessian's
iteration counts (167 against 165 with the base schedule) and at the 0.5%
rule is never ahead of the better of the two Hessians. Nothing gained.

### The passes on threads (`threaded_pass.jl`, `FILTERDDP_THREADED=1`)

One worker per period was a count read from a sequential run. The passes
now exist in that form, on Julia threads, loaded through
`BACKWARD_PASS_HOOK`, `AFFINE_DIRECTIONS_HOOK` and `AFFINE_ROLLOUT_HOOK`:

- backward, phase 1, all stages at once: derivatives, barrier terms, the
  stage matrix, the network factorization and its preparation for the battery
  block, the first network solve;
- phase 2, in sequence: what carries the value function of `t+1` (battery
  right-hand side, battery block, battery solution, value update);
- phase 3, all stages at once: the second network solve, the multiplier
  updates, the stored policy and the battery rows of the feedback;
- forward: the states of all periods from a recursion on those rows, then
  every period's policy solve, directions and ratio test at once, and every
  period's part of each trial point at once.

Each stage does the solver's own operations in the solver's order, so the
backward pass reproduces the sequential iterates to the last bit (ieee123
`T=6/24`, med2522 `T=6`, large10k `T=6`, both Hessians, with and without
reuse, 4 to 8 threads). With the threaded forward pass the accepted steps and
backtracks are the same and the objectives agree to `1e-16`; the state
recursion differs from the rollout in rounding. Outside the supported
configuration, or on a failed factorization, the solver's own pass runs.

Checked for correctness only. The user does not want parallel timings from
the 309 PC; this is the code for a machine with a core per period (Kamiak
account requested 2026-10-11, ticket KSD-7186). Still sequential: phase 2,
which for the exact Hessian is where the pairwise combination of Section 10
would go.

### Files

| | |
|---|---|
| `DDP4OPF.jl/src/reuse_pass.jl` | the repeated pass after a barrier update |
| `threaded_pass.jl` | the backward and forward passes with one task per period |
| `tree_kkt.jl` | `tree_kkt_feedback_rows` (battery rows of the feedback after the exact sweep) |
| `run_barrier_sweep.sh`, `barrier_sweep_summary.py` | barrier schedules and the hybrid, iteration counts |
| `clean_hp_table.py` | times at a stopping rule; `--primal=original`, `--threshold=`, `--tag=`, `--arms=` |
| `results/barrier_sweep/logs/`, `results/threaded_pass/` | the sweep; the first threaded runs |

## 12. The pairwise sweep inside the solve, and the collector in the count (2026-10-11 afternoon to night)

Branch `ddp-oct10-new-tests`. High performance plan throughout. Every run
below was followed by its OpenDSS check; none violates a voltage limit.

### `T=144`

Ipopt MA57 at the reporting rule, default barrier update / adaptive:
ieee123 3.62 s / not met (3.44 s to tolerance), med2522 119.5 / 110.2 s,
large10k 157.0 / 138.9 s. FilterDDP with the ordinary sweep, sequential /
counted: med2522 483.3 / 106.9 s (diagonal) and 511.1 / 70.1 s (exact);
ieee123 17.9 / 10.4 s and 10.3 / 5.1 s. The exact Hessian's iterations hardly
move from `T=96` (66) to `T=144` (69) on med2522; the diagonal's go from 49
to 65. No FilterDDP run on large10k `T=144`: the user wanted it only "maybe",
and the PC has 31.7 GB (Ipopt alone takes 13.6 GB there). `T=288` is not to
be run without asking.

### The sweep by pairs in a real solve (`paired_pass.jl`)

A pass is: step 1, every stage at once (`_tp_pre` of `threaded_pass.jl`, and
the stage as a function `G_t(a,b)` of the battery energies before and after
it); step 2, the rounds (`ceil(log2(T+1))` of them, in each of which every
function is combined with the first one of the block to its right); step 3,
every stage at once with the value function after it. `FILTERDDP_PAIRED_SWEEP=count`
does all of it in sequence on one BLAS thread and times every piece;
`=1` runs the three steps on Julia threads (correctness only). It runs
through the wrapper `results/paired_pass/paired_main.jl` and the launcher
`results/paired_pass/run_paired.sh <system> <T> <count|1> [tag]`; the driver
has no switch of its own yet.

Same iterations as the ordinary exact solve in every case run (twelve),
objectives equal to 15 digits, value functions within 1.2e-9 of the ordinary
sweep's where the two were run side by side (`PAIRED_CHECK=1`).

The first counted run (med2522 `T=96`, allocating version) gave 65.0 s
against 47.1 s for the ordinary exact run, although the rounds counted
0.037 s a pass against 0.41 s for the sweep. The rest of that section is why.

### Kept storage (`FILTERDDP_PAIRED_INPLACE=1`, the launcher's default)

Each stage's function lives in storage kept from pass to pass, and a
combination writes its result over its left operand with three `n_B x n_B`
blocks of scratch per worker. In a round the left operands are distinct and
none of them is a right operand, so this is safe on threads. Before, a pass
allocated three blocks per stage and six per combination (1.1 GB on med2522
`T=96`, 18 GB on large10k). Same solve as the allocating version.

### The collector

Measured, not assumed:

- It takes 17-22% of a run on the two larger feeders (74 of 421 s on med2522
  `T=96`, pairwise).
- A pass of med2522 `T=24` allocates 0.71 GiB in step 1 (29 MB a stage),
  nothing in the rounds (kept storage) and 0.11 GiB in step 3; the collector
  takes about 0.2 s of a 1.25 s pass, nearly all in step 1.
- Held off for the whole pass and run once, it takes 0.15 s and the stages
  get slower by as much (fresh memory instead of reused). Its cost follows
  the volume allocated, not how often it runs. No switch was kept for this.
- The allocation has no single source (`alloc_profile.jl` with
  `FILTERDDP_PAIRED_SWEEP` set): the largest site is 8.6% of the bytes; the
  tree solver's per-node arrays are about a third of the allocations.
  `tree_kkt_inplace.jl` (written, not wired, untested) would take the factor's
  share only.
- A pause lands in whichever stage is on the clock. In the forward pass of
  med2522 `T=144` single periods showed 0.16 to 0.80 s against a mean of 5 to
  7 ms, and the forward pass was 22.6 of 33.5 s counted.

So the timing was changed. Counting runs take each stage's time without the
collections that fell in it and print the collector's time separately:
`pre_gc_s`, `post_gc_s`, `pair_gc_s` on `FILTERDDP_PARSIM_BACKWARD` (and the
volumes on `FILTERDDP_PAIRED_GC`), and, with `FILTERDDP_PARSIM_GCNET=1`,
`gc_s` and `stages` on `FILTERDDP_PARSIM_FORWARD` (solver source:
`affine_linesearch.jl`). `parsim_from_log.py` and `clean_hp_table.py` then
count a worker per period as collecting its own share, `1/stages` of the
collection in the per-period work; `--gc=sequential` (table script) counts
all of it in sequence. Collections during sequential work stay sequential.

**The user's decision (2026-10-11): per-worker collection is the default
assumption of the count, and the paper says so.** It is the same premise as
the rest of the count: one worker per period, with its own memory. It also
means a real parallel run should be one process per period, not threads in
one program, where every pause stops all the threads.

Logs dated before this (the clean batch, the first pairwise batch tagged
`inplace_hp_tol1e5`) have the pauses inside the stage times: the slowest stage
of a pass carries its longest pause and the pauses in the other stages are
dropped. Their counts are upper bounds of the per-worker reading, not the
same thing.

### The ordinary sweep on the same timing (`FILTERDDP_PAIRED_ROUNDS=0`)

Steps 1 and 3 as in the pairwise pass and the ordinary sweep (`_tp_seq`)
between them, stage after stage. One BLAS thread like the rest of the pass,
which costs the sweep on large10k: 103 ms a step against 66 ms in the
solver's own exact pass, which lets the sweep use all cores. For large10k
the solver's own counts (clean batch: 9.4, 27.3, 51.6, 110.1 s) are the
ordinary sweep's best figures.

### Results (0.5% rule, per-system thresholds; the pairwise figures are the same at 0.1%)

Seconds. Ipopt: its faster variant. Ord. and pair.: counted, one worker per
period, per-worker collection, exact Hessian. In seq.: the pairwise count
with every collection in sequence. One proc.: the pairwise run as measured.

| | `T` | Ipopt | ord. | pair. | Ipopt/pair. | in seq. | one proc. |
|---|---|---|---|---|---|---|---|
| ieee123 | 6 | 0.10 | 0.09 | 0.07 | 1.4 | 0.07 | 0.24 |
| | 24 | 0.33 | 0.23 | 0.14 | 2.4 | 0.66 | 1.59 |
| | 96 | 1.68 | 1.30 | 0.27 | 6.2 | 3.09 | 7.43 |
| | 144 | 3.62 | 1.66 | 0.35 | 10.3 | 4.11 | 10.5 |
| med2522 | 6 | 2.58 | 3.42 | 2.61 | 0.99 | 4.81 | 13.2 |
| | 24 | 13.3 | 9.25 | 4.16 | 3.2 | 16.7 | 71.9 |
| | 96 | 67.8 | 42.2 | 8.43 | 8.0 | 81.0 | 421.4 |
| | 144 | 110.2 | 72.0 | 14.4 | 7.7 | 153.1 | 750.3 |
| large10k | 6 | 3.97 | 12.3 | 8.32 | 0.48 | 11.5 | 34.1 |
| | 24 | 21.9 | 40.8 | 12.4 | 1.8 | 24.0 | 173.2 |
| | 48 | 44.4 | 78.0 | 15.8 | 2.8 | 37.6 | 390.3 |
| | 96 | 101.4 | 175.1 | 22.6 | 4.5 | 70.2 | 973.2 |

ieee123 `T=24` was run to a tolerance of 1e-6 (`PAIRED_TOL`): at 1e-5 it ends
with infeasibility 3.7e-6, above its 1e-6 threshold, and never meets the rule.

Where a counted pass goes, per pass at `T=96`: med2522 0.04 s slowest stage
of step 1, 0.03 s the seven rounds, 0.01 s step 3, 0.02 s forward; large10k
0.13, 1.0, 0.12 and 0.1 s. On large10k the rounds are the largest part: a
combination of 1020 batteries takes 0.14 s on one core.

Run with:

    python ddp/examples/power_system/clean_hp_table.py --tag=pairwise_gc_hp_tol1e5,reuse_gc_hp_tol1e5,reuse_hp_tol1e5 --arms=exact,paired --gaps=5e-3,1e-3
    python ddp/examples/power_system/clean_hp_table.py --tag=ordinary_gc_hp_tol1e5 --arms=exact,paired --gaps=5e-3

(the `paired` column is the log `fddp_<system>_T<T>_paired_<tag>.log`; the
`exact` column is the clean batch's own exact run).

### What the count leaves out, and what is open

- It is a count from one machine: one core per period and per combination,
  each worker with its own memory, no time for passing data. A combination
  needs three `n_B x n_B` blocks of another period: 1.5 MB on med2522, 25 MB
  on large10k. Ipopt is one process.
- One process runs the pairwise arrangement slower than the ordinary sweep
  (973 against 351 s on large10k `T=96`): it does every combination in a row.
- The repeated passes after a barrier update are counted in sequence (0.8 s
  of the 8.4 s on med2522 `T=96`); they could be split like the others.
- Allocation: reducing what a pass allocates would shrink the collector under
  either way of charging it. No single source; not started.
- Not measured: any of this on a machine with a core per period (Kamiak,
  ticket KSD-7186, one process per period).

### Longer horizons: med2522 `T=288` and large10k `T=144` (2026-10-12, user's go)

At the 0.5% rule (pairwise figures the same at 0.1%), seconds:

| | Ipopt | DDP seq. (diag / exact) | ord., same timing | pair. | Ipopt/pair. | in seq. | one proc. |
|---|---|---|---|---|---|---|---|
| med2522 `T=288` | 255.2 | 1072.3 / 1248.8 | 203.4 | 26.4 | 9.7 | 495.7 | 2036 |
| large10k `T=144` | 138.9 | 401.2 / 632.8 | 351.9 | 80.6 | 1.7 | 192.9 | 2238 |

Ipopt: the adaptive variant in both (default 255.7 and 157.0 s). The solver's
own passes count 237.7 s (diagonal) and 173.7 s (exact) on med2522 `T=288`,
99.9 and 194.9 s on large10k `T=144`. Same iterations for the two sweeps (75
and 18), objectives equal to 15 digits, all ten dispatches pass OpenDSS.
Warm-up 12 iterations on med2522 `T=288` (its first barrier update is near
iteration 9), 8 on large10k.

**large10k `T=144` is short of memory on the 309 PC (31.7 GB).** A
three-iteration test (`chk_memprobe`) showed Julia at 17-20 GB with 1.5-5 GB
free. In the timed pairwise run the slowest period of step 3 takes 1.3-2.4 s
in 15 of 19 passes (0.2-0.4 s in the other four, 0.12 s at `T=96`) and the
last four passes have one combination of about 2 s against 0.14 s. The
80.6 s carries those stalls; about 45 s without them, which is an estimate.
The cause of the slow combinations in the last passes was not established
(memory, or very small numbers late in the solve). A clean figure needs a
rerun with the other applications closed, or a machine with more memory;
`T=288` on large10k cannot run here at all.

On med2522 `T=288` the collector takes about a quarter of the one-machine
run, which is why the in-sequence figure (495.7 s) is so far from the
per-worker one.

### Corrections to commit messages

- `66d38fb9` gives 137 s and 74 s for the pairwise counts with every
  collection in sequence at med2522 `T=144` and large10k `T=96`. The measured
  figures are 153.1 s and 70.2 s (table above).
- `71b96878` states that collector pauses land in the slowest stage of step 1.
  At that point it was an inference from the 12-fold gap between the slowest
  and the mean stage; it was measured afterwards (the forward-pass figures
  above, and the net-of-collection maxima: 0.04 s against 0.5 s a pass).

### Files

| | |
|---|---|
| `paired_pass.jl` | the pairwise pass; kept storage; `FILTERDDP_PAIRED_ROUNDS=0` for the ordinary sweep on the same timing |
| `results/paired_pass/run_paired.sh`, `paired_main.jl` | launcher (OpenDSS check included) and wrapper; `PAIRED_WARMUP`, `PAIRED_INPLACE`, `PAIRED_ROUNDS`, `PAIRED_TOL`, `PAIRED_CHECK`, `PAIRED_MAIN` |
| `DDP4OPF.jl/src/affine_linesearch.jl`, `backward_pass.jl`, `forward_pass.jl` | `FILTERDDP_PARSIM_GCNET=1`: forward per-period times without collections |
| `clean_hp_table.py`, `parsim_from_log.py` | read the collector fields; `--gc=sequential` |
| `alloc_profile.jl` | allocation profile, of the pairwise pass when `FILTERDDP_PAIRED_SWEEP` is set |
| `tree_kkt_inplace.jl` | the tree factor in kept storage: written, not wired, untested |
| logs `*_paired_pairwise_gc_hp_tol1e5.log`, `*_paired_ordinary_gc_hp_tol1e5.log` | the final runs; `*_paired_inplace_hp_tol1e5.log` the first batch (forward pass on the old timing); `*_paired_count_hp_tol1e5.log` the allocating version |
