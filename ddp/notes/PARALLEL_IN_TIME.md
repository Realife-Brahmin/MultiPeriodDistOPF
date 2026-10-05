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
