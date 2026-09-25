# MA57 and fill-reducing orderings for FilterDDP's stage KKT solves

Agenda item for R. Gupta, meeting of 2026-09-25: test MA57 as a sparse linear
solver, compare the reordering and fill-in of Ipopt's solver and FilterDDP's
UMFPACK, and test whether another ordering reduces fill or time. Measured
2026-09-24 on the 309 lab PC (i9-10900X, 31.7 GiB), all runs sequential.

Not repeated here, because they are already established: sequential UMFPACK
beating sequential MUMPS on captured matrices; this Ipopt build using
`MUMPS_seq` whatever the thread variables; reuse of UMFPACK's symbolic analysis
slowing the full ieee2522 `T=3` run; exact-versus-diagonal Hessian convergence
across horizons; the existing sparsity plots of the three systems.

## Answers

1. **Is MA57 available?** No. There is no HSL binary for Julia or Ipopt; the
   only licensed MA57 on the machine is inside MATLAB, which would not run
   headless. Nothing below is an MA57 result. (Section 1.)
2. **Does exploiting symmetric-indefinite structure beat UMFPACK's LU?** Not
   with the symmetric-indefinite solver we have. MUMPS LDLᵀ (`sym=2`, the same
   algorithm family as MA57 and Ipopt's own solver) is **4.5-11x slower** on
   factorization plus wide solve in all six cases, and 3.7-13x slower at
   factorization alone. UMFPACK's own symmetric strategy is 1.5-2.5x slower at
   large10k.
3. **Which ordering minimizes fill?** Minimum-degree orderings, everywhere.
   Among the orderings UMFPACK actually factors with, its default
   (COLAMD, unsymmetric strategy) is within 3% of the least-fill choice at
   large10k and has the least fill outright on the exact Hessian (1.45x
   `nnz(K)`). By ordering-only fill, AMD (UMFPACK symmetric or MUMPS) wins every
   case. Nested dissection (METIS, SCOTCH) never wins.
4. **Does lower fill reduce the wide-solve time?** Within UMFPACK at large10k,
   yes, monotonically. At med2522 not reliably, and across solvers not at all:
   MUMPS stores 1.7-3.1x more factor entries but its single-RHS solve is
   32-77x slower, so per-call overhead, not fill, decides.
5. **Does the best isolated gain survive a full FilterDDP run?** No. The only
   clear isolated gain (ieee123 exact, UMFPACK symmetric strategy, 0.79x)
   gives identical iterations and solutions, and a wall-time change inside
   run-to-run noise. It reverses at med2522 and large10k, so longer runs were
   not justified.
6. **Does diagonalizing the Hessian change the best ordering?** No. UMFPACK's
   default stays best or tied in both arms at med2522 and large10k.
   Diagonalization changes cost, not ranking: at large10k it halves
   factorization (0.31 -> 0.15 s) and cuts the wide solve by 40%
   (3.59 -> 2.16 s).
7. **Is solver replacement still important after diagonalization?** Little,
   for the factorization. At large10k with the diagonal Hessian, factorization
   is 7% of the stage's linear algebra, so even an infinitely fast factorization
   saves at most 7%. A replacement matters only if its **multi-RHS solve** is
   much faster than UMFPACK's.
8. **Is the (n_x+1)-column solve still the dominant linear-algebra cost?** Yes,
   more than before: 91% (diagonal) and 89% (exact) of ordering + factorization
   + solve at large10k, 77% and 61% at med2522. Only at ieee123 is it a
   minority.

## 1. MA57 availability

Evidence: `ddp/results/kkt_ordering/ma57_availability_julia.txt`, from
`ddp/examples/power_system/ma57_availability_probe.jl`.

| Route | Finding |
|---|---|
| Julia packages | `HSL.jl` and `HSL_jll` not installed. Solver stack in `envs/ddp2026`: Ipopt.jl 1.15.0 / Ipopt_jll 300.1400.1902 (Ipopt 3.14.19), MUMPS.jl 1.6.2 / MUMPS_seq_jll 500.900.100 (MUMPS 5.9), SuiteSparse_jll 7.8.3, METIS_jll 5.1.3, SCOTCH_jll 7.0.11, SPRAL_jll 2025.9.18 |
| HSL libraries | none of `libhsl`, `libcoinhsl`, `libma57`, `libma27`, `libhsl_ma57`, `libhsl_ma97` loadable |
| Ipopt `linear_solver` | `mumps` works. `ma27`, `ma57`, `ma77`, `ma86`, `ma97` all fail: `DYNAMIC_LIBRARY_FAILURE ... Error 126 while loading DLL libhsl.dll`. The build supports HSL; the binary is missing. `pardiso` fails the same way (`libpardiso.dll`); `pardisomkl` and `wsmp` are not compiled in |
| MATLAB | R2023a and R2024b ship `bin/win64/libmwma57.dll`, used by sparse `ldl`, `decomposition(A,'ldl')` and backslash on symmetric indefinite matrices. `matlab -batch` of `ma57_availability_probe_matlab.m` printed nothing in over 10 min with 10 s of CPU (consistent with a sign-in or licence prompt) and was stopped |
| MATLAB's DLL from Julia | not attempted: undocumented interface, and calling MathWorks' licensed copy outside MATLAB would sidestep that licence |

**Exact blocker:** no `libhsl.dll`. HSL issues libHSL free for academic use;
it provides `HSL_jll` for Julia and the `libhsl` Ipopt looks for, so both
`Ipopt(linear_solver="ma57")` and a direct MA57 benchmark would work with no
code change. MATLAB's copy could give a restricted test (MATLAB's defaults
only: no ordering choice, analysis not separable from factorization) once
MATLAB runs headless.

Side result: Ipopt's `spral` (SPRAL SSIDS, an open-source symmetric-indefinite
solver from the HSL group) fails with `info.flag = -53` unless
`OMP_CANCELLATION=TRUE` and `OMP_PROC_BIND=TRUE`, and then solves correctly.

## 2. Method

**Matrices.** `capture_ordering_kkt.sh` ran FilterDDP on the matched protocol
(periodic profile, `C_B = 1e-3`, soft terminal SOC, `T = 3`) to iteration 20 and
kept stage 1's `K` and its full `(n_x+1)`-column right-hand side, for both the
exact and the diagonal stage Hessian. Iteration 20 is mid-solve: the barrier
parameter has fallen and the active set has moved away from the cold start.
The captures are gitignored (up to 1.6 GB); the logs record how they were made.

| | n | nnz(K), diag / exact | RHS columns |
|---|---:|---:|---:|
| ieee123 | 1,353 | 5,269 / 8,072 | 52 |
| med2522 | 23,691 | 96,025 / 162,818 | 250 |
| large10k | 96,968 | 393,075 / 1,453,094 | 1,021 |

`K` is exactly symmetric in every capture, with a structurally zero (2,2) block.

**Configurations** (`kkt_ordering_benchmark.jl`):

- **UMFPACK** (unsymmetric LU). `UMFPACK_default` is exactly FilterDDP's
  `lu(K)`: Julia's default control, which differs from SuiteSparse's C default
  (ordering AMD/COLAMD rather than CHOLMOD, iterative refinement off). UMFPACK
  reported using the unsymmetric strategy with COLAMD on every matrix. Also:
  unsymmetric COLAMD, symmetric AMD, unsymmetric and symmetric METIS, CHOLMOD
  (AMD/COLAMD, METIS if fill is high), and BEST (tries several, keeps the least
  fill).
- **MUMPS, `sym=2`** (LDLᵀ with 1x1/2x2 pivoting). ICNTL(7) = AMD, AMF, QAMD,
  PORD, SCOTCH, METIS, automatic. By default MUMPS orders a graph compressed
  into 2x2 pivot candidates, so AMD, AMF and METIS were also run on the plain
  graph (ICNTL(12)=1), plus METIS on the compressed graph explicitly, plus
  Ipopt's own MUMPS option set (automatic ordering, scaling 77, permuting
  scaling 7, pivot tolerance 1e-6).

**Timing.** One BLAS, OpenMP and OpenBLAS thread for every solver. Per
configuration, one warm-up repetition is discarded (compilation) and medians are
taken over 21 (ieee123), 9 (med2522) or 3 (large10k) repetitions of:
ordering/symbolic analysis; numerical factorization; a one-column solve; and the
full `(n_x+1)`-column solve done exactly as FilterDDP does it (`ldiv!` for
UMFPACK, one multi-RHS solve for MUMPS). Every configuration factorizes the same
`K` and solves the same right-hand side; the relative residual
`||KX - B||_F / ||B||_F` is checked for both solves. MUMPS's `INFOG(1)` is
checked after every phase.

**Fill, defined.** *Stored fill*: UMFPACK `(nnz(L)+nnz(U))/nnz(K)`; MUMPS
`2 x INFOG(29) / nnz(K)`, i.e. the entries an LU with the same pattern would
hold. MUMPS stores dense frontal blocks, so its count includes explicit zeros.
*Ordering-only fill*: the pivot order applied symmetrically to `K`'s pattern,
with the exact Cholesky fill of an SPD surrogate (`|K|` plus a dominant
diagonal), as `2 nnz(L)/nnz(K)`. That measures the ordering alone, with no
pivoting and no supernode amalgamation. For UMFPACK it uses the column order,
which only approximates what an unsymmetric factorization does.

**Shared machine.** Before each case the runner waited for background load
below 1.5 cores, and logged everyone else's CPU during it. Median background
load was 0.27-0.61 cores in every case, the machine's idle floor.

## 3. Results

All 114 configurations solved, with the worst relative residual 2.1e-10. Every
MUMPS factorization reports exactly `nc` negative pivots and no null pivots:
the captured systems have the correct KKT inertia.

| System | Hessian | UMFPACK default: order / factor / wide solve (s) | wide share | best other UMFPACK (f+w vs default) | best MUMPS LDLᵀ (f+w vs default) | stored fill: UMFPACK default / MUMPS auto |
|---|---|---|---:|---|---|---|
| ieee123 | diag | 0.0006 / 0.0014 / 0.0009 | 31% | BEST 0.97x | Ipopt options 11.3x | 2.09 / 5.35 |
| ieee123 | exact | 0.0011 / 0.0038 / 0.0015 | 23% | symmetric AMD 0.79x | Ipopt options 4.6x | 2.77 / 6.61 |
| med2522 | diag | 0.011 / 0.023 / 0.112 | 77% | BEST 0.94x | SCOTCH 10.1x | 2.24 / 5.73 |
| med2522 | exact | 0.036 / 0.066 / 0.159 | 61% | unsym COLAMD 0.99x | AMD plain 5.9x | 2.65 / 4.63 |
| large10k | diag | 0.053 / 0.154 / 2.16 | 91% | unsym COLAMD 0.99x | SCOTCH 8.1x | 2.22 / 6.91 |
| large10k | exact | 0.148 / 0.31 / 3.59 | 89% | unsym COLAMD 1.05x | QAMD 4.5x | 1.45 / 2.73 |

"f+w" is factorization plus wide solve; "wide share" is the wide solve's share
of ordering + factorization + wide solve for the baseline. The per-configuration
tables are `ddp/results/kkt_ordering/SUMMARY_TABLES.md`, and every number is in
`kkt_ordering_benchmark.csv`.

**No alternative ordering helps at scale.** The apparent winners are ties
(0.94-1.05x), or pay for themselves in ordering time. BEST's 0.94x at
med2522 diag costs 0.24 s of ordering against a 0.14 s factor+solve, and
FilterDDP refactorizes every stage of every iteration, so that cost recurs. At
large10k the alternatives lose outright:

| large10k, exact | stored fill | wide solve (s) | f+w vs default |
|---|---:|---:|---:|
| UMFPACK default (COLAMD) | 1.45 | 3.59 | 1.00x |
| CHOLMOD | 1.67 | 4.33 | 1.20x |
| BEST | 1.78 | 4.48 | 1.25x (plus 11.6 s ordering) |
| symmetric AMD | 2.00 | 5.25 | 1.47x |
| symmetric METIS | 2.23 | 6.03 | 1.68x |
| unsymmetric METIS | 4.52 | 9.88 | 2.70x |

Here the wide-solve time rises with stored fill, monotonically. At med2522 exact
the relation breaks: unsymmetric METIS stores less (2.35 vs 2.65) yet solves
slower (0.179 vs 0.159 s).

**Why the symmetric-indefinite solver loses.** MUMPS's factorization is
3.7-13x slower than UMFPACK's, and its solves far more so: a single right-hand
side takes 150 ms at large10k against UMFPACK's 2.5 ms, and 34 ms against
0.5 ms at med2522 (32-77x across the six cases). Its stored factor is only
1.7-3.1x larger, so the gap is per-call solve overhead in this MUMPS build, not
arithmetic. MA57's solve phase is a different implementation, which is why this
result does not settle the MA57 question. MUMPS also delays pivots (2,200-10,800
of 97k at large10k exact, depending on ordering), which enlarges fronts. Ipopt's
option set (pivot tolerance 1e-6) delays fewer at med2522 (4 against 304 with
the diagonal Hessian) but leaves residuals ~1e-11 instead of ~1e-14.

**Why nested dissection never wins.** The figures show it. Every ordering
places the dense battery-curvature block (exact Hessian) and the battery
coupling rows last, and that trailing block holds most of the factor: 1.08M of
UMFPACK's 2.0M entries at large10k exact. The rest of the graph comes from
radial distribution feeders, which are tree-like, and minimum degree eliminates
a tree with almost no fill. Nested dissection instead cuts the graph with
separators, and those separators fill against the dense block (the scattered
off-diagonal fill in the METIS zoom panels). This is the opposite of
mesh-like systems, where nested dissection is the right choice. It is an
explanation consistent with the figures, not a separate experiment.

Figures (`ddp/results/kkt_ordering/figures/`), for med2522 and large10k in
both Hessian arms. Each shows the captured `K`, then for UMFPACK's default,
MUMPS's automatic (Ipopt-like) ordering and symmetric METIS: the reordered
matrix, its factor, and a zoom on the factor's trailing 4,000 rows and columns.
The MUMPS factor panels are the ordering-only Cholesky pattern, since MUMPS
keeps its factors internal.

## 4. Full FilterDDP runs

The gate for integrating an ordering was a meaningful isolated benefit. Only
one case cleared it: ieee123 exact with UMFPACK's symmetric strategy (0.79x on
f+w, residual 5.5e-11 -> 1.6e-13). It was tested end to end through a new
opt-in hook, `FILTERDDP_UMFPACK_STRATEGY` / `FILTERDDP_UMFPACK_ORDERING` in
`ddp/DDP4OPF.jl/src/backward_pass.jl`. Unset, FilterDDP still calls `lu(K)`
exactly as before. Runs were ieee123 `T=3`, both arms, baseline and candidate
alternated, two repeats each, with full per-iteration logging and
near-optimality against the matched Ipopt objective (2743.846) at primal `1e-6`
(`run_umfpack_strategy_fullrun.sh`, traces in
`ddp/results/kkt_ordering/fullrun/`):

| arm | variant | iterations | wall (s) | time to near-optimality (s) | objective |
|---|---|---:|---|---|---|
| exact | baseline | 47 | 20.75, 21.46 | 18.41, 18.67 | 2743.846660020481 |
| exact | symmetric | 47 | 20.64, 20.50 | 18.06, 18.32 | 2743.846660020481 |
| diag | baseline | 52 | 21.15, 20.69 | 18.54, 18.38 | 2743.846660021356 |
| diag | symmetric | 52 | 20.55, 20.39 | 18.17, 17.85 | 2743.846660021355 |

Solutions agree to `3.5e-15` in controls and `7e-18` in states, and every run
has max equality residual ~6e-7. The wall-time difference (~2%) is inside the
run-to-run spread. At ieee123 the summed factorization time is 0.8-1.6 s of an
~18 s solve and varies by 0.5 s between identical repeats. Because the isolated
benchmark shows the same strategy 1.2-2.5x *slower* at med2522 and large10k,
the med2522 `T=3`/`T=12` and large10k `T=3` runs were not launched: no
configuration showed an improvement worth carrying to scale.

## 5. Recommendation for the meeting

- **Keep UMFPACK with its default control.** At the sizes that matter, no
  ordering and no symmetric-indefinite solver we can run beats it.
- **Get the HSL academic licence** to close the MA57 question properly. It is
  the only blocker, and it needs no code change. Expect MA57 to matter only if
  its multi-RHS solve is much faster than UMFPACK's: after diagonalization,
  factorization is 7% of large10k stage linear algebra.
- **The lever is the (n_x+1)-column solve** (89-91% at large10k), not the
  ordering. Its columns feed only the value-function term
  `beta' B + omega' c_x = -R' K^{-1} R` (with `R` the stacked `[B; c_x]`
  right-hand side), which is the Schur complement of the bordered matrix
  `[K R; R' 0]`. A solver that returns Schur complements (MUMPS `ICNTL(19)`,
  for example) could form it inside one factorization instead of 1,021
  triangular solves. The right-hand side is also sparse (active battery rows and
  `c_x` only), which UMFPACK's column-by-column solve cannot exploit. **Both are
  untested proposals.**
- Nested dissection is the wrong family for these tree-like feeder KKT
  systems. The comparison with NREL's documented large sparse systems (agenda
  item 1.4) was not done tonight.

## Reproduce

```bash
bash ddp/examples/power_system/capture_ordering_kkt.sh 20          # ~20 min, captures (gitignored)
bash ddp/examples/power_system/run_kkt_ordering_benchmark.sh       # ~55 min, benchmark CSV
python ddp/examples/power_system/summarize_kkt_ordering.py > ddp/results/kkt_ordering/SUMMARY_TABLES.md
bash ddp/examples/power_system/run_umfpack_strategy_fullrun.sh ieee123C_1ph 3 2
julia --project=envs/ddp2026 ddp/examples/power_system/ma57_availability_probe.jl
```

Figures: `kkt_ordering_patterns.jl` then `plot_kkt_ordering_patterns.py`
(usage in each file's header).

A note on the earlier comparison, `ddp/results/agenda_pipeline/C_solvers_*.csv`:
several of its MUMPS rows have relative residuals near `1e+22`, i.e. failed
solves recorded as timings. This benchmark checks `INFOG(1)` and residuals on
every solve and had no failures. The earlier conclusion (UMFPACK faster than
MUMPS) stands on valid solves.
