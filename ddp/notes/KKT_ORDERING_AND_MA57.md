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

1. **Is MA57 available?** Now yes. At first there was no HSL binary for Julia
   or Ipopt (Section 1); MA57 and HSL_MA97 have since been built locally from
   the user's licensed HSL sources and measured directly (Section 9).
2. **Does exploiting symmetric-indefinite structure beat UMFPACK's LU?** No,
   and MA57 and MA97 themselves confirm it (Section 9): they factor 2.8-6.5x
   faster at large10k but their 1,021-column solve is 4.4-6.9x slower.
   Before HSL was available: MUMPS LDLᵀ (`sym=2`, the same
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

**Follow-up (Section 6): the solve is slow because UMFPACK solves one column at
a time, and that is fixable without a new solver.** A blocked solve over
UMFPACK's own factors (16 columns at a time) is 1.9-3.0x faster on one thread,
with identical numerics. In complete FilterDDP runs (diagonal Hessian) it cuts
time to near-optimality by 13-21% at med2522 (`T` = 3 to 24) and **23-30% at
large10k** (`T` = 3, 6; 30% in the paper's Table II configuration), with the
same iterations and objectives. This is a real saving, but it does not close
the gap to centralized Ipopt. These results are in the TPEC paper (Sections
IV-E and IV-F, Tables V and VI, Fig. 10).

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

- **Keep UMFPACK's factorization and default ordering, and replace only its
  solve** with the blocked multi-column solve (Section 6):
  `FILTERDDP_BLOCKED_SOLVE=16`, opt-in, same answers, 13-30% faster FilterDDP
  runs measured so far. At present it is opt-in; making it the default is a
  one-line change once it has run on a longer horizon.
- **Get the HSL academic licence** to close the MA57 question properly. It is
  the only blocker, and it needs no code change. MA57 now has to beat the
  blocked solve, not UMFPACK's column-by-column one. After diagonalization,
  factorization is 7% of large10k stage linear algebra, so MA57's case rests on
  its multi-RHS solve.
- **None of this closes the gap to Ipopt.** Even free linear algebra would
  leave FilterDDP roughly 12-18x slower at large10k. The rest is FilterDDP's
  per-stage derivatives, assembly and rollout, and ~1.7x more iterations.
- **Untested next steps on the same solve:**
  - *Parallel blocks.* Column blocks are independent and the extracted factors
    can be shared read-only, which UMFPACK's own solve could not do. That gave
    6.5-9.2x on the isolated large10k solve with 8 threads.
  - *Schur complement.* The columns only feed the value-function term
    `beta' B + omega' c_x = -R' K^{-1} R` (with `R` the stacked `[B; c_x]`
    right-hand side), which is the Schur complement of the bordered matrix
    `[K R; R' 0]`. A solver that returns Schur complements (MUMPS `ICNTL(19)`,
    for example) could form it inside one factorization.
- Nested dissection is the wrong family for these tree-like feeder KKT
  systems. (The comparison with NREL's documented ACOPF matrices, agenda item
  1.4, already exists in the paper's KKT sparsity appendix; the new ordering
  figure, Fig. 10, sits beside it.)

## 6. Blocked multi-column solve over UMFPACK's factors

**The question.** UMFPACK's solve accepts one right-hand side per call, so
`ldiv!(F, B)` walks the entire factor once per column: 1,021 times at large10k,
at one multiply-add per factor entry read. Does solving many columns per pass
help, on the same factors?

**Method** (`blocked_multirhs_solve_benchmark.jl`). The same UMFPACK factors
are used, `L U = (Rs .* K)[p, q]`, extracted with `F.L, F.U, F.p, F.q, F.Rs`.
Columns are solved `w` at a time with the block stored transposed, so each
factor entry becomes one contiguous update of length `w`. The arithmetic is
identical to UMFPACK's. Extraction time is charged to the blocked method, and
the baseline is exactly FilterDDP's `ldiv!`. Captures, repeats and load
logging are as in Section 2, sequential. An 8-thread variant (blocks spread
over threads) is reported separately.

**Isolated result** (`blocked_multirhs_solve.csv`; speedup includes extraction):

| System | Hessian | UMFPACK `ldiv!` (s) | blocked, 1 thread (best w) | blocked, 8 threads |
|---|---|---:|---|---|
| ieee123 | diag / exact | 0.0009 / 0.0015 | 1.27x / 1.61x | -- |
| med2522 | diag | 0.10-0.12 | 1.87x (w=16) | 5.4x |
| med2522 | exact | 0.18-0.20 | 2.64x (w=32) | 7.1x |
| large10k | diag | 2.1-2.5 | 2.77x (w=16) | 6.5x |
| large10k | exact | 3.6-3.7 | 3.02x (w=32) | 9.2x |

Solutions match UMFPACK's to `5e-14` relative or better, with the same
residuals. At `w=1` the kernel is slower than UMFPACK (0.38-0.64x), so the gain
is the blocking, not a better kernel. The optimum is 16-32 columns, where one
block stays in cache; at 256 or all 1,021 columns the gain disappears. Factor
extraction costs under 1 ms (ieee123) to 13-46 ms (large10k).

**Full FilterDDP runs** (`run_blocked_solve_fullrun.sh`, opt-in hook
`FILTERDDP_BLOCKED_SOLVE=16` in `ddp/DDP4OPF.jl/src/blocked_solve.jl`).
Matched protocol, full per-iteration logging, near-optimality at the per-system
primal thresholds, baseline and blocked alternated, median background load
0.26-0.35 cores in every run. Traces and logs are in
`ddp/results/kkt_ordering/fullrun_blocked/`.

| Case | iterations | time to near-optimality: baseline -> blocked | saving | wide solve inside the run |
|---|---:|---|---:|---|
| ieee123 T=3, diag (2 repeats) | 52 | 17.1, 17.6 -> 17.6, 17.8 s | none | 0.20 -> 0.12 s |
| ieee123 T=3, exact (2 repeats) | 47 | 17.9, 17.8 -> 17.4, 18.0 s | none | 0.28 -> 0.17 s |
| med2522 T=3, diag | 55 | 87.2 -> 75.8 s | **13%** | 24.6 -> 12.1 s |
| med2522 T=3, exact | 42 | 91.5 -> 75.5 s | **17%** | 28.3 -> 11.2 s |
| med2522 T=12, diag | 73 | 311.0 -> 247.8 s | **20%** | 132.7 -> 56.7 s |
| med2522 T=24, diag | 79 | 613.8 -> 486.8 s | **21%** | 276.8 -> 123.7 s |
| large10k T=3, diag | 99 | 1785.0 -> 1368.2 s | **23%** | 685.9 -> 277.3 s |
| large10k T=6, diag | 103 | 3230.5 -> 2364.6 s | **27%** | 1489.1 -> 595.1 s |
| large10k T=6, diag, with assembly rewrites (Table II configuration) | 103 | 2991.9 -> 2103.9 s | **30%** | 1502.5 -> 571.8 s |

Every pair has the same iteration count, the same objective to every printed
digit and the same final equality residual. The largest control difference is
`1.7e-12` against controls of magnitude ~1,000 at large10k, and states agree to
`4e-16`. The saving is what the solve's share predicts: the wide solve becomes
2.2-2.5x faster inside the run, and it was 38% (large10k `T=3`) and 46%
(large10k `T=6`) of the baseline run.

**Trend with horizon.** At med2522 the saving rises from 13% (`T=3`) to 20-21%
(`T=12`, `T=24`) and levels off, because the solve's share does too (45% at
`T=24`, and 44-45% at `T=48`/`96` in the race logs). At large10k it is still
rising (23% -> 27%), and the solve's share keeps growing with horizon (58% at
`T=24` in the race log), so roughly a third is expected at `T=24`. That is **not
yet measured**.

**Absolute times vs Table V.** These pairs run FilterDDP without the three
exact assembly rewrites (`FILTERDDP_DIRECT_DIAG_HESSIAN`,
`FILTERDDP_TRIPLET_SECOND_DERIVATIVES`, `FILTERDDP_CACHE_KKT_PATTERN`), which
Table V's large10k rows use. Both arms of every pair share the same settings,
so the savings are fair, but the large10k baselines here (e.g. 3230.5 s at
`T=6`) are slower than Table V (2330.7 s at `T=6`, same 103 iterations and
objective). The med2522 rows of Table V did not use the rewrites and agree with
the baselines here (613.8 s vs 587.6 s at `T=24`). (The paper numbers that
table "Table II"; this note calls it Table V after the agenda.)

Re-running large10k `T=6` **with** the three rewrites (both arms,
`RUN_TAG_SUFFIX=_rewrites`) gives 2991.9 -> 2103.9 s, a 30% saving: the
rewrites shrink the rest of each iteration, so the solve's share and the
saving grow. The rewrite baseline is still 28% slower than Table II's
2330.7 s.

**Threading, verified with a controlled pair.** The matched race and the
rewrite queue set no thread variables, so Table II's FilterDDP ran with
Julia's default of **10 BLAS threads**, while every pair above pins BLAS and
OpenMP to one. Re-running large10k `T=6` with the rewrites and default threads
(`PIN_BLAS_THREADS=0`, `RUN_TAG_SUFFIX=_rewrites_blasdefault`; the log records
`blas_threads=10`):

| large10k T=6, rewrites | baseline | blocked | saving |
|---|---:|---:|---:|
| 1 BLAS thread | 2991.9 s | 2103.9 s | 30% |
| 10 BLAS threads (Table II's setting) | 2545.6 s | 1782.9 s | 30% |

All four runs take 103 iterations to the same objective. The threads cut the
dense value-function update ("update" time 315 -> 168 s summed over the run)
and leave the UMFPACK solve nearly unchanged (1503 -> 1433 s), so pinning to
one thread makes FilterDDP 17% slower. The blocked solve does not use BLAS, so
its 30% saving holds either way. The remaining 9% gap to Table II's 2330.7 s at
identical settings is variation between sessions five days apart. The paper
now states this in its protocol paragraph: Ipopt's MUMPS build is sequential,
FilterDDP's Table II runs used ten BLAS threads.

A side observation, not investigated: in the blocked runs the summed
*factorization* time is often higher (e.g. 131 -> 188 s in the rewrite pair)
although the factorization code is unchanged. The likely cause is garbage
collection from the per-stage work arrays the blocked solve allocates, landing
inside the factorization timer. Preallocating that workspace may recover a few
percent.

What this changes in the answers above: Q2 and Q7 stand for *factorization*,
but the multi-RHS solve is a real lever. It can be pulled on the existing
UMFPACK factors, so a new solver is not needed to get it.

## Reproduce

```bash
bash ddp/examples/power_system/capture_ordering_kkt.sh 20          # ~20 min, captures (gitignored)
bash ddp/examples/power_system/run_kkt_ordering_benchmark.sh       # ~55 min, benchmark CSV
python ddp/examples/power_system/summarize_kkt_ordering.py > ddp/results/kkt_ordering/SUMMARY_TABLES.md
bash ddp/examples/power_system/run_umfpack_strategy_fullrun.sh ieee123C_1ph 3 2
julia --project=envs/ddp2026 ddp/examples/power_system/ma57_availability_probe.jl
bash ddp/examples/power_system/run_blocked_solve_benchmark.sh      # ~20 min, isolated blocked solve
bash ddp/examples/power_system/run_blocked_solve_fullrun.sh large10kC_1ph 3 diag 16 1
```

Figures: `kkt_ordering_patterns.jl` then `plot_kkt_ordering_patterns.py`
(usage in each file's header).

A note on the earlier comparison, `ddp/results/agenda_pipeline/C_solvers_*.csv`:
several of its MUMPS rows have relative residuals near `1e+22`, i.e. failed
solves recorded as timings. This benchmark checks `INFOG(1)` and residuals on
every solve and had no failures. The earlier conclusion (UMFPACK faster than
MUMPS) stands on valid solves.

## 7. Table II re-run with the blocked solve (2026-09-25)

Every FilterDDP cell of the paper's Table II (`tab:near_opt`) was re-run with
`FILTERDDP_BLOCKED_SOLVE=16` in Table II's configuration: diagonal Hessian,
the three exact assembly rewrites, factor-backed policy (as the race used),
Julia's default ten BLAS threads, per-system primal thresholds
(`run_table2_blocked.sh`; logs `fullrun_blocked/*_tableII_*`). Only the blocked
arm was run; the Ipopt column is unchanged.

| System | T | iterations | earlier FilterDDP (s) | blocked (s) | ratio to Ipopt |
|---|---:|---:|---:|---:|---:|
| ieee123 | 6 | 67 | 21.0 | 21.6 | 57.2x |
| ieee123 | 24 | 85 | 32.2 | 34.2 | 25.5x |
| ieee123 | 96 | 120 | 97.7 | 97.0 | 4.3x |
| med2522 | 6 | 69 | 161.9 | 136.7 | 19.4x |
| med2522 | 24 | 79 | 587.6 | 471.7 | 11.5x |
| med2522 | 96 | 94 | 2693.2 | 2139.9 | 10.0x |
| large10k | 6 | 103 | 2330.7 | 1782.9 | 36.8x |
| large10k | 24 | 96 | 7664.9 | 5332.9 | 21.6x |
| large10k | 48 | 81 | 15343.2 (see below) | 8804.5 | 15.8x |

Every run reaches near-optimality at the same iteration as its predecessor.
ieee123 is unchanged within noise; med2522 and large10k fall 16-30%
(large10k `T=48` 43% against a pre-rewrite baseline). The gap to Ipopt now
narrows with horizon on all three feeders.

**Threshold inconsistency found and fixed.** The earlier Table II entry for
large10k `T=48` (iteration 109, 20499.5 s) was the run's time at primal
`1e-6`, whereas large10k is reported at `1e-4`
(`matched_ipopt_race/NEAR_OPT_SUMMARY.txt` lists 81 / 15343.2 s at `1e-4`).
The blocked run qualifies at the same iteration 81 at `1e-4`. The paper's
Table II now reports iteration 81 and says so in the text.

## 8. Blocked solve on eight threads (2026-09-25)

The blocked solve now splits its column blocks across Julia threads
(`JULIA_NUM_THREADS`; `blocked_solve.jl`). Blocks write disjoint columns and
only read the factors, and the result is bit-identical to one thread (ieee123
`T=3`). Three Table II cells were re-run with eight threads, otherwise in Table
II's configuration (`RUN_TAG_SUFFIX=_tableII_jt8`, logs
`fullrun_blocked/*_tableII_jt8_*`, which record `julia_threads=8`):

| System | T | iterations | blocked, 1 thread (s) | blocked, 8 threads (s) | saving | ratio to Ipopt |
|---|---:|---:|---:|---:|---:|---:|
| med2522 | 24 | 79 | 471.7 | 381.0 | 19% | 9.2x (was 11.5x) |
| med2522 | 96 | 94 | 2139.9 | 1700.2 | 21% | 7.9x (was 10.0x) |
| large10k | 6 | 103 | 1782.9 | 1353.0 | 24% | 28.0x (was 36.8x) |

Every pair reaches near-optimality at the same iteration with the same
objective to all printed digits. Background load stayed below 1.2 other cores
throughout. Against the pre-blocked Table II times, the two changes together
take 35-42% off (e.g. large10k `T=6` 2330.7 -> 1353.0 s).

**Why eight.** It was not tuned. The machine has 10 physical cores (20
logical); eight left two for the OS, the load sampler and other work. On a
single stage solve at `w=16` it gives 2.3x (med2522) and 2.7x (large10k) over
one thread, far from linear. The thread count is also capped by the number of
16-column blocks: 4 at ieee123 (52 columns), 16 at med2522, 64 at large10k. A
2/4/6/10-thread sweep on the captured stage KKTs is added to
`run_blocked_solve_benchmark.sh` to settle the choice.

Not yet run at eight threads: ieee123 (all three horizons), med2522 `T=6`,
large10k `T=24` and `T=48`. Table II keeps the one-thread blocked times until
those are done, so that all its cells share one setting.

## 9. MA57 and HSL_MA97, measured (2026-09-25)

**Build.** From the user's licensed HSL packages (MA57 3.11.3, HSL_MA97
2.8.1), compiled on the 309 PC with `hsl/build_hsl_windows.sh` (MinGW
gfortran 8.3, `-O3 -march=haswell`, LP64 OpenBLAS32 from Julia's MUMPS
artifact). METIS is replaced by HSL's placeholder, so the automatic orderings
fall back to AMD. Sources and DLLs live in `Documents/hsl/`, outside the repo.
MA57 is called directly (`ma57id/ad/bd/cd`), MA97 through the C shim
`hsl/s97.c`. Every configuration factors the same stage-1 diagonal-Hessian
captures as Sections 3 and 6 and solves the same `(n_x+1)`-column RHS in one
call (`hsl_kkt_benchmark.jl`, `run_hsl_benchmark.sh`; CSV
`hsl_kkt_benchmark.csv`, per-case files in `hsl/`). All residuals are at or
below `6e-12`. Every MA57/MA97 factorization reports the same inertia (e.g.
42,303 negative pivots at large10k) and full rank.

**large10k (n = 96,968, 1,021 RHS), one thread unless marked:**

| Configuration | order (s) | factor (s) | 1 RHS (ms) | 1,021 RHS (s) | factor + wide (s) | stored entries | delayed pivots |
|---|---:|---:|---:|---:|---:|---:|---:|
| UMFPACK default (FilterDDP) | 0.069 | 0.146 | 2.5 | 2.10 | 2.24 | 870,685 (L+U) | -- |
| UMFPACK + blocked w=16 | 0.069 | 0.146 | 2.5 | 0.99 | 1.13 | same | -- |
| UMFPACK + blocked w=16, 8 threads | 0.056 | 0.134 | 2.5 | 0.32 | 0.45 | same | -- |
| MA57 defaults (AMD) | 0.015 | 0.051 | 4.4 | 14.43 | 14.48 | 2,565,910 | 22,832 |
| MA57 MA27 min. degree | 0.015 | 0.051 | 4.6 | 14.00 | 14.05 | 2,532,184 | 22,483 |
| MA57, Ipopt's settings (pivtol 1e-8, no scaling) | 0.016 | 0.022 | 4.3 | 13.38 | 13.40 | 2,414,760 | 19,159 |
| MA97 defaults (AMD) | 0.034 | 0.052 | 9.6 | 9.66 | 9.72 | 1,781,056 | 16,564 |
| MA97, no scaling | 0.036 | 0.025 | 9.5 | 9.29 | 9.31 | 1,804,234 | 17,262 |
| MA97 defaults, 8 OpenMP threads | 0.038 | 0.047 | 29.5 | 5.92 | 5.97 | same | same |
| MA97 no scaling, 8 threads | 0.033 | 0.025 | 25.1 | 5.47 | 5.50 | same | same |

**Reading.**

- **Factorization is where MA57 wins**: 0.022-0.051 s against UMFPACK's
  0.146 s, the fastest of every solver measured. But with the diagonal
  Hessian factorization is 7% of the stage's linear algebra (Q7), so this
  cannot pay for itself.
- **The solve is where it loses.** A pattern-only symmetric ordering (AMD)
  ignores that the constraint block of `K` has a zero diagonal. 17-23k of the
  96,968 pivots are delayed, and the factors store 2.0x (MA97) to 2.9x (MA57)
  as many entries as UMFPACK's L and U together, even though LDLᵀ stores only
  one triangle. UMFPACK's unsymmetric COLAMD ordering pivots numerically
  during elimination and gets a much sparser factor. Every RHS column
  traverses the whole factor, so the wide solve is 4.4-6.9x slower than
  UMFPACK column by column, and 9-15x slower than the blocked solve.
- **MA57's multi-RHS call is slower per column than its single-RHS call**
  (13-14 ms against 4.3-4.6 ms). Its best case, calling it one column at a
  time, extrapolates to about 4.4 s at large10k. That is still 2x UMFPACK's
  column-by-column solve and 4.4x the blocked one. This is an extrapolation,
  not a measurement.
- **MA97's OpenMP helps (9.3 -> 5.5 s) but not enough**: the 8-thread blocked
  UMFPACK solve takes 0.32 s, 17x less.
- Tuning barely matters: the ordering choice moves either solver's wide
  solve by under 4%. Ipopt's settings (smaller pivot threshold, fewer delays)
  are the best MA57 variant.

**med2522 and ieee123.** Same ranking at med2522. Factor plus wide solve takes:

| med2522 | factor + wide (s) |
|---|---:|
| UMFPACK | 0.146 |
| UMFPACK blocked | 0.085 |
| UMFPACK blocked, 8 threads | 0.052 |
| MA57 | 0.29-0.36 |
| MA97 | 0.34-0.35 |
| MA97, 8 threads | 0.31-0.33 |

At ieee123 everything is 2.0-3.3 ms per stage. MA57 with Ipopt's settings is
the fastest there (2.0 ms against 2.3 ms blocked), which is immaterial in
absolute terms.

**Caveat.** Without METIS, the nested-dissection orderings of MA57/MA97 are
untested. Section 3 found nested dissection never reduced fill on these
radial-feeder KKTs (METIS and SCOTCH through UMFPACK and MUMPS), so it is
unlikely to close a 4-7x solve gap. It is not measured here.

**Answer for the agenda.** MA57 is not worth adopting for FilterDDP's stage
solves. It would help a method dominated by factorization. FilterDDP with the
diagonal Hessian is dominated by the `(n_x+1)`-column solve, and there UMFPACK's
sparser unsymmetric factor, applied in blocks, is the fastest option measured.

Reproduce (after building the DLLs):

```bash
bash ddp/examples/power_system/hsl/build_hsl_windows.sh <hsl_src_dir> <build_dir> <libopenblas.dll>
HSL_LIB_DIR=<build_dir> bash ddp/examples/power_system/run_hsl_benchmark.sh
```
