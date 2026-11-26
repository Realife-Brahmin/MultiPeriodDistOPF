# KKT matrix updates: magnitudes, stale Jacobians, block-diagonal Jacobians, HSL in Ipopt

Agenda for the meeting with R. Gupta on 2026-10-02. Branch
`ddp-kkt-matrix-updates-sep27`. Every experiment is on the paper's current
formulation unless stated:
- `T = 6`, periodic profile
- per-system `C_B`
- soft terminal SOC
- diagonal stage Hessian with the exact assembly rewrites

The companion study of how the matrix changes between iterations is
`KKT_EVOLUTION.md`.

## 2a. Are many KKT entries tiny?

**Setup.**
- **Captures.** Stage-1 KKT systems of all three systems at iteration 5, at
  iteration 40, and near the Table II near-optimality iteration (ieee123 67,
  med2522 68, large10k 100). Script: `capture_kkt_magnitude.sh`; the captures
  are gitignored.
- **Analysis** (`kkt_magnitude_analysis.jl`, driven by
  `run_kkt_magnitude.sh`):
  - The magnitude of every stored entry, absolute and after symmetric scaling
    `r_ij = |K_ij| / sqrt(d_i d_j)`, `d_i = max_k |K_ik|`.
  - For thresholds `1e-12 .. 1e-2`, every off-diagonal entry below the
    threshold is dropped (the diagonal is always kept). The thresholded matrix
    is then factored and solved, and the analysis reports:
    - the entries dropped
    - the factor fill
    - how far the `(n_x+1)`-column solution moves
- **Outputs.**
  - `ddp/results/kkt_magnitude/magnitude_threshold.csv`
  - `magnitude_blocks.csv`
  - `figures/kkt_magnitude_<system>_iter<k>.png`, which shows the full matrix,
    the matrix thresholded at `1e-4` and `1e-2`, the entries dropped at `1e-2`,
    and the magnitude histograms.

**What changes is the "primal" solution.** That is the control rows, from which
FilterDDP takes its feedforward and feedback gains. The multiplier rows are
far more sensitive and are reported in the CSV.

| threshold (absolute) | entries dropped | factor size | change in the control directions |
|---:|---|---|---|
| `1e-8` | 0.0-0.8% | -0.1% to -0.6% | <= 9e-5 (mostly <= 1e-6) |
| `1e-6` | 1.7-5.5% | -1.5% to -5% | ieee123: up to **13%** (iteration 40); med2522 <= 0.07%; large10k <= 0.04% |
| `1e-4` | 3.9-16% | -2% to -20% | ieee123 13%, med2522 up to **41%** (iteration 40), large10k <= 0.2% |
| `1e-2` | 16-37% | -19% to -55% | 2-43% everywhere |

**Where the small entries are.** They are the line-impedance coefficients:
- in the voltage-drop rows, `2r·P`, `2x·Q` and `(r²+x²)·ℓ`. At large10k 40% of
  those entries are below `1e-4` and 60% below `1e-2`.
- the `r·ℓ` and `x·ℓ` loss terms of the power-balance rows.
- a few Hessian diagonal entries sitting at the `1e-8` floor or with a tiny
  barrier term.

The SOCP and battery-energy rows are `O(1)`. Line impedances are small in
per-unit; that does not make them unimportant.

**Scaling is no guide.** Dropping entries that are small *relative to their
row and column scale* makes the matrix **singular** in late iterations, even
at a threshold of `1e-8` (ieee123 and med2522). Some constraint rows are tied
to their variables only through entries that look negligible next to a
barrier-dominated diagonal of up to `1e10`. Those entries are structurally
essential.

**Answer.**
- **Harmless:** only entries below about `1e-8` in absolute value can be
  dropped harmlessly, and that is under 1% of the matrix, so the sparsity plot
  does not change.
- **Useful but instance-dependent:** thresholds that visibly thin the matrix
  (`1e-4` .. `1e-2`) change the Newton directions by 13-43% on ieee123 and
  med2522.
- **large10k is the exception.** There `|K_ij| < 1e-3` removes 26% of the
  entries and 45% of the factor while moving the directions by only 2%. That
  would make an approximation worth a full run, but it cannot be assumed on
  other feeders.

## 3. Can the KKT refresh be pinpointed? (stale Jacobian)

**Motivation.** The evolution study showed which parts of the matrix move. The
linear constraint rows never change, so they never need refreshing. Only the
barrier diagonal and the SOCP rows of the Jacobian move. The barrier diagonal
moves by orders of magnitude and must be refreshed. So the pinpointed
question is whether the SOCP Jacobian rows can be refreshed less often.

**Mechanism.** `FILTERDDP_STALE_JACOBIAN_PERIOD = p` (`backward_pass.jl`)
builds the KKT matrix from each stage's Jacobian of the last iteration that
was a multiple of `p`. Everything else stays current: Hessian, barrier terms,
gradient `cu'φ`, residuals. The result is an inexact-Newton step.

**Runs.** Table II configuration, per-system `C_B`, `run_stale_jacobian.sh`,
logs `fullrun_blocked/*_stalejac<p>_cbsys_*`. `p = 1` is the Table II run.

| system T=6 | p = 1 (Table II) | p = 2 | p = 5 | p = 10 |
|---|---|---|---|---|
| ieee123 | 67 it / 20.7 s | 83 it / 23.1 s | 97 it / 23.3 s | 113 it / 24.5 s |
| med2522 | 68 it / 102.0 s | 66 it / 101.6 s | **fails** at iteration 4 (line search) | **fails** at iteration 4 |

**Answer.**
- **Short lags are tolerated.** A Jacobian one iteration old costs nothing on
  med2522 and 24% more iterations on ieee123.
- **Longer lags hurt.** On ieee123 they cost 45-69% more iterations. On
  med2522 they fail at once: in the first iterations the iterates move far
  and the stale linearization steers the step out of the filter's
  acceptance.
- **Nothing is gained either way.** The barrier diagonal changes every
  iteration, so the matrix must be refactored every iteration regardless, and
  the Jacobian evaluation that a stale Jacobian saves is a small share of a
  backward pass. So there is no time saving to trade the extra iterations
  against: every stale run above is slower than the fresh one or equal to it.
- **Refresh can be pinpointed to the SOCP rows plus the diagonal**, which is
  what DDP4OPF already does in effect: the kkt-pattern cache re-uses the
  structure and rewrites only the values. Skipping even that is not worth it.

## 2b. Can the Jacobian be block-diagonal?

**Why not "diagonal" like the Hessian.** The Jacobian is rectangular: rows are
constraints, columns are variables. Every constraint (a power balance, a
voltage drop) necessarily involves several variables, so there is no "own"
entry per row to keep.

There is also a more basic difference:
- **Approximating the Hessian** changes only the step's curvature. The
  constraints stay exactly linearized, so every step still heads for
  feasibility.
- **Approximating the Jacobian** changes the linearized constraints
  themselves, so feasibility progress is no longer guaranteed.

The meaningful analogue is block-diagonal. Two readings are possible:
- **By area**, tested here.
- **By quantity**, as in fast-decoupled load flow: active vs reactive. Not
  tested yet; to confirm with R. Gupta which he meant.

**Test.** `kkt_block_jacobian.jl` and `run_kkt_block_jacobian.sh`, on the
stage-1 captures at iteration 40 and near-optimality.
- **Areas.** The radial feeder is cut into k subtree areas. Every variable
  and constraint row is assigned to an area, and only the Jacobian entries
  linking two areas are dropped. These sit at the cut lines: the parent bus's
  balance-row entry for the outgoing flow, and the cut line's voltage-drop and
  SOCP entries for the upstream voltage.
- **Effect.** With the diagonal Hessian the KKT matrix then separates into k
  independent blocks.
- **Measured.** How far the Newton step moves, and whether block-Jacobi
  iterative refinement (`x += M \ (b - K x)`) recovers the exact step, i.e.
  whether the blocks at least work as a preconditioner.
- **large10k areas.** The greedy cut groups nothing at the root, so at
  large10k it only produces a real split at its natural ~100 areas
  (`k = 128` gives 103). Smaller k leave one area, and those rows are empty.

Change in the control directions / refinement steps to a `1e-10` residual:

| system | areas | entries dropped | iteration 40 | near-optimality |
|---|---:|---:|---|---|
| ieee123 | 2 | 8 (0.15%) | 44% / 6 steps | 69% / exact already |
| ieee123 | 29 | 224 (4.3%) | 61% / 56 steps | 95% / exact already |
| med2522 | 2 | 8 (0.008%) | 81% / 27 steps | 80% / 1 step |
| med2522 | 7 | 48 | 95% / 70 steps | 97% / **diverges** |
| med2522 | 14-58 | 104-456 | **diverges** | **diverges** |
| large10k | 103 (natural areas) | 816 (0.2%) | 90% / **diverges** | 64% / 49 steps |

("exact already" means the first block solve already met the residual
tolerance: the dropped entries hardly affected the feedforward column. The
feedback-gain columns, which carry most of the change above, still move.)

**Answer: no.**
- **As a replacement.** Dropping as few as 8 of 96,025 entries on med2522
  changes the feedback gains by 80%. Those gains are what FilterDDP
  propagates backward, so the boundary couplings carry the whole
  intertemporal and spatial sensitivity. The feedforward direction alone is
  more robust (0.3-44%).
- **As a preconditioner.** Block-Jacobi recovers the exact step only with few
  areas, and it diverges for 8 or more areas on med2522 and for large10k's
  103 natural areas mid-solve. Even where it converges, each refinement step
  costs one block solve and one multiplication by K. The 6-70 steps would
  have to be paid back by solving the areas in parallel, and at large10k a
  single solve is already the dominant cost.

## 1. MA57 / MA97 in centralized Ipopt

**Setup.** On the FilterDDP side this was settled earlier: on the stage KKT,
MA57 and MA97 lose to UMFPACK's blocked solve (`KKT_ORDERING_AND_MA57.md`,
Sections 9-10). Here the same locally built HSL libraries are loaded into
centralized Ipopt through `hsllib` (`run_ipopt_hsl.sh`,
`IPOPT_EXTRA_OPTIONS` in `centralized_ipopt_matched.jl`). The runs cover the
nine Table II instances at the per-system `C_B`, with
`linear_system_scaling = none` for every solver (MC19 is not in the licensed
packages) and default BLAS threads, as in Table II. Logs are in
`ddp/results/ipopt_hsl/logs/`.

| system | T | MUMPS (s) | MA57 (s) | MA97, 1 thread (s) |
|---|---:|---:|---:|---:|
| ieee123 | 6 / 24 / 96 | 0.39 / 1.47 / 23.2 | 0.19 / 0.74 / 3.91 | 0.25 / 0.86 / 4.49 |
| med2522 | 6 / 24 / 96 | 7.2 / 41.1 / 220.0 | 4.45 / 24.4 / 127.3 | 4.69 / 25.7 / 140.6 |
| large10k | 6 / 24 / 48 | 51.2 / 283.0 / 623.7 | 27.6 / 144.7 / 440.7 | 28.5 / 132.4 / 272.0 |

**Findings.**
- **Same solutions.** Objectives agree to about `1e-15`, with the same
  iteration counts except ieee123 `T=96` (44 iterations with HSL, 97 with
  MUMPS).
- **HSL is faster.** MA57 is 1.4-5.8x faster than MUMPS, and MA97 on one
  thread is fastest at large10k `T=24` and `T=48`.
- **MA97 on 8 OpenMP threads is slower: 1.5-3.3x slower than on one thread.**
  It is *not* a clash with OpenBLAS's threads, as first suspected. A full
  re-run with OpenBLAS pinned to one thread (`BLAS1=1 run_ipopt_hsl.sh`, logs
  tagged `_blas1`) reproduces every time within about 10% at small sizes and
  2% at large10k, including the slow 8-thread MA97 (large10k `T=48`: 899.8 s
  pinned, 882.9 s unpinned).
  - The likely cause is that the fronts are tiny (at most about 40 rows), so
    OpenMP task overhead and synchronization outweigh the work.
  - This matches the stage-KKT benchmark only partly: there, 8 threads did
    speed up MA97's many-right-hand-side solve (9.7 -> 5.9 s). Ipopt
    factorizes every iteration and solves one right-hand side.

**Paper (user's decisions, 2026-09-27 and 2026-09-29).** Separate tables:
Table II keeps Ipopt with MUMPS, and Table III compares the same FilterDDP
runs with Ipopt using the faster of MA57 and MA97 (one thread) in each row.
MA97 is used at large10k `T=24` and `T=48`, MA57 everywhere else. There
FilterDDP is 14-110x slower (ieee123 110 / 43 / 24x, med2522 23 / 15 / 14x,
large10k 44 / 30 / 23x), against 4-64x with MUMPS.

### Longer horizons and threads (2026-09-29/30)

**Longer horizons.** MA57 against MA97 on one thread, per-system `C_B`,
`run_ipopt_hsl.sh` with `CELLS`/`SOLVERS`. Seconds; the last column is peak
resident memory (GiB) for MA57 / MA97.

| system | T | MA57 | MA97 | faster | memory |
|---|---:|---:|---:|---|---|
| ieee123 | 192 / 288 | 8.0 / 13.2 | 8.6 / 14.2 | MA57, ~7% | 1.0 / 1.0 at T=288 |
| med2522 | 144 / 192 / 288 | 212.7 / 276.6 / 446.1 | 239.8 / 322.7 / 502.4 | MA57, 11-17% | 3.8-6.9 / 4.1-7.3 |
| large10k | 96 | 1871 | 926 | MA97 2.0x | 10.6 / 9.6 |
| large10k | 144 | 5244 | 1739 | MA97 3.0x | 16.9 / 14.3 |
| large10k | 192 | 6860 | 1433 | MA97 4.8x | 19.3 / 19.4 |

- **Same solutions.** Every pair has identical iteration counts and
  objectives.
- **Why MA97 pulls ahead at large10k.** Its lead grows with the horizon
  (1.1x at T=24, 1.6x at 48, then 2.0 / 3.0 / 4.8x). The factorization is
  what differs: at T=192 MA57 spends 6367 s factorizing, against 788 s for
  MA97, while MA97's back-solves are slower (309 vs 196 s).
- **Probable explanation (not measured).** The fronts are set by the battery
  count, about 1,020 at large10k and 250 at med2522, and MA97's dense-block
  factorization pays off only on the wide ones. On the smaller feeders MA57's
  faster back-solves win at every horizon.

**Threads.** MA97 on 1, 4 and 8 OpenMP threads, with OpenBLAS pinned to one
thread:

| cell | 1 thread | 4 threads | 8 threads |
|---|---:|---:|---:|
| med2522 T=96 | 142.7 | 422.0 | 444.6 |
| large10k T=24 | 129.0 | 368.9 | 406.0 |
| large10k T=48 | 272.9 | 798.1 | 899.8 |
| large10k T=96 | 926.5 | 2829.1 | -- |

Four threads are already about 3x slower, nearly as slow as eight. So the
cost does not scale with the number of threads: entering MA97's parallel
code path is what costs.
- **Likely causes.**
  - Its task scheduling and synchronization over thousands of tiny tree
    nodes.
  - A long, thin elimination tree along the time chain.
  - Possibly the MinGW OpenMP runtime.
- **Consequence.** For Ipopt on these instances MA97 should run on one
  thread. MA57 has no threads of its own.

### med2522 to the memory limit (2026-09-30)

**Setup.** MA57 and MA97 on one thread (BLAS pinned), per-system `C_B`, run
by `run_ipopt_oom_sweep.sh`. A watchdog kills any Julia process above 27 GB
of private memory; the PC has 32 GB. Seconds / peak resident memory in GiB:

| T | MA57 | MA97 |
|---:|---|---|
| 384 | 589.8 / 9.0 | 691.2 / 9.3 |
| 576 | 917.6 / 12.9 | 1067.8 / 13.9 |
| 768 | 1263.7 / 15.9 | 1500.1 / 17.3 |
| 1152 | 2045.5 / 19.9 | **out of memory** (killed at 28.9 GB) |
| 1536 | **out of memory** (killed at 29.4 GB) | -- |

- **Same solutions.** Every successful pair reaches the same objective; MA97
  takes one extra iteration at T=768.
- **MA57 stays ahead.** It is 15-19% faster than MA97 at every horizon here,
  and its time grows about linearly: roughly 1.5-1.8 s per period.
- **MA57 goes further.** Its last solved horizon is T=1152, against MA97's
  T=768.
- **The failures are sudden jumps, not smooth growth.** MA57 peaked at
  19.9 GiB at T=1152, then passed 27 GB at T=1536. Both solvers grow their
  factor storage in large steps when their forecast falls short.
- **The limit is ours, not the machine's.** The 27 GB watchdog sits below
  the 32 GB of RAM, so both limits may be conservative by one step in T.
- **No like-for-like MUMPS comparison.** The earlier MUMPS escalation on
  med2522 (Table IV of the paper) ran out of memory between T=576 and 768,
  but on the original instance family (different `C_B` and price sampling),
  so it is not directly comparable.
