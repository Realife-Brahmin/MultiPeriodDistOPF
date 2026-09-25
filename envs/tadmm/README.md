# tADMM: temporal ADMM for multi-period OPF

**Temporal ADMM (tADMM)** splits the multi-period OPF along time into small
per-period branch-flow OPFs, and coordinates them through an augmented-Lagrangian
consensus on the battery states. This solves horizons that the monolithic
formulation cannot. A monolithic branch-flow (BF) solver is included as the
baseline.

> Companion code for the IAS-Trans paper *"Scalable Multi-Period Optimal Power Flow
> in Distribution Systems via Temporal Decomposition."* The repository state behind
> the paper is archived in the
> [`260731-ias-trans-tadmm`](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/releases/tag/260731-ias-trans-tadmm)
> release.

## Install

From the repository root:

```bash
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
```

Key dependencies: JuMP, Ipopt (default solver), OpenDSSDirect, Plots. Gurobi is
optional and needs a license.

## Quickstart

The run scripts activate `envs/tadmm` themselves, so after instantiating:

```bash
julia envs/tadmm/root_level/run_bf.jl                           # monolithic BF baseline
JULIA_NUM_THREADS=16 julia envs/tadmm/root_level/run_tadmm.jl   # tADMM, parallel over periods
```

Configure through [`root_level/config.jl`](root_level/config.jl), or override
per run with environment variables:

```bash
# system and horizon
SYSTEM_OVERRIDE=ieee123C_1ph T_OVERRIDE=24 julia envs/tadmm/root_level/run_tadmm.jl
# tADMM penalty and tolerance
RHO_OVERRIDE=15000 EPS_PRI_OVERRIDE=1e-4 julia envs/tadmm/root_level/run_tadmm.jl
```

Outputs go to `envs/tadmm/processedData/<system>_T<T>/`, which is gitignored and
regenerated on each run. The curated results behind the paper are in
[`results/`](../../results/).

> **Why the scripts are in `root_level/`.** tADMM is parked while other work is
> active at the repository root, so its entry points live in
> [`root_level/`](root_level/). They resolve their own paths and run the same from
> there or from the repository root.

## Layout

```text
envs/tadmm/
  root_level/config.jl          configuration (edit this, then run)
  root_level/run_bf.jl          monolithic branch-flow baseline
  root_level/run_tadmm.jl       temporal ADMM entry point
  root_level/run_rho_sweep.jl   penalty (rho0) tuning sweep
  root_level/tadmm_socp.jl      single-file interactive runner (VS Code)
  parse_opendss.jl              OpenDSS feeder parser
  opendss_validator.jl          checks a solution against an OpenDSS power flow
  solution_validator.jl         constraint-violation checks
  logger.jl, Plotter.jl         run logs and figures
  tex/tadmm_formulation.pdf     the formulation, written out
```

## Reproducing the paper's results

Every reported `(system, T)` cell (winning penalty `rho0`, BF time, tADMM
near-optimal time, speedup) is in [`results/summary.csv`](../../results/summary.csv).
[`results/README.md`](../../results/README.md) explains the layout and how the
trajectory figures are produced.

**Paper name to system ID.** `summary.csv` and the paper use short names; the code
takes the feeder directory name:

| Paper | `SYSTEM_OVERRIDE` |
|-------|-------------------|
| ieee123  | `ieee123C_1ph`  |
| med2522  | `ieee2522C_1ph` |
| large10k | `large10kC_1ph` |

**Start here.** The cheapest cell that shows a real speedup takes well under a
minute:

```bash
SYSTEM_OVERRIDE=ieee2522C_1ph T_OVERRIDE=6 RHO_OVERRIDE=4000 \
  JULIA_NUM_THREADS=16 julia envs/tadmm/root_level/run_tadmm.jl     # ~30 s, expect ~1.27x over BF
```

Reproduce any other cell by passing its `rho0_winner` from `summary.csv` as
`RHO_OVERRIDE`, with `T_OVERRIDE` set to that row's `T`.

**Know the cost before you start.** Times are from a 16-thread run. tADMM
parallelises over periods, so fewer threads are proportionally slower.

| Cell | tADMM | BF baseline | Note |
|------|-------|-------------|------|
| `ieee123C_1ph`, any `T` | seconds | seconds | solution-quality check; tADMM is *slower* here by design |
| `ieee2522C_1ph`, `T=6..48` | 0.5–9 min | 0.6–12 min | good middle ground |
| `ieee2522C_1ph`, `T=144` | ~23 min | ~3.5 h | the 9.31x headline |
| `large10kC_1ph`, `T=48` | ~1.6 h | **fails after ~4.9 h** | see below |

**The `large10k, T=48` BF run is *supposed* to fail.** That is the paper's
headline result: the monolithic solve exhausts memory during MUMPS factorisation
(~10.5 GB) and never converges, while tADMM completes. If `run_bf.jl` fails on
that cell, nothing is misconfigured; the reported "BF time" is the wall-clock time
at failure. Budget ~16 GB of RAM to observe it.

To regenerate a cell's full penalty sweep rather than a single run:

```bash
SYSTEM_OVERRIDE=ieee2522C_1ph T_OVERRIDE=144 julia envs/tadmm/root_level/run_rho_sweep.jl
```
