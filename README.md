# MultiPeriodDistOPF

Scalable **multi-period optimal power flow** (MPOPF) for distribution feeders with
battery storage and solar PV. The headline method is **temporal ADMM (tADMM)**, which
decomposes the multi-period problem along time into small per-period branch-flow OPFs and
coordinates them through an augmented-Lagrangian consensus — solving horizons that the
monolithic formulation cannot. A monolithic Branch-Flow (BF) solver is included as the
baseline.

> Companion code for the IAS-Trans paper *"Scalable Multi-Period Optimal Power Flow in
> Distribution Systems via Temporal Decomposition."*

## Install

Requires **Julia ≥ 1.10**. The project environment lives in `envs/tadmm/`:

```bash
git clone https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
```

Key dependencies: JuMP, Ipopt (default solver), OpenDSSDirect, Plots. Gurobi is optional
(needs a license). See [`docs/first-time-installation-guide.md`](docs/first-time-installation-guide.md)
for a step-by-step VS Code walkthrough.

> **Tip for a fast clone:** this repo has a long history. To grab just the latest snapshot,
> `git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git`
> (or `--filter=blob:none` for a blobless partial clone).

## Quickstart

The run scripts auto-activate `envs/tadmm`, so after instantiating you can just:

```bash
julia run_bf.jl                       # monolithic Branch-Flow baseline
JULIA_NUM_THREADS=16 julia run_tadmm.jl   # temporal ADMM (parallel over periods)
```

Configure via [`config.jl`](config.jl), or override per run with environment variables:

```bash
# system + horizon
SYSTEM_OVERRIDE=ieee123C_1ph T_OVERRIDE=24 julia run_tadmm.jl
# tADMM penalty / tolerances
RHO_OVERRIDE=15000 EPS_PRI_OVERRIDE=1e-4 julia run_tadmm.jl
```

Outputs are written to `envs/tadmm/processedData/<system>_T<T>/` (gitignored — regenerated
on each run). The curated final-winner results behind the paper live in [`results/`](results/).

## Test systems

Balanced single-phase OpenDSS feeders in [`rawData/`](rawData/):

| Name | Buses | Notes |
|------|-------|-------|
| `ieee123C_1ph`  | 128    | IEEE 123-node, small (solution-quality reference) |
| `ieee2552C_1ph` | 2,552  | synthetic medium feeder |
| `large10kC_1ph` | 10,321 | synthetic large feeder (primary scalability benchmark) |

## Repository layout

```text
config.jl, run_bf.jl, run_tadmm.jl   core entry points (edit config, then run)
run_rho_sweep.jl                     penalty (rho0) tuning sweep
tadmm_socp.jl                        single-file interactive runner (VS Code)
envs/tadmm/                          Julia project + parser/validators/logger/plotter
rawData/                             OpenDSS feeder models
results/                            curated final-winner results (summary.csv + trajectories)
plots/                               figure-generation scripts
scripts/                             helpers; scripts/probes/ = one-off experiment scripts
docs/                                installation guide
```

## Reproducing paper results

Winning `rho0` per `(system, T)` is listed in [`results/summary.csv`](results/summary.csv).
To regenerate a cell's full penalty sweep:

```bash
SYSTEM_OVERRIDE=ieee2552C_1ph T_OVERRIDE=144 julia run_rho_sweep.jl
```

## Notes

- **Windows "filename too long" error?** Run, in an elevated Git Bash:
  `git config --global core.longpaths true`
  ([reference](https://stackoverflow.com/questions/22575662/filename-too-long-in-git-for-windows)).
