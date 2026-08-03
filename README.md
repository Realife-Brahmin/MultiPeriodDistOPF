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
julia envs/tadmm/root_level/run_bf.jl                       # monolithic Branch-Flow baseline
JULIA_NUM_THREADS=16 julia envs/tadmm/root_level/run_tadmm.jl   # temporal ADMM (parallel over periods)
```

> **Where the run scripts live.** They sit in
> [`envs/tadmm/root_level/`](envs/tadmm/root_level/) while tADMM is parked and other work is
> active at the repo root — the same convention as `envs/ddp/root_level/` and
> `envs/multi_poi/root_level/`. They resolve their own paths and run correctly from either
> `envs/tadmm/root_level/` or the repo root, so moving them back up needs no code change.
> The state of the repo behind the IAS-Trans paper, with these scripts at the root, is
> archived in the [`260731-ias-trans-tadmm`](../../releases/tag/260731-ias-trans-tadmm) release.

Configure via [`config.jl`](envs/tadmm/root_level/config.jl), or override per run with
environment variables:

```bash
# system + horizon
SYSTEM_OVERRIDE=ieee123C_1ph T_OVERRIDE=24 julia envs/tadmm/root_level/run_tadmm.jl
# tADMM penalty / tolerances
RHO_OVERRIDE=15000 EPS_PRI_OVERRIDE=1e-4 julia envs/tadmm/root_level/run_tadmm.jl
```

Outputs are written to `envs/tadmm/processedData/<system>_T<T>/` (gitignored — regenerated
on each run). The curated final-winner results behind the paper live in [`results/`](results/).

## Test systems

Balanced single-phase OpenDSS feeders in [`rawData/`](rawData/):

| Name | Buses | Notes |
|------|-------|-------|
| `ieee123C_1ph`  | 128    | IEEE 123-node, small (solution-quality reference) |
| `ieee2552C_1ph` | 2,522  | synthetic medium feeder |
| `large10kC_1ph` | 10,321 | synthetic large feeder (primary scalability benchmark) |

## Repository layout

```text
envs/tadmm/                          Julia project + parser/validators/logger/plotter
  root_level/config.jl               configuration (edit this, then run)
  root_level/run_bf.jl               monolithic Branch-Flow baseline
  root_level/run_tadmm.jl            temporal ADMM entry point
  root_level/run_rho_sweep.jl        penalty (rho0) tuning sweep
  root_level/tadmm_socp.jl           single-file interactive runner (VS Code)
envs/ddp/, envs/multi_poi/           other methods in this repo (not used by the paper)
rawData/                             OpenDSS feeder models
results/                            curated final-winner results (summary.csv + trajectories)
plots/                               figure-generation scripts
scripts/                             helpers; scripts/probes/ = one-off experiment scripts
docs/                                installation guide
```

## Reproducing paper results

Every reported `(system, T)` cell — winning penalty `rho0`, BF time, tADMM near-optimal
time, speedup — is in [`results/summary.csv`](results/summary.csv);
[`results/README.md`](results/README.md) explains the layout and how the trajectory
figures are produced.

**Paper name → system ID.** `summary.csv` and the paper use the short names; the code
takes the feeder directory name:

| Paper | `SYSTEM_OVERRIDE` |
|-------|-------------------|
| ieee123  | `ieee123C_1ph`  |
| med2522  | `ieee2552C_1ph` |
| large10k | `large10kC_1ph` |

**Start here.** The cheapest cell that shows a real speedup takes well under a minute:

```bash
SYSTEM_OVERRIDE=ieee2552C_1ph T_OVERRIDE=6 RHO_OVERRIDE=4000 \
  JULIA_NUM_THREADS=16 julia envs/tadmm/root_level/run_tadmm.jl     # ~30 s, expect ~1.27x over BF
```

Reproduce any other cell by taking its `rho0_winner` from `summary.csv` and passing it
as `RHO_OVERRIDE`, with `T_OVERRIDE` set to that row's `T`.

**Know the cost before you start.** Times below are from a 16-thread run; tADMM parallelises
over periods, so fewer threads will be proportionally slower.

| Cell | tADMM | BF baseline | Note |
|------|-------|-------------|------|
| `ieee123C_1ph`, any `T` | seconds | seconds | solution-quality check; tADMM is *slower* here by design |
| `ieee2552C_1ph`, `T=6..48` | 0.5–9 min | 0.6–12 min | good middle ground |
| `ieee2552C_1ph`, `T=144` | ~23 min | ~3.5 h | the 9.31x headline |
| `large10kC_1ph`, `T=48` | ~1.6 h | **fails after ~4.9 h** | see below |

**The `large10k, T=48` BF run is *supposed* to fail.** That is the paper's headline
result: the monolithic solve exhausts memory during MUMPS factorisation (~10.5 GB) and
never converges, while tADMM completes. If you run `run_bf.jl` on that cell expecting a
solution, you have not misconfigured anything — the reported "BF time" is the wall-clock
at failure. Budget ~16 GB RAM to observe it.

To regenerate a cell's full penalty sweep rather than a single run:

```bash
SYSTEM_OVERRIDE=ieee2552C_1ph T_OVERRIDE=144 julia envs/tadmm/root_level/run_rho_sweep.jl
```

## Notes

- **Windows "filename too long" error?** Run, in an elevated Git Bash:
  `git config --global core.longpaths true`
  ([reference](https://stackoverflow.com/questions/22575662/filename-too-long-in-git-for-windows)).
