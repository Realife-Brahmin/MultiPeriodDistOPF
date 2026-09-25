# MultiPeriodDistOPF

Julia research code for **multi-period optimal power flow (MPOPF) on distribution
feeders** with batteries and solar PV. Batteries tie the periods together through
their state of charge, so the horizon cannot be solved one period at a time. This
repository holds three lines of work on that problem. Each has its own code, Julia
environment, branches and README.

| Workstream | What it is | Status | Start here |
|---|---|---|---|
| **DDP** | Differential dynamic programming with **DDP4OPF**, our fork of FilterDDP.jl, benchmarked against centralized Ipopt on identical instances | **Active, main focus** | [ddp/README.md](ddp/README.md) |
| **tADMM** | Temporal ADMM: splits the horizon into per-period OPFs coordinated by an augmented-Lagrangian consensus | Companion code for the IAS-Trans paper; parked | [envs/tadmm/README.md](envs/tadmm/README.md) |
| **MSOPF** | Multi-source OPF: feeders supplied by several substations | Active | [MSOPF section below](#msopf-multi-source-opf) |

All three read the same OpenDSS feeder models in [`rawData/`](rawData/).

## Getting started

Requires **Julia 1.10 or newer** (the DDP work runs on 1.12).

```bash
git clone https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

Then instantiate the environment for the workstream you need. The workstream
READMEs give the exact commands:

| Environment | Used for |
|---|---|
| `envs/ddp2026` | all DDP work: DDP4OPF, JuMP, Ipopt |
| `envs/tadmm` | tADMM, exporting feeder instances for DDP, and the MSOPF angle-sweep drivers |
| `envs/multi_poi` | the other MSOPF scripts |

Ipopt is the default solver everywhere. Gurobi is optional and needs a license.
For a step-by-step VS Code setup, see
[`docs/first-time-installation-guide.md`](docs/first-time-installation-guide.md).

> **Fast clone.** The history is long. `git clone --depth 1 <url>` fetches only the
> latest snapshot.
>
> **Windows "filename too long" error?** Run `git config --global core.longpaths true`
> in an elevated Git Bash
> ([reference](https://stackoverflow.com/questions/22575662/filename-too-long-in-git-for-windows)).

## Test systems

Balanced single-phase OpenDSS feeders in [`rawData/`](rawData/). Diagrams of the
main three are in [`assets/networks/`](assets/networks/).

| Directory | Paper name | Size | Used by |
|---|---|---|---|
| `ieee123C_1ph` | ieee123 | 128 buses | DDP, tADMM |
| `ieee2522C_1ph` | med2522 | 2,522 buses | DDP, tADMM |
| `large10kC_1ph` | large10k | 10,321 buses | DDP, tADMM |
| `ieee123_5poi_1ph` | | IEEE 123-node with five substations | MSOPF |
| `small2poi_1ph` | | two-substation test feeder | MSOPF |

## MSOPF: multi-source OPF

Code in [`envs/multi_poi/`](envs/multi_poi/), which is being renamed
**MultiSourceOPF**. Earlier multi-source multi-period OPF results are archived in
the [PESGM 2026 release](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/releases/tag/251117-pesgm-multi-source-mpopf).
Current work studies how the voltage-angle difference between two substations
drives power flow, and how batteries widen the range of angles in which no
substation back-feeds:

- [`full_angle_pf.jl`](envs/multi_poi/full_angle_pf.jl): AC power flow and OPF
  with exact bus-voltage angles, verified point by point against OpenDSS.
- [`root_level/two_source_angle_sweep.jl`](envs/multi_poi/root_level/two_source_angle_sweep.jl):
  sweeps the angle difference between two substations.
- [`root_level/battery_angle_sweep.jl`](envs/multi_poi/root_level/battery_angle_sweep.jl):
  the same sweep with battery dispatch.

MSOPF does not have its own README yet. Until it does, the header comment of
each driver documents its usage.

## Repository layout

```text
ddp/                DDP workstream: solver, drivers, results, notes
envs/
  ddp2026/          Julia environment for all DDP work
  tadmm/            tADMM code and environment
  multi_poi/        MSOPF code and environment
  ddp/              the original first-order DDP scheme (2025), kept as read-only reference
rawData/            OpenDSS feeder models
assets/networks/    feeder diagrams
results/            curated tADMM paper results
plots/              tADMM figure scripts
scripts/            helpers and run queues; scripts/probes/ holds one-off experiments
docs/               installation guide
```

## Reference papers

Each `resources/` folder (`ddp/resources/`, `envs/ddp/resources/`,
`envs/multi_poi/resources/`) has a `MANIFEST.txt` listing the papers it relies on,
with source links. Other authors' PDFs are not committed. Fetch the openly
available ones with:

```bash
bash scripts/fetch_resources.sh
```

## Branches and contributor notes

`master` holds merged work. New work goes on a dated branch named for its
workstream, e.g. `ddp-hessian-rewrites-sep20` or `msopf-initial-testing-sep17`.

[`CLAUDE.md`](CLAUDE.md) and [`AGENTS.md`](AGENTS.md) record standing rules and
context for AI coding agents, and are worth reading before producing DDP numbers.
They include the settings every reported comparison must share.
