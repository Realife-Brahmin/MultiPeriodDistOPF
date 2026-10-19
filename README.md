# MultiPeriodDistOPF

Julia research code for optimal power flow on distribution feeders with
batteries and solar PV. The repository contains three workstreams:

- **SPAR-DDP** — Sparse Parallel Differential Dynamic Programming for
  multi-period OPF. Active focus: [spar_ddp/README.md](spar_ddp/README.md).
- **tADMM** — temporal ADMM for MPOPF. Retained companion research code:
  [envs/tadmm/README.md](envs/tadmm/README.md). Earlier results remain available
  for reference; fresh testing and potentially further research are deferred.
- **MSOPF** — multi-source OPF:
  [envs/multi_poi/README.md](envs/multi_poi/README.md).

## Start here

Use Julia 1.12 for the current DDP workflow. Clone the repository, then
instantiate the environment belonging to the workstream you need:

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
julia --project=envs/spar_ddp -e 'using Pkg; Pkg.instantiate()'
```

SPAR-DDP uses `envs/spar_ddp` for the solver and centralized JuMP/Ipopt
reference. Exporting instances from OpenDSS uses `envs/tadmm`. Current MSOPF
drivers also use `envs/tadmm`; `envs/multi_poi` holds an older environment
for its branch-flow scripts. See the workstream READMEs before running.

There is no root Julia project or universal main.jl entry point. Gurobi is
optional and requires a license. HSL MA57/MA97 require separately granted
access and the appropriate library configuration; the recorded Kamiak
comparison used those backends. See the
[installation guide](docs/first-time-installation-guide.md).

## Repository map

```text
spar_ddp/                 DDP4OPF solver, MPOPF drivers, evidence and research notes
envs/spar_ddp/            Julia environment for SPAR-DDP and centralized references
envs/tadmm/               tADMM code, OpenDSS exporter and shared dependencies
envs/tadmm/plots/         historical tADMM plotting tools
envs/multi_poi/           multi-source OPF code and older environment
rawData/                  OpenDSS feeder models and input data
figures/networks/         editable network SVG diagrams and vector PDF exports
results/                  retained tADMM result collection
archives/first_order_ddp/ retired original first-order DDP and later diagnostics
scratch/                  exploratory helpers, old run queues and one-off probes
tools/                    maintained repository utilities
docs/                     installation and repository-maintenance guidance
```

The active workstream and its environment deliberately have separate homes:
`spar_ddp/` contains code/evidence; `envs/spar_ddp/` declares dependencies.
The old first-order DDP is archived, not a second supported environment.
Personal editor workspaces are local and ignored.

## Test systems and network illustrations

OpenDSS inputs live in [rawData/](rawData/). The main MPOPF systems are:

- `ieee123C_1ph` → **ieee123**, 128 modeled buses.
- `ieee2522C_1ph` → **med2522**, 2,522 buses.
- `large10kC_1ph` → **large10k**, 10,321 buses.

Multi-source studies also use `ieee123_5poi_1ph` and `small2poi_1ph`.
Raw IDs are stable; display names need not match the original IEEE feeder
node count. PV/storage placement and scheduling inputs are study choices.

[Network figures and attribution](figures/networks/README.md) explains each
SVG/PDF, its underlying test-system sources, and reuse. Please credit
Aryan Ritwajeet Jha and contributors, link to this repository, identify the
figure and commit used, and indicate modifications; cite the underlying
network publication separately. Retain applicable license notices. The
citation request adds no restrictions to the [MIT license](LICENSE).

## Current SPAR-DDP evidence

The [completed October 16 campaign](spar_ddp/results/proper_logging_campaign_oct16/README.md)
records 96 validated runs across three networks and four horizons, with two
repeats of each of four configurations. The current main-paper criterion is
0.1% reference-objective agreement together with original-model feasibility:
1e-6 for ieee123 and 1e-5 for med2522/large10k. Runs continued to full
convergence; crossing times were scored afterward. Standalone diagnostics
and older instance families are documented separately.

The manuscript is maintained in the separate
[SPAR-DDP paper repository](https://github.com/Realife-Brahmin/TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF).
It is in preparation; no new paper DOI or software release is implied.
When using DDP4OPF, retain its [NOTICE](spar_ddp/DDP4OPF.jl/NOTICE.md) and cite
the original FilterDDP method and implementation as described in its README.

## Historical paths and resources

The [October 18 path map](docs/repository-layout-oct18.md) explains renamed
folders. Recorded logs, CSVs, manifests and frozen campaign harnesses retain
their original contents, including old paths. Use the evidence commit for
exact reproduction; old harness launch commands are not current entry points.

Reference-paper `resources/MANIFEST.txt` files record sources. Download
missing papers with `bash tools/fetch_resources.sh`; existing files remain
untouched. Third-party source licenses and citations remain applicable.

## Contributions

`master` contains merged code. Work on a dated workstream branch, use focused
atomic commits, and open a PR. See [AGENTS.md](AGENTS.md) for current shared
guidance. `scratch/` is reviewed at research milestones; it is not a queue
of automatically authorized jobs. No fresh tADMM or DDP runs are required
merely to reorganize the repository.
