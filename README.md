# MultiPeriodDistOPF

Research code and OpenDSS test networks for optimal power flow on distribution
feeders with **solar PV, batteries and multiple substations**.

**[Research lines](#three-research-lines)** ·
**[Test networks](#test-networks)** ·
**[Getting started](#getting-started)** ·
**[Citation](#citation)**

## Three research lines

The repository holds three independent lines of work on the same branch-flow
model. Each has its own folder, Julia environment and README.

| | What it does | Start here |
|---|---|---|
| **SPAR-DDP** | Multi-period OPF by sparse, parallel Differential Dynamic Programming: battery energy is the only state linking periods, and each period's network is eliminated locally. | [problems/mpopf/spar_ddp](problems/mpopf/spar_ddp/README.md) |
| **tADMM** | Multi-period OPF by temporal ADMM: single-period OPFs solved in parallel, coordinated with whole-horizon battery problems through consensus on the battery powers. | [problems/mpopf/tadmm](problems/mpopf/tadmm/README.md) |
| **MSOPF** | Multi-source OPF: feeders supplied by several substations with fixed voltage angles, and the dispatch that keeps them feasible over a day. | [problems/msopf](problems/msopf/README.md) |

SPAR-DDP and tADMM are two ways to solve the same multi-period problem (MPOPF);
MSOPF asks a different question on related networks. [problems/](problems/README.md)
lists all three with their environments.

## Test networks

The OpenDSS decks are usable on their own, **without Julia or any optimization
license**. Small to large:
[ieee123](test_networks/ieee123C_1ph),
[med2522](test_networks/ieee2522C_1ph),
[large10k](test_networks/large10kC_1ph), the
[multi-source ieee123](test_networks/ieee123_5poi_1ph), and converted EPRI,
IEEE 8500/9500 and SMART-DS cases.

They are research adaptations, not untouched official feeder models. The
[network catalogue](test_networks/README.md) gives their sizes, provenance,
shared files and how to run them in OpenDSS.

## Getting started

[Download ZIP](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/archive/refs/heads/master.zip),
or clone:

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

Each research line is a separate Julia project (Julia 1.12). From the repository
root, instantiate the one you want:

```bash
julia --project=problems/mpopf/spar_ddp -e 'using Pkg; Pkg.instantiate()'
julia --project=problems/mpopf/tadmm    -e 'using Pkg; Pkg.instantiate()'
julia --project=problems/msopf          -e 'using Pkg; Pkg.instantiate()'
```

then follow that folder's README. There is no single root `main.jl`.

- **Solvers.** Everything runs with open-source Ipopt and its bundled MUMPS.
  The HSL linear solvers (MA57/MA97), Gurobi and MOSEK are optional and only
  speed up or extend the comparisons.
- **Scale.** OpenDSS power flows and the small ieee123 cases run on a
  workstation. The large networks at long horizons were run on a university
  computing cluster and need substantial memory and many cores; the recorded
  job scripts are a starting point, not a guarantee for every machine.
- **Setup help.** The [setup guide](docs/first-time-installation-guide.md) covers
  Windows and Bash helpers and optional solver downloads.

## Citation

Please cite this repository when you use its **code, test networks or editable
figures**:

> Aryan Ritwajeet Jha and contributors, *MultiPeriodDistOPF*,
> https://github.com/Realife-Brahmin/MultiPeriodDistOPF.

Use GitHub's **"Cite this repository"** button or [CITATION.cff](CITATION.cff).
Please also cite the [original network sources](test_networks/README.md#cite-the-networks-you-use)
and the method you build on. SPAR-DDP's solver is a fork of FilterDDP.jl and keeps
its [attribution and citation guidance](problems/mpopf/spar_ddp/SPARDDP.jl/NOTICE.md).
The citation request adds no restriction to the [MIT license](LICENSE);
third-party notices remain in force.

## More

- [Editable network illustrations](figures/networks/README.md)
- [SPAR-DDP studies and reproducibility evidence](problems/mpopf/spar_ddp/experiments/README.md)
- [MSOPF interactive region viewer](https://realife-brahmin.github.io/MultiPeriodDistOPF/four_regions.html)
- [Contributing and repository maintenance](docs/repository-maintenance.md)

Questions, reproduction reports and suggestions are welcome as
[GitHub issues](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/issues).
If you try a case, please give its name, the commit and the command you ran.
