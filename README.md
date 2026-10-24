# MultiPeriodDistOPF

Explore optimal power flow on distribution networks with **solar PV, battery
storage and multiple substations**. This repository provides Julia research
implementations and OpenDSS test cases you can use on their own.

**[Browse the test networks](test_networks/README.md)** ·
**[Try the optimization code](#try-the-optimization-code)** ·
**[How to cite](#citation)**

## Just want a network to experiment with?

Start with the [IEEE123 single-phase adaptation](rawData/ieee123C_1ph/), or
try [med2522](rawData/ieee2522C_1ph/) and
[large10k](rawData/large10kC_1ph/). There are also
[multi-source examples](rawData/ieee123_5poi_1ph/) and converted EPRI,
IEEE 8500/9500 and SMART-DS cases.

These are research adaptations, not untouched official feeder models.
The [network catalogue](test_networks/README.md) explains their sizes,
provenance, shared files and how to run them in OpenDSS. **You do not need
Julia or an optimization license just to use the OpenDSS decks.**

[Download ZIP](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/archive/refs/heads/master.zip),
or clone:

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

## Try the optimization code

Choose the method matching your question:

- **[SPAR-DDP](spar_ddp/README.md)** — sparse parallel Differential Dynamic
  Programming for multi-period OPF, exploiting radial-network and temporal structure.
- **[tADMM](envs/tadmm/README.md)** — temporal ADMM decomposition, with a
  monolithic branch-flow baseline and reproduction instructions.
- **[MSOPF](envs/multi_poi/README.md)** — multi-source OPF with explicit source
  voltage angles and multi-period battery scheduling.

For SPAR-DDP, use Julia 1.12.4+ within the 1.12 series. From the repository
root, instantiate its solver and OpenDSS-export environments (Bash/Git Bash):

```bash
julia --project=envs/spar_ddp -e 'using Pkg; Pkg.instantiate()'
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
```

Then follow the [small-instance walkthrough](spar_ddp/README.md#environments-and-entry-points).
Each method has its own environment and modeling assumptions; there is no
single root `main.jl`. Gurobi and HSL backends require their respective licenses
or access; read the chosen method's guide before selecting a solver.
See the [installation guide](docs/first-time-installation-guide.md) for setup.

## Citation

If you use the **code, OpenDSS network adaptations or figures** in research,
please cite this repository and identify the network and commit you used:

> Aryan Ritwajeet Jha and contributors, *MultiPeriodDistOPF*,
> https://github.com/Realife-Brahmin/MultiPeriodDistOPF, commit `<commit>`.

Use GitHub's **“Cite this repository”** option or [CITATION.cff](CITATION.cff).
Also cite the [original network sources](test_networks/README.md#cite-the-networks-you-use)
and the method you build upon. SPAR-DDP's DDP4OPF fork retains the upstream
[FilterDDP attribution and citation guidance](spar_ddp/DDP4OPF.jl/NOTICE.md).
This citation request does not add restrictions to the [MIT license](LICENSE);
applicable third-party notices remain in force.

## Explore further

- [Editable network illustrations](figures/networks/README.md).
- [SPAR-DDP studies and reproducibility evidence](spar_ddp/experiments/README.md).
- [MSOPF interactive region viewer](https://realife-brahmin.github.io/MultiPeriodDistOPF/four_regions.html).
- [Contributing and repository maintenance](docs/repository-maintenance.md).

Questions, reproduction reports and suggestions are welcome through
[GitHub issues](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/issues).
If you try a case, please include its name, commit and the command you ran.
