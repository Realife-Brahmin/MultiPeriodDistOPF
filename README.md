# MultiPeriodDistOPF

Explore optimal power flow on distribution networks with **solar PV, battery
storage and multiple substations**. This repository provides Julia research
implementations and OpenDSS test cases you can use on their own.

**[Browse the test networks](test_networks/README.md)** ·
**[Try the optimization code](#try-the-optimization-code)** ·
**[How to cite](#citation)**

## Just want a network to experiment with?

Start with the [IEEE123 single-phase adaptation](test_networks/ieee123C_1ph), or
try [med2522](test_networks/ieee2522C_1ph) and
[large10k](test_networks/large10kC_1ph). There are also
[multi-source examples](test_networks/ieee123_5poi_1ph) and converted EPRI,
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

Browse [implementations by problem](problems/README.md):

- **MPOPF — multi-period scheduling**

  - **[SPAR-DDP](problems/mpopf/spar_ddp/README.md)** — sparse parallel Differential Dynamic
    Programming for multi-period OPF, exploiting radial-network and temporal structure.
  - **[tADMM](problems/mpopf/tadmm/README.md)** — temporal ADMM decomposition, with a
    monolithic branch-flow baseline and reproduction instructions.
- **[MSOPF](problems/msopf/README.md)** — multi-source OPF with explicit source
  voltage angles and multi-period battery scheduling.

For SPAR-DDP, use Julia 1.12.4+ within the 1.12 series. From the repository
root, instantiate its environment, including OpenDSS (Bash/Git Bash):

```bash
julia --project=problems/mpopf/spar_ddp -e 'using Pkg; Pkg.instantiate()'
```

Then follow the [small-instance walkthrough](problems/mpopf/spar_ddp/README.md#environments-and-entry-points).
Each method has its own environment and modeling assumptions; there is no
single root `main.jl`. **SPAR-DDP needs no proprietary solver license.**
The centralized Ipopt reference can use its bundled open-source MUMPS backend.
MA57/MA97 (HSL backends for Ipopt) and Gurobi (an alternative centralized
optimizer) are optional ways to make the centralized comparison more competitive,
not requirements for running DDP or the basic Ipopt reference.
See the [installation guide](docs/first-time-installation-guide.md) for setup.

### Where to start

OpenDSS power flows and the small IEEE123 import/export checks are suitable
workstation starting points. Larger networks and long-horizon optimization may
need substantial RAM, worker/thread configuration and cluster job setup.
The latest large-case benchmarks were run on Kamiak; the recorded job scripts
are not a promise of out-of-the-box execution on every machine.

## Citation

Please cite this repository when using its **code, test networks or
editable Inkscape figures**:

> Aryan Ritwajeet Jha and contributors, *MultiPeriodDistOPF*,
> https://github.com/Realife-Brahmin/MultiPeriodDistOPF.

Use GitHub's **“Cite this repository”** option or [CITATION.cff](CITATION.cff).
Also cite the [original network sources](test_networks/README.md#cite-the-networks-you-use)
and the method you build upon. SPAR-DDP's SPARDDP fork retains the upstream
[FilterDDP attribution and citation guidance](problems/mpopf/spar_ddp/SPARDDP.jl/NOTICE.md).
This citation request does not add restrictions to the [MIT license](LICENSE);
applicable third-party notices remain in force.

## Explore further

- [Editable network illustrations](figures/networks/README.md).
- [SPAR-DDP studies and reproducibility evidence](problems/mpopf/spar_ddp/experiments/README.md).
- [MSOPF interactive region viewer](https://realife-brahmin.github.io/MultiPeriodDistOPF/four_regions.html).
- [Contributing and repository maintenance](docs/repository-maintenance.md).

Questions, reproduction reports and suggestions are welcome through
[GitHub issues](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/issues).
If you try a case, please include its name, commit and the command you ran.
