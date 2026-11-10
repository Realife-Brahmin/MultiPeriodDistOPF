# Multi-source, multi-period distribution OPF

This is the focused MSOPF reference snapshot for the PESGM 2027 study. It contains the model, adapted OpenDSS networks, recorded Julia environments, and archived simulation results for distribution networks supplied simultaneously by multiple substations with prescribed upstream voltage phasors.

The study maps sampled feasible operating points and infeasibility bounds with and without battery operation. It provides a mathematical model for evaluating this operating option; it does not establish general superiority to radial reconfiguration.

## Explore the results first

- [Interactive feasibility and certificate maps](https://realife-brahmin.github.io/MultiPeriodDistOPF/farkas_maps.html): IEEE123, med3980, and sf6937.
- [Paper-result index](problems/msopf/RESULTS.md): exact archived run folders and the meaning of each phase.
- [Release and validation notes](problems/msopf/RELEASE.md): provenance, certificate-check coverage, and interpretation limits.
- [Recorded environments and hardware](problems/msopf/environments/kamiak/README.md): the sf6937 Kamiak campaign is the reference configuration.
- [Software credits](problems/msopf/SOFTWARE_CREDITS.md): modelling, numerical solvers, simulation engines, and plotting tools.

No Julia installation is needed to inspect the archived CSVs or HTML maps. The included [HTML snapshot](docs/farkas_maps.html) contains its numerical data and loads Plotly.js from a CDN; opening it needs internet access. The live site can change, so retain the HTML from this version for a fixed reference. [Earlier four-source polytopes](docs/four_regions.html) are supporting historical material, not the current certificate maps.

## Obtain this snapshot

Download **Source code (zip)** from the [MSOPF release](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/releases/tag/msopf-pesgm2027-v1.0.0), or select the `msopf-pesgm2027-release` branch on GitHub. This branch presents only the MSOPF workstream, its shared network inputs, and supporting download/export tools. The full research repository, including the separate DDP workstreams, remains on [`master`](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/tree/master).

## Use the code

Start with the [MSOPF guide](problems/msopf/README.md). The reference environments record Julia **1.12.7**, JuMP **1.30.0**, MathOptInterface **1.50.0**, Ipopt **1.14.1**, Clarabel **0.11.1**, and OpenDSSDirect **0.9.9**. Whole-day certificates used MKL Pardiso; AC solves used an externally supplied MA57 library. Installation and licence requirements are explained in the environment notes. The original manifests are preserved rather than re-resolved.

The largest reference tasks requested **80 GB RAM and 8 CPUs each**. Those resources concern recomputation, not browsing the results. No new numerical solves were performed when assembling this release.

## Repository contents

- `problems/msopf/`: equations, solvers, campaign drivers, archived evidence, and release documentation.
- `test_networks/`: shared profiles, adapted network decks, conversion tools, and upstream provenance.
- `docs/`: interactive result viewers.
- `tools/`: optional public-input download and network-export helpers.

Historical results retain their original paths and labels. Use the paper-result index to select the rated-current runs matching this study. External solver binaries and downloaded upstream datasets are subject to their own availability and licences.

## Attribution

Please cite the software using [CITATION.cff](CITATION.cff), retain the [MIT licence](LICENSE), and credit the [upstream networks](test_networks/README.md) and [software packages](problems/msopf/SOFTWARE_CREDITS.md). This snapshot accompanies research in preparation and does not imply paper acceptance.
