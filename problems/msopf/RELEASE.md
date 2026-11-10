# MSOPF reference state — PESGM 2027 study

This reference state packages the rated-current, multi-period MSOPF study with two to four substations, PV reactive dispatch and batteries returning to their initial charge. It preserves source code and archived evidence; it does not establish a globally certified continuous operating region or comparative superiority to radial reconfiguration.

Start with the [interactive Farkas maps](https://realife-brahmin.github.io/MultiPeriodDistOPF/farkas_maps.html), the [recorded environments](environments/kamiak/README.md), the [paper-result index](RESULTS.md), and the [software credits](SOFTWARE_CREDITS.md). The [HTML snapshot](../../docs/farkas_maps.html) is included in the repository. The live site can change after a release; use the released snapshot when citing a particular version.

## Validation coverage

All **78/78 sf6937 whole-day certificates pass exact-arithmetic checks**: 59 during the original solves and 19 in post-run VD audits. The [six follow-up audits](results/kamiak/farkas_day148b_nov03/reference_config/vd_audits_nov09/README.md) complete directions 2, 11, 14, 16, 17 and 19 at the original published CD radii. Direction 11 is certified by a negative exact objective upper bound despite a slow-progress solver status; the other five follow-up checks use Farkas certificates. Per-direction coverage is recorded in [cd_exact_status.csv](results/kamiak/farkas_day148b_nov03/reference_config/cd_exact_status.csv). The original report/archive retains the earlier 72/78 snapshot as historical provenance; this release supersedes that coverage statement.

The follow-up audits do not change any numerical operating point, idle certificate or validated battery-reach result. Curves, shells and polytope displays interpolate or predict sampled behaviour as labelled; they are not certificates between tested directions. A failed local solve is not an infeasibility proof.

## Canonical hardware and numerical stack

The reference is the sf6937 campaign on **Kamiak**, using Julia **1.12.7**. AC phases used Ipopt 1.14.1 with a private MA57 3.11.3 build; idle certificates used Clarabel 0.11.1 with QDLDL. Whole-day CD jobs used Clarabel 0.11.1 with Pardiso.jl 1.1.2 and MKL_jll 2025.2.0+0. JuMP was 1.30.0, MathOptInterface 1.50.0, OpenDSSDirect 0.9.9 and the recorded DSS C-API engine 0.14.3. Exact manifests and external-library qualifications are in the environment README.

The 78 whole-day tasks requested **8 CPUs and 80 GB RAM each** on shared nodes: 74 tasks on Intel Xeon 6760P and four on Xeon Gold 6438Y+. Slurm batch-step peak RSS was about **43.7 GB median / 45.7 GB maximum**; wall time per direction, including repeated relaxations, was **2.71 h median / 9.03 h maximum**. These timings are observations from that campaign, not resource or runtime guarantees for another machine. Per-relaxation application RSS is a different measurement (43.5 / 46.4 GB median/max).

One Julia process ran per task. JULIA_NUM_THREADS, OPENBLAS_NUM_THREADS and OMP_NUM_THREADS were 1; MKL_NUM_THREADS was 8. Actual Pardiso thread count was not logged. The observed CPU utilization was about 2.2 busy cores median, so resource allocation must not be treated as measured parallel speedup.

## Provenance and known limits

The checkout moved among commits 335a30e3, 4f29d098, 5712c808 and 865e26d1 while arrays ran. The reference audit found the physical model code and sf6937 deck unchanged across those commits; changes concerned certificate checks, audits and derived bounds. Uncommitted edits at execution time cannot be ruled out. The ten idle-edge directions run on the 309 PC are explicitly distinguished in the campaign README.

Slurm 17.11 did not preserve each job's environment or submit line. Settings were reconstructed from archived launcher metadata and recovered submit commands. Effective Pardiso thread count, exact runtime BLAS binding and compute-node system OpenBLAS version remain unknown. See the [canonical configuration](results/kamiak/farkas_day148b_nov03/reference_config/canonical_configuration.json) for machine-readable package, resource and coverage metadata.

The HTML viewer loads **Plotly.js 2.27.0** from jsDelivr and optional Google Fonts. Opening that snapshot therefore needs access to those resources. Original network sources, software credits and licences remain applicable; proprietary solver binaries are not included.

The current `root_level/farkas_viewer.py` builder uses Python standard-library modules; it does not require Julia or a new solve. Supporting prediction/plotting scripts have their own NumPy/Matplotlib dependencies. Their historical Python package versions were not captured by the canonical Julia manifests and are not presented as pinned campaign software.
