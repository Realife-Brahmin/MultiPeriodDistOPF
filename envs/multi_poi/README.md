# MSOPF (multi-source OPF)

Code for the multi-source OPF workstream: distribution feeders fed by several
substations whose voltage angles are fixed upstream. The paper is the separate
`PESGM2027_Multi-Source-Multi-Period-OPF` repository; its Section II documents
`solve_full_angle` in [`full_angle_pf.jl`](full_angle_pf.jl).

## Setup on any machine

The current drivers run in **`envs/tadmm`**, not `envs/multi_poi`.
`envs/tadmm` has a **tracked `Manifest.toml`**, so every machine gets identical
package versions:

| Package | Version |
|---|---|
| JuMP | 1.30.0 |
| Ipopt | 1.14.1 |
| OpenDSSDirect | 0.9.9 |
| Gurobi | 1.9.2 |

- **Julia 1.12.x.** The Manifest was resolved on 1.12.5; any 1.12 patch release works.
- **Instantiate once per clone.** This installs the exact versions; it never resolves new ones:
  ```bash
  julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
  ```
  Run Julia from **Git Bash**, not PowerShell. PowerShell mangles the quotes in `-e '...'`.
- **OpenDSS needs no separate install.** OpenDSSDirect ships the engine as a Julia artifact.
- **Gurobi is optional.** Its binary also comes as an artifact, and only a license is needed. Without one,
  `gurobi_usable()` skips the Gurobi cross-checks with a message; Ipopt does everything else.
- **Harmless error on Julia 1.12.** On load, OpenDSSDirect 0.9.9 prints
  `ERROR: Method overwriting is not permitted during Module precompilation`, plus some
  `WARNING: Constructor ... extended` lines. The package then loads without its
  precompiled cache and works, just a few seconds slower. Ignore it.
- **Checked from a fresh clone.** On 2026-09-28 the branch was cloned into an empty
  folder, instantiated, and run: solve and OpenDSS agree to 1.6e-5 kW and 5e-10 pu.
- **Leave the Manifest alone.** Don't `Pkg.update`, `Pkg.resolve` or delete `envs/tadmm/Manifest.toml`
  to fix an error. Pkg also strips comments from `Project.toml`, so put notes in READMEs.

`envs/multi_poi/Project.toml` belongs to the older branch-flow scripts
(`root_level/multi_poi_mpopf.jl` and similar). It has no tracked Manifest.
`README_HOME_PC.md` is a November 2025 note about that environment.

## Drivers

Run from the repo root with `julia --project=envs/tadmm <script> [args]`.

| Script | What it does | Output |
|---|---|---|
| `root_level/two_source_angle_sweep.jl small2poi \| ieee123 [subsA subsB]` | power flow over the substation angle difference, checked point by point against OpenDSS | `processedData/` (gitignored) |
| `root_level/battery_angle_sweep.jl small2poi \| ieee123` | the no-backflow window with and without battery dispatch | `processedData/` (gitignored) |
| `root_level/source_impedance_study.jl` | no-backflow windows versus the source impedance, with and without batteries and the 0.95–1.05 pu band | `results/source_impedance/` (**tracked**: `windows.csv`, and `table_windows.tex` for the paper) |
| `check_voltage_bases.jl` | checks each deck's `VoltageBases` against its Vsource setpoints | stdout |

## Decks

- `rawData/ieee123_5poi_1ph` has five substations; `rawData/small2poi_1ph` has two. Both are 20 kV line-to-line.
- Since 2026-09-28, every substation has a realistic source impedance, Zs = 0.032 + j0.799 Ω.
  - |Zs| = U²/S_k, using the typical 20 kV short-circuit level of 0.50 GVA
    (Traupmann & Kienberger 2020, Table 19).
  - X/R = 25, inside the 15–40 range of IEEE C37.010-2016, Table 24.
- `solve_full_angle(...; ideal_sources = true)` recovers the old stiff sources.
