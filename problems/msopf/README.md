# MSOPF (multi-source OPF)

Code for the multi-source OPF workstream: distribution feeders fed by several
substations whose voltage angles are fixed upstream. The implementation starts
with `solve_full_angle` in [`full_angle_pf.jl`](full_angle_pf.jl).

## Setup on any machine

The local **Project.toml** is MSOPF's current Julia environment, beside its
models, scripts and results. Run with `--project=problems/msopf` from the
repository root. The older branch-flow environment is preserved in
`archive/legacy-environment/`; it is not the main environment.
See [environment notes](ENVIRONMENT.md). Historical reference versions were:

| Package | Version |
|---|---|
| JuMP | 1.30.0 |
| Ipopt | 1.14.1 |
| OpenDSSDirect | 0.9.9 |
| Gurobi | 1.9.2 |

- **Julia 1.12.x.** Historical environments used 1.12.5.
- **Instantiate once per clone.** Exact historical reproduction needs the recorded
  campaign Manifest; a fresh environment without one resolves package versions:
  ```bash
  julia --project=problems/msopf -e 'using Pkg; Pkg.instantiate()'
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
- **Leave the Manifest alone.** Don't `Pkg.update`, `Pkg.resolve` or delete `problems/msopf/Manifest.toml`
  to fix an error. Pkg also strips comments from `Project.toml`, so put notes in READMEs.

`problems/msopf/Project.toml` belongs to the older branch-flow scripts
(`root_level/multi_poi_mpopf.jl` and similar). It has no tracked Manifest.
`README_HOME_PC.md` is a November 2025 note about that environment.

## Drivers

Run from the repo root with `julia --project=problems/msopf <script> [args]`.

| Script | What it does | Output |
|---|---|---|
| `root_level/two_source_angle_sweep.jl small2poi \| ieee123 [subsA subsB]` | power flow over the substation angle difference, checked point by point against OpenDSS | `processedData/` (gitignored) |
| `root_level/battery_angle_sweep.jl small2poi \| ieee123` | the no-backflow window with and without battery dispatch | `processedData/` (gitignored) |
| `root_level/source_impedance_study.jl` | no-backflow windows versus the source impedance, with and without batteries and the 0.95–1.05 pu band | `results/source_impedance/` (**tracked**: `windows.csv`, and `table_windows.tex` for the paper) |
| `root_level/multiperiod_study.jl [stage ...]` | the same model over a whole day, with battery state of charge: how far a constant angle difference can go. Stages: `main horizon zs duration pairs three five bound tables`, plus `ieee2522`, `large10k` and `large10kd` (large10k with loads, PV and batteries divided by 20), which run only when named | `results/multiperiod/` (**tracked**: CSVs, figures, and `tables.tex` for the paper) |
| `check_voltage_bases.jl` | checks each deck's `VoltageBases` against its Vsource setpoints | stdout |

## Multi-period model

[`multiperiod_full_angle.jl`](multiperiod_full_angle.jl) repeats the single-period equations in every
period and links the periods through each battery's energy. Its conventions follow the DDP work:

- **Profiles:** load, price and PV shapes from `problems/mpopf/tadmm/root_level/config.jl`, with periodic sampling.
- **Batteries:** one power variable each, lossless, 30–95% state of charge from 62.5%.
- **Costs:** `C_B` comes from the sizing rule in
  `problems/mpopf/spar_ddp/examples/power_system/terminal_soc_penalty.jl`, applied to this deck's fleet.
- **Terminal SOC is hard: every battery ends the day where it started** (`terminal = :cyclic`, the
  default). This is an MSOPF-only exception to the DDP rule of a soft terminal SOC, decided by the user on
  2026-10-11 because feasibility is the question here. The soft rule (`terminal = :soft`, with `gamma`
  from the same sizing file) is kept as a diagnostic: it lets batteries end the day fuller, so its
  feasible window grows without bound with battery size, while the hard rule levels off.
- **PV:** real power is fixed at `lambda_PV * Pmpp`; reactive power is a decision within the inverter rating.

Three things to know before running sweeps with it (each cost hours on 2026-10-11):

- **Do not ask Ipopt to prove a point infeasible.** That takes 500+ iterations, about 30 s at T = 24.
  To find a window edge, maximize the angle directly (`objective = :ray`): one solve per edge. To test
  one point, use `objective = :violation`, which adds slacks and is feasible exactly when the minimum
  violation is zero.
- **Each solve must release its Ipopt problem** (`release_solver!`, called inside `solve_multiperiod`).
  Otherwise about 22 MB of solver memory per solve piles up outside Julia's heap, and a few hundred
  solves swap a 16 GB machine to a standstill.
- **On Windows, stopping a background job does not stop its Julia process.** Stop it by process ID.

On the 10k-bus decks (2026-10-13):

- **Use HSL MA57 in Ipopt** through `IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=<dir>/libma57.dll;linear_system_scaling=none"`
  (the DDP scripts' convention; the library is built locally, not in the repository).
- **Gurobi's nonconvex solver is not an alternative there.** `solve_multiperiod(...; solver = :gurobi)` builds the
  same model, but on large10kD at T = 2 it found no feasible point in 15 minutes (8 threads), with or without
  batteries, where Ipopt with MA57 solves the same import-margin problem in under 30 s.
- **Do not find window edges with `objective = :ray` there.** It stopped early, ran out of iterations or ended
  locally infeasible. `window_secant` (a secant on the import margin) needs six or seven solves per window.

## Decks

- `test_networks/ieee123_5poi_1ph` has five substations; `test_networks/small2poi_1ph` has two. Both are 20 kV line-to-line.
- Since 2026-09-28, every substation has a realistic source impedance, Zs = 0.032 + j0.799 Ω.
  - |Zs| = U²/S_k, using the typical 20 kV short-circuit level of 0.50 GVA
    (Traupmann & Kienberger 2020, Table 19).
  - X/R = 25, inside the 15–40 range of IEEE C37.010-2016, Table 24.
- `solve_full_angle(...; ideal_sources = true)` recovers the old stiff sources.

## Customer-bus test-system studies (2026-10-25)

`root_level/multiperiod_study.jl` now recognizes these four-source stages:

- `four_ieee9500xf`: the service-transformer IEEE9500 deck, five feeder ties closed, and the fourth source at `l2804256` selected in `choose_fourth_9500.jl`. Its 2,702 primary line definitions are identical to the primary-only conversion; the added secondary branches do not change the primary series effective impedances or the site choice.
- `four_smartds_p1u_typeiii`: `subs1`–`subs4`, switching read from `Master_4sub_typeIII.dss`.
- `four_smartds_p1u_typeii`: `subs1`, `subs2`, `subs4`, `subs5`, switching read from `Master_4sub_typeII.dss`.

For decks containing `svc_*` lines, `setup` selects their receiving buses as customers.
The selection is passed through `problem` to every day/hour/edge solve, and the same
selection is used to report voltage extrema and classify voltage-limited edges.
Primary-only decks retain their existing band on non-source buses.

Run the local preparation, serially, from Git Bash at the repository root:

```bash
bash problems/msopf/root_level/run_test_edition_sensitivities.sh
```

This runs IEEE9500xf and P1U Type III at T=24 with PV reactive dispatch and the
existing DER fleets. `MULTI_SENS_ONLY=1` computes the baseline plus six central-difference
solves (three moving angles, +/-0.1 degree). Every successful solve is replayed for
all 24 hours in OpenDSS. Replay tolerances: 1e-3 kW/kvar/battery-kW and 1e-7 pu;
OpenDSS must converge and keep the loads constant-PQ. A failed solve/replay stops the
stage before a sensitivity prediction is published.

The sensitivity CSV contains hourly imports and derivatives. The predicted-reach CSV
contains 78 directions for both idle batteries (worst hour) and cyclic batteries
(daily mean). **These are estimates based on idle-dispatch sensitivities, not measured
feasibility boundaries or global certificates.** Voltage limits and PV redispatch can
change the boundaries. A region should only be called a simplex after its geometry
and the independent edge searches support that statement.

Each case writes a new directory under `results/test_networks/sensitivity_oct25/`, with
raw stdout/stderr, replay CSV, environment/solver/thread settings, source and deck hashes,
Git revision and source patch, wall time, and exit code. Existing run directories are
never overwritten; set `MSOPF_SENS_OUT` to a new root for a rerun. Each solve has a
600-second wall cap. The local run uses one Julia/BLAS/OpenMP thread and lean JuMP mode;
these preparation timings are not matched benchmark timings.

The regression `tests/customer_band.jl` solves a primary bus above 1.05 pu with a
customer inside the band, compares customer-only versus legacy margin diagnostics,
and replays both dispatches in OpenDSS. It also checks the legacy primary-only rule,
variant source selections, and the identical IEEE9500 primary impedances:

```bash
julia --startup-file=no --project=problems/msopf problems/msopf/tests/customer_band.jl
```

The single-substation K1, Ckt24 and 8500xf decks remain covered by
`test_edition_smoke.jl`; four-source stages do not invent sites for those systems.
The full 78-direction edge campaign is a separate, much larger run. Use
`MULTI_DIRSET=fib78 MULTI_VMARGIN=20000 MULTI_SAMPLES=0` with the desired new stage,
without `MULTI_SENS_ONLY`, after checking cluster access and the shared queue.

Validate a finished preparation directory and regenerate its scientific figures without
running any optimization:

```bash
python problems/msopf/root_level/check_test_edition_sensitivities.py <run-directory>
python problems/msopf/root_level/plot_test_edition_sensitivities.py <run-directory>
```

The validator requires a successful process exit, all seven converged OpenDSS replays,
complete hourly sensitivities, and 156 consistent directional predictions. The plotter
uses NumPy and Matplotlib, verifies the four intersections of the daily-mean import
planes, and exports PNG plus editable SVG. These figures show predictions only;
measured voltage cuts and nonlinear corrections await the edge campaign.
