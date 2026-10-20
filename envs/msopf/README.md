# envs/msopf

The Julia environment of the MSOPF (multi-source OPF) workstream: the drivers in
`envs/multi_poi/` run as `julia --project=envs/msopf <script>`.

- **Created 2026-10-19 as an exact copy** of `envs/tadmm/Project.toml` and
  `envs/tadmm/Manifest.toml`, where the MSOPF drivers had run until then (user's
  request: MSOPF should have its own environment). `Pkg.instantiate()` leaves the
  copied Manifest byte-identical, so every committed MSOPF result before that date
  was produced with these same package versions.
- **Pinned:** JuMP 1.30.0, Ipopt 1.14.1, OpenDSSDirect 0.9.9, Gurobi 1.9.2 (Julia 1.12).
- **Some dependencies are not used by MSOPF** (Plots, Revise, Crayons, LaTeXStrings,
  Parameters, the tADMM driver's own needs). They are kept on purpose: removing them
  would make Pkg resolve the environment again and could change the versions above.
  Prune them only in a deliberate, separately checked change.
- Do not `Pkg.update` or `Pkg.resolve` it to fix an error. Pkg strips comments from
  `Project.toml`, so notes belong here.
- Checked on Windows (309 PC) and Linux (Kamiak) with
  `envs/multi_poi/root_level/portability_check.jl`, which recomputes committed IEEE123
  results; see `envs/multi_poi/README.md`.
