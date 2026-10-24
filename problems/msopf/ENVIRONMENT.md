# problems/msopf

The Julia environment of the MSOPF (multi-source OPF) workstream: the drivers in
`problems/msopf/` run as `julia --project=problems/msopf <script>`.

- **Created 2026-10-19 as an exact copy** of `problems/mpopf/tadmm/Project.toml` and
  `problems/mpopf/tadmm/Manifest.toml`, where the MSOPF drivers had run until then (user's
  request: MSOPF should have its own environment). `Pkg.instantiate()` leaves the
  copied Manifest byte-identical, so every committed MSOPF result before that date
  was produced with these same package versions.
- **Pinned:** JuMP 1.30.0, Ipopt 1.14.1, OpenDSSDirect 0.9.9, Gurobi 1.9.2 (Julia 1.12).
- **Some dependencies are not used by MSOPF** (Plots, Revise, Crayons, LaTeXStrings,
  Parameters, the tADMM driver's own needs). They are kept on purpose: removing them
  would make Pkg resolve the environment again and could change the versions above.
  Prune them only in a deliberate, separately checked change.
- For exact historical reproduction use its recorded snapshot/Manifest, not a
  fresh resolve. Do not `Pkg.update` a pinned campaign to fix an error. Pkg strips comments from
  `Project.toml`, so notes belong here.
- Checked on Windows (309 PC) and Linux (Kamiak) with
  `problems/msopf/root_level/portability_check.jl`, which recomputes committed IEEE123
  results; see `problems/msopf/README.md`.
