# Restored July convergence view

The user requested the original single-run objective/r-norm/s-norm/rho plot.
Style recovered from the July 31, 2026 release:
`50b8901c324ba9eefd0d672c6e484e0ededc7a57`,
`envs/tadmm/Plotter.jl::plot_tadmm_ldf_convergence`.

Four stacked panels share iteration on the horizontal axis. The original
Computer Modern font and metric colors are restored: objective dodgerblue
(#1E90FF), primal darkgreen (#006400), dual darkorange2 (#EE7600), and rho
purple (#800080). Exact RGB values were checked against Julia Colors' named
color definitions. Circle/square/diamond/hexagon markers match the old plot;
markers are sparse for readability, while all raw samples remain in the lines.

The solid green r-norm is the historical owned-period RMS, and the orange
s-norm is unique dual RMS. Dotted green additionally shows every-copy RMS.
These residuals are not substituted with original-model feasibility error.
The blue curve is actual objective value, with a dashed centralized reference,
not a relative-cost-gap plot. Rho is shown as the value actually used for each
iteration. Red lines show the actual saved run targets; the gray 1e-3 guide is
explicitly not the stopping threshold of these tighter runs. End-of-run status
is taken from the saved assessment, including the failed-to-converge arm.

## Current fixed-penalty runs

- [T48 PNG](T48-fixed48.png), [PDF](T48-fixed48.pdf),
  [exact plotted values](T48-fixed48.csv), [provenance/quality flags](T48-fixed48.json).
- [T96 PNG](T96-fixed48.png), [PDF](T96-fixed48.pdf),
  [exact plotted values](T96-fixed48.csv), [provenance/quality flags](T96-fixed48.json).

## Earlier controller for comparison

- [T48 PNG](T48-July-controller.png), [PDF](T48-July-controller.pdf).
- [T96 starting rho240000 PNG](T96-July-controller-240k.png),
  [PDF](T96-July-controller-240k.pdf).
- [T96 starting rho60000 PNG](T96-July-controller-60k.png),
  [PDF](T96-July-controller-60k.pdf). This arm hit the iteration budget.

All figures use the current completed runs, not historical numerical results.
No solves were launched, and no source trace, threshold or quality flag was
changed. CSV export is checked for exact numeric round-trip equality with the
selected raw columns. Each JSON records input/generator/style hashes, column
mapping, actual flags and style provenance. Figures were visually checked.

Run from the repository root, for example:

```powershell
python problems/mpopf/tadmm/review_pilot/plot_july_convergence.py `
  problems/mpopf/tadmm/review_pilot/results/targeted_nov09/results/run_T96_primal_dual_31657825 `
  problems/mpopf/tadmm/review_pilot/results/july_view_nov09/T96-fixed48
```

For multi-start cases, supply `--rho 60000` or `--rho 240000`. The code refuses
to silently choose an arm. The all-copy RMS and strict physical-check flags
remain available in the exports/audits; an RMS stopping flag is not a claim
that the SOC bounds are exact.
