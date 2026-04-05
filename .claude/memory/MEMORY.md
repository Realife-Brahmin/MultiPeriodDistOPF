# MultiPeriodDistOPF Project Memory

## Project Overview
- [Project basics](project_overview.md) — Julia MPOPF, BF vs tADMM, Ipopt only, 16 threads
- [User preferences](user_preferences.md) — CSVs over JLS, atomic commits, nonlinear-only PhD focus

## Key Results
- [Results table](results_table.md) — all BF vs tADMM comparisons across systems and T values
- [Convergence analysis](convergence_analysis.md) — adaptive ρ diagnosis, vanilla ADMM sweep results

## Systems
- [large10k_1ph](system_large10k.md) — 10321 buses, 102 batteries, original system
- [large10kB_1ph](system_large10kB.md) — 10321 buses, 1020 batteries (small, 80-120 kVA)
- [large10kC_1ph](system_large10kC.md) — 10321 buses, 1020 batteries (5x scaled), T=4 and T=6 done

## Current Work
- [Active investigation](active_investigation_mar31.md) — k=1 warm-start experiment, loose-tol SOCP pre-solve running
- [Known issues](known_issues.md) — git long filenames, Plotter.jl NaN, config in mid-experiment state
