# Exact rewrites around the diagonal-Hessian FilterDDP path

## Result

Three opt-in implementation rewrites preserve the existing diagonal-Hessian
algorithm while avoiding work whose result is discarded:

1. `FILTERDDP_DIRECT_DIAG_HESSIAN=1` evaluates
   `diag(fu' * Vxx * fu)` directly instead of constructing the full battery
   block and deleting its off-diagonal entries.
2. `FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1` assembles the handwritten SOCP
   second derivatives from one set of sparse triplets instead of repeatedly
   inserting into a CSC matrix.
3. `FILTERDDP_CACHE_KKT_PATTERN=1` constructs each stage KKT sparsity pattern
   once and overwrites its numerical values on later backward passes. It does
   **not** reuse a numerical factorization.

Each switch was tested by itself on IEEE123 and IEEE2522 at `T=3`. Every
recorded objective, primal infeasibility, dual infeasibility and
complementarity value matched the baseline through the approved stopping
iteration. The combined implementation also preserves those trajectories and
the IEEE2522 `T=24` trajectory.

The strongest wall-clock result is IEEE2522 `T=24`: `755.189 s -> 591.613 s`
to the same iteration 87, a **21.66% reduction**. The final objective is
`8877.4115983`, maximum equality violation `9.279e-7`, dual infeasibility
`4.866e-7`, and complementarity `2.571e-10`.

The combined large10k `T=6` run reached the matched near-optimal criterion at
iteration 103 in `2330.671 s`: objective `3178533.3170` (relative gap
`1.177e-5`), maximum equality violation `8.362e-5`, zero dynamics and bound
violations, dual infeasibility `2.327e-3`, and complementarity `1.318e-5`.
There is no pre-rewrite run on this exact soft-terminal matched instance, so a
large10k wall-clock speedup must not be invented from an older instance.

## Where the saving comes from

On IEEE2522 `T=24`, cumulative second-derivative callback time through the
target iteration falls from `44.625 s` to `4.947 s`; derivative time overall
falls from `105.999 s` to `49.030 s`. KKT assembly falls from `54.289 s` to
`30.516 s`. Factorization is unchanged (`97.537 s` versus `98.751 s`), as it
should be, because the numerical KKT matrix and its factorization are not
approximated. The dominant multi-RHS solve remains dominant (`290.726 s`
versus `236.953 s` in these single runs); some of that apparent difference is
ordinary timing noise because the solved systems are identical.

The warm IEEE2522 allocation probe confirms the intended mechanical savings.
At a representative stage, derivative allocation falls `10.702 -> 9.159 MiB`,
algebra allocation `16.388 -> 9.867 MiB`, and KKT assembly allocation
`5.668 -> 1.041 MiB`. Factor, solve and update allocations are unchanged.

## Decision point

The rewrites are safe and worthwhile. They do not change the algorithm and
produce a material medium-system improvement. A matched large10k `T=24`
rerun is justified if a publishable large-system before/after number is worth
roughly another two hours of machine time. It has deliberately **not** been
started automatically: the completed large10k `T=6` run lacks an exact
pre-rewrite baseline, while the already completed large10k `T=24` baseline is
the proper comparison target.

Raw timing logs and complete iteration traces are under
`ddp/results/hessian_rewrites/`; tabulated comparisons are in
`runtime_comparison.csv` and `warm_allocation_comparison.csv`.
