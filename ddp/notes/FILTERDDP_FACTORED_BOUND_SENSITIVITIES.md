# Factoring redundant bound sensitivities out of FilterDDP

## Result

The two dense bound-dual sensitivity matrices stored by every FilterDDP stage
are algebraically redundant. Removing them cuts the retained large10k update
rule from 1606.982 MiB to 757.012 MiB per stage without changing the algorithm
or its numerical result.

FilterDDP originally stored

```math
\zeta_l=-\operatorname{diag}(\Sigma_L)\beta,\qquad
\zeta_u= \operatorname{diag}(\Sigma_U)\beta.
```

The forward pass used these matrices only through `ζl * δx` and
`ζu * δx`. Since it already computes `βδx = β * δx` for the control
feedback, the identical dual updates are

```math
\zeta_l\delta x=-\Sigma_L\odot(\beta\delta x),\qquad
\zeta_u\delta x= \Sigma_U\odot(\beta\delta x).
```

The update rule therefore stores the two length-`nu` vectors `Σ_L` and `Σ_U`
instead of two `nu x nx` matrices. This changes storage from `2*nu*nx` to
`2*nu` floating-point values and adds no matrix solve.

## Validation

The complete IEEE123C `T = 3` regression remained identical to the committed
baseline: 48 iterations, objective `2808.924770488853`, and maximum equality
residual `6.227e-8`.

A bounded one-iteration large10kC `T = 3` memory run measured:

- `beta`: 425.402 MiB;
- `omega`: 329.201 MiB;
- factored bound data: 0.834 MiB, down from 850.804 MiB;
- complete retained update rule: 757.012 MiB, down from 1606.982 MiB.

That is a 99.90% reduction in the bound-sensitivity component and a 52.89%
reduction in the full retained update rule. Because one update rule is retained
per stage, the saving scales linearly with horizon: about 2.49 GiB at `T = 3`
and 79.68 GiB at `T = 96` for large10k dimensions.

This does not eliminate `beta` itself. The audit confirms that full `beta` is
the state-feedback policy used by the forward rollout, while full `omega` is
used to update equality multipliers. Eliminating either requires a different
forward-policy representation, recomputation, or retained factorization. The
bound matrices were uniquely attractive because they were exact row-scaled
copies of `beta`.

## Full convergence regressions

Four strict `1e-7` end-to-end regressions were run with the factored update.
The iteration counts and final solutions are unchanged. On large10k `T = 3`,
all 116 printed rows (initial row plus 115 iterations) have identical
objectives, primal/complementarity infeasibilities, step sizes, and line-search
backtracks. The only printed trajectory difference is `1e-9` in dual
infeasibility at iteration 105, consistent with floating-point operation
ordering and far below either stopping tolerance.

| System | Horizon | Historical | Factored | Runtime change | Strict iterations | First practical `1e-6` iteration |
|---|---:|---:|---:|---:|---:|---:|
| IEEE123C | 3 | 13.856 s | 15.635 s | +12.84% | 48 | 44 |
| IEEE2522C | 3 | 233.098 s | 245.175 s | +5.18% | 56 | 52 |
| large10kC | 3 | 7123.782 s | 6660.908 s | -6.50% | 115 | 111 |
| IEEE2522C | 12 | 1902.103 s | 1183.039 s | -37.80% | 79 | 76 |

The small cases show that this is not a universal per-iteration speedup:
startup noise and unchanged KKT work dominate when little memory is removed.
The benefit appears when retained horizon storage becomes material. The
large10k `T = 3` run saves about 2.49 GiB across three stages and is 6.5%
faster. IEEE2522 `T = 12` saves about 609 MiB and is 37.8% faster in this
single-run comparison. Because wall times are one historical run versus one
new run, the exact memory reduction and identical convergence are stronger
claims than the precise speed percentages; repeated timings would be needed
for confidence intervals.

Raw traces and the comparison row set are stored as
`factored_bounds_*_trace.csv` and
`factored_bound_sensitivity_regression.csv` in the network results folder.

## Reproduction

After applying `dynamic_network_scaling.patch`, apply
`ddp/patches/factor_bound_sensitivities.patch` to the upstream FilterDDP clone.
The raw before/after row is stored in
`ddp/results/network_filterddp/factored_bound_sensitivity_memory.csv`.
