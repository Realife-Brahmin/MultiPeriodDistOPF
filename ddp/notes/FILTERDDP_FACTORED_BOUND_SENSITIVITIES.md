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

## Reproduction

After applying `dynamic_network_scaling.patch`, apply
`ddp/patches/factor_bound_sensitivities.patch` to the upstream FilterDDP clone.
The raw before/after row is stored in
`ddp/results/network_filterddp/factored_bound_sensitivity_memory.csv`.
