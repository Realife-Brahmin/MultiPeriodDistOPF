# Cached network constraint Jacobian

## Why the full sensitivity maps remain

After the first seven exact improvements, source tracing showed that `beta`
and `omega` are not retained accidentally. During the backward recursion they
enter `beta' * B`, `omega' * cx`, `beta' * Qu`, and `omega' * c`. During every
line-search rollout they are applied to the state displacement as
`beta * delta_x` and `omega * delta_x`. Because `delta_x` changes along the
rollout and between trial step sizes, replacing the complete maps requires a
new representation of their action or additional KKT solves; simply deleting
columns would change the algorithm.

## Post-optimization profile

A bounded large10kC `T = 3` profile therefore re-ranked the remaining work.
At the two warm stages, approximate per-stage times were:

- first-order derivatives and matrix construction: `2.24 s`;
- local second-order algebra: `1.27 s`;
- sparse KKT assembly and factorization: `0.76 s`;
- the 1021-column KKT solve: `3.35 s`;
- policy/value update: `0.52 s`.

Second-derivative callbacks themselves remained negligible at about `0.016 s`.
Splitting the first-order work showed that almost all of it was reconstruction
of the sparse control-constraint Jacobian `cu`: approximately `1.87--1.91 s`
per warm large10k stage. Constraint values required about `0.09--0.10 s`, the
state Jacobian about `0.0001 s`, and dynamics Jacobians about `0.004 s`.

## Exact network-specific change

For the BFM-SOCP model, `cu` has a fixed sparsity pattern. All entries are
constant except four derivatives per branch in the nonlinear SOC equality.
The network driver now constructs one sparse `cu` matrix for each stage during
model building. Each callback updates only those four iterate-dependent values
per branch and returns the same stage-owned matrix. The equations, derivatives,
KKT system, and FilterDDP iteration are unchanged.

On warm large10k stages, `cu` callback time fell to `0.002--0.003 s`, about a
99.9% reduction. A bounded one-iteration solve fell from `37.497 s` to
`29.987 s`; model building increased from `1.455 s` to `12.165 s` because the
three large sparse structures are created upfront. That one-time cost is
amortized after only a few stages/iterations.

## Full validation

- IEEE123C `T = 3` retained 48 iterations and a byte-for-byte identical
  49-row trace.
- IEEE2522C `T = 12` retained 79 iterations and a byte-for-byte identical
  80-row trace. Practical `1e-6` convergence remains at iteration 76.
- IEEE2522C solve time fell from `701.245 s` to `577.190 s`, a further 17.69%
  single-run reduction. Model building increased from `1.415 s` to `3.340 s`,
  so the end-to-end gain remains essentially the same.
- Relative to the original sparse IEEE2522C runtime of `1902.103 s`, the full
  sequence now gives a 69.66% single-run solve-time reduction.

Timings are individual runs rather than statistical estimates. The next
fundamental target remains the many-right-hand-side sensitivity solve and the
dense `beta`/`omega` maps, not second-derivative evaluation.
