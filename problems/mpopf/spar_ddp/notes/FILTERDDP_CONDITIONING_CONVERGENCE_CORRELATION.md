# Does worsening KKT conditioning cause FilterDDP's slow convergence?

## Method

Stage-1 KKT diagnostics sampled every ten outer iterations were aligned with
the stored strict `T=3` convergence traces. The comparison uses barrier level,
UMFPACK pivot ratio, isolated factor and representative-solve time, relative
solve residual, accepted step size, and recorded line-search backtracks.

This is a small deterministic trajectory study, not a statistical timing
experiment. Isolated timings include normal run-to-run noise. The pivot ratio
is a conditioning proxy, not a formal condition-number estimate. IEEE123 also
has the independently computed equilibrated spectral conditions reported in
`FILTERDDP_KKT_NUMERICAL_DIFFICULTY.md`.

## Results

The barrier and pivot ratio are very strongly associated. Pearson correlation
between `log10(barrier)` and `log10(pivot ratio)` is `0.980` for ieee123 and
`0.988` for ieee2522. As the barrier is reduced, the smallest LU pivot becomes
many orders of magnitude smaller relative to the largest pivot.

This numerical deterioration does not produce a corresponding runtime
explosion:

- ieee123 factor time stays within `0.312--0.319 s`, and the representative
  solve within `0.080--0.089 s`, while the equilibrated condition grows from
  `5.77e4` to roughly `5.38e15`;
- ieee2522 factor time grows from `0.290 s` initially to `0.373 s` at iteration
  50, while representative solve time grows from `0.252 s` to `0.385 s`.
  This 29%/53% increase is measurable, but tiny relative to the roughly
  16-order pivot deterioration.

The stored linear residuals remain at most `2.95e-10` for ieee123 and
`3.20e-11` for the sampled ieee2522 stage-1 systems. UMFPACK therefore remains
accurate on these accepted trajectories.

Most importantly, late ill-conditioning is not associated with filter
rejection in these runs. IEEE123 and ieee2522 record zero backtracks at every
sample; after the initial point their accepted step is consistently `0.5`.
The large10k trace shows the opposite timing from the conditioning hypothesis:
its one-to-four backtracks occur predominantly during iterations 10--75 while
the barrier remains `1`. Once the barrier decreases at iteration 85, accepted
steps are `0.5` and recorded backtracks are zero through convergence.

## Conclusion

Barrier reduction clearly drives the KKT system toward numerical singularity,
and this increases late solve cost somewhat on ieee2522. It is a legitimate
robustness and solver-scaling concern. The tested evidence does **not** show it
causing late filter rejection or the dominant wall-clock burden.

The principal cost remains structural: at every outer iteration FilterDDP
factorizes `T` stage matrices, solves each against `nx+1` right-hand sides, and
propagates the resulting full sensitivity maps. Conditioning is a secondary
multiplier on that workload, not its source.
