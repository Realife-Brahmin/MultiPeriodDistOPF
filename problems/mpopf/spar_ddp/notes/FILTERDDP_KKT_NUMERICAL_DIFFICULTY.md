# FilterDDP KKT numerical difficulty across convergence

## Question

Does FilterDDP become slow because its stagewise coefficient matrix becomes
indefinite, nearly singular, or expensive to factor as the network and horizon
grow?

The matrix is the constrained backward-pass KKT system. It is symmetric and
indefinite **by construction**: a well-posed saddle-point system should have
`nu` positive and `nc` negative eigenvalues. Indefiniteness alone is therefore
not a defect. The more relevant failure modes are a zero/near-zero eigenvalue,
extreme coefficient scaling, unstable pivots, and sparse-factor fill-in.

## Diagnostic

`analyze_kkt_numerics.jl` reports KKT dimension and density, UMFPACK LU fill,
the ratio of the smallest to largest absolute diagonal pivot in `U`, and the
relative residual of a representative solve. For ieee123 only, where a dense
eigendecomposition is affordable, it first equilibrates the matrix by a
diagonal congruence scaling and then reports inertia and spectral condition.
The congruence scaling preserves inertia.

Full cold-start `T=3` solves were sampled every ten iterations for ieee123 and
ieee2522. A one-iteration large10k run establishes its initial matrix facts.
The snapshots are diagnostic temporaries; the compact measurements are stored
in `kkt_numerics_over_iterations.csv`.

## Findings

1. **The intended indefiniteness is present.** At ieee123 iteration 0 the
   equilibrated stage-1 matrix has inertia `(791 positive, 562 negative, 0
   near-zero)`, exactly `(nu,nc,0)`. Thus “indefinite” does not explain a
   failure; it is the expected KKT structure.

2. **The matrix approaches numerical singularity as the barrier shrinks.** At
   ieee123 stage 1, the equilibrated spectral condition grows from `5.77e4`
   at iteration 0 to `3.24e10` at iteration 20 and roughly `5.38e15` at
   iteration 40. One eigenvalue becomes numerically unresolved by iteration
   30. The LU pivot ratio simultaneously falls from `2.72e-4` to `1.17e-20`.

3. **The same pivot collapse occurs at network scale.** For ieee2522 stage 1,
   the pivot ratio falls from `1.45e-6` at iteration 0 to `2.51e-22` at
   iteration 50. This is strong evidence of late-iteration ill-conditioning,
   even though a full inertia/condition calculation is deliberately avoided
   for the `23695 x 23695` matrix. The initial large10k stage-1 pivot ratio is
   `1.52e-5` for its `96968 x 96968` KKT matrix.

4. **UMFPACK remains numerically successful in these runs.** All sampled
   representative solves have relative residual at most `1.35e-9`; most are
   `1e-11` or smaller. No backward regularization was invoked in the existing
   periodic traces. The matrix is therefore difficult but not failing on the
   tested trajectories.

5. **Fill-in and per-solve cost do not explode with the conditioning.** The
   ieee2522 stage-1 LU fill ratio rises only from `2.22` to `2.94`. Isolated
   factor/backsolve timings rise modestly, not by the many orders of magnitude
   seen in the pivot ratio. These timings are diagnostic single measurements,
   not benchmark estimates.

## Interpretation

There are two distinct scaling problems:

- **Numerical robustness near convergence:** barrier and bound terms create
  extreme scaling and a nearly singular KKT system. Better scaling,
  regularization, or an inertia-aware symmetric-indefinite solver may reduce
  fragility and possibly avoid late iterations or rejected steps.
- **Total computational work:** FilterDDP still performs one stagewise KKT
  factorization and a many-right-hand-side sensitivity solve at every stage of
  every outer iteration. Earlier profiling shows that, on large10k, the dense
  sensitivity solve and value/policy propagation dominate while sparse LU
  factorization itself is only a small fraction. Conditioning alone therefore
  does not explain the current wall time.

The next useful test is to correlate barrier value, pivot/conditioning proxy,
factor time, multi-RHS time, and filter acceptance for every sampled iteration.
That can determine whether the late ill-conditioning causes extra iterations
or line-search rejection, rather than merely coexisting with convergence.
