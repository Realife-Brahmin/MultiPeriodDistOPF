# Factor-backed FilterDDP policy actions

## Exact reformulation

FilterDDP normally retains the complete dense feedback maps `beta` and
`omega` at every horizon stage.  A rollout uses them only through
`beta * delta_x` and `omega * delta_x`.  The same actions can be recovered
from the already-factorized stage KKT system:

```text
K * [beta * delta_x; omega * delta_x]
    = -[B * delta_x; cx * delta_x].
```

For the guarded network path, `B` is already represented by its exact active
`nx x nx` block and `cx` is sparse.  The optional factor-backed policy stores
the UMFPACK factor, that compact block, `cx`, and two work vectors instead of
the dense maps.  The backward pass still computes all `nx+1` KKT columns,
because the full second-order value update remains necessary; the change is
to what survives for subsequent line-search rollouts.

The mode is enabled with `FILTERDDP_FACTOR_BACKED_POLICY=1`.  Other paths keep
the original dense maps.

## Action validation

On a captured IEEE123 stage, four state directions reproduced the stored map
actions with maximum relative errors `1.60e-14` for `beta` and `1.49e-13` for
`omega`.  On a captured large10k stage, the corresponding errors were
`2.44e-15` and `2.33e-14`; the maximum KKT residual was `6.71e-14`.

For the large10k stage, the two dense maps occupy `754.603 MiB`, while Julia's
reported size for the retained factor is `23.763 MiB`.  A factor-backed action
took `0.00862 s`, versus `0.03748 s` to read and multiply both dense maps in
the isolated benchmark.

## Full regressions and tradeoff

- IEEE123 `T = 3` preserved the complete 49-row trace byte-for-byte and
  converged in 48 iterations.
- A one-iteration large10k `T = 3` profile reduced the retained stage update
  from `757.012 MiB` to `34.888 MiB` (95.39%).  Solver construction fell from
  `1.538 s` to `0.476 s`.  The bounded solve was `28.109 s`, versus `26.589 s`
  for the dense-map active-row version.
- IEEE2522 `T = 12` preserved the complete 80-row trace byte-for-byte and
  converged in 79 iterations.  An initial run while the PC was in use took
  `507.278 s`, versus `481.518 s` for an earlier dense-map run.  A subsequent
  sequential idle-machine comparison put dense and factor-backed modes under
  the same conditions: `523.357 s` dense versus `474.091 s` factor-backed,
  making factor-backed 9.41% faster.  The mode avoids preallocating
  `542.336 MiB` of dense beta/omega maps across this horizon, although retained
  native factor memory means this number is not the net process-RAM reduction.

This is therefore an exact low-memory mode with promising runtime behavior,
not merely a memory fallback.  The controlled evidence supports making it the
preferred network mode, while retaining the dense path as a compatibility
reference.  It does not reduce the backward pass's `nx+1` sensitivity solves;
reducing those columns remains the deeper algorithmic target.
