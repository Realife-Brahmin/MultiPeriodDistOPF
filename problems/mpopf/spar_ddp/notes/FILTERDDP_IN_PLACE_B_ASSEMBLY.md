# In-place mixed-derivative assembly

## What was wasteful

After the typed-constraint fix, a warm large10k stage still allocated about
`2277.949 MiB` while forming the local second-order matrices. The largest
source was the mixed state/control matrix

`B = lux + ux_tmp * fx + fux + cux`.

Here `ux_tmp * fx` is a dense `nu x nx` product, while `lux`, `fux`, and
`cux` are sparse matrices of the same shape. Julia's chained `+` expressions
created a fresh dense matrix for each intermediate sum. Those copies changed
no mathematics; they only moved the same large array repeatedly.

## Exact change

The backward pass now forms the dense product once and adds the three sparse
terms into that owned result:

```julia
B = ux_tmp * fx
B .+= lux
B .+= fux
B .+= cux
```

The last addition remains conditional on equality constraints. No matrix entry,
derivative, elimination, or factorization is approximated.

## Validation

- The strict IEEE123C `T = 3` run retained 48 iterations and produced a
  byte-for-byte identical 49-row CSV trace, including the initial row.
- The strict IEEE2522C `T = 12` run retained 79 iterations and produced a
  byte-for-byte identical 80-row trace. Practical `1e-6` convergence remains
  at iteration 76. Runtime fell from `737.043 s` to `701.245 s`, a further
  4.86% single-run reduction.
- In warm large10kC stages, total algebra allocation fell from `2277.949 MiB`
  to `1001.743 MiB` (56.02%).
- The `B` phase fell from `1276.207 MiB` to `425.402 MiB` (66.67%), and the
  later constraint-algebra phase fell from `450.971 MiB` to `25.569 MiB`
  (94.33%).

These are allocation-traffic figures, not simultaneous resident-memory peaks.
The remaining large allocations are principally the one necessary dense
`nu x nx` product and the Hessian-related stage algebra. Relative to the
original sparse IEEE2522C `T = 12` runtime of `1902.103 s`, the complete exact
optimization sequence now gives a 63.13% single-run reduction. Timings are not
statistical estimates.
