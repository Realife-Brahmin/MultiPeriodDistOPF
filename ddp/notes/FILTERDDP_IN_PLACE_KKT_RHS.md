# In-place sparse-KKT right-hand-side solve

## Change

The sparse FilterDDP backward pass previously evaluated `F \\ rhs`. Julia
therefore retained the dense `rhs` and allocated a second equally sized dense
matrix for the solution. The new path calls `ldiv!(F, rhs)` and treats the
overwritten right-hand side as the solution. This changes only buffer ownership,
not the KKT equations or factorization.

KKT-capture mode is deliberately exceptional: it copies the unsolved right-hand
side before `ldiv!` so captured systems remain reproducible. Ordinary runs make
no such copy.

## Validation

- A full IEEE123C, `T = 3`, strict regression followed exactly the established
  48-iteration trajectory and ended at objective `2808.924770488853`, primal
  infeasibility `6.227371345322e-8`, dual infeasibility
  `2.303870314635e-8`, and complementarity `1.624842969030e-9`.
- A captured system retained the original right-hand side and reproduced its
  stored solution with residual `3.1904825916612145e-10`.
- A bounded one-iteration large10kC, `T = 3`, memory probe reduced warm-stage
  solve allocation from `1510.686 MiB` to `755.343 MiB`: exactly one
  `755.343 MiB` dense matrix, or 50%. The update-rule allocation remains
  `3525.163 MiB`, as expected.

The one-iteration probe took 55.949 s, but it is not a runtime benchmark. The
established result is the exact removal of the second dense KKT buffer. Full
timing should only be repeated if a publication-quality runtime number is
needed.

Raw memory comparisons are in
`ddp/results/network_filterddp/in_place_kkt_rhs_memory.csv`.
