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
- A full strict IEEE2522C, `T = 12`, cold-start run preserved the complete
  79-iteration trace exactly. It finished in `1092.507 s`, versus
  `1145.392 s` for the immediately preceding no-copy implementation: a
  `52.885 s` or 4.62% reduction. Relative to the factored-only run
  (`1183.039 s`), the combined reduction is 7.65%. The final objective was
  `8512.516760782939`, with primal infeasibility `2.290052331922e-8`, dual
  infeasibility `8.669033484674e-8`, and complementarity
  `6.803853214214e-10`. All three practical `1e-6` criteria were first met at
  iteration 76.

The runtime percentages are single-run measurements, not statistical
estimates. The exact results are the identical trajectory and removal of the
second dense KKT buffer.

Raw comparisons and the full trace are in
`ddp/results/network_filterddp/in_place_kkt_rhs_memory.csv`,
`ddp/results/network_filterddp/in_place_kkt_rhs_runtime.csv`, and
`ddp/results/network_filterddp/in_place_kkt_rhs_ieee2522C_1ph_T12_trace.csv`.
