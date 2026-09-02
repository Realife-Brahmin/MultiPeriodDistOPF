# KKT right-hand-side block benchmark

## Question

The post-optimization large10k profile identified the `nx + 1 = 1021`
right-hand-side sparse KKT solve as the largest remaining warm-stage operation.
This experiment asks whether solving those columns in smaller blocks improves
cache behavior or reduces UMFPACK solve time.

## Method

One exact large10kC `T = 3`, stage-1 system was captured after the eight exact
implementation improvements:

- KKT dimensions: `96968 x 96968`;
- right-hand side: `96968 x 1021`;
- factorization: one shared UMFPACK LU;
- block widths: 1, 4, 16, 32, 64, 128, 256, 512, and 1021;
- three repeats, reporting median solve time.

The benchmark times `ldiv!` separately from copying each source block. Every
layout produced the same solution norm. Two representative columns had relative
residual `5.08e-15`.

## Result

The current all-at-once layout is fastest: `3.127 s`. Other block widths took
`3.192--3.384 s`, except the narrow cases, which were still about `3.24--3.25 s`.
Width 128 was the best smaller-block result at `3.192 s`, 2.08% slower than the
full solve. Therefore RHS blocking is rejected as a runtime optimization.

Total solve allocation is approximately `755.34 MiB` for every layout. This is
expected: FilterDDP still computes all 1021 columns. Blocking changes the peak
size of a temporary workspace, not the total data processed. If the solved
blocks were copied directly into `beta` and `omega`, width 128 could reduce the
temporary workspace from about `755 MiB` to about `95 MiB`, at a small measured
solve-time penalty. That optional memory tradeoff is not enabled because it
would complicate the production path without improving the present runtime
objective.

## Consequence

UMFPACK is already handling the dense multi-RHS solve efficiently. The next
speed improvement cannot come from merely changing column batch size. It must
reduce the number of sensitivity columns, exploit additional structure in the
right-hand sides/solutions, or change how the required `beta` and `omega`
actions are represented while preserving both the backward recursion and
line-search rollout.
