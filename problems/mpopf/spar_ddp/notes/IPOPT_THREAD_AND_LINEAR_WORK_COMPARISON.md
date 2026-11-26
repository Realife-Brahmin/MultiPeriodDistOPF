# IPOPT threading and FilterDDP linear-work comparison

## Installed IPOPT/MUMPS capability

The Julia Ipopt artifact on the lab PC loads `MUMPS_seq_jll`. It is therefore
the sequential MUMPS build. Setting `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
and the visible BLAS thread count to 1, 2, 4, or 8 does not turn MUMPS into a
parallel factorizer. A true parallel-MUMPS experiment would require a different
MPI/OpenMP-enabled binary and execution setup.

Matched centralized diagnostic models were run sequentially at nominal thread
counts 1, 2, 4, and 8. Every setting preserved the objective and iteration
count. IEEE2522 `T=12` was fastest at one thread in overall wall time
(`17.525 s`; the nominal four-thread result was `17.960 s`). Large10k `T=3`
also degraded monotonically from `18.441 s` at one thread to `19.796 s` at
eight threads. The linear backsolve did not improve: large10k increased from
`4.978 s` to `5.403 s`. These are capability probes, not repeated statistical
benchmarks, but they decisively show that environment thread settings do not
parallelize the installed MUMPS.

## What FilterDDP repeats

At one stage FilterDDP factorizes its KKT coefficient once, then solves the
same factor against one feedforward column plus `nx` state-sensitivity columns.
It does **not** refactor the coefficient `nx` times. For large10k this is one
`96968 x 96968` sparse LU followed by a dense `96968 x 1021` solve.

The captured large10k stage benchmark measured:

- one UMFPACK factorization: `0.697 s`;
- all 1021 RHS columns in one call: `3.127 s`;
- the same columns in 1021 one-column blocks: `3.249 s`;
- explicit solution copying: about `0.198--0.225 s`.

Thus batching is already effective and result collation is secondary. The
material cost is obtaining all 1021 required directions and propagating them
into the policy/value recursion. At `T=3`, coefficient factorization plus the
wide solve alone is approximately `3*(0.697+3.127)=11.47 s` per outer
FilterDDP iteration before derivative construction, value updates, rollout,
filter evaluation, and other work.

For comparison, the centralized diagnostic large10k `T=3` model has 133035
variables, 95949 equalities, and 30960 inequalities. Ipopt completed 44 Newton
iterations in `18.441 s`; its complete primal-dual linear-system work totalled
`13.758 s`, or about `0.313 s` per Newton iteration. This is not an
apples-to-apples matrix solve--Ipopt factorizes one global space-time KKT matrix
while FilterDDP factors stagewise matrices and requests a wide sensitivity
map--but it exposes the algorithmic contrast. FilterDDP's three stagewise
factor-plus-wide-solve operations alone cost roughly 37 times one complete
Ipopt primal-dual linear-system step in this captured comparison.

## Interpretation

The evidence rejects the simple explanation that “UMFPACK is slow while MUMPS
is fast.” Single-process MUMPS was already slower on captured FilterDDP stage
systems, and the installed Ipopt MUMPS is sequential. Ipopt's advantage is
mainly that it computes one global Newton direction per iteration using a
mature sparse primal-dual formulation. FilterDDP deliberately computes the
response to every temporally coupled state direction at every stage so it can
construct and propagate its value-function model.

The next experiment should correlate the late-iteration conditioning evidence
with factor time, wide-RHS time, accepted step length, and filter rejection.
That will distinguish whether poor conditioning merely accompanies barrier
convergence or actively causes extra work.
