# Where FilterDDP's Hessian-related work becomes bulky

This note addresses the questions raised for the 10 September 2026 meeting.
It separates derivative evaluation, construction of the backward-pass matrices,
factorization, the many-right-hand-side solve, and the value-function update.
The timings are diagnostic measurements, not new convergence runs.

## Short answer

For the network models, **evaluating the Hessian callbacks is not the dominant
cost**. The main burden comes after those derivatives have been evaluated:

1. FilterDDP solves one stage KKT system against `nx + 1` right-hand sides, not
   one direction vector. The extra columns produce every control sensitivity
   with respect to the temporally coupled state.
2. The resulting `beta`, `omega`, bound-dual sensitivities, and `Vxx` are dense.
3. The update of the dense value Hessian, including products such as
   `beta' * B`, is the largest measured backward-pass component.

Thus "Hessian-related computation" is bulky mainly because FilterDDP carries
and propagates **dense second-order sensitivity information**, not because the
individual second-derivative formulas are unusually expensive to evaluate.

## What is differentiated, and how?

FilterDDP's generic constructors use `Symbolics.gradient`,
`Symbolics.jacobian`, and `Symbolics.hessian` in
`src/ocp/{objective,dynamics,constraints}.jl`. Those derivatives are generated
symbolically and compiled; they are not finite differences.

The network experiment bypasses that generic symbolic construction and
provides handwritten analytic sparse callbacks in
`ddp/examples/power_system/ieee123c_filterddp.jl`:

- objective second derivative `luu`: lines 60--65;
- affine battery dynamics and zero dynamics Hessians: lines 32--46;
- constraint Jacobians and weighted Hessian `cuu`: lines 168--240.

The centralized JuMP/Ipopt model also uses exact derivatives. Its affine and
quadratic expressions are differentiated by JuMP/MathOptInterface; Ipopt uses
the supplied exact sparse Lagrangian Hessian. No limited-memory Hessian option
is enabled in the recorded centralized run.

## Warm backward-pass timing

The diagnostic timer was placed around six disjoint sections of each stage.
The figures below sum the second backward pass, avoiding first-use compilation
as far as practical.

| Component | IEEE2522C, T=3 | Share | large10kC, T=3 | Share |
|---|---:|---:|---:|---:|
| All derivative/model callbacks | 0.445 s | 11.5% | 11.254 s | 23.0% |
| Matrix algebra forming `C`, `Hhat`, and `B` | 0.314 s | 8.1% | 4.466 s | 9.1% |
| Sparse KKT assembly | 0.493 s | 12.7% | 1.691 s | 3.5% |
| Sparse LU factorization | 0.263 s | 6.8% | 1.651 s | 3.4% |
| Solve for all `nx+1` right-hand sides | 0.791 s | 20.4% | 10.032 s | 20.5% |
| Policy and dense value-function update | 1.585 s | 40.9% | 19.807 s | 40.5% |

The IEEE2522 rerun further split its warmed callback time into approximately
0.395 s for first-order/model evaluation and only **0.046 s for
second-derivative callbacks**. This makes the distinction concrete: evaluating
the second derivatives is cheap here; using the dense second-order information
is expensive.

The timers are intentionally coarse and do not constitute a full profiler.
They nevertheless locate the order-of-magnitude bottleneck reliably. The
large10k run has `nx = 1020`, so its KKT right-hand side has 1021 columns and
its `1020 x 1020` `Vxx` is fully dense (1,040,400 nonzeros).

## Does translating equations into matrices take substantial time?

Not by itself. On large10k, sparse KKT assembly was about 3.5% of the measured
warm backward pass. The earlier algebra that combines derivatives with `Vxx`
was another 9.1%. Both matter, but neither explains the runtime alone. The
dense many-right-hand-side solve and value update together consumed about 61%.

Model construction is also a separate one-time cost. It should not be confused
with rebuilding and solving the stage systems during every FilterDDP iteration.

## UMFPACK versus MUMPS

Captured stage-1 KKT systems were solved with the existing Julia sparse LU
(UMFPACK) and MUMPS 1.6.2. Median direct timings were:

| Captured system | UMFPACK total | MUMPS unsymmetric | MUMPS symmetric |
|---|---:|---:|---:|
| IEEE123C, `1353 x 1353`, 52 RHS | 0.0065 s | 0.0594 s | 0.0905 s |
| IEEE2522C, `23695 x 23695`, 251 RHS | 0.2842 s | 0.5246 s | 0.4867 s |

All residuals were numerically small. In these single-process tests MUMPS was
slower, although its relative gap narrowed on IEEE2522C. This does not test a
carefully configured multi-process MUMPS run. More importantly, factorization
was only about 3.4% of the measured large10k backward-pass time. Replacing
UMFPACK can therefore offer only a limited total speedup unless the dense RHS
solve and value recursion are changed as well.

The benchmark is reproducible with
`ddp/examples/power_system/benchmark_captured_kkt.jl`; compact result CSVs are
stored beside the network results. The captured matrices themselves are local
diagnostic artifacts because the IEEE2522 capture is about 50 MB.

## Practical research direction

The most promising question is not simply "which sparse factorizer is faster?"
It is whether the MPOPF structure permits FilterDDP to obtain the necessary
backward information without explicitly solving for and retaining all `nx`
sensitivities. Candidate directions include selected right-hand sides,
low-rank or structured representations of `Vxx`, matrix-free products, and
using one or more stagewise Ipopt solves as an oracle for the local optimum and
the particular sensitivity information the recursion actually needs.

## Allocation and retained-memory profile

A second diagnostic pass measured both Julia allocation traffic and the
storage retained by the principal arrays. Allocation traffic is cumulative
memory created during an operation; it is **not** the same as simultaneous
peak RAM. Retained sizes are direct object-storage measurements. UMFPACK's
native internal factor storage is not fully visible through Julia object-size
inspection, so no unsupported factor-memory claim is made.

On warmed IEEE2522 stages:

- sparse `K`: 1.7--2.7 MiB;
- dense RHS and its solution: 45.4 MiB each;
- `beta`: 25.5 MiB; `omega`: 19.7 MiB;
- two bound-sensitivity matrices: 51.0 MiB;
- complete stored update rule: 96.5 MiB per stage;
- temporary allocation traffic while building the update: about 386 MiB per
  stage.

On warmed large10k stages:

- sparse `K`: 7.1--22.9 MiB;
- dense RHS and its solution: 755.3 MiB each;
- `beta`: 425.4 MiB; `omega`: 329.2 MiB;
- two bound-sensitivity matrices: 850.8 MiB;
- `Vxx`: only 7.9 MiB by itself;
- complete stored update rule: 1607.0 MiB per stage;
- temporary allocation traffic while building the update: about 6408 MiB per
  stage.

This sharpens the earlier conclusion. `Vxx` is important because it causes
dense information to spread, but storing `Vxx` itself is not the main RAM
consumer. The dominant retained storage consists of the tall `nu x nx` and
`nc x nx` sensitivity matrices, especially the bound sensitivities. The raw
profile rows are in
`ddp/results/network_filterddp/backward_pass_memory_profile_T3.csv`.
