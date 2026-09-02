# What a stagewise Ipopt oracle would have to provide

The proposed experiment is to let Ipopt solve a single-period network OPF and
retain FilterDDP's backward convergence mechanism. This note defines the
minimum interface before implementing another solver path.

## What Ipopt gives directly

For a fixed entering battery-energy state, an ordinary JuMP/Ipopt stage solve
can return:

- the optimal time-local network and battery variables;
- equality-constraint multipliers;
- lower- and upper-bound multipliers;
- the exact sparse constraint Jacobian and Lagrangian Hessian used by the NLP.

These outputs can replace the task of locating the stage optimum and evaluating
its local first- and second-order equations. They do not, by themselves,
replace FilterDDP's complete backward message.

## What FilterDDP additionally needs

The preceding stage needs to know how the current stage's optimum and duals
change when each component of the entering battery state changes. In the
current implementation those derivatives are the dense matrices `beta`,
`omega`, and the two bound-dual sensitivities. They are obtained by solving
the stage KKT system against `nx` sensitivity columns in addition to the one
feedforward column. Those sensitivities then produce the value gradient and
value Hessian passed backward.

Standard JuMP/Ipopt output does not expose that full parametric sensitivity
map. Merely replacing the local nonlinear solve with Ipopt would therefore
still leave one of the following tasks:

1. solve the linearized KKT system for all required state directions;
2. obtain equivalent sensitivities through an Ipopt sensitivity facility;
3. prove that MPOPF needs only selected columns or combinations of columns;
4. represent the map in a structured, low-rank, or matrix-free form.

## Smallest useful experiment

The first experiment should use IEEE123C at `T = 3` and one fixed backward
stage:

1. solve the stage NLP with JuMP/Ipopt at the same entering state and barrier
   point used by FilterDDP;
2. compare the primal solution and all multipliers;
3. extract the exact Jacobian and Lagrangian Hessian;
4. reconstruct the same stage KKT coefficient;
5. compute selected sensitivity columns for one battery, one feeder region,
   and the aggregate all-battery direction;
6. compare those columns with FilterDDP's `beta` and `omega`.

This test answers whether Ipopt can serve as a trustworthy local oracle before
attempting a full hybrid algorithm. It also reveals whether selected MPOPF
directions contain enough structure to avoid the complete dense `nx`-column
solve. A full-horizon hybrid should not be implemented until this equivalence
test passes.

## IEEE123 stage-1 result

The experiment was implemented in
`ddp/examples/power_system/stagewise_ipopt_oracle.jl` and run on stage 1 of the
converged IEEE123C, `T = 3` trajectory. The Ipopt stage received the same fixed
entering battery state and the same incoming quadratic future-cost message
(`future_Vx`, `future_Vxx`) captured from FilterDDP. This is essential: an
ordinary isolated one-period OPF without that message is a different problem.

For three state directions, complete Ipopt stage solves at `x +/- h*d` were
used to estimate the parametric sensitivities independently. Results were
stable for `h = 1e-4`, `1e-5`, and `1e-6`. At `h = 1e-6`:

| Direction | Relative error in `beta*d` | Relative error in `omega*d` |
|---|---:|---:|
| First battery | 0.994% | 0.165% |
| Middle battery | 2.665% | 0.040% |
| Unit aggregate battery direction | 2.742% | 0.049% |

The perturbed Ipopt solves satisfy their stage equalities to about `2.8e-12`.
JuMP's equality-dual sign is opposite to FilterDDP's convention and is mapped
accordingly. The raw equality-multiplier levels should not be compared: the
slack transcription is highly dual-degenerate/ill-scaled, with FilterDDP's
captured multiplier infinity norm around `3.9e6` versus Ipopt's `1.2e3`.
Despite that non-uniqueness in multiplier levels, their directional derivative
`omega*d` agrees closely.

The captured FilterDDP point precedes its final forward update. Its Ipopt stage
re-optimization moves by `1.69e-3` in infinity norm, while FilterDDP's predicted
feedforward step is `1.35e-3`; the displacement differs from that prediction by
`3.37e-4`. This explains why exact primal equality is not expected in this
snapshot.

## Interpretation

The equivalence test passes in the useful sense: when supplied with the same
future-cost message, Ipopt reproduces FilterDDP's local control and equality-
dual sensitivities to a few percent or better. But finite-differencing complete
Ipopt solves is only a validation method. Recovering every column this way
would require `2*nx` additional nonlinear solves per stage—102 on IEEE123 and
2040 on large10k—which is plainly not the replacement algorithm.

The next implementation question is therefore narrow: can Ipopt's linearized
KKT system be reused to return selected sensitivity products, or can MPOPF be
shown to need only a structured subset? The experiment supports an Ipopt-based
local oracle, but it also confirms that obtaining the backward sensitivity
message remains the central computational problem.

The follow-up investigation is recorded in
`IPOPT_SENSITIVITY_REUSE_AND_MAP_COMPRESSION.md`: ordinary Ipopt.jl does not
expose factorization reuse, the bundled sIPOPT library is a viable integration
target, and a simple low-rank/selected-column approximation of `beta` is not
accurate enough on IEEE123.
