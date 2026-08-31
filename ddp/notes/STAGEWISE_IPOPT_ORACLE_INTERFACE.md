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
