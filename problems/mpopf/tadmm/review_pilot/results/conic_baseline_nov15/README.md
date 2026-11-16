# Centralized conic baseline (2026-11-15)

Reviewer items AE-5 and R3-1: the centralized baseline is Ipopt, a general
nonlinear solver, applied to a second-order cone program. `conic_baseline.jl`
builds the same centralized problem with the cone stated explicitly (rotated
second-order cone, every original bound kept) and solves it with Clarabel
(open source) and Gurobi's barrier method, next to Ipopt/MA57 on the same
machine. One thread each unless stated; each solver gets a warm-up solve and is
timed on a freshly built model. Environment: `conic_env/Project.toml`.

## 309 PC (Julia 1.12.6, Gurobi 13.0.3, Clarabel 0.11.1)

Wall time of the solve call, seconds; iterations in brackets.

| case | Ipopt/MA57 | Clarabel | Gurobi barrier |
|---|---|---|---|
| ieee123 T24 | 0.46 (25) | 0.77 (20) | 0.61 (18) |
| ieee123 T96 | 2.42 (29) | 6.60 (28), "almost optimal" | 3.63 (24) |
| med2522 T24 | 17.7 (49) | 28.2 (34), "almost optimal" | 43.5 (51), suboptimal termination |
| med2522 T24, Gurobi with 8 threads | | | 37.2 (51), suboptimal termination |

Objectives agree with Ipopt's to between 2e-10 and 3e-7 relative. On med2522
Clarabel and Gurobi stop short of their 1e-8 tolerances; their dispatches
violate the original model by 1e-6 and 7e-8, Ipopt's by 8e-16.

Neither conic solver is faster than Ipopt/MA57 on these instances: 1.3 to 2.7
times slower on one thread, and Gurobi gains little from eight. The Ipopt
baseline is therefore not what makes centralized look slow. Kamiak runs of
Ipopt/MA57 and Clarabel on the test nodes, including large10k, were queued when
this was written (`~/kamiak_jobs/tadmm_sc_nov15_v12_conic`).

Dispatch files (`*_dispatch.jls`) are not committed.
