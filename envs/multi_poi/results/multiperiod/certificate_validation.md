# Certificate validation under exploratory current limits, IEEE 123 hour 12 (2026-10-22)

Driver `root_level/certificate_validation.jl` (envs/msopf_sdp, Clarabel 0.11.1; `CV_EXTRA=1..4` for the
follow-up cells). Library changes in `feasibility_certificate.jl` (`dual_bound` with the corrected PSD
pairing and per-term detail, `tighten_boxes!`, restoration mode) and `multiperiod_full_angle.jl`
(`margin_band_hard`). Logs `certificate_validation{,_extra,_extra2,_extra3,_extra4}.log`; CSVs
`certificate_validation_{cases,contributors,restoration,envelope_runs,envelope,envelope_v2}.csv`.
The failed first attempt (wrong PSD pairing) is kept as `certificate_validation_halved_pairing_partial*`.

Fixed throughout: hour 12, batteries idle, PV reactive power dispatched, substation 1 the reference.
A test point is `delta3 = r cos(theta)`, `delta4 = r sin(theta)` in the plane of the two angle
differences: `theta` is a direction on the window map, `r` the size. Inside and outside are 98% and
102% of the hour's boundary `r*` (292.5 deg: `r* = 0.118008` deg; 157.5 deg: `r* = 0.100075` deg).
Exact margins `L`: about `+2.73` / `-2.73` kW (292.5) and `+2.87` / `-2.87` kW (157.5).

## Terms

* **Validated feasible:** an exact AC dispatch (Ipopt) passing `validate_dispatch`: equations <= 1e-9,
  band, no-backflow, and every added current limit.
* **Validated infeasible under the limits:** the independent check `dual_bound` turns the solver's
  multipliers into a bound that holds for any multipliers. A Farkas bound > 0, or an objective
  bound < 0 (no-backflow margin) or > 0 (restoration `rho`), proves the *relaxation* infeasible, so the
  exact model **with the added limits** is infeasible. It says nothing about the model without them.
* **Solver-only:** Clarabel's status, the primal objective, or currents and losses taken from an
  infeasible or failed solve. None of these are used as results.

## Comparison (current envelope, v2)

| Case | Point | Exact | Relaxation status (solver) | Independently checked | Outcome |
|---|---|---|---|---|---|
| A: sources+lines 1.2x, no nodal additions, 292.5 | inside | L = +2.727, valid | optimal, U = 4.56 kW | U <= 4.561 kW | feasible |
| | outside | L = -2.728 | NUMERICAL_ERROR | margin bound U = -1.5e14 kW from that solve's multipliers; restoration `rho >= 0.891` kW | **infeasible under the limits** (see caveat) |
| B: sources+lines 1.5x, both nodal additions, 292.5 | inside | +2.727, valid | optimal, 5.09 kW | U <= 5.103 kW | feasible |
| | outside | -2.728 | infeasible | Farkas +653 (value 6.60e5, box term 6.59e5); `rho >= 0.353` kW | **infeasible under the limits** |
| A at 157.5 | inside | +2.866, valid | optimal, 4.35 kW | U <= 4.348 kW | feasible |
| | outside | -2.865 | infeasible | Farkas +6,919 (value 3.31e5); `rho >= 1.381` kW | **infeasible under the limits** |
| B at 157.5 | inside | +2.866, valid | optimal, 4.48 kW | U <= 4.489 kW | feasible |
| | outside | -2.865 | infeasible | Farkas +140.5 (value 8.76e3); `rho >= 1.240` kW | **infeasible under the limits** |

The exact restoration problem gives `rho = 2.7275` kW (292.5) and `2.8646` kW (157.5), validated: the
cheapest exact dispatch backflows by that much. The relaxation's proven lower bounds are 13-48% of it.

Each inside dispatch, lifted into every relaxation (rank one), satisfies all constraints to <= 4.3e-9
and lies inside every tightened box (violation 0). This was checked after every reformulation.

Earlier envelope versions gave the same classification:

| Envelope | A outside, 292.5 | B outside, 292.5 | A / B outside, 157.5 |
|---|---|---|---|
| frozen (attainable_currents) | Farkas +601.5; `rho >= 0.943` | Farkas +8,913; `rho >= 0.615` | not run |
| v1 (first challenge) | Farkas +601.5; `rho >= 0.908` | Farkas +9,306; `rho >= 0.522` | U = -6.5e16 (NUMERICAL_ERROR multipliers), `rho >= 1.392` / Farkas +138.3, `rho >= 1.330` |
| v2 (second challenge) | as in the table above | as above | as above |

**Caveat, A at 292.5 on v2:** Clarabel stopped with NUMERICAL_ERROR. The bound is still valid, since
`dual_bound` holds for any multipliers, but those multipliers come from a failed solve. The restoration
problem confirms the result independently (optimal solve, `rho >= 0.891` kW > 0). Primal values from
that solve (1,240 kW) are not a solution.

**Numerical uncertainty:** each bound already subtracts a rounding allowance on its accumulated terms
(<= 1e-9 kW). Tight-box and loose-box bounds differ by <= 0.02 kW (restoration) and agree to four
digits (Farkas). The Farkas margins are thin relative to the certificate's size: bound/value = 0.1%
(B, 292.5) to 2.1% (A, 157.5). They are positive but not robust to large perturbations of the ray.
The restoration lower bounds (0.35-1.38 kW against box terms of 4-18 kW) are the more readable
evidence.

## What the verifier needed: a convention bug, not loose boxes

* The first attempt halved the off-diagonal PSD dual. It treated MOI's `PSDConeTriangle` dual as
  scaled-svec, whereas it is the upper triangle of the dual matrix `Z`, paired with the primal by the
  Frobenius product (weight 2 off the diagonal).
* That bound was valid but very loose: the toy `min x, [1 x; x 1/4] PSD` gave -5.5 against the
  optimum -0.5, and the two toy Farkas cases gave -28.9. Every real certificate failed (box terms
  1e8-1e10 against values 1e5-1e6). This was not a property of the problem.
* `conventions()` now checks the pairing first on toys. The optimum is -0.5 (dual bound -0.500000),
  and both Farkas toys give +1.0. The raw triangle dual matches JuMP's dual matrix. Multiplier signs
  are projected onto the dual cones, and the violation is reported (0 in every validated case).
* With the pairing fixed, **tightening the boxes barely matters**: A and B certify with the loose
  boxes too, and the restoration bounds move by <= 0.02 kW.
* The box tightening used:
  - source discs `|V_s - E_s| <= |Z_s| I_s,max`;
  - path discs by the triangle inequality along the tree (`|V_i - V_j| <= |z_l| I_l,max`);
  - interval products for every moment variable.

  These follow from the restricted model for every feasible point, not from sampled dispatches.
  Scaling and PSD null-direction removal were not needed. No second conic solver was available
  (none installed), so there is no cross-solver check.

**Largest box terms** (`certificate_validation_contributors.csv`): term = |residual coefficient| x
half-width of the variable's box, where the residual is the uncancelled part of the multipliers.

* **Farkas, A and B outside:** almost all of the box term sits on PV reactive powers: `Q_PV[pv67]`
  (box +-0.108 pu; residual coefficient 3.0e5 in A, 8.2e6 in B), `pv116` (+-0.015), then `pv3`, `pv26`,
  `pv35` and `pv18` (+-0.031).
  - Physical meaning: the ray weights the reactive-power balance at those PV buses without fully
    cancelling the inverter variables, and the inverter capability `sqrt(S^2 - P^2)` absorbs the rest.
    The certificate pays for not knowing each inverter's reactive output.
  - By constraint kind, the certificate value is a near-cancellation of the inequality rows (+3.7e8)
    against the PSD constant entries (-3.7e8). The net 1.9e5 survives a box term of 1.89e5.
* **Margin and restoration bounds:** the largest terms are the source-terminal voltages `e[1s]`,
  `e[3s]`, `e[4s]` (boxes about 1.044-1.056 pu from the source discs, coefficients 1e-7 to 3e-6), then
  `Wee` at a source and `Q_PV[pv67]`. All are below 0.01 kW.

**Caps carrying the rays** (share of the line-cap multipliers' total weight):

* A: l115 15%, l110 11%, l116 9%, l112 8%, l61 5%, l114, l118, l2 (about 3% each). These are the
  lines into substation 4's POI, the fictitious-loss lines of the earlier experiment.
* B: l71 23%, l67 8%, l100 and l107 (5% each), then l125, l76, l120, l12.

## Challenging the empirical limits

No current limits were imposed in this search. Source angles 3 and 4 were free in [-5, +5] deg
(subs1 = 0), and every source and every line was maximized from several feasible starts: cold, the
inside dispatch, the earlier current maximizers, and cross-seeded from each other's maxima. Every
validated dispatch records every element's current (`certificate_validation_envelope_runs.csv`,
1,836 dispatches).

* **First challenge:** 7 elements exceeded the frozen envelope, 3 of them by more than 1.2x. l120,
  l121 and l125 reached 29.97-30.00 A against 24.83 A (x1.207). Envelope v1 =
  `certificate_validation_envelope.csv`.
* **Second challenge**, on the 11 influential caps not maximized before (l107, l100, l71, l61, l120,
  l76, l125, l67, l12, l114, l112):
  - 5 elements exceeded v1: l125 29.97 -> 35.81 A, l120 29.98 -> 35.82, l121 30.00 -> 35.83,
    l71 30.18 -> 32.40, l61 44.12 -> 44.81.
  - l120 and l125 came within 0.5% of 1.2x of v1.
  - None exceeded 1.2x of v1. Envelope v2 = `certificate_validation_envelope_v2.csv`.
* Every comparison was rerun on each envelope version, with the same classification (tables above).

**These are empirical maxima of a nonconvex problem, not upper bounds.** Each round of search raised
some lines by up to 20%, so the envelope is not settled. The 1.2x cushion was nearly used up on two
lines in a single round. Failing to find a larger current is evidence, not proof.

## Sensitivity: cushion, sources versus lines

Outside point, 292.5 deg, envelope v1.

| Limits | Base relaxation | With zero-injection and device-current blocks |
|---|---|---|
| sources only 1.2x | margin 158.4 kW, unresolved | 0.27 kW, unresolved |
| sources only 1.5x | 195.3 kW, unresolved | 2.33 kW, unresolved |
| sources only 2.0x | 257.5 kW, unresolved | 6.70 kW, unresolved |
| sources + lines 1.2x | **certified** (Farkas +601.5; `rho >= 0.908`) | not needed |
| sources + lines 1.5x | 0.49 kW, unresolved | **certified** (Farkas +9,306; `rho >= 0.522`) |
| sources + lines 2.0x | 3.54 kW, unresolved | 1.17 kW, unresolved |

* Source limits alone never certify, at any cushion. Restoration lower bounds are <= 0 in every
  unresolved case.
* Line limits do the work. Their effect fades between 1.2x and 2x: at 2x nothing certifies, even with
  the nodal blocks.

## Recommendation

1. **Validated feasible AC dispatches** are solid: every inside point, at both directions and on every
   envelope.
2. **Validated infeasibility holds only for the model with the empirical line and source limits
   added**, at 1.2x (base relaxation) or 1.5x (with nodal blocks). Four outside points are certified,
   each confirmed twice (Farkas or margin bound, plus restoration `rho > 0`). One certificate relies
   on multipliers from a NUMERICAL_ERROR solve, which the restoration problem backs.
3. **Solver-only indications**: Clarabel statuses, and the 2x and sources-only cases, which stay
   unresolved.
4. **Empirical maxima versus proven bounds:** the limits are cushioned empirical maxima that rose by
   up to 20% per search round. They are not proven bounds, so these certificates do **not** prove the
   original model (without limits) infeasible. They become a statement about the real feeder only
   once real ampacities replace the empirical envelope (the user is collecting them), or once a
   proven current bound is available.
5. Next steps, if the user approves:
   - Repeat with the real ampacities when they arrive (the same driver takes a Dict by line).
   - Otherwise, keep the 1.5x configuration with the nodal blocks, which has more cushion over the
     envelope.
   - Do not move to the cyclic-battery day until the envelope stops moving.
