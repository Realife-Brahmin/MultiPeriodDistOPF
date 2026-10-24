# Exact AC dispatch against an SDP outer relaxation: proof of concept (2026-10-21)

Question: at a prescribed angle vector, does a daily AC dispatch exist? Lower bound L: the clean
maximum-margin problem (margin alone, `margin_reg = 0`, `margin_band_kw = 20000`) solved locally
by Ipopt and re-validated against the original constraints. Upper bound U: a sparse SDP relaxation
of the same problem at fixed angles (`feasibility_certificate.jl`), solved by Clarabel 0.11.1
(`envs/msopf_sdp`). Feasible if L >= 0, infeasible if U < -0.01 kW, unresolved otherwise.

Case: IEEE 123-node, sources 1, 3, 4, T = 24, PV reactive power dispatched, hard cyclic SOC.
Directions in (d3 - d1, d4 - d1): 292.5 deg (substation-4 face) and 157.5 deg (substation-3
face); reaches from `simplex_check.csv`. Files: `feasibility_poc_ieee123_T24{,_ibounds,_amp5}.csv`
(all columns) and `.log`; driver `root_level/feasibility_poc.jl poc123` (FEAS_CUTS=ibounds,
FEAS_AMP=5).

## Results

| Batteries | Query | d3, d4 (deg) | L (kW) | Original-model violation | U, relaxation (kW) | U, with 5 pu ampacity (kW) | Class |
|---|---|---|---|---|---|---|---|
| idle | equal | 0, 0 | +77.64 | 4e-10 | 1898.6 | 487.7 | feasible |
| idle | 98% (292.5) | +0.044, -0.107 | +2.73 | 9e-11 | 1899.9 | 487.7 | feasible |
| idle | 102% (292.5) | +0.046, -0.111 | -2.73 | 9e-11 (equations) | 1907.6 | 487.7 | unresolved |
| idle | beyond prediction (292.5) | +0.039, -0.094 | +18.43 | 9e-11 | 1904.5 | 487.7 | feasible |
| idle | 98% (157.5) | -0.091, +0.038 | +2.87 | 4e-14 | 1895.9 | 487.7 | feasible |
| idle | 102% (157.5) | -0.094, +0.039 | -2.86 | 4e-14 (equations) | 1903.6 | 487.7 | unresolved |
| dispatched | equal | 0, 0 | +82.13 | 4e-10 | 1853.5 | 629.3 | feasible |
| dispatched | 98% (292.5) | +0.091, -0.219 | +5.57 | 2e-12 | 1901.6 | 629.3 | feasible |
| dispatched | 102% (292.5) | +0.094, -0.228 | -5.57 | 2e-8 (equations) | 1898.7 | 629.3 | unresolved |
| dispatched | beyond prediction (292.5) | +0.084, -0.203 | +25.45 | 7e-8 | 1905.0 | 629.3 | feasible |
| dispatched | 98% (157.5) | -0.176, +0.073 | +5.57 | 3e-10 | 1887.2 | 629.3 | feasible |
| dispatched | 102% (157.5) | -0.183, +0.076 | -5.57 | 5e-11 (equations) | 1891.5 | 629.3 | unresolved |

U here is the solver's primal = dual objective (gap 2e-3 to 4e-2 kW; primal residuals 1e-13 to
4e-8; primal PSD blocks' smallest eigenvalue >= -2e-8; dual PSD and sign violations 0 after
projection). That is a numerical certificate of the relaxation optimum, within those tolerances.
The safeguarded bound built from the projected multipliers and the variables' box bounds
(valid for any multipliers up to rounding) is 5e6 to 1e8 kW: the box term (residual x box) is
huge because the current and current-product boxes are wide. Either way U >> 0, so no point
is certified infeasible, and the classification does not depend on which U is used.

Feasible points: every L >= 0 dispatch was re-validated (equations to <= 7e-8 pu, no band,
no-backflow, reactive, battery or SOC violation), and the equal-angle and 98% days replay in
OpenDSS to 5e-5 kW, 8e-5 kvar, 6e-10 pu. At the 102% points the equations hold (violation is
the backflow / band only, as the negative margin says).

Times (per query): Ipopt 1-38 s; relaxation build about 1 s (5 s for the first), Clarabel 10-45 s.

## Why the relaxation cannot certify here

* Its optimum is about 1.9 MW at every angle tested (0.5 MW with ampacity), almost independent
  of the angle: it is limited by what it can fake, not by the physics of the angle.
* It books fictitious line losses: at equal angles 4,838 kW per period beyond what its own
  currents imply (sum_l R_l (|I|^2 moment - |I moment|^2)). Losses must be imported, so every
  source's import rises and no-backflow holds for free. Relaxations of power flow can only
  over-estimate losses, which is exactly the direction that helps P_Subs >= 0.
* It also lets |V|^2 exceed the square of the voltage it reports (by up to 0.08), so the band
  room in the margin sits near its maximum (the margin's band ceiling is 2,000 kW).
* Moment blocks are far from rank one (lambda2/lambda1 0.85 to 1.0); recovered states violate the
  original model by 12-20 pu.
* Strengthening 1, valid derived current bounds in the solver: 1.80 MW. The bounds come from the
  tree (each line carries at most the injections on either side), but between sources they reach
  530 pu, leaving 1.9 GW of loss capacity.
* Strengthening 2, line ampacity 5 pu (an assumed value, now an option of the model): 488 /
  629 kW. It caps the fictitious losses at sum_l R_l I_max^2, about 1.1 MW per period; the
  relaxation uses about that much. Better, but the boundary margins are a few kW.
* Validity check: the nonlinear solution, lifted into the relaxation, satisfies it to 4e-10.

## Two bracket rows (2705-bus feeder, 1.03 pu, dispatched, directions 47 and 65)

`check_bracket_rows.csv`: margins at 100 / 100.1 / 100.2 / 100.5% of the recorded reach,
warm-started, are +1.36 / -0.19 / -1.32 / -4.60 kW (47) and +0.76 / -0.23 / -1.13 / -4.67 kW
(65); cold at 100 / 100.1%: +0.83 / -0.06 kW and +1.11 / -0.12 kW. Both are slightly interior
endpoints (the edge lies within 0.1%, about 0.007 deg, beyond them); `binding = none` is the
1 kW labelling threshold. The secant gave up because locally optimal margins at the same angle
differ by about 0.5 kW with the start, as large as a secant step's change near the edge.

## Literature

* Chiang and Jiang, "Feasible Region of Optimal Power Flow: Characterization and Applications",
  IEEE Trans. Power Syst. 33(1):236-244, 2018: the feasible region equals the union of the stable
  equilibrium manifolds of a quotient gradient system (QGS) built from the constraints; the QGS is
  completely stable, so its trajectories converge to feasible components. A constructive way to
  find feasible points and other components; no certificate that a region is empty.
* Li, Chiang, Wang, QGS-based Trust-Tech for distinct feasible solutions in multiple/all
  components of ACOPF (IEEE Trans. Power Syst., 2025; repository litengmu/QGS-TT, MATLAB with
  MATPOWER and IPOPT, single-period cases): finds solutions across components; the repository
  states no completeness guarantee and notes extensive sampling.
* Molzahn, "Computing the Feasible Spaces of Optimal Power Flow Problems", arXiv:1608.00598:
  discretizes inequality constraints, finds all power flow solutions per point by numerical
  polynomial homotopy continuation, and prunes points with a hierarchy of convex (moment)
  relaxations; provably complete within the discretization, shown on 5- and 9-bus systems.
* For a controlled multi-period problem like ours, the QGS / Trust-Tech side applies (equality
  and inequality constraints of any size, with slack variables): it can find more feasible
  components, i.e. better L. Certifying that no dispatch exists needs an outer bound, which is
  where the homotopy-plus-relaxation route stops scaling (a 133-bus, 24-period system is far
  beyond all-solutions homotopy) and where the plain SDP relaxation is weak for no-backflow.

## Follow-up (2026-10-22): two valid strengthenings, one hour

Hour 12 (the limiting hour of the 292.5-deg idle edge), batteries idle, no ampacity, relaxation on
m >= 0. Points 2% inside (validated feasible, L = +2.727 kW) and 2% outside (local margin -2.727 kW)
the single-hour boundary r* = 0.118008 deg. `feasibility_hour1_ieee123.csv` (default Clarabel
settings) and `_reg1e-6.csv` (static regularization 1e-6, 600 iterations); driver mode `hour1`.

* `:kcl0`, current balance at zero-injection banded buses (J = 0, and its lifted form M a = 0 on a
  bus moment block over the voltage and the incident line currents): no measurable change.
* `:jbound`, |J| <= Smax / 0.95 at banded device buses, on the lifted |J|^2: fictitious losses -8%.
* Both: -10% (about 4.9 MW per period left), relaxed margin about 1.75-1.79 MW instead of 1.94 MW.
* The validated inside solution lifts into every variant (<= 4.3e-9). The outside point stays
  relaxation-feasible in every variant (no infeasibility status, so no Farkas ray to validate):
  unresolved. `:jbound` and both did not solve cleanly (slow progress / numerical error); their
  returned points violate the constraints by <= 2e-4 pu with about 1.8 MW of margin.
* Where the slack sits: the largest |V|^2 gaps are at source terminal 3s and its POI bus 120; 60-78%
  of the fictitious losses are on five lines, mostly L115 and L116 into bus 116 (substation 4's POI).
  A source current is (E - V_s)/z_s with |z_s| = 0.006 pu, so a voltage-moment gap of about 0.07 at a
  source terminal permits a current second moment of about 2,000 pu^2, circulating between sources.
  Neither addition touches the source currents; bounding <|E - V_s|^2> at the source terminals from
  the fixed source phasors and the linear network response is the natural next step.
