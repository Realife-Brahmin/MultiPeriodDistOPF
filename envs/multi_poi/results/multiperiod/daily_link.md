# Dispatch-dependent network response in the daily relaxation (2026-10-22)

**Setup.** IEEE 123, sources 1, 3 and 4, T = 24, PV reactive power dispatched, batteries dispatched
(lossless), hard cyclic SOC. Restoration objective min rho with P_s,t >= -rho, band and device
limits hard, rho <= rho_exact + 1 kW. No empirical caps or ampacity.

**Baseline.** The corner-enclosure relaxation of `daily_refined.jl`: derived per-hour source-current
bounds, both nodal strengthenings, scaled source rows.

**Results are pointwise.** The daily reaches r* (0.241513 deg at 292.5 deg, 0.194598 deg at 157.5
deg) are locally measured.

## Answer

* **Keeping the dependence of the source currents on the battery dispatch improves certification
  substantially.**
  - At 1.02 r*, the checked lower bound rises from 0.333 to **2.460 kW** at 292.5 deg, and from
    unresolved (-0.461) to **2.234 kW** at 157.5 deg. The exact rho is 5.57 kW at both.
  - At 1.10 r*, it rises from 22.59 / 21.23 to **24.77 / 24.57 kW** (exact 27.86 / 27.84).
  - Every positive bound is confirmed by an exact rational re-check of the exported certificate.
* **The remaining weakness is mainly the size of the link's error allowance eps.**
  - The bound falls steeply and smoothly as eps grows: at 157.5 deg it goes 2.23, 1.85, 1.28, 0.52
    and -0.48 kW for eps x 1, 1.1, 1.25, 1.5 and 2.
  - Validated real dispatches use at most 0.47 of eps.
  - The relaxation exploits the rest through moment freedom: the lifted source-current second moments
    stay 0.62-0.64 pu^2 above the squares of their first moments, and fictitious line losses stay at
    96-134 kWh a day.
  - A valid but tighter eps is the direct lever. Larger moment blocks are not yet needed to explain
    the gap.

## Derivation, checked against the model and the code

* Device current: J_b = conj(S_b) / conj(V_b), where S_b is the device power at bus b. In the model
  (`relax_margin`, `validate_dispatch`) that is:
  - PV real power at the hour's PV factor (fixed);
  - minus the load at the hour's load factor (fixed);
  - plus the battery variables P_B, and j times the PV reactive-power variables Q_PV.

  Batteries are lossless, so S_b is affine in the model's variables.
* From the elimination (`current_bounds.jl`, derivation in `source_current_bounds.md`):
  I_s = a_s + sum_b c_sb J_b exactly, with a_s = (E_s - V0_t(s)) / z_s and c_sb = -Z_t(s),b / z_s.
* With |V_b - V0_b| <= rho_b and |V_b| >= vmin_b > 0, write
  J_b - conj(S_b)/conj(V0_b) = conj(S_b) conj(V0_b - V_b) / (conj(V_b) conj(V0_b)). Hence
  |J_b - conj(S_b)/conj(V0_b)| <= Smax_b rho_b / (vmin_b |V0_b|), and summing over b gives

      | I_s - a_s - sum_b (c_sb / conj(V0_b)) conj(S_b) | <= eps_s = sum_b |c_sb| Smax_b rho_b / (vmin_b |V0_b|).

* Every quantity is the one `current_bounds(...; corners = true)` already proves, per hour:
  - rho_b: the final disc radius of its iteration;
  - vmin_b = max(0.95, |V0_b| - rho_b): computed in the last iteration from the previous rho_b,
    hence smaller than with the final one, still valid;
  - Smax_b: the largest |S_b| over the bus's power rectangle, attained at a corner.

  These bound every state allowed by the band, the device limits and the network equations, over
  every allowed battery and PV dispatch. Nothing is conditioned on a schedule.
* **Nonzero voltages.**
  - V_b != 0 at every device bus, because every device bus is banded (|V| >= 0.95);
    `device_rectangles` raises an error if a device sits at a source terminal, and none does.
  - V0_b != 0: checked at every query and hour. The smallest |V0_b| at a device bus is 1.0495 pu, and
    the code raises an error below 0.5.
  - The terminal-voltage lemma (`source_current_bounds.md`) gives V != 0 at the terminals, which the
    elimination needs.
* **In the relaxation.**
  - I_s = (E_s - V_s) / z_s is affine in the terminal's first moments (e, f).
  - conj(S_b) is affine in P_B and Q_PV.
  - So each (source, hour) gets one second-order cone [eps; Re D; Im D], scaled by |z_s| to give
    unit-size coefficients: 72 cones per day.
* **Numerical allowance on the cone data.** eps is multiplied by 1 + 1e-9 and increased by
  1e-12 * (|a_s| + sum_b |k_sb| Smax_b + 3 |1/z_s|). That covers rounding a_s, k_sb = c_sb / conj(V0_b)
  (from the 256-bit Z) and 1/z_s to Float64.
* **Verifier.** `dual_bound` now projects second-order-cone multipliers onto the cone (self-dual), and
  `export_certificate` writes them. The exact verifier projects them too, then confirms membership in
  rational arithmetic, raising the head if needed. In every case all 72 cone multipliers were already
  in the cone, so none was adjusted.

## Comparison at the 1.02 r* queries

| Query | Case | Clarabel | rho relaxed (kW) | **Checked lower bound (kW)** | Exact re-check (kW) | Box term (kW) | Residual | Fictitious losses (kWh) | Source moment gap (pu^2) | eps (A) | Build / solve / verify (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 292.5 deg, 1.02 | baseline | NUMERICAL_ERROR | 0.657 | 0.3327 | 0.332714 | 45.6 | 4.4e-7 | 270 | 0.817 | | 1.0 / 26.5 / 1.6 |
| | **dispatch link** | NUMERICAL_ERROR | 2.750 | **2.4598** | 2.459787 | 30.3 | 4.9e-7 | 96 | 0.637 | 0.406-0.710 | 2.6 / 27.2 / 1.5 |
| 157.5 deg, 1.02 | baseline | NUMERICAL_ERROR | 0.018 | -0.4614 (unresolved) | | 1.99 | 2.6e-7 | 296 | 1.007 | | 1.3 / 25.7 / 1.5 |
| | **dispatch link** | NUMERICAL_ERROR | 2.493 | **2.2340** | 2.234008 | 7.58 | 4.5e-7 | 134 | 0.621 | 0.407-0.710 | 1.7 / 27.3 / 1.6 |

eps values in per unit (solver base, I_base = 86.603 A): 4.69e-3 to 8.20e-3. Every value is in
`daily_link_eps.csv`, by source and hour.

Costs:

* The derived bounds take 13 s per angle vector (24 hours, corner enclosure) and supply eps at no
  extra cost.
* The cones add about 1 s to the build. Solve and verification are unchanged.
* The exact re-check of a 24-hour certificate takes 7-8 s.

Notes:

* The checked bounds are valid for any multipliers, NUMERICAL_ERROR solves included. Primal values
  from those solves are not solutions, and the fictitious-loss and moment-gap figures from them are
  indicative only.
* "Source moment gap" is the largest, over sources and hours, of the lifted |I_s|^2 minus
  |first-moment I_s|^2.

## Controls at 0.98 r*

| Query | Case | Checked bound (kW) | Lifted witness: largest violation, all constraints including the new cones | Box violation |
|---|---|---|---|---|
| 292.5 deg | baseline / link | -0.3329 / -0.0736 | 9.9e-7 (a pre-existing reactive-balance row, at the witness's solve tolerance) | 0 |
| 157.5 deg | baseline / link | -0.0634 / -0.0546 | 4.3e-9 | 0 |

Both controls stay unresolved, as they must.

The cones were also evaluated directly on every validated exact dispatch at the six daily queries,
margin and restoration witnesses alike (`daily_link_check.log`):

* |I_s - a_s - sum k conj(S_b)| / eps_s is at most **0.468** at 292.5 deg and **0.436** at 157.5 deg.
* No real state comes near its allowance, which supports the validity of eps and shows it is about 2x
  conservative on these states.

## Stability of the thin 292.5 deg, 1.02 r* certificate (baseline, no link)

| Variation | Clarabel | Checked lower bound (kW) | Exact re-check (kW) |
|---|---|---|---|
| as before (static regularization 1e-7) | NUMERICAL_ERROR | 0.3327 | 0.332714 |
| static regularization 1e-8 | NUMERICAL_ERROR | 0.3704 | 0.370387 |
| static regularization 1e-6 | ALMOST_OPTIMAL | 0.4427 | 0.442711 |
| equilibration off | NUMERICAL_ERROR | 0.3713 | 0.371335 |

The positive certificate persists under every variation, at 0.33-0.44 kW, exactly re-checked each
time, whatever the solver status. With the link it is 2.46 kW.

## Regression: 1.10 r* with the link (the formulation changed)

| Query | Before (corner baseline) | With the link | Exact re-check | Exact rho |
|---|---|---|---|---|
| 292.5 deg, 1.10 | 22.586 kW | **24.770 kW** (ALMOST_OPTIMAL) | 24.770108 | 27.864 |
| 157.5 deg, 1.10 | 21.229 kW | **24.575 kW** (ALMOST_OPTIMAL) | 24.574958 | 27.838 |

## Is it the error allowance or the moment freedom?

At 157.5 deg, 1.02 r* (`daily_link_eps_scan.csv`, `daily_link_robust.csv`):

| eps scale | Valid | Clarabel | Checked lower bound (kW) | Fictitious losses (kWh) | Source moment gap (pu^2) |
|---|---|---|---|---|---|
| 2.0 | yes | NUMERICAL_ERROR | -0.357 | 292 | 1.008 |
| 2.0 (regularization 1e-6) | yes | ALMOST_OPTIMAL | -0.484 | 289 | 0.997 |
| 1.5 | yes | ALMOST_OPTIMAL | 0.523 | 206 | 0.962 |
| 1.25 | yes | NUMERICAL_ERROR | 1.278 | 169 | 0.846 |
| 1.1 | yes | NUMERICAL_ERROR | 1.852 | 148 | 0.712 |
| **1.0** | yes | NUMERICAL_ERROR | **2.234** | 134 | 0.621 |
| 1.0 (regularization 1e-6) | yes | ALMOST_OPTIMAL | 2.210 | 131 | 0.620 |
| 1.0 (regularization 1e-8) | yes | NUMERICAL_ERROR | 2.156 | 132 | 0.626 |
| 0.1 | **no** (diagnostic) | NUMERICAL_ERROR | meaningless: the relaxation is nearly infeasible | 0.6 | |
| 0 | **no** (diagnostic) | ALMOST_INFEASIBLE | | | |

Reading the scan:

* **The bound is continuous and steep in eps.** It is not a numerical accident: at eps x 1 three
  regularizations agree to within 0.08 kW.
* **The relaxation uses the whole allowance.** Fictitious losses and the source second-moment gap grow
  with eps, and the bound falls with them.
* **eps x 0.1 and x 0 cannot be used,** even as a target. They exclude real states, which need up to
  0.47 of eps, so the relaxation becomes nearly infeasible and the bound is meaningless.
* **The link is first-order.** It constrains only first moments (e, f at the terminals, P_B, Q_PV).
  The second moments it does not touch, the lifted |I_s|^2 and the line-current moments behind the
  fictitious losses, are where the relaxation spends the allowance.
* **The gap is driven by the allowance, not by moment freedom independent of it.** Within eps x 1,
  moment freedom still costs about 3.3 kW of the 5.57 kW. But shrinking eps, not adding moment
  blocks, is what the data says moves the bound.

## Where this leaves the method

* **Valid ways to tighten eps (not done here).**
  - Bound |S_b| hour by hour from the dispatch, rather than by its rectangle maximum Smax_b.
  - Use per-hour rho_b, which current_bounds already computes but which could be tightened by
    iterating with the link itself.
  - Replace the product bound Smax_b rho_b by a joint bound on conj(S_b) conj(V0_b - V_b), whose worst
    cases rarely coincide: the real states use under half of eps.
* **If eps cannot be tightened,** a second-order link would bound the lifted |I_s - A|^2 through
  moments of (V, P_B, Q_PV) products. That means moment blocks coupling battery and voltage variables:
  a larger relaxation, to consider only after the eps route.
* **Not done, as instructed:** no empirical caps, no paper edits, no larger feeder, no angular sweep.
  Conclusions are pointwise. Nothing here concerns global AC optimality: the reaches and the
  fixed-schedule experiments are local.
* **Verification levels, kept separate.**
  - The exact re-checks establish the certificate bounds exactly, given the exported data.
  - The derived bounds, eps and the verifier boxes are numerical constructions (Float64 with stated
    allowances; the Z matrix in 256-bit arithmetic).

## Reproduction

From the repository root, in envs/msopf_sdp:

    unset NO_COLOR FORCE_COLOR
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_link.jl        # comparison, stability, regression (about 10 min)
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_link_eps.jl    # eps x 2, 1, 0.1, 0
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_link_check.jl  # witnesses vs eps; eps x 1.1-1.5; regularizations
    gunzip -k envs/multi_poi/results/multiperiod/certificates/daily_link_<case>.txt.gz
    julia envs/multi_poi/root_level/verify_certificate.jl envs/multi_poi/results/multiperiod/certificates/daily_link_<case>.txt

Results: `daily_link{.csv,_eps.csv,_eps_scan.csv,_robust.csv,.log,_eps.log,_check.log}`, and
`certificates/daily_link_*.txt.gz` (8 certificates, gzipped, about 5 MB each).
