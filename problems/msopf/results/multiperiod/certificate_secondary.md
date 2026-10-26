# Certificates on the 2705-bus feeder with service transformers (2026-10-26)

The 123-node certificate pilot (`certificate_validation.md`) moved to `ieee9500xf_1ph`: three substations,
the five ties closed, the 0.95-1.05 band at customer buses only, batteries idle, PV reactive power
dispatched, the limiting hour alone. Driver `root_level/certificate_secondary_pilot.jl` (environment
`problems/msopf/environments/sdp`, Clarabel 0.11.1); logs `certificate_secondary*.log`, CSVs
`certificate_secondary_{hours,rstar,envelope,cases*}.csv`.

## Setup

* **Hour.** Along 135 deg, the T = 24 idle edge of the 2026-10-25 campaign is limited by hour 13.
* **Hour boundary.** At hour 13 alone, the boundary along each of the 16 directions was found by
  bisecting the sign of the exact, re-validated margin. It matches the campaign's T = 24 reach to
  within 0.004 deg in every direction, as it must with idle batteries.
* **Current envelope.** Each line's and source's largest current over 17 frozen, validated dispatches:
  equal angles and 98% of the hour boundary in every direction. The limits are this envelope itself
  (1.0x, no cushion), with a floor for lines that carry almost nothing.
* **What is certified.** These are empirical limits, not ratings: a certificate proves the model
  **with these limits** infeasible, nothing about the model without them.

## Results (hour 13, idle batteries)

| Direction | Point | Exact margin L (kW) | Bound U, nodal cuts, floor 1e-4 pu | Bound U, nodal cuts, floor 0.01 pu | Class |
|---|---|---|---|---|---|
| 135 deg | inside, 98% | +20.37 | +23.99 | +29.67 | feasible |
| 135 deg | outside, 102% | -20.34 | **-16.57** | **-11.22** | **infeasible under the limits** |
| 315 deg | inside, 98% | +19.98 | +26.77 | +32.10 | feasible |
| 315 deg | outside, 102% | -19.95 | **-13.10** | **-7.71** | **infeasible under the limits** |

* Inside points are validated AC dispatches; their OpenDSS replays agree within 1e-7 kW.
* Without the nodal cuts (`:kcl0`, `:jbound`), nothing certifies: U = +89 to +100 kW outside (floor
  1e-4), +179 to +191 kW (floor 0.01).
* With them, the bound is within 4-12 kW of the exact margin, a tight relaxation on a 3,980-bus meshed
  feeder.
* Clarabel ends with NUMERICAL_ERROR or SLOW_PROGRESS on most of these solves. The bounds are computed by
  `dual_bound` from the multipliers it returns, valid whatever the solver status; a better linear solver
  should tighten them, not change their validity.

## Changes that made it run

* **Customer-only band** (`feasibility_certificate.jl`). The relaxation and `validate_dispatch` band only
  `kw.band_buses`. Other buses get `|V| >= 0.5 pu` (V_FLOOR), added to the certified model and checked on
  every dispatch, so that injection currents stay bounded.
* **Loops.** The derived line-current bounds assumed a tree. A line on a loop (closed ties) got bound zero,
  which would have made certificates invalid; it now gets the total bound.
* **Near-ideal sources.** The substation sources of this deck have |z_s| = 1.9e-10 pu, so the relaxation
  carried 5e9 coefficients and stalled at SLOW_PROGRESS with bounds near 1e20 kW. A source with
  |z_s| < 1e-6 pu is now written as ideal inside the relaxation: terminal voltage fixed at E, power free in
  the bus balance with the no-backflow margin, and its current limit as the exact cone |S| <= |E| I_max.
  Coefficients now span 1e-6 to 12 (`certificate_scaling_check.log`). The certified model has ideal sources;
  the difference is below |z_s| |I_s|, about 1e-9 pu.

## Next

* MA57 or MA97 as Clarabel's linear solver (its HSL extension), to clear the numerical-error statuses.
* More directions, then the San Francisco deck, then the cyclic-battery day.
