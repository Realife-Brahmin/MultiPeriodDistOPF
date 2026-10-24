# Certificate audit and the daily pilot (2026-10-22)

IEEE 123, sources 1, 3 and 4, PV reactive power dispatched. No empirical current caps or assumed
ampacity anywhere. Everything is pointwise: each certificate covers one angle vector. The AC reaches r*
are locally measured.

## Summary

1. **The two hour-12 certificates at 1.02 r* survive an independent exact check.**
   * Exported certificate data, recomputed in exact rational arithmetic with no solver, gives
     1.090440037 kW (292.5 deg) and 1.862403089 kW (157.5 deg). The Float64 checks gave
     1.090440 and 1.862403, so the Float64 rounding error was 1.5e-7 and 2.3e-7 kW, inside the
     allowance.
   * Every PSD multiplier is exactly PSD as returned (exact LDL').
   * The derived source-current bounds recomputed in 256-bit arithmetic lie 1.00e-9 (relative) below
     the values used, exactly the allowance.
   * A numerically better but equivalent form (source-current rows scaled by |z_s|^2) raises the
     checked bounds to 2.33 and 1.96 kW.
2. **Batteries: valid per-hour source-current bounds exist.**
   * With the sector enclosure they are 73-106 A (up to 2.5x the idle-battery bounds). With the
     corner enclosure, 56-97 A.
   * They are loosest at night. Nine buses hold load, PV and battery together, and there the device
     power set contains zero.
3. **The daily pilot (T = 24, batteries dispatched, hard cyclic SOC) remains useful, with a wider
   gap.**
   * 1.10 r* is certified in both directions, with checked lower bounds of 22.59 and 21.23 kW
     (exact rho 27.86 and 27.84 kW).
   * 1.02 r* is certified at 292.5 deg (0.333 kW, a thin margin from a NUMERICAL_ERROR solve) and
     unresolved at 157.5 deg.
   * The 0.98 r* points are unresolved, as they must be, and their witnesses lift cleanly.
   * Diagnosis: the loss comes from the wider device-current sets, not from temporal coupling or
     conditioning.
4. **Cost per angle vector:**
   * Bounds: 3-13 s.
   * Relaxation build: about 1 s.
   * Solve: 17-33 s (6,288 PSD blocks, 88,587 variables).
   * Float check: about 1.5 s.
   * Optional exact check: about 7.5 s.
   * The network response (5 s) is computed once per network.

## Reproduction

From the repository root, in envs/msopf_sdp (JuMP 1.30, Ipopt 1.14.1, Clarabel 0.11.1):

    unset NO_COLOR FORCE_COLOR
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/certificate_audit.jl   # part 1 (about 35 min, mostly the 256-bit bounds)
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_pilot.jl         # parts 2-3, sector enclosure
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_diagnose.jl      # loss diagnosis
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/daily_refined.jl       # corner enclosure
    julia envs/multi_poi/root_level/verify_certificate.jl envs/multi_poi/results/multiperiod/certificates/<file>.txt   # no solver, no JuMP

Code:

| File | What |
|---|---|
| `current_bounds.jl` | bounds for device power rectangles; sector and corner enclosures; generic in the number type |
| `feasibility_certificate.jl` | per-hour source bounds, row scaling, `export_certificate` |
| `root_level/verify_certificate.jl` | exact verifier |
| `root_level/certificate_audit.jl` | part 1 |
| `root_level/daily_pilot.jl` | parts 2-3, sector enclosure |
| `root_level/daily_diagnose.jl` | loss diagnosis |
| `root_level/daily_refined.jl` | corner enclosure |

Results:

| File | What |
|---|---|
| `certificate_audit{.csv,_bounds.csv,.log}` | part 1 |
| `daily_pilot_{reach,bounds,queries}.csv`, `daily_pilot.log` | the pilot |
| `daily_diagnose{.csv,.log}` | the diagnosis |
| `daily_refined{.csv,.log}` | corner enclosure |
| `certificates/*.txt` | exported certificate data: 10 hour-12 files, 5 daily |

## 1. The two near-boundary hour-12 certificates

Settings as committed: hour 12, batteries idle, both nodal strengthenings, derived source-current
bounds (sector enclosure), restoration relaxation with rho <= rho_exact + 1 kW.

### Dependency chain, audited

1. **Terminal-voltage lemma.**
   * Uses the model's box |e|, |f| <= 1.5, the POI line z = 7.5e-9 (1 + j) pu, and purely capacitive
     shunts (real part exactly 0 in the data).
   * Re-derived; it is unchanged by this round.
2. **Network elimination.**
   * Re(Y) is positive definite: 256-bit Cholesky, smallest pivot 16.9.
   * Z = Y^-1 in 256 bits, |YZ - I| = 1.7e-71.
3. **Device-current sectors and their iteration.**
   * The disc radius around the no-load voltage goes 0.006727 -> 0.003794 -> 0.003773 pu (256-bit),
     the same as Float64.
   * Each step uses only containments proved by the previous one.
4. **Support functions and the direction-grid correction.**
   * The grid spacing is 0.1 deg, and the 1/cos(0.05 deg) factor covers the discretization.
   * Recomputed in 256-bit arithmetic on the same grid, the bounds come out at 0.99999999899 to
     0.99999999900 of the values used. The 1e-9 allowance is what separates them, so the Float64
     error itself is about 1e-11.
5. **Source-current constraints and verifier boxes.**
   * The bounds enter the model as constants, sa^2 in the lifted |I_s|^2 row.
   * The boxes come from `tighten_boxes!`: source and path discs, then interval products.
   * Widening every box by 1e-12 relative (more than the Float64 rounding of a few interval products)
     changes no bound in the ninth digit.
6. **Projected multipliers, residual correction, final bound.**
   * Recomputed exactly from the exported data (below).

### Exact recomputation

`export_certificate` writes the objective, every variable's box, every scalar row with its dual value,
and every PSD block with its dual triangle and its entries' terms. Numbers are in shortest round-trip
form, so they parse back exactly; each hour-12 file is about 0.87 MB.

`verify_certificate.jl` reads such a file with Base and LinearAlgebra only. It recomputes, in
Rational{BigInt}:

* the sign-projected scalar multipliers;
* the PSD multipliers, after an exact LDL' test: the raw dual triangle is used if it is PSD, otherwise
  the Float64 eigenvalue clip shifted by tau*I until the exact LDL' succeeds;
* the Frobenius pairing, the residual coefficients and the box minimum.

The result is an exact lower bound, given the file's data.

| Point | Form | Clarabel | Float64 check (kW) | Allowance (kW) | **Exact (kW)** | Exact minus Float64 | Raw PSD duals PSD | Residual | Solve (s) |
|---|---|---|---|---|---|---|---|---|---|
| 292.5 deg, 1.02 r* | committed | ALMOST_OPTIMAL | 1.090440 | 1.5e-7 | **1.090440037** | +1.45e-7 | 262 / 262 | 3.6e-4 | 6.3 (with compilation) |
| | rows scaled by \|z_s\|^2 | OPTIMAL | 2.330506 | 3.2e-7 | **2.330505835** | +3.2e-7 | 262 / 262 | 1.5e-7 | 0.7 |
| | tight tolerances (1e-9) | ALMOST_OPTIMAL | 1.367856 | 1.8e-7 | **1.367856419** | +1.8e-7 | 262 / 262 | 5.7e-4 | 0.5 |
| | scaled + tight | ALMOST_OPTIMAL | 2.286625 | 3.2e-7 | **2.286624828** | +3.2e-7 | 262 / 262 | 1.9e-6 | 0.6 |
| 157.5 deg, 1.02 r* | committed | OPTIMAL | 1.862403 | 2.3e-7 | **1.862403089** | +2.3e-7 | 262 / 262 | 8.6e-6 | 0.8 |
| | rows scaled | ALMOST_OPTIMAL | 1.912365 | 2.6e-7 | **1.912365267** | +2.6e-7 | 262 / 262 | 3.8e-7 | 0.8 |
| | tight tolerances | ALMOST_OPTIMAL | 1.130138 | 1.8e-7 | **1.130138211** | +1.8e-7 | 262 / 262 | 3.2e-4 | 0.9 |
| | scaled + tight | ALMOST_OPTIMAL | 1.957611 | 2.4e-7 | **1.957611328** | +2.4e-7 | 262 / 262 | 4.0e-7 | 0.8 |

* **The rounding allowance covers the Float64 error.** "Exact minus Float64" equals the allowance
  that was subtracted, to two digits. The actual Float64 accumulation error is far smaller than the
  allowance: the exact value sits at the Float64 value plus the allowance.
* **No PSD shift was needed.** Every raw PSD dual was exactly PSD, so the eigenvalue-clip path, whose
  rounding the Float64 allowance did not formally cover, was never exercised.
* **Box rounding is negligible.** Widening every box by 1e-12 relative leaves every bound unchanged
  to nine digits.
* **The scaling matters.** The lifted source-current row has coefficients 1/|z_s|^2 = 2.8e4. Scaling
  it by |z_s|^2 lets Clarabel reach OPTIMAL at 292.5 deg (residual 1.5e-7 against 3.6e-4) and more
  than doubles the checked bound (2.33 against 1.09 kW). Tighter tolerances alone make things worse.
  The daily runs use the scaled form.
* **No other conic solver is available.** Gurobi is installed but has no semidefinite cones. None of
  SCS, Mosek, COSMO or Hypatia is installed, so there is no cross-solver check.

### Controls at 0.98 r*

The validated dispatch, lifted into the strengthened relaxation with rho = 0, gives:

| Control | Largest constraint violation | Largest violation of a tightened verifier box | Source currents vs. derived bounds |
|---|---|---|---|
| 292.5 deg | 4.3e-9 (a reactive balance row, the witness's own tolerance) | 0 | at most 0.982 |
| 157.5 deg | 4.3e-9 | 0 | at most 0.800 |

The restoration bounds are -0.0133 and -0.0045 kW (exact): unresolved, as they must be.

### Level of verification achieved

* **Achieved:** an exact-arithmetic recomputation of the Lagrangian lower bound from the exported
  data. That covers the projection, the PSD check and the coefficient accumulation, with an exact PSD
  proof for every block.
* **Achieved:** a 256-bit recomputation of the derived source-current bounds, which agrees with the
  values used to within the 1e-9 allowance.
* **Not formal:**
  - the box construction is in Float64 (shown insensitive at 1e-12);
  - the exporter and verifier mirror the code they check (the floating results agree with the exact
    ones to within the allowance);
  - the model data are taken as given;
  - no interval arithmetic or proof assistant was used.

  These are independently checked numerical certificates, with the final bound exact given the data.

## 2. Source-current bounds with dispatched batteries

Battery model as in the original:

* lossless: B_t = B_{t-1} - P_B,t dt, with dt = 1 h;
* |P_B,t| <= rating;
* SOC within [0.30, 0.95] x rated energy, B_0 = 0.625;
* cyclic B_24 = B_0, hard.

### Device power sets

At bus b and hour t, every allowed device power lies in the rectangle
[P0 - Pb, P0 + Pb] x [Q0 - qd, Q0 + qd]:

* P0: PV real power minus load, both fixed;
* Q0: minus the load's reactive power;
* qd: the PV reactive range sqrt(S_rated^2 - P_PV^2);
* Pb: the battery rating.

**SOC coupling adds no tightening at T = 24.** Every battery can run at full rating for at least
1.30 h in either direction from B_0 (energy-to-power ratio 4 h; 0.325 of capacity on each side). So
no hour's range shrinks below [-P_rated, P_rated]. Nothing is conditioned on a schedule.

### Two enclosures of J_b = conj(S_b) V_b / |V_b|^2

* **Sector** (`daily_pilot.jl`).
  - If the rectangle misses zero, arg S lies in the arc it subtends (< 180 deg), and J lies in a
    sector, as for idle batteries.
  - If the rectangle contains zero, arg S is unrestricted, and J lies in the full disc of radius
    Smax / Vmin.
* **Corner** (`daily_refined.jl`, `corners = true`).
  - For a fixed V, J is conj(rectangle) rotated by arg V and scaled by 1/|V|. Its support function in
    direction psi is max_i |c_i| cos(psi - arg V + arg c_i) / |V| over the four corners c_i.
  - Over arg V within asin(rho / |V0|) of arg V0 and |V| >= Vmin, each corner term is at most a
    sector term in closed form: |c_i| / Vmin times cos of the angular excess beyond that window,
    floored at 0.
  - This holds for every rectangle, including those that contain zero.
* **Check of the corner enclosure.** At hour 12 with batteries idle, the corner bounds equal or
  tighten the sector bounds (subs4: 31.4 against 33.1 A). They stay above every current the earlier
  falsification search found, the largest ratio being 0.9966.
* A two-wedge refinement was also tried, for batteries whose power is purely real through zero. It
  applies to no bus of this deck: the nine zero-containing buses all have load, PV and battery
  together.

### Bounds by source and hour (1.02 r*; every hour is in `daily_pilot_bounds.csv` and the logs)

| Direction, source | Idle batteries | Dispatched, sector | Dispatched, corner | Loosest hours | Largest current in the witness |
|---|---|---|---|---|---|
| 292.5, subs1 | 35-67 A | 87-99 A | 67-90 A | night (largest at hour 6) | 57 A |
| 292.5, subs3 | 45-73 | 91-104 | 77-95 | night (hour 6) | 67 |
| 292.5, subs4 | 41-51 | 73-81 | 57-69 | daytime (hours 8-17; largest at 11) | 42 |
| 157.5, subs1 | 41-73 | 92-106 | 73-96 | night (hour 6) | 51 |
| 157.5, subs3 | 44-50 | 81-89 | 61-73 | daytime (hours 10-11) | 31 |
| 157.5, subs4 | 36-64 | 77-89 | 61-80 | night (hour 6) | 55 |

Why the bounds widen:

* The nine buses with load, PV and battery contain zero power in every hour. With the sector
  enclosure, their device currents are unrestricted in phase.
* At night the PV real power is zero, so the full inverter rating is available as reactive power.
  Together with the battery rating, the rectangles there are largest.
* Bounding each device separately adds up 13 battery ratings with no cancellation. The corner
  enclosure keeps the rectangle's shape and recovers 10-30%.
* The witness currents use at most 0.80 of the corner bounds and 0.70 of the sector bounds.

## 3. The daily pilot

Daily reach, secant on the validated exact margin, the best of a cold and a warm start at each iterate
(`daily_pilot_reach.csv`):

| Direction | Daily r* | Hour-12 r* | Validated margin there |
|---|---|---|---|
| 292.5 deg | **0.241513 deg** | 0.118 deg | -9e-5 kW |
| 157.5 deg | **0.194598 deg** | 0.100 deg | -6e-9 kW |

Batteries roughly double the reach over the day, consistent with the earlier studies.

**Restoration relaxation.** The full coupled 24-hour relaxation:

* battery power and SOC variables kept, with their dynamics and the cyclic constraint;
* band and device limits hard;
* min rho with P_s,t >= -rho for every source and hour;
* rho <= rho_exact + 1 kW, which keeps every rho = 0 candidate;
* both nodal strengthenings;
* the derived bounds as lifted source-current rows, scaled.

The relaxation has 6,288 PSD blocks and 88,587 variables.

### Witnesses

All witnesses below were validated:

| Query | Angles (d3, d4 deg) | Exact margin L (kW) | Exact restoration rho (kW) | Equations | SOC residual | Cyclic residual |
|---|---|---|---|---|---|---|
| 292.5, 0.98 | +0.0906, -0.2187 | +5.573 | 0 | 6e-11 | 1e-16 | 0 |
| 292.5, 1.02 | +0.0943, -0.2276 | -5.572 | 5.573 | 5e-10 | 1e-16 | 0 |
| 292.5, 1.10 | +0.1017, -0.2454 | -27.863 | 27.864 | 3e-10 | 1e-16 | 0 |
| 157.5, 0.98 | -0.1762, +0.0730 | +5.568 | 0 | 5e-10 | 1e-16 | 0 |
| 157.5, 1.02 | -0.1834, +0.0760 | -5.568 | 5.568 | 3e-11 | 1e-16 | 0 |
| 157.5, 1.10 | -0.1978, +0.0819 | -27.838 | 27.838 | 5e-11 | 1e-16 | 0 |

### Restoration test results

| Query | Enclosure | Relaxation status | rho relaxed (kW) | **Checked lower bound (kW)** | Exact re-check (kW) | Box term / rounding (kW) | Residual | Lifted witness (constraints / boxes) | Bounds / build / solve / verify (s) | Outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| 292.5, 0.98 | sector | ALMOST_OPTIMAL | -0.0000 | -0.0905 | | 0.088 / 3e-8 | 1e-8 | 4.3e-9 / 0 | 4.6 / 6.3 / 29.2 / 4.1 | unresolved (feasible control) |
| 292.5, 1.02 | sector | ALMOST_OPTIMAL | 0.0054 | -0.2165 | | 1.39 / 5e-8 | 4e-8 | | 2.9 / 1.1 / 24.2 / 1.5 | unresolved |
| 292.5, 1.10 | sector | ALMOST_OPTIMAL | 20.784 | **20.554** | 20.553979 | 33.7 / 4e-7 | 2e-7 | | 2.8 / 1.2 / 24.9 / 1.6 | certified |
| 157.5, 0.98 | sector | ALMOST_OPTIMAL | -0.0000 | -0.0568 | | 0.051 / 3e-8 | 7e-9 | 4.3e-9 / 0 | 2.9 / 1.1 / 16.5 / 1.4 | unresolved (feasible control) |
| 157.5, 1.02 | sector | ALMOST_OPTIMAL | 0.0002 | -0.0852 | | 0.108 / 4e-8 | 4e-9 | | 2.8 / 1.2 / 21.3 / 1.4 | unresolved |
| 157.5, 1.10 | sector | ALMOST_OPTIMAL | 18.932 | **18.662** | 18.662317 | 7.11 / 3e-7 | 3e-7 | | 2.9 / 1.3 / 22.2 / 1.4 | certified |
| 292.5, 0.98 | corner | NUMERICAL_ERROR | -0.0000 | -0.3329 | | 0.325 / 4e-8 | 6e-8 | 9.9e-7 / 0 | 13.4 / 6.4 / 31.2 / 3.8 | unresolved (feasible control) |
| 292.5, 1.02 | corner | NUMERICAL_ERROR | 0.657 | **0.3327** | 0.332714 | 45.6 / 3e-7 | 4e-7 | | 13.0 / 1.0 / 26.4 / 1.4 | certified (thin; multipliers from a failed solve) |
| 292.5, 1.10 | corner | ALMOST_OPTIMAL | 22.834 | **22.586** | 22.585695 | 47.5 / 3e-7 | 3e-7 | | 13.0 / 1.1 / 26.5 / 1.4 | certified |
| 157.5, 0.98 | corner | ALMOST_OPTIMAL | -0.0000 | -0.0634 | | 0.056 / 3e-8 | 8e-9 | 4.3e-9 / 0 | 13.3 / 1.2 / 16.8 / 1.3 | unresolved (feasible control) |
| 157.5, 1.02 | corner | NUMERICAL_ERROR | 0.0182 | -0.4614 | | 1.99 / 9e-8 | 3e-7 | | 13.1 / 1.2 / 25.9 / 1.6 | unresolved |
| 157.5, 1.10 | corner | NUMERICAL_ERROR | 21.548 | **21.229** | 21.229019 | 7.19 / 3e-7 | 4e-7 | | 13.0 / 1.3 / 22.2 / 1.4 | certified |

Notes on these results:

* Every dual multiplier lies in its cone: PSD and sign violations are 0. The exact re-checks shifted
  no PSD block (0 of 6,288).
* Where Clarabel stopped with NUMERICAL_ERROR, the checked bound is still valid: it holds for any
  multipliers. But the multipliers come from a failed solve, and primal values from those solves are
  not solutions.
* The 292.5 deg, 1.02 r* certificate (0.333 kW) is positive but thin: about 0.7% of its box
  correction.
* The first build of each run includes compilation (6.3 s).

### Diagnosis at 1.02 r*

`daily_diagnose.csv`: the same angle, same relaxation, sector enclosure.

| Case | 292.5: checked bound (kW) | Fictitious losses | 157.5: checked bound (kW) | Fictitious losses | Valid? |
|---|---|---|---|---|---|
| V: dispatched-battery bounds, batteries free | -0.39 | 269 kWh | -0.09 | 258 kWh | yes |
| D1: idle-battery bounds, batteries free | 6.19 | 50 kWh | 2.30 | 167 kWh | **no**: 6.19 exceeds the exact 5.57, so these bounds cut off the true optimum |
| D2: dispatched bounds, batteries fixed to the witness schedule | -0.46 | 280 kWh | -0.14 | 252 kWh | for the fixed-schedule model |
| D3: bounds for that fixed schedule, batteries fixed | 3.68 | 75 kWh | 2.20 | 147 kWh | for the fixed-schedule model |

**The main loss is the wider device-current sets.**

* With them, the relaxation books about 260 kWh a day of fictitious line losses, and those absorb the
  5.6 kW of backflow.
* Fixing the battery schedule does not help (D2), so the loss is not from temporal coupling or
  battery flexibility inside the relaxation. The SOC dynamics are exact linear constraints in the
  relaxation anyway.
* Not numerical either: the residuals are at most 4e-7 and the bounds were checked exactly.
* With current sets that match a fixed schedule, the per-hour network relaxation is tight enough to
  certify (D3).
* The corner enclosure closes part of the gap: it certifies 292.5 deg at 1.02 r*, but not 157.5 deg.

## 4. Cost

| Item | Time | Reused |
|---|---|---|
| Network response Z (256-bit) | 5.1-5.4 s | once per network (and per load model) |
| Source-current bounds, 24 hours | 2.8-4.6 s (sector); 13 s (corner) | per angle vector |
| Hour-12 bounds in 256-bit, for checking only | 465-515 s per point | not needed for classification |
| Relaxation build | 1.0-1.3 s (6.3 s with compilation) | per angle vector |
| Relaxation solve, T = 24 | 16.5-33 s | per angle vector |
| Float check (dual_bound) | 1.3-4.1 s | per angle vector |
| Exact re-check from the file | 7-8 s (T = 24); 0.2-0.4 s (hour 12) | optional, per certificate |

Classifying one more daily angle vector costs about 25-50 s with the sector bounds, or 40-60 s with
the corner bounds, plus the exact re-check if wanted.

## Where this stands

* **Hour 12:** both near-boundary certificates hold under an exact recomputation from exported data,
  and their derived bounds hold in 256-bit arithmetic. The scaled form gives larger margins.
* **Daily, with batteries and cyclic SOC:** usefulness is kept at a wider gap. Points 10% beyond the
  daily reach are certified in both directions, with checked bounds at 68-81% of the exact restoration
  value. At 2% beyond, one direction is certified, with a thin margin.
* **What limits the daily case:** the device-current enclosures at buses where load, PV and battery
  can together take zero power. Each such bus enters every source's bound at its worst corner.
  Tighter valid enclosures there are the lever. Candidates (not done):
  - an enclosure per bus that tracks the joint |S| along the rectangle's edges, not just its corners;
  - a per-hour OBBT on the daily relaxation, if its numerical issue (the 1/|z_s|^2 amplification,
    now reduced by the row scaling) allows.
* Held, as instructed: the larger feeder, any broader sweep, and paper edits.
