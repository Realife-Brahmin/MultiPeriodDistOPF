# Source-current bounds from the original model (2026-10-22)

**Question.** Can current bounds derived from the original model replace the empirical caps, and still give
useful infeasibility certificates?

**Answer, for these queries: yes.**

* Source-current bounds derived analytically from the network equations, the band and the device limits
  were tested at 10 queries: IEEE 123, hour 12, batteries idle, PV reactive power dispatched,
  directions 292.5 and 157.5 deg.
* The bounds are 31-39 A per source, below even the empirical envelope (46-56 A), and less than half
  the 1.5x caps (69-84 A).
* With the two nodal strengthenings, every tested point at or beyond 1.02 r* is certified (8 of 8).
  The checked lower bound on the restoration objective is 1.09-34.7 kW. Both 0.98 r* points stay
  unresolved, as they must: they are feasible.
* No ampacity, line limit or empirical cap enters the model, the bounds or the certificates.

Files:

| File | What |
|---|---|
| `root_level/source_current_bounds.jl` | bounds, OBBT and restoration tests |
| `root_level/source_current_falsify.jl` | falsification search |
| `source_current_bounds{.log,_bounds.csv,_queries.csv}` | bounds and query results |
| `source_current_falsify.log`, `source_current_bounds_falsify.csv` | falsification results |

Reproduce, from the repository root (envs/msopf_sdp: JuMP 1.30, Ipopt 1.14.1, Clarabel 0.11.1):

    unset NO_COLOR FORCE_COLOR
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/source_current_bounds.jl
    MULTIPERIOD_SOLVE_LOG=0 julia --project=envs/msopf_sdp envs/multi_poi/root_level/source_current_falsify.jl

## The model the bounds hold for

The original single-hour model at fixed source angles. Its parts:

* Power-form KCL at every bus, line series and shunt admittances, and sources as fixed EMFs E_s behind
  z_s = 0.00024 + j0.0060 pu.
* Constant-power loads (hour-12 load factor 0.9873), PV real power fixed (factor 1.0), PV reactive power
  in [-q_max, q_max] with q_max = sqrt(S_rated^2 - P^2), and batteries idle.
* The 0.95-1.05 pu band at every bus except the three source terminals.
* The model's box |e|, |f| <= 1.5 at every bus.

**No-backflow is not used by the analytical bounds.** They hold for every state satisfying the above,
whatever the source imports.

Units: S_base 1,000 kVA per phase, V_base 11.547 kV, I_base 86.603 A (unchanged).

## Derivation

**1. Every source terminal voltage is nonzero.** This lets power-form KCL be divided by V.

* Suppose V_t = 0 at the terminal t of source s. The terminal's POI line has z = 7.5e-9 (1 + j) pu, and
  its far bus is banded (|V| >= 0.95). The line therefore dissipates at least
  0.95^2 * r / |z|^2 = 6.0e7 pu of real power.
* Summing power-form KCL over all buses, total real losses plus loads equal source imports plus PV.
  Line shunts are purely capacitive here, with real part exactly 0 in the data.
* Each source imports at most |V_s| |I_s| <= 2.12 * (1.05 + 2.12) / 0.006 = 1,120 pu, using the box.
  Three sources give at most 3,360 pu, far below 6e7. So V_t != 0.
* All other buses are banded, so V != 0 everywhere. Power-form KCL S_b = V_b conj(J_b) is then
  equivalent to current-form KCL, with device current J_b = conj(S_b / V_b) (zero at buses without
  devices).

**2. Eliminating the network.**

* Current-form KCL and the source equations I_s = (E_s - V_t(s)) / z_s give
  Y V = J + sum_s (E_s / z_s) e_t(s), where Y is the bus admittance matrix of lines, half-shunts and the
  source series admittances.
* **Invertibility.** Re(Y) = Laplacian(Re 1/z_l) + diag(Re 1/z_s), plus the shunts' real parts (0).
  Every line has r > 0 (smallest 7.5e-9 pu), so the tree Laplacian is positive semidefinite with only
  the constant vector in its null space. Re(1/z_s) > 0 at the terminals removes that. So Re(Y) is
  positive definite (Cholesky in 256-bit arithmetic, smallest pivot 16.9).
* Y is complex symmetric, so x^H Y x = x^H Re(Y) x + j x^H Im(Y) x, whose real part is > 0 for x != 0.
  Y is therefore nonsingular.
* Z = Y^-1 is computed in BigFloat (256-bit), with |YZ - I| = 1.7e-71. Then
  V = V0 + Z J, with V0 = Z sum_s (E_s / z_s) e_t(s) the no-load voltage, and
  **I_s = a_s + sum_b c_sb J_b**, with a_s = (E_s - V0_t(s)) / z_s and c_sb = -Z_t(s),b / z_s.
* a_s is the current at zero device injection: 3-15 A at these angles (column `no_load_term_A`).

**3. Bounds on J_b.**

* S_b lies on the segment P_b + j(Q0_b + q), q in [-q_max,b, q_max,b]. Hence |S_b| <= Smax_b, and
  arg S_b lies in the arc the segment subtends: under 180 deg, since P_b != 0 at every device bus here.
* |J_b| = |S_b| / |V_b| <= Smax_b / Vmin_b, the existing `:jbound` bound with Vmin = 0.95.
* Sector refinement: V_b lies in the disc of radius rho_b around V0_b. A first radius comes from the
  termwise bound, rho_b = sum_k |Z_bk| Jmax_k.
* Hence |V_b| >= max(0.95, |V0_b| - rho_b), and arg V_b is within asin(rho_b / |V0_b|) of arg V0_b.
* Since J_b = conj(S_b) V_b / |V_b|^2, arg J_b = arg V_b - arg S_b. J_b therefore lies in a circular
  sector: known centre, half-width asin(rho_b / |V0_b|) plus half the S arc, radius
  Smax_b / max(0.95, |V0_b| - rho_b).
* With those sectors, rho_b is recomputed through the support function (below), and the step repeats.
* Each step uses only containments proved by the previous one, so every radius stays valid. The disc
  radius falls from 0.0067 to 0.0038 pu in one step and then stays put.

**4. Three source-current bounds:**

* **naive:** KVL and the box alone, (|E_s| + 2.12) / |z_s| = 528 pu = 45,772 A. This is what the
  relaxation already had.
* **termwise:** the triangle inequality on the eliminated relation, |a_s| + sum_b |c_sb| Smax_b / 0.95.
  It gives 59-73 A.
* **sector, combined:** the linear relation kept whole.
  - Each c_sb J_b lies in the rotated and scaled sector, and a_s is fixed. So I_s lies in the Minkowski
    sum, and for every direction u, Re(conj(u) I_s) <= Re(conj(u) a_s) + sum_b h_b(u), where h_b is
    the support function of the sector's hull.
  - On a grid of N = 3,600 directions (spacing 0.1 deg), every vector w has a grid direction within
    0.05 deg, so |w| <= max_k H(u_k) / cos(0.05 deg). That makes the bound rigorous, not a sample.
  - Result: 31-39 A.
  - For contrast, the same sectors bounded device by device (triangle inequality per device) give
    54-67 A.

So the gain comes from keeping the combined relation. Loads and PV inject currents of nearly opposite
phase at noon, and the termwise bound adds them instead of letting them cancel.

**Numerical status.**

* Analytical validity: steps 1-4 hold for every state of the model above. The direction-grid
  discretization is covered by the 1/cos factor.
* Numerical verification:
  - Z is computed in 256-bit arithmetic and rounded to Float64.
  - The support sums run in Float64 with a 1e-9 relative allowance. That is far above the Float64
    summation error of a few hundred terms (about 1e-13 relative), but no interval arithmetic was used.
* The derived bounds are therefore independently checked numerical bounds, not exact-arithmetic
  proofs.

## Optimization-based tightening (OBBT): no gain

Each source's lifted |I_s|^2 was maximized over the restoration relaxation, starting from the analytic
bounds:

* band hard;
* P_s >= -rho_max;
* both nodal strengthenings;
* upper bound checked by `dual_bound`.

The checked upper bounds were 36.8-43 A, above the analytic 31-39 A, so no bound tightened and one
round ran per query.

The cause is numerical: the lifted |I_s|^2 is a voltage-moment expression divided by |z_s|^2 = 3.6e-5,
so a 1e-6 error in the moments moves it by about 3%. The solver's own maxima even exceed the imposed
cap by up to 8%: they are not solutions and are not used. **The operative bounds are the analytic sector
bounds.** OBBT and its dependency chain (round r uses the bounds of round r-1, starting from the
analytic ones) are kept in the code and the log.

## Results

Bounds, at 1.02 r* (every query in `source_current_bounds_bounds.csv`; they vary by at most 2.7 A across
the queries):

| Source | Envelope v2 (empirical) | 1.5x cap (empirical) | naive | termwise | sectors, separately | **sectors, combined** | No-load term | Validated dispatch |
|---|---|---|---|---|---|---|---|---|
| 292.5 deg, subs1 | 56.2 A | 84.3 A | 45,772 A | 68.0 A | 62.0 A | **35.63 A = 0.4115 pu** | 3.0 A | 35.0 A |
| subs3 | 51.9 | 77.8 | 45,772 | 67.5 | 62.2 | **37.37 A = 0.4315 pu** | 10.0 | 34.8 |
| subs4 | 45.9 | 68.9 | 45,772 | 63.9 | 59.1 | **33.61 A = 0.3881 pu** | 12.1 | 28.6 |
| 157.5 deg, subs1 | 56.2 | 84.3 | 45,772 | 70.2 | 64.1 | **37.82 A = 0.4367 pu** | 5.2 | 29.0 |
| subs3 | 51.9 | 77.8 | 45,772 | 69.8 | 64.5 | **34.09 A = 0.3936 pu** | 12.2 | 23.6 |
| subs4 | 45.9 | 68.9 | 45,772 | 59.5 | 54.7 | **31.16 A = 0.3598 pu** | 7.6 | 24.9 |

Validity domain of every derived bound: the stated model at that one fixed angle vector. The bounds use
no no-backflow, no sampled dispatch and no empirical cap.

Restoration test with the derived bounds and both nodal strengthenings (restoration relaxation: band
hard, minimize rho with P_s >= -rho, rho <= rho_max = exact rho + 1 kW):

| theta | Point | Exact witness: L (kW), rho (kW), equations | Relaxation status | rho relaxed (kW) | **Checked lower bound on rho (kW)** | Box term / rounding (kW) | Primal residual | Solve (s) | Outcome |
|---|---|---|---|---|---|---|---|---|---|
| 292.5 | 0.98 r* | +2.727, 0, 9e-11 | OPTIMAL | -0.0001 | -0.0133 | 0.028 / 5e-9 | 1.6e-7 | 6.6 (first, compilation) | feasible (validated); unresolved, as it must be |
| | 1.02 r* | -2.728, 2.728, 9e-11 | ALMOST_OPTIMAL | 1.408 | **1.090** | 14.0 / 1.5e-7 | 3.6e-4 | 0.6 | certified |
| | 1.04 r* | -5.455, 5.455 | ALMOST_OPTIMAL | 4.498 | **4.358** | 17.1 / 1.8e-7 | 2.2e-4 | 0.6 | certified |
| | 1.10 r* | -13.635, 13.635 | ALMOST_OPTIMAL | 13.121 | **13.109** | 17.4 / 2.5e-7 | 2.4e-5 | 0.8 | certified |
| | 1.25 r* | -34.086, 34.086 | OPTIMAL | 33.556 | **33.548** | 17.5 / 2.6e-7 | 1.5e-5 | 0.8 | certified |
| 157.5 | 0.98 r* | +2.866, 0, 4e-14 | OPTIMAL | -0.0001 | -0.0045 | 0.009 / 1e-9 | 6.3e-8 | 0.5 | feasible (validated); unresolved |
| | 1.02 r* | -2.865, 2.865 | OPTIMAL | 1.867 | **1.862** | 4.37 / 2.3e-7 | 8.6e-6 | 0.9 | certified |
| | 1.045 r* | -6.446, 6.446 | ALMOST_OPTIMAL | 5.395 | **5.380** | 4.38 / 2.2e-7 | 2.7e-5 | 0.8 | certified |
| | 1.10 r* | -14.325, 14.325 | OPTIMAL | 13.262 | **13.323** | 4.42 / 2.4e-7 | 5.1e-6 | 0.8 | certified |
| | 1.25 r* | -35.813, 35.813 | OPTIMAL | 34.749 | **34.739** | 4.37 / 2.5e-7 | 1.9e-5 | 0.9 | certified |

Notes:

* Every analytic bound takes 0.5 s (1.9 s the first time). The BigFloat Z takes 5.4 s once.
* The dual multipliers lie inside their cones in every case (PSD and sign violations 0).
* The exact restoration rho is from a validated local solve: a witness, an upper bound on the true
  minimum.
* The checked lower bounds never exceed it, which is consistent.
* At 157.5 deg, 1.10 r*, the checked bound (13.323) exceeds the solver's primal objective (13.262). The
  solver's primal point is slightly infeasible (residual 5e-6, amplified), so its objective is not a
  bound. The dual-checked number is.
* The validated dispatch column in the bounds table is the exact margin (= restoration) dispatch at that
  angle: feasible at 0.98 r*, backflowing beyond.
* Certificates hold for the individual tested angles only. They do not show that every point beyond
  1.02 r* is infeasible.
* The AC boundary r* remains locally measured.

## Falsification search (evidence, not proof)

At every query angle, each source current was maximized in the exact model with the angles fixed:

* over the bound's domain (no-backflow dropped), and over the original model (no-backflow kept);
* from the margin dispatch, the restoration dispatch, a cold start, and each other maximizer;
* every dispatch validated.

A validated current above a derived bound would show a bound wrong. See `source_current_falsify.log`
for the per-query maxima; the summary is filled in below.

**Result: no bound violated.** There were 2,418 validated maximizations: 930 over the bound's domain
(93 per query) and 1,488 over the original model at the two 0.98 r* points.

* The largest ratio of a validated current to its derived bound is **0.9966**. Each query's maximum
  lies 0.13-0.14 A below the subs1 bound and 0.4-4.6 A below the subs3 and subs4 bounds:

  | Query | subs1 found / bound (A) | subs3 | subs4 |
  |---|---|---|---|
  | 292.5, 0.98 | 35.52 / 35.66 | 35.08 / 37.05 | 30.12 / 33.13 |
  | 292.5, 1.02 | 35.50 / 35.63 | 35.21 / 37.37 | 30.28 / 33.61 |
  | 292.5, 1.04 | 35.49 / 35.62 | 35.28 / 37.52 | 30.37 / 33.85 |
  | 292.5, 1.10 | 35.45 / 35.59 | 35.48 / 38.00 | 30.63 / 34.56 |
  | 292.5, 1.25 | 35.37 / 35.50 | 36.03 / 39.19 | 31.33 / 36.36 |
  | 157.5, 0.98 | 37.62 / 37.75 | 33.70 / 33.82 | 29.70 / 30.91 |
  | 157.5, 1.02 | 37.69 / 37.82 | 33.82 / 34.09 | 29.78 / 31.16 |
  | 157.5, 1.045 | 37.73 / 37.87 | 33.89 / 34.39 | 29.84 / 31.32 |
  | 157.5, 1.10 | 37.83 / 37.97 | 34.06 / 35.06 | 29.96 / 31.66 |
  | 157.5, 1.25 | 38.11 / 38.24 | 34.58 / 36.88 | 30.32 / 32.60 |

  So the derived bounds are close to sharp: within 0.4% for the largest source, and within 0.4-16% for
  the others.
* With no-backflow kept, validated dispatches exist only at the 0.98 r* points. At all eight outside
  points, no maximization returned a no-backflow-feasible dispatch. That is consistent with the
  certificates, though it is not evidence of its own.
* This is a local search. It can miss a larger maximum, so it does not prove the bounds; the derivation
  above does. A global check is possible: Gurobi's spatial branch-and-bound on the single-hour
  maximization would give a proven maximum with a gap. It was not run.

## Where this leaves the question

* For the source-only formulation, at these angles, the empirical-cap assumption is removed. Bounds
  derived from the original equations alone (network elimination, the band and device limits, plus one
  lemma that excludes zero terminal voltage via the model's box) are tighter than the empirical caps.
* With the nodal strengthenings, they certify no-backflow impossible at every tested point from
  1.02 r* outward, with checked lower bounds of 1.09-34.7 kW on the restoration objective.
* The certificates are independently checked floating-point certificates.
* What made this work: at fixed angles the source currents are nearly determined by the network and the
  fixed loads. Only PV reactive power and the voltage moves them.
* The empirical caps were looser because they spanned all hours and free angles.
* Not done, as instructed: no larger feeder, no 24-hour cyclic-battery day, no paper edits.
* With batteries dispatched, J_b at battery buses gains a bounded real-power range. The same derivation
  applies, with wider sectors.
