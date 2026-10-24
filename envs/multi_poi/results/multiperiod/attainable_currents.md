# Attainable currents in the exact AC model, and empirical limits in the relaxation (2026-10-22)

Driver `root_level/attainable_currents.jl` (envs/msopf_sdp); results `attainable_currents_runs.csv`
(every maximization), `attainable_currents_lines.csv` (largest current per element over every
validated schedule), `attainable_currents_limits.csv` (exploratory limits), log `attainable_currents.log`.

## Units

`ieee123_5poi_1ph` is a balanced single-phase equivalent of a 20 kV (line-to-line) feeder, solved per
phase: buses at 11.547 kV line-to-neutral; every load is a third of the three-phase load (S3P is 13.33 kW;
the deck totals 1,163 kW per phase, a third of the IEEE 123's 3.49 MW). Quantities are per phase.

| Base | Value |
|---|---|
| S_base | 1,000 kVA per phase (3,000 kVA three-phase) |
| V_base | 11.547 kV line-to-neutral (20.00 kV line-to-line) |
| I_base | S_base / V_base = 86.603 A (phase current) |
| Z_base | V_base^2 / S_base = 133.33 ohm (the 0.80 ohm source impedance is 0.006 pu) |

Check: replaying the largest-current dispatch in OpenDSS (imports match to 1.8e-5 kW), line L115 carries
27.66 A in OpenDSS against 0.3193 pu x 86.603 A = 27.65 A in the model; L116 28.22 against 28.21 A.

A demand base would read better: in hour 12 the feeder's apparent demand is 1,310.8 kVA per phase, so
S_rep = 1,310.8 kVA gives I_rep = 113.52 A (= 1.3108 pu of the solver base), and total demand is 1 pu.
Below, currents are given in A, in solver pu, and in demand-base pu (a reporting conversion only).

## Largest validated currents, hour 12 alone

Hour 12 (load factor 0.9873, PV factor 1.0), batteries idle, PV reactive power dispatched, every model
constraint hard, no current limit. Substation 1 is the reference (0 deg); substations 3 and 4 free in
[-5, +5] deg (never active). Starts: equal angles (cold), the hour's inside/outside boundary dispatches
(only feasible ones used), and 95%-of-reach dispatches at 0, 90, 225 deg. Local maxima of a nonconvex
problem: the largest found, not upper bounds.

| Maximized | Largest found | Other local maxima | Angles d3, d4 (deg) | Sources at that point, P + jQ (kVA) | Losses | Voltages (band; source terminals) | Equations | Solve |
|---|---|---|---|---|---|---|---|---|
| source 1 | 44.0 A = 0.508 pu = 0.388 demand pu | 38.2 A | -0.167, -0.197 | S1 379.3 + j373.2 (532 kVA, 44.0 A); S3 0.0 + j386.3 (31.9 A); S4 0.0 + j359.5 (29.7 A) | 0.78 kW | 1.0464-1.0479; 1.0478-1.0479 | 7e-10 | <= 0.03 s (3.7 s first) |
| source 3 | 41.2 A = 0.476 pu = 0.363 | 35.0 A | +0.159, -0.057 | S1 0.0 + j427.6; S3 379.3 + j324.3 (499 kVA, 41.2 A); S4 0.0 + j367.1 | 0.79 kW | 1.0463-1.0481; 1.0476-1.0481 | 7e-10 | 0.03 s |
| source 4 | 36.8 A = 0.425 pu = 0.325 | 30.3 A | -0.027, +0.220 | S1 0.0 + j461.9; S3 0.0 + j423.5; S4 379.7 + j234.1 (446 kVA, 36.8 A) | 1.12 kW | 1.0462-1.0486; 1.0474-1.0486 | 7e-10 | 0.03 s |
| line L1, L2 (source 1's path) | 44.0 A = 0.508 pu | 38.2 A | as source 1 | as source 1 | 0.78 kW | | 7e-10 | 0.03 s |
| line L119 (source 3's path) | 41.3 A = 0.476 pu | 35.0 A | as source 3 | as source 3 | | | 7e-10 | 0.03 s |
| line L116 (into source 4's POI) | 37.7 A = 0.436 pu | 29.8 A | -0.027, +0.217 | S4 379.7 + j214.7 | 1.12 kW | | 5e-10 | 0.04 s |
| line L115 | 36.5 A = 0.422 pu | 29.2 A | -0.027, +0.217 | as L116 | | | 5e-10 | 0.04 s |
| line L110 | 34.1 A = 0.394 pu | 28.3 A | -0.027, +0.217 | as L116 | | | 4e-10 | 0.04 s |
| line L46 | 1.2 A = 0.014 pu | none | -0.027, +0.220 | | | | 7e-13 | 0.03 s |
| line L111 | 0.01 A | none | | | | | 3e-12 | 0.02 s |

* At every maximizer one source carries almost all the real power (about 379 kW) and the other two sit
  exactly at no-backflow (0 kW). All three export 230-460 kvar, which the PV inverters absorb: the large
  currents are reactive circulation, not real power. Losses stay near 1 kW. Voltages and source-terminal
  voltages stay inside 1.046-1.049 pu. Nothing unexpectedly large appeared.
* Start dependence: two local maxima per element (for example 38.2 and 44.0 A for source 1).
* The SDP experiment's fictitious-loss lines (L115, L116, L110, L111, L46) carry at most 37.7 A.
* Over every validated schedule, including the 24-hour ones whose other hours carry more load and no PV,
  the largest currents are 56.3 A (L1, L2, source 1), 51.9 A (L119, source 3), 50.4 A (L54) and 45.9 A
  (source 4); 24 of the 132 lines exceed 40 A at some point.

For scale: the relaxation's fictitious current at a source terminal reaches a second moment of about
2,000 pu^2, about 45 pu = 3,900 A, roughly 100 times the largest current any validated dispatch carries.

## Exploratory limits (empirical, cushioned)

Limits are 1.5 or 1.2 times the largest current observed for each source and each line (over every
validated schedule), with no line below 0.05 pu (4.3 A). They are empirical: they keep every observed
dispatch, which does not show they keep every feasible dispatch. No equipment rating is known (the deck
carries only OpenDSS's default 400 A), and no proven upper bound is available. Each limit is applied to
the exact model (|I|^2 <= I_max^2) and to the relaxation's squared-current moments (for a source,
|E - V_s|^2 / |z_s|^2 lifted). Hour 12, the points 2% inside and outside the hour's boundary.

| Limits | Point | Exact L (kW) | Relaxation, no additions: status, optimum, fictitious losses | Relaxation with zero-injection and device-current blocks | Outcome |
|---|---|---|---|---|---|
| none (previous test) | inside / outside | +2.727 / -2.727 | optimal, 1,937 / 1,948 kW, 5,385 / 5,373 kW | about 1,746 / 1,786 kW, about 4,860 / 4,976 kW | feasible / unresolved |
| sources 1.5x (84.3, 77.8, 68.9 A) | inside | +2.727 | optimal, 195.3 kW, 728 kW | optimal, 7.8 kW, 12.9 kW | feasible |
| | outside | -2.728 | optimal, 195.6 kW, 728 kW | optimal, 2.3 kW, 12.9 kW | unresolved |
| sources and lines 1.5x | inside | +2.727 | optimal, 5.9 kW, 3.2 kW | optimal, 4.8 kW, 2.0 kW | feasible |
| | outside | -2.728 | optimal, 0.45 kW, 3.2 kW | solver: infeasible | unresolved |
| sources and lines 1.2x | inside | +2.727 | optimal, 4.5 kW, 1.8 kW | optimal, 4.0 kW, 1.3 kW | feasible |
| | outside | -2.728 | solver: infeasible | solver: infeasible | unresolved |

* The validated inside dispatch lifts into every limited relaxation (<= 4.3e-9) and satisfies every limit.
* Fictitious losses fall from about 5.4 MW to 728 kW with source limits alone, to 2-13 kW with source
  limits plus the bus blocks or line limits, and the relaxed margin at the outside point from 1.9 MW to
  0.45 kW (or infeasibility). The empirical current envelope is what the relaxation was missing.
* No point is certified infeasible. Where Clarabel reports infeasibility, its ray does not survive the
  independent check: for sources and lines 1.2x, no additions, the ray's certificate value is 5.0e6 against
  a box term of 1.9e10, and tighter infeasibility tolerances (1e-10, 1e-12) give 1.8e7 against 6.8e10.
  A numerical indication of infeasibility under empirical limits, not a certificate; and even a validated
  one would hold only for the model with those empirical limits added.
