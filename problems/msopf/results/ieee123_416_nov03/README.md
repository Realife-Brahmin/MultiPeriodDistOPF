# ieee123 at 4.16 kV (2026-11-03, 309 PC)

`ieee123_5poi_1ph` at its own 4.16 kV (until today it ran at 20 kV with the same line impedances), source impedance
the feeder's substation transformer (5000 kVA, 1 + j8 %: 0.0346 + j0.2769 ohm), PV and batteries from `ieee123C_1ph`
(51 units each, 506.6 kW, 44% of peak load). T = 24, PV reactive power dispatched, exact certificate check on.
`run.sh <label> <ratings> <script:phase:dir> ...`, with `SYS=ieee123_s13|ieee123_s14`.

Edges (deg; A no batteries, AD batteries over the day) and certified edges (C, CD; ratio to the edge):

| pair | ratings | A + / - | C + / - | AD + / - | CD + / - |
|---|---|---|---|---|---|
| S1-S3 (bus 120, trunk) | none | 4.118 / 3.487 | | 5.769 / 4.631 | |
| S1-S3 | conductor tiers | 4.118 / 3.487 | 4.742 (1.152) / 3.999 (1.147) | 5.769 / 4.631 | 6.454 (1.119) / 5.224 (1.128) |
| S1-S4 (bus 116, 1-ph lateral) | none | 7.509 / 5.356 | | 10.602 / 7.105 | |
| S1-S4 | conductor tiers | equal angles infeasible (hour 6) | | | |

Conductor tiers (`uniform_ratings.jl tiers`): 336 ACSR 530 A (3-phase), 353 A (2-phase), 1/0 ACSR 76.7 A (1-phase,
230 A x 1/3), 1/0 AA cable 195 A. Substation 4 connects at the end of a 1-phase 1/0 lateral: without ratings its
connection carries 271 A (406 A with batteries) against 76.7 A, so S1-S4 has no feasible dispatch under the ratings.
S1-S3's largest current is 332 A (441 A with batteries, line l2 at the head), under the 530 A trunk: the ratings do
not move its OPF edges, only tighten its certificates.

`old_fleet_blanket686A/`: the earlier pass on the 5poi fleet (PV 770 kW, batteries 1318 kW) with every line at
686 A (1.2 x the largest current): S1-S4 edges 4.52 / 3.68 (10.03 / 7.33 with batteries), certified at 1.27-1.38x.

## Tighter certificates (2026-11-03, user: 15% too loose)

The relaxation beyond the edge booked its fictitious losses (51.6 kW at 1.10x on +subs3, `relaxation_slack.jl`) on the
lines between S1 and S3, with |I|^2 moments up to the 530 A rating against ~0.5 pu first moments. `FB_DERIVED=1`
bounds every line's and source's current by the model's own response at the point's angles (current_bounds.jl,
lines = true): limit = min(rating, derived). `lsec_only/`: the lsec cut alone (12.8% / 14.2%).
`derived_bounds/` (cuts kcl0, jbound, lsec + derived bounds): idle certified at 1.019x / 1.014x the edge, with
batteries over the day 1.091x / 1.078x (batteries at +-rating in every hour leave the trunk and source bounds near
their ratings). All certificates exact-checked.
