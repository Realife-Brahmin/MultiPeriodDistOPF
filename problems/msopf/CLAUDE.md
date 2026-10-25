# MSOPF (multi-source multi-period OPF): context for Claude

Claude Code loads this file for work under `problems/msopf/`. It holds only MSOPF; the
repository-wide rules (test decks, OpenDSS validation, commits, Kamiak) are in the root
`.claude/CLAUDE.md`, and SPAR-DDP has its own notes (`problems/mpopf/spar_ddp/CLAUDE.md`, `.agents/AGENTS.md`). MSOPF and
SPAR-DDP share only the test decks in `test_networks/` and the OpenDSS checks; their rules differ
(e.g. terminal SOC is hard here, soft in DDP and tADMM). Keep this file updated when a session
establishes something a later session on any machine would need.

## MSOPF quick facts (2026-09-28)

- Setup, drivers and decks: [`README.md`](README.md) in this folder.
  The full-angle drivers run in **`problems/msopf`** (tracked Manifest), MSOPF's own
  environment since 2026-10-19: an exact copy of `problems/mpopf/tadmm`'s pinned files, where
  they ran before, so earlier results used the same versions. Not `problems/msopf`
  (older branch-flow scripts, no Manifest).
- The paper is the separate **`PESGM2027_Multi-Source-Multi-Period-OPF`** repo. Its
  Section II documents `solve_full_angle` in `problems/msopf/full_angle_pf.jl`, so
  keep the two in sync. Do not confuse it with the rejected PESGM2026 repo.
- Decks `ieee123_5poi_1ph` and `small2poi_1ph` carry a realistic source impedance,
  Zs = 0.032 + j0.799 ohm (cited in the deck comments). `ideal_sources = true` gives
  the old stiff sources.
- Battery model: one variable `P_B` per battery, lossless, with an optional
  `C_B * P_B^2` cycling cost (the user's standard). Feasibility means the network
  equations, no-backflow (`P_Subs >= 0`, the key constraint) and the 0.95-1.05 pu
  band on buses not set by a substation.
- **Multi-period MSOPF (2026-10-11):** `problems/msopf/multiperiod_full_angle.jl`
  and `root_level/multiperiod_study.jl`; results in
  `problems/msopf/results/multiperiod/`. User decisions: the angle difference is
  constant over the day; **terminal SOC is HARD in MSOPF, `B[T] == B[0]`** (an
  MSOPF-only exception to the soft terminal SOC rule, which still holds for DDP
  and tADMM: feasibility is the question here, and a soft terminal lets batteries
  end the day fuller, so its window grows without bound with battery size); `C_B`
  by the DDP sizing rule on this deck's fleet; PV real power fixed, PV reactive
  power an optimization variable. **No new grid assets, ever** (no phase shifters,
  soft open points or added storage as a "solution"): only dispatch of what
  exists. **PV reactive power dispatched (`pv = "qfree"`) is the only PV configuration
  of interest** (user, 2026-10-20): no unity-power-factor runs, panels or comparisons in
  new work or the paper. New decks get our own PV and batteries, never the model's
  built-in DERs: about 10% of buses each, partly overlapping, each fleet 30-50% of the
  load (e.g. `test_networks/ieee9500_1ph`), all written by `convert_to_1ph.jl` (switches and
  regulators as 1e-3 + j1e-3 ohm lines, never merged buses; checked by `check_1ph.jl`).
  Provenance (user, 2026-10-20): `ieee2522C_1ph` is the single-phase, primary-only version
  of the IEEE 8500 (so `ieee8500_1ph` is a second conversion of the same feeder);
  `large10kC_1ph` is purely synthetic. All existing decks are primary-only (one voltage
  level, loads on primary buses, no transformers).
  Public OpenDSS systems to convert from are in
  `test_networks/sources/` (`bash tools/fetch_test_systems.sh`). OpenDSS only (no CYME, CIM).
  **Commit and push MSOPF work as it is done**, code and paper alike.
  Headline: over a 24 h day batteries roughly double the feasible angle window
  (subs3+4, PV on: 0.51 -> 1.05 deg), about half of what hour-by-hour solves
  suggest. What limits it: a constant angle difference drives a constant
  circulating power (826 kW/deg here = V^2/|Z_loop|), batteries only move energy
  between hours, so the sustainable difference is about Pbar*|Z|/V^2 rad (the
  constrained substation's mean import), and more storage or power beyond what
  flattens the day adds nothing. Method notes that cost hours are in
  `problems/msopf/README.md` (never ask Ipopt to prove infeasibility; release
  each Ipopt problem; stop Julia by PID).

### MSOPF handoff (written 2026-10-12 on the Surface, for a session on any machine)

Branch `msopf-oct21-scoping-study` since 2026-10-21 (`msopf-oct09-further-investigations` was merged into master in PR #174 and deleted); re-resolve from the remote first.
Everything below is one driver, `problems/msopf/root_level/multiperiod_study.jl`,
run from the repo root as
`julia --project=problems/msopf problems/msopf/root_level/multiperiod_study.jl <stage ...>`;
results are tracked in `problems/msopf/results/multiperiod/`.

- **Done, ieee123 (five substations):** stages `main horizon zs duration pairs
  three five bound estimate tables`, all under the hard terminal SOC.
- **Done, the limit is explained and checked:** stage `estimate`. With batteries
  idle, the imports at equal angles and their slope per degree give the window:
  the worst hour without batteries, the daily mean with them. It matches the
  measured windows to 0.001 deg when only batteries are dispatched (PV off or
  unity power factor). With PV reactive power free the measured upper edge is
  0.07 deg wider, because reactive power moves 20-40 kW between the substations.
- **Done, several substations: the window is a simplex (stage `simplex`,
  2026-10-19, `simplex_check*.csv`).** Each substation's import is close to linear
  in the angles (the angle-to-import sensitivities, an N x (N-1) matrix), so each
  no-backflow limit is a half-space and N substations give a simplex in the N-1
  angle differences (three: a triangle). Predicted from one idle solve plus
  central differences, against measured windows on ieee123 at T = 24:
  - PV at unity power factor: without batteries exact in all 26 directions
    (three-substation rays and five-substation cuts, within 1e-4 deg); with
    batteries (daily mean) exact in 24 of 26 (within 4e-4 deg). The two
    exceptions both raise subs4: there the 1.05 pu upper voltage limit binds first,
    with batteries at full power, a face the backflow-only prediction lacks
    (probed: beyond the measured edge the violation is overvoltage, not backflow).
  - PV reactive power free (the committed windows): the prediction is
    conservative everywhere, exact on the subs1 and subs3 faces and 0.04-0.05 deg
    short on the subs4 face, which reactive dispatch pushes outward.
- **Done, ieee2522 (stage `ieee2522`, finished 2026-10-12):** a second
  substation is added in memory at buses 947, 996, 1120, 1216 (along the main
  trunk) and 546 (farthest bus), with Zs = 0.311 ohm, T = 12. Rows are in
  `large_windows.csv`, merged by (system, bus, T), so a run can be split with
  `LARGE_AT=bus[,bus]` and `LARGE_T=n`; the table is in `tables.tex`. All five
  cost-optimal days replay in OpenDSS to 2e-5 kW, with no failed edge checks.
  - Windows are several degrees wide (ieee123: about half a degree) and widen
    with distance: with batteries, [-8.3, +3.7] deg at bus 947 (3.4 ohm away),
    [-12.7, +10.7] at bus 1120 (6.1 ohm). Batteries widen them about 1.3x.
  - **The estimate only holds while no-backflow is the binding limit.** At 947,
    996 and 1120 both edges sit at zero import and the daily-mean estimate is
    4-19% narrow. At 1216 (10.1 ohm) and 546 (14.7 ohm) the upper edge is set by
    the voltage band with hundreds of kW of import to spare ([-19.2, +3.4] and
    [<= -20, +2.2] deg), and the estimate is far off (+23.8 against +3.4).
  - The whole run took 75 minutes on the Surface, of which about 10 were CPU:
    the rest was swapping. With free memory expect about 2 minutes per bus.
- **Done, 10k buses (309 PC, 2026-10-12/13): use stage `large10kd`, not
  `large10k`.** User decisions: the C deck is not a realistic feeder (1032 MW,
  and R/X = 7 on every section), so MSOPF uses **large10kD** = large10kC with
  loads, PV and batteries divided by 20 and line sections at R/X = 1 (R kept,
  X raised to 0.0175 ohm), built in memory by the stage; T = 24 is the standard
  horizon; feasibility tolerance 1e-5 as in DDP; Ipopt with HSL MA57. Run as
  `IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=<hsl>/libma57.dll;linear_system_scaling=none" MULTIPERIOD_LEAN=1 LARGE_T=24 julia ... multiperiod_study.jl large10kd`.
  - Second substation at bus 1120, 2.475 ohm from the first; **angle-to-import
    sensitivity** (the user's name for kW of import per degree) 855 kW/deg.
    Windows, no batteries / batteries over the day: T = 24 [-0.236, >= +20] /
    [-0.277, >= +20] deg; T = 48 -0.233 / -0.274; T = 96 -0.231 / -0.273. The
    lower edge is no-backflow and the import estimate gives it to 0.008 deg;
    there is no upper edge inside the 20 deg limit (voltage 1.033-1.050 pu), so
    only line ratings, not modeled, would close it. OpenDSS replays to 3e-4 kW.
  - Cost: T = 24 1318 s and 5.9 GB; T = 48 3567 s, 12.7 GB; T = 96 6745 s,
    16.7 GB (with `MULTIPERIOD_LEAN=1`; the PC has about 16 GB free).
  - How edges are found there: `window_secant`, a secant on the import margin
    (`objective = :margin`), six solves per window. **Do not use the
    maximize-angle solve on 10k buses** (it stopped early, hit the iteration
    limit or went infeasible), nor Gurobi (no feasible point in 15 minutes).
  - On large10kC itself nothing binds inside +-20 deg, with or without
    batteries (reactive circulation, and its losses raise the import on both
    sides); not committed as a result.
  - Still easy to get wrong: the sources stay stiff; never call
    `source_separation` on 10k buses (dense; the stage uses `separation_sparse`).
  - Not written to `test_networks/` yet: large10kD exists only in the stage. The R/X
    review (`rx_ratio_report.md`, `rx_ratio_refs.bib`) was not on the 309 PC.
- **Done 2026-10-20 (Kamiak), three substations on the large decks, T = 24, PV q dispatched**
  (`three_large_windows.csv`, figures `fig_three_<label>`, `plot_three_large.jl`, hourly
  sensitivities `three_sens_*.csv`; evidence `results/kamiak/three_oct19`, `three9500_oct20b`,
  `three9500xf_oct20`):
  - **ieee9500** (its own three substations, the five ties between feeders closed): nested
    triangles; all 32 edges no-backflow, 3-11% beyond the sensitivity prediction (conservative).
  - **ieee2522** (947, 1120): prediction good on backflow faces; 67.5-135 deg cut by 1.05 pu.
  - **large10kD** (1120, 1060): on the prediction (idle), 1.3-2% inside (cyclic); a thin wedge.
  - **Secondaries** (`ieee9500xf_1ph`, service transformers modeled): edges move outward
    1.2-1.8%, still no-backflow, so primary-only decks are conservative by under 2%.
  - ieee2522 at 1120 at T = 12/24/96: the window narrows 1.3-3% from 24 to 96 (8.7 h job).
  - Edges now by `reach_secant` (import-margin secant, about 4 solves per point instead of 11 for
    the bracket, warm-started from the neighbouring direction; `MULTIPERIOD_WARM`).
  - SMART-DS P1U solves (5,064 buses, 70-86 s per day) but has no ties between its five
    substations below 69 kV: no window until the user decides how they connect.
- **Done 2026-10-20, four substations (stages `four_ieee123a/b`, `four_ieee2522`):** the region
  is the predicted tetrahedron on ieee123 (planar faces, normals within ~1 deg, the predicted
  substation binds 298/312 edges at `MULTI_DIRSET=fib78`); substation 4's face is pushed out by PV
  reactive dispatch; on ieee2522 the 996 face bends where voltage binds (78-direction rerun was
  still running on Kamiak, jobs 31328984-98). Results `multi_windows*.csv`, `check_four.jl`.
  Interactive viewer: `results/multiperiod/viewer/four_regions.html` (README there; public link
  via raw.githack; also a private claude.ai artifact). Rebuild with `four_viewer_data.jl fib78`.
- **Done 2026-10-21, four substations on the 2705-bus feeder** (`ieee9500_1ph`, stage `four_ieee9500`:
  its three substations, ties closed, plus one at l2804256 chosen by `choose_fourth_9500.jl`; evidence
  `results/kamiak/four9500_oct20`). Idle: 66 no-backflow edges all at or beyond the simplex; 12 voltage-bound
  (1.05 pu upper limit, next to the 1.05 pu sources), up to 65% short. Substation 4's face is planar and
  33% beyond the prediction. In the paper and to collaborators it is the **2705-bus feeder** (cite Bobo 2021
  for the per-phase convention, Sadnan2021 for the 2522-bus lineage). Viewer on GitHub Pages:
  https://realife-brahmin.github.io/MultiPeriodDistOPF/ (served from `docs/`; copy the rebuilt page
  there).
- **`MULTI_VMARGIN=20000`** (2026-10-21): margin solves count voltage room, so voltage-bound edges take the
  secant in 4-6 solves instead of 40-55 slow ones or a capped hour. Use it for every new edge search.
- **Solver hygiene (2026-10-20):** solves past an edge are 10x slower than feasible ones and a
  cold margin solve can take over an hour, so searches use warm-started margin solves, a per-solve
  time cap (8x the reference solve, `MULTIPERIOD_CAP`), and park capped solves in
  `parked_solves.csv`; every solve logs one `[solve]` line.
- **Paper (PESGM2027) as of 2026-10-21:** Section II-E on the angle-to-import sensitivity;
  Fig. 3 two-source windows, Fig. 4 three regions (123 / 2705 / 2522), Fig. 5 four-source faces
  seen straight on (123 / 2705, rendered from `four_faces_*.csv` of `four_face_data.jl`), Table I
  cross-system with two-, three- and four-source rows; battery case named "dispatched" (the
  cyclic constraint keeps its name); one-line captions; 7 pages (not trimmed yet).
- **Paper: edit only when asked.** Since 2026-10-12 the user also has ChatGPT
  editing the PESGM2027 paper; it rewrote the draft as prose with booktabs
  tables (one file per table, `tables/mp_*.tex`) and vector figures rendered by
  `scripts/render_figures.py`. Always pull first, match that style, and do not
  regenerate its tables or figures from this repo's `tables.tex` or PNGs.
  The paper has the ieee123 results and, since commit `c160f83` (added at the
  user's request), the ieee2522 subsection and `tables/mp_ieee2522.tex`. It
  has nothing on large10k.
- **Still open:** per-line ampacity ratings (not modeled; the user is collecting
  them), and the substation voltage setpoint (1.05 pu sits on the band edge).
- **Test systems, 2026-10-24 (branch `msopf-oct24-test-system-edition`; user decisions):**
  - **Secondaries in, triplex out.** New studies use the `*xf_1ph` decks: each service transformer is a line
    to its secondary bus, where its loads sit; the 120/240 V triplex lines are not modeled. DERs stay on
    primary buses (the utility's network), as in the primary-only deck of the same system.
  - **The voltage band applies only at the secondary (customer) buses** (ANSI C84.1 is a service-point limit);
    the primary has no band. The lower limit may be set below 0.95 for these decks, chosen from first runs
    (at nominal load with no DERs they reach 0.85-0.96 pu). Implemented through `band_buses` in the solver; the study driver now selects the customer buses automatically (2026-10-25).
  - **Up to four substations** (two, three, four); no five-source cases.
  - SMART-DS P1U (`smartds_p1u_xf_1ph`, 6,937 buses): its normally-open ties are published as disabled
    switches, now read by `convert_to_1ph.jl` (39 ties; 6 between substations, chain hs0-hs2-hs1-hs3). Use
    `Master_4sub_typeIII.dss` (hs0 on one transformer; the usual problem) or `Master_4sub_typeII.dss` (hs0's two
    transformers, Type II, hs3 out). The four substations also exist alone (`smartds_p1u_hs0..3_1ph`).
  - Other decks with secondaries: `ieee9500xf_1ph` (about 3,980 buses, 3 substations), `ieee8500xf_1ph` (3,698;
    one substation, -5 to -6% import mismatch inherited from the primary conversion), `epri_ckt24xf_1ph`
    (2,111; 34.5/13.2 kV), `epri_k1xf_1ph` (962). Single-substation decks need added sources for MSOPF.
  - The restructure (PR #178) untracked `problems/msopf/Manifest.toml`; it is identical to
    `problems/mpopf/tadmm/Manifest.toml` (same Project.toml): copy it locally until it is tracked again.


- **Continuation, 2026-10-25:** the study driver now propagates customer-only `band_buses` through every problem and voltage diagnostic. New stages: `four_ieee9500xf`, `four_smartds_p1u_typeiii`, `four_smartds_p1u_typeii`. Both requested T=24 sensitivity runs (9500xf and P1U Type III) completed, seven solves each, all 336 hourly OpenDSS comparisons passing (worst 5.2e-5 kW, 1.35e-9 pu). Independent export validator and PNG/SVG predicted-region figures are included. Evidence and pinned environment: `results/test_networks/sensitivity_oct25/README.md`. **Predictions only**; the 78-direction measured edge campaign still awaits Kamiak login/queue check. Type II wiring is tested but its sensitivity campaign is not run. The regression is `tests/customer_band.jl`; the serial launcher is `root_level/run_test_edition_sensitivities.sh` (fresh output root required for reruns).

### 2026-10-25: secondary-feeder OPF evidence and paper update

- New measured study: `root_level/secondary_8500_study.jl`, launcher
  `run_secondary_8500_screen.sh`, input checker `check_secondary_8500_inputs.py`,
  validator/figure/viewer builder `publish_secondary_study.py`. Evidence and exact
  scientific qualifications: `results/secondary_8500_oct25/README.md`.
- Matched new IEEE8500 primary/secondary conversions share primary impedances,
  DERs and aggregate loads. Customer 0.95–1.05 pu limits vs primary attachment
  proxies; no primary transit band. Added source sites m1047401/m1108258/m1142880,
  with original substation delivery equivalent. This is not the old 2522 conversion.
- Full T24 2-source search has one verified positive cyclic bracket and one
  unresolved negative direction. T24 3-source has eight idle and five cyclic
  sign brackets, three cyclic timeouts. Never convert unresolved zero saved radii
  into boundaries. Some positive-margin stopping brackets are wide; retain them.
- Full T24 equal-angle checks exist for both battery cases, primary/secondary
  2/3/4 sources and SF4. SF uses native source sites and our constructed Type III
  ties, not a published operating configuration. All successful reportable OPFs
  pass strict OpenDSS replay and original SOC/device checks.
- T6 is explicitly exploratory: legacy profile resampling changes daily PV energy,
  so T6 vs T24 is not a resolution comparison. SF4 completed twelve T6 axis
  searches; six axes cannot establish tetrahedron faces or volume. New sparse
  viewer tabs suppress face fits and disclose their coverage and prediction caps.
- Paper updated and pushed in PESGM2027_Multi-Source-Multi-Period-OPF main:
  initial 9e0a123, full-day SF baseline and legibility follow-up 3f7ec8b. New text
  is in Section IV-F, vector figures `three_source_T24.pdf` and `sf_axis_T6.pdf`,
  exact copied inputs in `figure-data/secondary_oct25/`. main.pdf compiled and
  pages inspected. Direct Overleaf sync was not verified (no browser access).
- #180 remains the broad test-system/modeling PR. #181 is deliberately a separate
  docs-only interactive-viewer publication PR; do not merge all of #180 just to
  publish figures. Inherited untracked `envs/` outputs are unrelated and untouched.
- Final publication: #181 merged to master as adb28724 (docs only). Both new
  four-source T6 screens are present, with 12 SF and 6 IEEE8500 cyclic axis points.
  Final validation: 15 complete cases, 320 successful OPFs replayed; the viewer's
  seven tabs and face views passed 192 finite-coordinate rendering checks.
  Research archive commit 696c4512; paper evidence follow-up a709042. The broader
  modeling PR #180 remains open, separate from published documentation.

### 2026-10-25: comparison and encoding correction

User correctly identified inconsistent new figure styling and corrupted viewer labels.
Repeated default-codepage read/write cycles caused mojibake; the generator now reads
UTF-8 explicitly, HTML declares UTF-8, original readable template labels are restored,
and source colors match the paper palette. Two rebuilds are byte-stable; live viewer
checked for corrected labels/palette after docs-only PR #182, merge 008e980d.
Paper commit 9cc09cc replaces the new diagnostic plots with two-source comparisons,
full-prediction plus zoom three-source panels, and four-source in-plane face views.
It imports the existing paper vector-renderer helpers and palette. The original
connected eight-ray polygon was misleading and is superseded, not a measured contour.
Legacy 2522 and new secondary conversion differ in source sites and other assumptions;
these are presentation comparisons, not controlled service-transformer ablations.
Four-source new cases remain sparse T6 screens, not matched to legacy T24/78 rays.
No numerical CSVs changed and no new OPFs were launched for this correction.
