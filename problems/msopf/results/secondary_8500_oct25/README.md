# New secondary-feeder OPF evidence, 25 October 2026

## Scope and model

`secondary_8500_study.jl` runs full-angle JuMP/Ipopt OPF with MA57, constant
source EMF angles over a 24-hour day, fixed PV real output, free inverter reactive
power, existing batteries and hard cyclic SOC. No Farkas or SDP study is included.

The IEEE8500 primary and secondary conversions have identical primary lines and
DER definitions; secondary-load aggregation differs by at most 6e-5 kVA. The
voltage contract is 0.95–1.05 pu at explicit secondary customers, or at their primary
attachment proxies in the matched primary-only control. Primary transit buses
have no band. This is a new matched conversion, not exactly the older 2522 case.
`input_check.json` records the checks and the separated source siting rule. Added
source branches use the same original substation transformer equivalent,
0.0759998+j0.877025 ohm. EMFs are stiff and fixed at 1.05 pu; that delivery
impedance remains in the OPF and OpenDSS models.

San Francisco uses the 6937-bus SMART-DS P1U conversion with 1873 delivery buses,
original source locations and a constructed four-source Type III topology:
original normally-open ties closed, fifth source at a shared substation disabled.
The configuration is our study variant, not a published operating condition.

## Reportable evidence

- `full_T24/secondary2`: positive cyclic endpoint 9.52828 degrees, negative-margin
  endpoint 9.57297 degrees. The negative-direction search timed out: **not a
  complete interval**, and its zero saved reach is **not** a physical boundary.
- `full_T24/secondary3`: eight idle and five cyclic directions sign-bracketed;
  three cyclic directions unresolved. Positive delta2 idle 3.76928–3.77309 degrees,
  versus 77.07298 predicted; cyclic 9.94889–9.99289 degrees. These voltage cuts
  occur with positive imports. Bracket widths are retained; a stopping criterion
  based on small positive margin does not guarantee a narrow radial bracket.
- `baseline_T24`: full-day equal-angle checks for primary/secondary 2,3,4 sources,
  and SF4. Secondary 3/4 and SF4 yield admissible idle/cyclic local solutions;
  secondary2 needs batteries at the origin. Primary proxies yield negative local
  margins. **Local failures/negative margins do not prove global infeasibility.**
- `quick_T6/sf4`: all 12 axis searches completed and replayed. Idle reaches
  7.176–18.290 degrees; cyclic 9.248–22.834 degrees. Each tested cyclic axis extends
  17–40% beyond idle. Six axes cannot determine face planarity, volume, or a full
  tetrahedral boundary.
- `quick_T6/secondary2`: cyclic interval [-24.46848,+8.16398] degrees; idle origin
  negative. Other quick cases retain their completion status in metadata.

T6 means six four-hour periods over the day, using the existing profile sampler.
**This sampler changes daily PV energy with T. T6 and T24 numerical windows must
not be mixed or interpreted as a time-resolution convergence comparison.**
Import-only predictions are separately recomputed for the matching horizon.
Connecting plot segments only guide the eye and do not certify feasibility
between measured points. All limits are local OPF/search results, not certificates.

## Checks and failures

`publish_secondary_study.py` validates successful-solver rows in every completed
case, rejects replay failures, and generates figures/viewer JSON. OpenDSS source
P/Q and battery dispatch differences <1e-3 kW/kvar; all bus voltage differences
<1e-7 pu. Original terminal SOC, SOC limits, lossless energy dynamics, battery power
limits, and inverter Q capability are checked independently of margin slacks.
Negative-margin trial points are network/device-valid but violate relaxed study
voltage or backflow limits: they are not reported as feasible endpoints.

Initial exploratory runs and deliberately interrupted runs remain in the root
and `refined/`. One preliminary T24 candidate failed the strict voltage replay
threshold, so feasibility accuracy was tightened (rather than replay tolerance
relaxed). `quick_T6/secondary3` later stopped on a negative-margin candidate with
1.15e-7 pu replay disagreement. That incomplete run is excluded from publication.
Solver timeouts are retained, not relabeled infeasible. Metadata wall times can
include computer suspension and must not be used as speed benchmarks.

## Reproduce and provenance

Run from repository root with the MSOPF pinned environment. The exact resolved
environment is preserved as `Manifest.snapshot.toml`; MA57 is locally licensed
and not included. Every case has metadata, source hashes/patch, stdout, stderr,
solver summaries, sensitivity inputs (where applicable), baselines and edge CSVs.
`sf_deck_hashes_postrun.json` supplements the launcher (which originally hashed
only the two IEEE8500 decks); these SF hashes are explicitly post-run snapshots.
No deck files were changed during these runs.

```
SECONDARY_T=24 SECONDARY_ROOT=... bash problems/msopf/root_level/run_secondary_8500_screen.sh secondary2 secondary3
SECONDARY_BASE_ONLY=1 SECONDARY_T=24 SECONDARY_ROOT=... bash problems/msopf/root_level/run_secondary_8500_screen.sh sf4
python problems/msopf/root_level/publish_secondary_study.py
```

The generator preserves the five historical viewer datasets and appends only
completed, validated new four-source cases. Sparse-screen faces are not fitted.
`validated_publication.json` summarizes the actual accepted checks. The manuscript
is in the separate PESGM2027_Multi-Source-Multi-Period-OPF repository. GitHub push
is verifiable independently of Overleaf; an Overleaf sync is not asserted here.

Endpoint category labels use proximity thresholds (minimum import <=1 kW, otherwise
customer voltage within 1e-4 pu of a limit); they are not active-set certificates.
A backflow-limited source identity was not archived by the screen. The viewer
therefore reports category only, never substitutes the predicted source identity.

## Figure comparison correction
The paper now uses scripts/render_secondary_comparison.py in its own repository,
which imports the existing vector renderer and palette. The original three_source_T24
and sf_axis_T6 plots are retained as superseded diagnostics. Do not reuse connected
eight-point polylines as measured contours. Legacy 2522 versus the new conversion
changes source siting and modeling assumptions: this is not a transformer-only ablation.
New four-face plots have only six axis samples per feasible case, with counts shown;
face projection suppresses normal offsets and does not prove tetrahedron membership.
