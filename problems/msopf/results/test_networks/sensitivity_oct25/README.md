# Customer-secondary four-source preparation, 2026-10-25

Continuation of HANDOFF_2026-10-24.md. Both local sensitivity preparations completed successfully. These outputs are local affine predictions; the full measured edge campaign has not been run.

## Model and verification

T = 24, periodic load/price/PV profiles, existing 40% DER fleets, PV reactive dispatch enabled. Batteries are idle for the central-difference solves, with each source angle fixed. No-backflow is released while measuring derivatives; the 0.95-1.05 pu customer band is enforced. Receiving buses of `svc_*` lines define customers; primary buses have no band. The hard cyclic battery rule remains the day-model convention.

Each case uses one equal-angle solve and +/-0.1 degree for each of three moving sources (subs1 reference). Every successful dispatch is replayed for all 24 hours in OpenDSS. Both runs exited 0. The independent Python validator checked all seven replays, all 96 hourly sensitivity rows, and all 156 predicted reaches per case (78 directions x two estimates). In total: 14 day solves and 336 hourly OpenDSS comparisons.

Idle estimates use the worst hour; cyclic-battery estimates use daily-mean imports and derivatives. The latter four import planes form a tetrahedron in this approximation; the plotter checks all four plane intersections against every plane. This geometric fact does not establish the true OPF region: voltage limits, reactive redispatch, nonlinear response and limited battery capacity may move its edges.

## Completed runs

### four_ieee9500xf

Elapsed process time: 272 s, including Julia startup, assembly, all solves and replays. Worst replay discrepancies: P 1.302e-05 kW; Q 1.089e-05 kvar; V 5.312e-10 pu; battery power 0.000e+00 kW.

- [Validation](four_ieee9500xf/validation.json), [raw log](four_ieee9500xf/run.log), [stderr](four_ieee9500xf/run.err), [metadata/hashes](four_ieee9500xf/metadata.txt).

- [Prediction figure](four_ieee9500xf/predicted_region.png), [editable SVG](four_ieee9500xf/predicted_region.svg), [plane-intersection vertices](four_ieee9500xf/predicted_daily_vertices.csv).

### four_smartds_p1u_typeiii

Elapsed process time: 471 s, including Julia startup, assembly, all solves and replays. Worst replay discrepancies: P 5.182e-05 kW; Q 5.352e-05 kvar; V 1.345e-09 pu; battery power 0.000e+00 kW.

- [Validation](four_smartds_p1u_typeiii/validation.json), [raw log](four_smartds_p1u_typeiii/run.log), [stderr](four_smartds_p1u_typeiii/run.err), [metadata/hashes](four_smartds_p1u_typeiii/metadata.txt).

- [Prediction figure](four_smartds_p1u_typeiii/predicted_region.png), [editable SVG](four_smartds_p1u_typeiii/predicted_region.svg), [plane-intersection vertices](four_smartds_p1u_typeiii/predicted_daily_vertices.csv).

## Provenance and reproduction

These are preparation runs on etrl309-aryan, Julia 1.12.7, Ipopt with its default linear solver, lean JuMP mode and one Julia/BLAS/OpenMP thread. They are not matched SPAR-DDP/tADMM benchmark results. Every solve log contains status, iterations, build time and optimizer time. Metadata records wall time, exit code, source/deck hashes and environment. The launcher caps each solve at 600 seconds.

The first run began at ef081491 plus the saved source.patch; its numerical driver was then committed as 8732ec28. The second run began at 8732ec28. Both metadata files record the same numerical-driver hash. Manifest.snapshot.toml preserves the actual resolved Julia environment, whose hash is recorded in each metadata file; use it as the local problems/msopf/Manifest.toml when reproducing. OpenDSSDirect 0.9.9 prints its known Julia 1.12 precompilation warnings; raw stderr is retained.

Run root_level/run_test_edition_sensitivities.sh from the repository root to reproduce into a fresh MSOPF_SENS_OUT. Use check_test_edition_sensitivities.py and plot_test_edition_sensitivities.py on each completed run directory to validate and render without solving again. See the MSOPF README for the full commands.

customer_band.log records the passing customer-only-band regression, its two OpenDSS replays, legacy-band propagation, both P1U switching configurations, and identical IEEE9500 primary line definitions. No new raw deck was written.

## Remaining work

- Full measured 78-direction idle/cyclic edge searches on Kamiak, after renewed login and a shared-queue check; use MULTI_VMARGIN=20000.

- P1U Type II is wired and its source/switching configuration is regression-checked, but has not had a new sensitivity campaign here.

- K1, Ckt24 and 8500xf retain the handoff smoke evidence; no extra substation sites were chosen for them.

- Numerical step-size robustness and measured voltage cuts remain to be assessed before using the prediction as a paper result.
