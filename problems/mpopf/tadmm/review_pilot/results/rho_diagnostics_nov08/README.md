# High-rho phase controller and diagnostic profiling

The high-rho start quickly enforces consensus: at iteration 5, owned RMS is
1.18707e-7 and original-model violation is 5.55996e-6. Cost is still 6.63%
above centralized, so this is a feasible early point, not a near-optimal result.
The subsequent task is to lower the penalty while improving cost and retaining
consistency. Restoring the original watchdog produces a historical exit at
**44 iterations / 39.6377 s core**, versus **48 / 51.0920 s** fixed-rho: a 22.4%
lower observed core time and 9.85x lower endpoint violation. Its cost gap is
0.2703%, so this does not establish a speedup at matched 0.1%-cost/original
feasibility. The distinction between historical exit and physical NO is retained.

## Controlled assessments

All runs use the same modified large10k T24 periodic input, July battery weight
8.013987384856764e-8, initial primal starts, free terminal SOC, corrected losses,
proved screening, exact Ipopt locals with MA57, and 24 pinned period workers
plus a pinned coordinator on the Gold 6438Y+ CPU class. The selected stopping
rule is owned RMS <=1e-4 AND (unique dual RMS <=1e-2 OR five-step objective
range/max <1e-3). This remains the selected monolithic historical stop variant;
the July split driver's both-residual exit differs. All-copy residuals and
original feasibility remain separately checked at every iteration.

- Fixed rho=48 with typed diagnostics, job 31652207: historical exit 48,
  core 51.0920 s, total loop 56.9572 s, endpoint violation 8.53850e-3,
  relative cost gap 1.60705e-5. Its complete numerical prefix and initial
  starts match the earlier fixed-rho trajectory.
- High rho0=15000, primal-first phase then residual balancing with decrease
  factor 5, job 31652208: exit 54, core 50.4490 s, loop 57.5365 s,
  violation 2.70224e-3, cost gap about 4.60e-4. This is 3.16 times less
  endpoint violation at similar core time, not a demonstrated timing speedup.
- Decrease factor 2, job 31652227: exit 86, core 78.8587 s, loop 88.5281 s,
  violation 4.89942e-3, cost gap 2.03003e-4. The gentler decrease is slower
  and less feasible than factor 5 in this observation.
- Factor 5 plus the original enabled watchdog, job 31652269: exit 44,
  core 39.6377 s, loop 46.3324 s, violation 8.67149e-4, cost gap 2.70306e-3
  (0.2703%). The numerical prefix matches the no-watchdog factor-five run
  through iteration 40; the watchdog then raises rho from 120 to 300.

The phase controller uses owned/unique RMS, balance factor 5, increase factor
2, update interval 5, rho bounds [1,1e6], and rescales all scaled dual copies
by rho_old/rho_new on each change. High rho0=15000 comes from the paper's T24
large10k appendix. Phase updates are recovered from run_tadmm.jl at release
50b8901c. This restores that mechanism, not the old network, profiles, solver
tolerances or later relaxed residual thresholds.

The watchdog is also recovered from that split source: increment a counter
when owned RMS >2e-4, reset it below 5e-5, retain it in between, and multiply
rho by 2.5 after 15 qualifying high-residual counts. It was absent from the
first phase-only assessments and is explicitly enabled only in job 31652269.
After restoration its effect is independently checked against the unchanged
40-step prefix. No unsupported timing or convergence claim is inferred from
the earlier incomplete controller trials.

None of those four historical exits passes the separate physical NO rule
(0.1% cost AND original maximum violation <=1e-5). All endpoint OpenDSS load
flows converge, with voltage differences 1.48e-10 pu. Voltage agreement alone
does not validate temporal storage consistency. No cone-inequality violation
was present in the first two endpoints; their maxima are original-bound errors.

The factor-five and watchdog runs have local post-first Ipopt median 10 and
maximum 18, versus 12/max 24 for the matched 48-step fixed-rho control. This
supports fewer/shorter local solves in the observed controller outcome, not
a need for the 30-step cap on this large-system case.

## Solve times versus validation times

**Core solve time excludes physical validation and full-dispatch retrieval.**
Native Ipopt solve times are separately logged. The runtime bars in the
comparison figure use core time; total loop times are secondary observations.
Every assessment has three discarded compile steps and a fresh model/state
reset. These are single exploratory observations, not full-warm repeat medians.

On the benchmark Gold host, seven alternating post-compilation evaluations of
the saved endpoint show the independent physical checker taking median
0.92440 s old versus 0.009754 s typed on large (94.8x), and 0.21850 versus
0.002258 s on medium (96.8x). Both implementations match every reported metric
on the saved endpoints and ten deliberately perturbed infeasible cases per
system. No checks or tolerances were removed. The implementation packs fixed
network data into typed numeric arrays once and retains every physical check.

The fixed-48 control still exits at 48 with identical values, while total loop
time falls from the prior 110.4559 s to 56.9572 s. Core remains 50.5359 versus
51.0920 s. The validation improvement is not counted as a solve-time speedup.
New component timers show full-dispatch retrieval 2.9317 s and physical checks
0.5228 s over the 48 iterations, both outside core time. The remaining core
communication/coordination accounting is discussed in the
[period-distribution and execution report](../kamiak_execution_insights/README.md).

## Evidence

- [Convergence and core-time comparison](analysis/rho-diagnostics.png),
  [PDF](analysis/rho-diagnostics.pdf), [summary](analysis/summary.json).
- [Watchdog prefix check](analysis/watchdog-prefix-audit.json) and
  [Gold-host diagnostic microprofile](diagnostic-profile.toml).
- Raw traces, stage iterations/statuses, model/input/source hashes, start
  hashes, hardware/affinity, endpoint OpenDSS logs, profile measurements and
  gate outputs are preserved in the archives. Serialized dispatches and
  checkpoints remain on Kamiak in the respective campaign directories.
- Initial unconstrained gate 31652189 passed diagnostic equivalence but
  failed MA57 initialization with INVALID_OPTION on cn186. Dependent jobs
  31652190/31652191 never ran and were cancelled. The identical frozen source
  passed the complete model/local-solve gate on Gold node cn167 in 31652196.
  Subsequent campaign gates were constrained to the benchmark CPU class.
  These statuses are retained separately; the failed gate is not a rho trial.

Archive SHA256:

- primal-first-31652207-31652208.tgz:
  b331f5468f2156171553d8309a0198de00a4a099d065819b0f99e6af3501c0bc
- gentle-rho-31652227.tgz:
  556cf5224a167050e254a16b9aa02381648090b3c28766af9dd9c3bdbcdc46b9
- diagnostic-gate-history.tgz:
  c963ca417d4b221d9af8a18d841b4eae9dff181d2ade575e5ce6d7193a94487e
- watchdog-rho-31652269.tgz:
  b31dbf00919a9022a49c3316aa3af564a70443920e591c35114650d25ce06426

Frozen numerical sources: c7492d23 (fixed/phase pair), 68c809c3 (factor 2).
The standalone diagnostic optimization was frozen initially in 39bbd6f4.
The watchdog assessment uses 2e772a34. All four assessments completed with
exit 0; their local-status/source/controller/endpoint checks also pass on the
downloaded artifacts. All assessment jobs have ended.
