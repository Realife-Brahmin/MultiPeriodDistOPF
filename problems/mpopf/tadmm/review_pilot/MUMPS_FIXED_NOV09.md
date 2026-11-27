# Fresh matched MUMPS comparison, November 9

User requested fresh MUMPS tests before pursuing the proposed optimized battery
schedule. That schedule is not implemented or used in this experiment.

Frozen numerical source: 4a074c167403ba9decb31a578a8df4de1e47f0c0, copied
unchanged from the completed MA57 targeted campaign. Input and execution-source
hashes are checked before launch. Remote campaign:
`/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_mumps_fixed`.

Compare both methods using MUMPS, with the same corrected model, July C_B,
periodic profiles, idle-battery/load-flow initialization, free terminal SOC,
exact local solves, primal-dual reuse, batched updates, single-thread policy,
Gold hardware, full warmup/reset, stopping criteria and validation as the
completed fixed-48 MA57 experiments. One worker per period. The rho schedule
starts at 120000 for T48 and 240000 for T96, divides by five every five
iterations to a floor of 48, then stays fixed. No solver-specific rho tuning.

Native repeated-solve gate: 31660025. T48 run: 31660026; energy audit: 31660027.
T96 run: 31660028; energy audit: 31660029. T96 depends on successful T48
completion including its collection/OpenDSS validator. No completed numerical
result is claimed by this submission record. Preserve unsuccessful outcomes.

Report core and loop times, initialization/build and warmup separately; all
diagnostics and OpenDSS checks remain outside core. The RMS/dual endpoint and
original-model physical feasibility are separate flags. Loose milestone times
must retain integrated battery-energy errors. These are exploratory single
trajectories, not repeated median paper benchmarks.

The older legacy-controller MUMPS T48 test took 635.6532 s tADMM core versus
179.6677 s centralized. It did not establish a MUMPS advantage. This fresh
test checks whether the improved rho schedule changes that outcome.

Reproduction launcher: `setup-mumps-fixed-nov09.sh` at the repository root.
