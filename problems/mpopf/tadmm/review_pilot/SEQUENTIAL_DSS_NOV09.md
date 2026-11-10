# Restoring chronological MPOPF validation

The former shared `opendss_check.jl` disables Storage and substitutes constant
power injections in snapshot mode. Its voltage/network checks are useful but
do not validate sequential battery execution. Prior “OpenDSS passed” labels
must be read as network-snapshot agreement only; original logs remain intact.

Historical sources actually inspected:

- MultiPeriodDistOPF92a6cf804b829f451394ea2215b06d5792aa6212 (2024-11-06):
  direct Julia storage commands began obeying dispatch under `DispMode=External`.
- a974207318ab23e363d36c00cd5de4ea0e97cf57 (2024-11-08): Follow -> External.
- 28ed0c233841872455237687a2872054d8596a63, `src/openDSSValidator.jl`:
  sequential live-storage solves, terminal SOC, bus voltages, source P/Q and
  line losses versus OPF, with exported discrepancy variables.
- MultiPeriod-DistOPF-Benchmark94298878310bbba9440a431630631093f6aed579,
  `functions/OneTimeScripts/validationAgainstOpenDSS_MATLAB/validateOPF_m.m`:
  daily mode, live storage/controllers/monitors and network discrepancy arrays.
  It is provenance, not a command template for the current engine.

Use the installed DSS C-API through OpenDSSDirect.jl. Explicit external kW
commands, SolveSnap then exactly one FinishTimeStep per interval, initialized
energy only once, and asserted clock/energy continuity. Never reset Storage
energy to OPF consensus copies between time steps. Do not substitute loads for
batteries. Compare requested/actual dispatch and before/after stored energy.
Detect native reserve/capacity clamping through energy-balance error even if
the stored energy remains numerically bounded. Read actual terminal powers,
not just the Storage kW setting.

Native unit cases cover charge/discharge signs and carried energy, reserve
clamp plus subsequent refusal, capacity clamp versus operating95% limit,
95% efficiencies, and quarter-hour increments. A full centralized replay is
the integration control; capped endpoints are tested without reoptimization.
CSV tables cover every battery, bus and branch at every time step, plus period
totals. Source/input/solution/table hashes and all failed checks are preserved.

Two explicitly different battery assumptions: `model` sets charge/discharge
efficiency100% to match the current lossless optimizer; `dss` retains the DSS
file95% values. The latter is an efficiency-mismatch stress test, not the same
optimization model. Both retain the documented network matching assumptions
(source1.05pu, line capacitance0, constant-PQ voltage range and PV injections).
No inverter efficiency curve or nonzero idle loss is silently ignored.

Sequential schedule pass requires network limits/agreement, actual dispatch
tracking, no energy-accounting clamp error, requested and actual SOC within
operating bounds and agreement with OPF energy states, at declared tolerances.
Power/energy differences are normalized by the instance base with1e-5 target;
voltage agreement1e-4pu. Raw discrepancy tables permit stricter re-scoring.
This validator is passive and outside optimization solve-time benchmarks.

Remote campaign: `/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss`.
Reproduce with repository-root `setup-sequential-dss-nov09.sh`.
