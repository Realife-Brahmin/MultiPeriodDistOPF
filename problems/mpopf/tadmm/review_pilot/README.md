# Corrected tADMM reviewer-response pilot — 309, November 7, 2026

This is a fresh experimental harness. Historical July drivers and result rows
are preserved. It implements the same overlapping temporal consensus ADMM:
each period owns one branch-flow OPF, local battery-energy copies span that
period and its immediate neighbors, their consensus is averaged, and scaled
duals are updated for every copy.

## Agreed problem

- Use the current shared test decks and SPAR-DDP exporter with periodic profiles.
- Use system-specific C_B from `terminal_soc_penalty.jl`: 0.001 IEEE123,
  0.0019 med2522, 0.000024 large10k. These weights are fixed across horizons.
- **No terminal-SOC penalty and no terminal equality** (explicit user decision).
  Initial SOC, SOC bounds and lossless battery dynamics remain enforced.
- Include incoming-line real/reactive losses r*ell and x*ell, voltage-drop loss
  terms, nonnegative current square, and the second-order-cone inequality.
- Preserve the current SPAR-DDP benchmark's root-device convention: root balances
  contain only source import and outgoing flows. Root-connected devices are
  excluded from the OpenDSS reinjection check. The med2522 root resource is
  therefore not an active network injection. This is an explicit matched-model
  convention, not a claim that all raw-deck resources are represented.
- Preserve the shared benchmark's absence of branch ampacity upper bounds.
  Do not describe a successful check as demonstrating thermal-limit compliance.

Existing soft-terminal SPAR-DDP timings are NOT matched comparators for this
free-terminal experiment. The new centralized solve here is the comparator.
Corrected large10k is also a different physical instance from July's deck.

## Verified July defects and differences

`260731-ias-trans-tadmm:run_tadmm.jl` omits r*ell and x*ell from the nodal
balances. Current master includes those terms, but its stopping residual still
checks only diagonal SOC copies. This harness checks all copies. Its dual
residual repeats consensus changes with the correct copy multiplicity.

The old phase-specific rho heuristics are replaced by standard residual balancing
every ten iterations, with scaled-dual rescaling when rho changes. rho starts at
3000*sqrt(T/24), before any tuning. This is an algorithm-implementation change
and must be disclosed in a July-to-current comparison.

Battery-power variables exist only where local dynamics need them; unconstrained
powers for distant periods are removed. Dynamics at t and t+1 constrain the
three-copy local window. At consensus, the network-owned battery power at every
period obeys the original full-horizon energy dynamics. T=1 and T=2 windows are
well-defined. No terminal condition is introduced by an end window.

## Environment and solver layers

tADMM does not call FilterDDP/SPARDDP. It uses JuMP and Ipopt. The small
`Project.toml` and pinned `Manifest.toml` avoid reinstantiating the historical
plotting/interactive environment or requiring Gurobi. The DDP environment and
defaults are unchanged. This pilot uses Julia 1.12.6 on ETRL309-ARYAN.

MUMPS is the initial Ipopt linear-system backend. `TADMM_LINEAR_SOLVER=ma57`
or `ma97` and `TADMM_IPOPT_OPTIONS` select alternatives for BOTH the centralized
reference and local OPFs. The local HSL builds are under
`C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build`. Include
`hsllib=<DLL path>;linear_system_scaling=none` when selecting them.
Library existence alone is not evidence that a backend loaded successfully.

A native conic backend is a separate future comparison. This harness currently
passes the cone as a smooth quadratic inequality to Ipopt, a general NLP solver.

## Execution and checks

From the repository root in PowerShell:

```powershell
& ./problems/mpopf/tadmm/review_pilot/run_case.ps1 -SystemName ieee123C_1ph -Horizon 6 -RunLabel pilot
```

The machine-local wrapper uses 309's existing Julia executable and depot. Each
new label gets a fresh result directory; existing evidence is never overwritten.
Source snapshots and deck hashes are frozen before execution.

Centralized and tADMM each receive a full warm-up, then restart from their
original battery-idle power-flow starts for the measured solve. Local OPFs reuse
their previous primal iterates within ADMM. No centralized optimum initializes
ADMM. Ipopt must actually terminate successfully; iteration-limit subproblems
are rejected. Warm-up/timed ADMM numerical traces must agree.

Default convergence requires maximum all-copy consensus error <=2.5e-7,
maximum dual residual <=1e-4, and original-model feasibility <=5e-7. The latter
also integrates dispatched battery powers from initial SOC and checks resulting
energy bounds, so small consensus errors cannot hide accumulated SOC drift.
An iteration cap is a failure, not convergence. Full traces allow later scoring
at looser objective/feasibility criteria without stopping the solve early.

`validate.py` requires objective agreement with the centralized reference within
1e-5 relative, independent original-model feasibility, matching warm-up/timed
trajectories, completed OpenDSS solves, no voltage violations at the shared
checker's 1e-4 p.u. limit tolerance, maximum voltage difference <=1e-4 p.u., and
OpenDSS/optimizer objective agreement within 1e-4 relative.

The shared OpenDSS checker reads raw DSS elements independently, reinjects only
the optimized device dispatch, and uses the same 1.05 p.u. source and constant
power device assumptions. Its `OPENDSS_ZERO_LINE_C=1` removes DSS default line
capacitance, because the branch-flow model has no shunt capacitance. This is a
check of the stated benchmark physics, not every native DSS modeling feature.
SOCP equality gaps and fictitious losses are recorded separately; convergence
of the relaxed problem alone is never called AC feasibility.

## Timing interpretation

Everything runs sequentially on 309: one Julia thread, one BLAS thread, one
OpenMP thread. No period-worker parallel simulation is performed.

- `central_solver_s`: Ipopt-reported centralized solve time, after warm-up.
- `central_wall_s`: wall time around its optimize call.
- `tadmm_algorithm_s`: measured sequential stage work plus coordination,
  excluding model builds, warm-up, physical diagnostic evaluation and CSV output.
- `tadmm_wall_with_diagnostics_s`: actual measured ADMM loop elapsed time,
  including per-iteration diagnostics/output, excluding build/warm-up/OpenDSS.
- `ideal_parallel_estimate_s`: sum of each iteration's maximum stage elapsed
  time plus measured coordination. This assumes one worker per period and zero
  communications/contention overhead. It is NOT measured parallel time.

Setup, builds, warm-up and validation are reported separately. Process samples
record working set, private bytes and CPU time every 0.5 s. The peak includes
setup/warm-up and resident sequential models; it is not parallel aggregate RAM.
These are exploratory single-repeat pilots, not publication-ready speedup rows.
