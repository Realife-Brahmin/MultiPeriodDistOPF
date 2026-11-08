# DDP research transferred to the tADMM assessment

Audit date: November 7, 2026. DDP baseline: origin/master 93738f83 (November 2).
This audit concerns the active SPAR-DDP implementation, not the retired
first-order DDP code. The fixed-root redundancy was already documented by DDP;
the earlier tADMM port missed that existing result.

## Acceptance target and baseline

Use the current main-paper NO rule: absolute relative objective gap <=0.001
AND maximum original-model violation <=1e-6 (ieee123) or <=1e-5 (med2522,
large10k). The user requests the first qualifying iteration **below 100**.
Continue final benchmark runs to strict convergence and keep that count separate.

The repaired, frozen C4 T24 MA57 baseline reaches NO at 66 iterations on
med2522 (job 31645746), versus 389 strict iterations; large10k reaches NO at
168 (job 31645790), versus 267 strict iterations. Both full runs passed
independent OpenDSS validation. The historical NO points were scored from their
traces; only the final dispatches were saved and OpenDSS-checked in C4.
New runs save the actual first NO dispatch and independently check it too.

At large10k iteration 100, original feasibility was 1.49484e-4 while relative
objective error was approximately 1.52e-6. Thus physical feasibility, rather
than objective accuracy, delays NO. One bounded candidate starts rho at 384
instead of 48 (the old controller increased 48 to 96 after iteration 10).
The relative controller and overlapping-window algorithm remain unchanged.
The completed 100-step diagnostics supply no clean benchmark timing or
full-convergence claim. The first three candidates all missed NO: initial
384 with the existing relative controller returned to 96 and ended at
1.396e-4 original violation; fixed 384 ended at 3.977e-4; fixed 48 improved
to 2.327e-5. All stage statuses and source checks passed. No NO dispatch was
claimed for these failures. See `results/ddp_transfer_NO` for the complete
preserved traces, including the failed gate-harness initialization.

The residual-balance ratio is now an explicit, recorded hyperparameter;
its default remains 10. A separate bounded trial uses ratio 2, starting at
48, to test earlier adaptation. This changes the penalty-controller setting,
not the consensus averaging or scaled-dual update. Success still requires
the unchanged original-feasibility and objective thresholds and an actual
saved NO dispatch passing OpenDSS.

## Exact, transferable changes

- **Redundant bounds:** `spar_ddp/notes/CONSTRAINT_SCREENING.md` and
  `spar_ddp/examples/power_system/voltage_screening.jl` prove that ell >= 0
  follows from P^2+Q^2 <= v*ell with positive voltage. tADMM now removes those
  lower bounds, retains every nonroot positive lower-voltage bound, and retains
  the fixed positive root equality. Root voltage bounds were already removed
  in C4. All original bounds remain in the independent physical checker.
- **Data-only voltage screening:** reuse the shared helper directly, with the
  same across-horizon union as DDP. No screen comes from a previous OPF optimum.
  The October 20 factorial study
  `spar_ddp/experiments/screening_ablation/results_oct20/FINDINGS.md` found that
  ell removal consistently helped, but voltage-upper removal slowed med2522.
  Therefore default flags are root+ell on med2522 and root+upper+ell on
  ieee123/large10k. Import-bound screening is available only with its positive
  data-derived proof; its mixed negligible speed benefit does not justify
  making it default. No lower-voltage screening is enabled.
- **Initialization and inputs:** shared `loadflow_start.jl`,
  `control_layout.jl`, periodic native exported profiles, system C_B and
  incoming r*ell/x*ell losses already match the corrected DDP formulation.
  Idle-battery, loss-aware power-flow starts are retained. The user explicitly
  overrides DDP's soft terminal penalty: terminal penalty and equality are
  absent here. Existing soft-terminal DDP timings are consequently not a
  matched final comparison.
- **Reuse without changing the method:** local JuMP models already persist
  between ADMM steps and reuse previous primal starts. Physical diagnostics
  and dispatch extraction now cache topology, control indices and membership
  maps instead of rebuilding them every iteration. Regression probes compare
  every metric against frozen C4, including deliberately violated original
  current, import, voltage, battery-power and energy bounds.
- **Accurate local solves:** retain exact sparse polynomial derivatives,
  unrelaxed bounds and explicit absolute complementarity tolerance 1e-8.
  `PLATEAU_AUDIT.md` records the controlled multiplier/termination diagnosis.
  Keep MA57, one numerical thread per actor, and the existing ordering/scaling.

## DDP optimizations that do not transfer as Ipopt settings

The notes index, CHANGES, experiment index and the following research were
checked: FILTERDDP_HESSIAN_EXACT_REWRITES, FILTERDDP_CACHED_CONSTRAINT_JACOBIAN,
FILTERDDP_TYPED_CONSTRAINT_VECTOR, FILTERDDP_ACTIVE_B_ROWS,
KKT_ORDERING_AND_MA57, KKT_MATRIX_UPDATES,
IPOPT_SENSITIVITY_REUSE_AND_MAP_COMPRESSION, the September 14-15 findings,
October tree ablation, battery-count, value-grouping and compact-coupling studies.

- Tree elimination, battery Schur complements, wide sensitivity RHS solves,
  cached KKT sparsity, factored policies, value recursion, parallel scans and
  affine DDP line searches operate inside SPAR-DDP. tADMM uses Ipopt for each
  local NLP; claiming these as enabled tADMM optimizations would be incorrect.
  JuMP/Ipopt supplies sparse derivatives and its own factorization workflow.
  Public Ipopt.jl does not expose reusable final factors for the proposed DDP
  sensitivity path. No stale numeric factors or ABI-dependent bridge is used.
- DDP barrier schedules and inner stopping rules are part of its outer method.
  Copying them into tADMM could destroy the accurate proximal solves required
  by ADMM. The earlier inaccurate local termination is evidence against that.
- Battery-count changes alter the physical problem. Value grouping introduced
  material accuracy/iteration costs; compact coupling failed its cost gate.
  Neither is adopted to manufacture a lower tADMM iteration count. Root-device
  population and all frozen exported input arrays remain unchanged.

## Validation and evidence boundaries

`validate_ddp_transfer.jl` checks the shared screen, original-metric equivalence
against both completed C4 systems, preserved positive voltage, removed current
bounds, and an actual MA57 local solve on each system. Its success is required
before the 10k calibration can run. `validate_no_calibration.py` checks every
stage status, complete iteration/period coverage, source hashes, the exact NO
crossing and OpenDSS at the saved NO point. A failed target is recorded as such.

The calibration changes exact constraint screening and initial rho together;
it can assess the resulting trajectory but cannot isolate their timing effects.
Retain immutable source/input hashes, all raw traces and earlier failed runs.
Only a subsequently warmed, fully converged and independently validated run
can supply a new benchmark timing row. No MA97/MUMPS sweep is authorized.
