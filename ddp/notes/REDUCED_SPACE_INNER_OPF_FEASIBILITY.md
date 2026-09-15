# Reduced-space MPOPF: is the network eliminable by a single-period inner OPF?

Phase-1 structural diagnostic for the proposed decomposition

```
outer : B^{t-1}, P_B^t, battery dynamics, battery bounds
inner : given fixed P_B^t, solve every algebraic network quantity
        (P_Subs, Q_Subs, branch P/Q, v, ell, DER reactive, SOCP slacks)
```

**Scope: IEEE123 (`ieee123C_1ph`), `T = 24` (hourly), one system only.** IEEE2522
and large10k are deliberately not run yet; cost estimates are at the end.
Nothing here modifies the production FilterDDP solver, and no full-horizon
FilterDDP run was launched.

Code (all opt-in, nothing imported by the production path):

| file | role |
|---|---|
| `ddp/examples/power_system/inner_network_opf.jl` | reusable inner-OPF driver + feasibility-restoration diagnostic |
| `ddp/examples/power_system/run_inner_opf_survey.jl` | dispatch-pattern feasibility survey |
| `ddp/examples/power_system/probe_reduced_value_function.jl` | `Phi_t` probe |

Raw results in `ddp/results/reduced_space/`.

## Modelling decision, stated rather than hidden

The primary inner problem is **network-only**. The energy-slack row
`B^{t-1} - dt*P_B^t - Bmin - s_E = 0`, `s_E in [0, Bmax-Bmin]`, is excluded from
it. This is not a relaxation of the network problem: with both `B^{t-1}` and
`P_B^t` fixed that row has no optimisation freedom at all, it merely evaluates
whether the dispatch respects the state-of-charge box, which the decomposition
assigns to the outer layer.

This mattered a great deal. Of the 111 surveyed dispatches, only **5 satisfy the
SOC box** at the entering energy taken from the centralized trajectory, while
**all 111 are network-feasible**. Folding the two together would have reported
106 spurious "infeasibilities" that are battery-box violations, not network
violations, and would have answered the wrong question. The SOC check is
recorded per row as `soc_box_ok`, and the energy row is available in the driver
via `include_energy_row=true`.

## Validation of the inner model

Two independent checks, both passed:

1. **Against the centralized reference.** Fixing `P_B^t` to the centralized
   optimum reproduces the centralized single-period network solution: `P_Subs`
   agrees to `5.4e-09`, `6.6e-10`, `7.8e-09` pu at the high/medium/low hours,
   and all bus voltages to `<= 1.2e-07` pu.
2. **Against FilterDDP's own constraint callback.** For every one of the 111
   surveyed solves, the recovered control vector was fed to FilterDDP's analytic
   stage-constraint function (a different code path from the JuMP model).
   Worst residual over the whole survey: **`7.2e-12`**.

A third, independent check against a captured FilterDDP stage
(`ieee123C_1ph_T3_stage1_oracle_capture.jls`, `T=3`, stage 1, `mu = 1e-8`): that
iterate is essentially network-feasible (`max|c| = 6.2e-08`). Fixing `P_B` to
FilterDDP's own captured dispatch and solving the inner OPF reproduces
FilterDDP's network solution to `5.8e-09` (`P_Subs`), `6.1e-07` (voltages) and
`8.5e-09` (branch real power).

## Feasibility survey

111 inner solves: 3 demand levels x (13 structured patterns + 24 reproducible
random dispatches). Demand levels chosen by **net** load (load minus PV):
low `t=13` (0.6291 pu), medium `t=22` (0.9335 pu), high `t=6` (1.1092 pu).
Battery groups are by electrical depth from the substation (IEEE123 is a single
feeder: the substation has one child, so depth, not feeder identity, is the
meaningful grouping axis; batteries span depth 4-24, median 13).

Patterns: `zero`, `centralized_optimal`, `all_max_charge`, `all_max_discharge`,
`deep/shallow_max_charge/discharge`, two `opposing_*` group directions,
`struct_max_voltage_drop` / `struct_max_voltage_rise` (deepest quartile only),
`struct_alternating_circulation`, and `random` seeds 1001-1024 drawn uniformly
throughout the battery-power box.

**Result: 111 / 111 feasible. Zero infeasible cases.** The restoration
diagnostic was therefore never required (it is implemented and ready in
`restore_feasibility`, but reporting it here would be fabricating a result).

Observed envelopes across all 111:

| quantity | range |
|---|---|
| `P_Subs` | 0.12448 .. 1.66462 pu |
| bus voltage | 1.00549 .. 1.05000 pu (limits are 0.95 / 1.05) |
| real loss | 0.00200 .. 0.04873 pu |
| IPOPT iterations | 29 .. 36 |
| solve time | 0.037 .. 1.390 s (mean 0.066 s) |

**Voltage limits never bind.** Minimum voltage over every dispatch tested is
1.0055 pu against a 0.95 floor; the only "active" upper-voltage bus is the
substation itself, whose 1.05 pu is a fixed equality, not a limit. The IEEE123
transcription also carries **no branch ampacity constraint at all** (`ell` has a
lower bound of 0 and no upper bound in `build_model`), so "branch-current
utilisation" has no rating to be measured against; `ell_max` / `imag_max` are
recorded as raw magnitudes instead, and the absence of a thermal limit is itself
worth flagging for the formulation.

**The one binding resource is DER reactive support.** `|q_norm|` reaches exactly
1.0 in every single solve, with 33-51 of the 51 DERs pinned at their reactive
limit. Reactive capability is fully consumed because it reduces losses and hence
`P_Subs`, which is the only priced quantity.

## Reverse export

`P_Subs >= 0` with no upper bound, so excessive discharge could in principle be
infeasible. Scanned all 24 hours at `all_max_discharge`:

**It cannot bind on this system.** Closest approach is `P_Subs = 0.12448 pu`
(124 kW) at `t=13`. The reason is structural, not numerical: total battery power
is 0.5066 pu against a minimum net load of 0.6291 pu, so full simultaneous
discharge still leaves ~0.12 pu to import, before losses (which push `P_Subs`
further up, never down). Reverse export would require a battery fleet ~24%
larger, or a lighter minimum net load, than this system has.

## The reduced stage-value function

`Phi_t(P_B^t) = min_y { c^t * S_base * dt * P_Subs^t : network constraints }`.
Probed at two bases (`zero`, `centralized_optimal`) x 3 hours x 6 directions
(3 single-battery: shallowest / median-depth / deepest; `aggregate_uniform`;
`group_deep`; `group_shallow`) x 5 perturbation sizes `h in {1e-4 .. 1e-2}` pu,
central differences. 180 probes, 137 usable; the other 43 are correctly detected
as leaving the battery box (at `t=6` the centralized optimum is *at* max
discharge, so every direction exits).

- **Smooth.** First derivatives are stable across two decades of `h`: relative
  spread `8.5e-10` to `7.3e-06`.
- **Locally convex along every direction tested.** Second differences are
  positive in all 137 usable probes, from 4.72 (single shallow battery) to
  398.8 (aggregate), and are themselves stable in `h` (spread `<= 5e-05`).
- **Physically sensible gradients.** `dPhi/dP_B` at `t=13` is `-132.52` USD/pu
  for the shallowest battery and `-136.00` for the deepest, against a raw energy
  price of `131.8` USD/pu: the marginal value of discharge is the price *plus*
  the avoided losses, and it is larger for electrically distant batteries. This
  is ordinary LMP structure.
- **Active-set changes do occur** (at `t=22`, the aggregate and group directions
  change which DER reactive limits are active) **without destroying local
  smoothness** at the probed scale.

### The gradient is analytic, not a finite-difference artefact

The dual of the real-power balance row at a battery bus **is** `dPhi_t/dP_B`.
Checked against central differences over all 137 usable probes: median relative
error `6.0e-09`, worst `7.3e-06`. The outer layer therefore never has to
finite-difference the inner solve -- one inner solve returns its own exact
gradient.

## Answers to the eight questions

1. **Can network variables be eliminated exactly through a single-period inner
   OPF?** On this evidence, yes. Given `P_B^t`, the inner solve reproduces the
   centralized and the FilterDDP network solutions to `1e-9`-`1e-7`, and the
   network block carries no inter-period coupling of its own.
2. **Is the battery-power box empirically fully recourse-feasible?** For IEEE123
   at `T=24`: yes, 111/111 including all box corners and 24 random interior
   points. **This is empirical coverage, not certification** -- 24 random draws
   plus 13 structured directions do not prove feasibility of an uncountable box.
   A certificate would need either a constructive argument or a global feasibility
   test, neither of which was attempted.
3. **If not, which combinations and constraints define its boundary?** No
   boundary was reached. The nearest thing to a limiting resource is DER reactive
   capability, which saturates everywhere but never causes infeasibility.
4. **Is `P_Subs` effectively determined by fixed battery dispatch and network
   losses?** Yes. `P_Subs = net_load - sum(P_B) + loss` holds identically, and
   the loss term varies only 0.002-0.049 pu across the entire survey. `P_Subs`
   is a dependent quantity, not an independent outer decision.
5. **Which genuine network actuators must remain as outer decisions?** None were
   identified for IEEE123. The only network actuator is DER reactive power, and
   it is a pure within-period quantity with no inter-period state -- the inner
   solve sets it optimally. Nothing in the network block couples periods.
6. **Does `Phi_t` appear smooth and convex enough for a battery-only DDP
   method?** Yes, locally: smooth, stable under `h` refinement, positive
   curvature along all 6 directions at all 3 hours. Caveat: convexity is
   established *along probed directions at probed points*, not globally, and
   active-set changes are present.
7. **What must the inner solve return to reproduce FilterDDP's backward
   recursion?** At minimum `Phi_t` and `dPhi_t/dP_B` (the battery-bus balance
   duals, free from the solve). A second-order outer method also needs
   `d2Phi_t/dP_B2`; that was obtained here only by finite differences, and
   getting it analytically would require a sensitivity solve against the inner
   KKT system -- exactly the cost structure the reduced space is meant to avoid,
   and the open question for phase 2.
8. **Exact reformulation or approximation?** Exact, *conditional on* the inner
   problem attaining its optimum and on the projected feasible set containing
   the battery box. Both hold empirically here. The reduction is a genuine
   variable elimination, not a relaxation: no constraint was dropped, the
   eliminated variables are recovered exactly, and the duals match.

## What this does not establish

- One system, one horizon, one profile. `T=24` hourly only.
- Empirical coverage is not a feasibility certificate (question 2).
- Convexity is local and directional, not global.
- The second-order information a DDP outer loop wants is not yet available
  analytically.
- No claim is made about whether the resulting method would be *faster*. This
  phase is a structural diagnostic; it says the decomposition is well-posed on
  IEEE123, not that it pays.

## Cost estimate for the larger systems

Measured here: 111 survey solves in 7.3 s wall (mean 0.066 s), 180 probe solves,
24 reverse-export solves -- the whole IEEE123 phase is well under a minute of
solver time and peaks at ~2 MiB RSS growth per solve.

Scaling by the measured per-period centralized cost ratio (ieee2522 ~ 20x
ieee123, large10k ~ 180x ieee123 per period):

| system | ~315 solves | note |
|---|---|---|
| ieee2522C_1ph | **~10-15 min** | recommended next, cheap |
| large10kC_1ph | **~1.5-2 h** | run only after ieee2522 confirms the workflow |

Both are feasible on an idle machine. The instruction's gate is respected: this
report stops before either.
