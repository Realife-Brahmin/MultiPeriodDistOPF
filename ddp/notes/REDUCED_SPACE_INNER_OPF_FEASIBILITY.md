# Reduced-space MPOPF: is the network eliminable by a single-period inner OPF?

Structural diagnostic for the proposed decomposition

```
outer : B^{t-1}, P_B^t, battery dynamics, battery bounds
inner : given fixed P_B^t, solve every algebraic network quantity
        (P_Subs, Q_Subs, branch P/Q, v, ell, DER reactive, SOCP slacks)
```

**Scope: `ieee123C_1ph` and `ieee2522C_1ph`, `T = 24` (hourly). large10k not yet
run.** Nothing here modifies the production FilterDDP solver, and no
full-horizon FilterDDP run was launched.

## Headline

**The elimination is exact on both systems. The claim that the battery-power box
is fully recourse-feasible is FALSE — it holds on IEEE123 and fails on
IEEE2522.** IEEE123 alone could not have detected this: that network is so
lightly constrained that no limit can activate (minimum voltage 1.0055 pu
against a 0.95 floor, no ampacity constraints in the model at all, reverse
export structurally impossible). IEEE2522 is genuinely voltage-constrained and
rejects 8 of the same 111 dispatches, all by undervoltage under charging.

Consequence for the proposal: an outer layer holding **only** battery energy,
dispatch, dynamics and bounds is **not sufficient**. It would propose dispatches
the network cannot serve — including one with essentially zero *aggregate*
battery power.

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

## IEEE123 feasibility survey: the box IS fully recourse-feasible here

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

## IEEE2522: the box is NOT fully recourse-feasible

Same 111-solve protocol, same three demand levels (`t=13/22/6` by net load).
250 batteries, 1.3301 pu total battery power, mean solve 1.78 s, 47 MiB peak
per solve, whole survey 197 s. Validation residual against FilterDDP's own
callback: worst `6.7e-11` over the 103 feasible solves.

**8 of 111 infeasible. Every one diagnosed as undervoltage** by the restoration
model; the reverse-export slack sits at zero in all 8.

| pattern | low | medium | high |
|---|---|---|---|
| `all_max_charge` | INFEASIBLE | INFEASIBLE | INFEASIBLE |
| `deep_max_charge` | INFEASIBLE | 0.9500 | INFEASIBLE |
| `shallow_max_charge` | 0.9635 | 0.9500 | INFEASIBLE |
| `struct_max_voltage_drop` | 0.9674 | 0.9558 | INFEASIBLE |
| `opposing_deep_chg_shallow_dis` | 0.9608 | 0.9512 | INFEASIBLE |
| `all_max_discharge` | 1.0301 | 1.0097 | 0.9952 |
| `zero` | 0.9908 | 0.9718 | 0.9557 |
| `centralized_optimal` | 0.9908 | 0.9949 | 0.9918 |

(cells are `vmin` in pu where feasible; floor is 0.95)

Three things this establishes:

1. **The boundary is undervoltage under charging**, and it tightens with demand:
   one pattern fails at low demand, one at medium, five at high. Restoration
   depths reach `7.4e-2` in `v` units at `t=6` `all_max_charge`, i.e. about
   0.91 pu voltage against a 0.95 floor.
2. **Feasibility is not a function of aggregate battery power.**
   `opposing_deep_chg_shallow_dis` has `sum(P_B) = -0.0170 pu` -- essentially
   zero net battery power -- and is still infeasible at high demand. The
   projected set `F_t` is genuinely multidimensional; it cannot be summarised as
   a tightened bound on total power.
3. **All 72 random interior draws were feasible** at all three hours. The
   infeasible region is confined to the charging corners and faces, not
   scattered through the interior.

`centralized_optimal` stays comfortably inside (0.9908-0.9949), so the
optimizer never approaches the boundary -- which is exactly why a study that
only ever looked at optimal trajectories would never have found this.

## A defect in the existing transcription, found in passing

The root power-balance rows in `ieee123c_filterddp.jl` (lines 118 and 130) are
`ps - sum(P_out)` and `qs - sum(Q_out)`: **no `pb`, `pD`, `pL` or `qD` term**,
and the per-bus loops that follow cover only `Nm1set`. Any resource sitting *on
the substation bus* is therefore invisible to the power balance.

- `ieee123C_1ph`, `large10kC_1ph`: nothing at the root bus. Unaffected.
- `ieee2522C_1ph`: has **both** a battery (bus 1, 0.00426 pu = 4.26 kW) **and**
  a DER/PV (`S_D_R` = 0.00511 pu = 5.1 kVA) at the root. Both are inert: their
  variables exist, carry bounds, appear in the energy-slack row and in the
  `C_B*P_B^2` cost, but have no physical effect.

**This has not corrupted anything published.** The centralized tADMM reference
shares the convention: the inner OPF reproduces the centralized ieee2522
solution to `2.3e-08`, and the centralized optimizer leaves that battery at
exactly `0.000000` at all 24 hours -- which is what an inert variable carrying a
positive quadratic cost does. Both models agree; the gap is representational.
At 0.320% of each fleet it is numerically negligible.

It did, however, break two things in this study's own probe code, both now
fixed: `lambda_bal` returned `NaN` for a battery with no balance row (correct
value is exactly `0.0` -- an inert battery genuinely has zero derivative), and
`NaN * 0` then poisoned every directional derivative; and the shallowest-battery
probe direction selected that inert battery, giving an identically flat probe.
Probe directions now exclude root-bus batteries and report how many were
excluded. Follow-up on whether to fix the transcription itself is tracked
separately -- it is deliberately **not** changed here, since altering it would
move published ieee2522 numbers.

## Reverse export

`P_Subs >= 0` with no upper bound, so excessive discharge could in principle be
infeasible. Scanned all 24 hours at `all_max_discharge`:

**It cannot bind on either system, for a structural reason.**

| system | total battery | min net load | closest `P_Subs` to zero |
|---|---|---|---|
| ieee123C_1ph | 0.5066 pu | 0.6291 pu | 0.12448 pu at `t=13` |
| ieee2522C_1ph | 1.3301 pu | (higher still) | 0.59872 pu at `t=13` |

Full simultaneous discharge still leaves a positive import in both cases, before
losses (which push `P_Subs` further up, never down). On IEEE123 reverse export
would need a battery fleet ~24% larger than the system has. **Not one of the 8
IEEE2522 infeasibilities is an export violation** -- that slack is zero in all
of them. On the three study systems the `P_Subs >= 0` floor is simply not the
active limit; undervoltage is.

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
2. **Is the battery-power box empirically fully recourse-feasible?** **No.** It
   is on IEEE123 (111/111) and it is **not** on IEEE2522 (103/111; 8 rejected by
   undervoltage). `F_t` is a strict subset of the box on a genuinely
   voltage-constrained network. Note also that even the IEEE123 "yes" is
   empirical coverage, not certification.
3. **If not, which combinations and constraints define its boundary?**
   **Undervoltage under charging**, on IEEE2522. The boundary tightens with
   demand (1 pattern fails at low demand, 1 at medium, 5 at high) and depends on
   the *spatial* distribution of charging, not only its total: the
   `opposing_deep_chg_shallow_dis` pattern has `sum(P_B) = -0.017 pu` and still
   fails at high demand. `F_t` is therefore genuinely multidimensional and cannot
   be expressed as a tightened bound on aggregate battery power. Reverse export
   is never the binding limit on any of the three study systems.
4. **Is `P_Subs` effectively determined by fixed battery dispatch and network
   losses?** Yes. `P_Subs = net_load - sum(P_B) + loss` holds identically, and
   the loss term varies only 0.002-0.049 pu across the entire survey. `P_Subs`
   is a dependent quantity, not an independent outer decision.
5. **Which genuine network actuators must remain as outer decisions?** None.
   The only network actuator is DER reactive power, it carries no inter-period
   state, and the inner solve already drives it to its limits (`|q_norm| = 1`
   throughout on IEEE123, heavily active on IEEE2522). Promoting it to the outer
   layer would buy nothing, because there is no unused reactive headroom left for
   an outer layer to exploit. What must be added to the outer problem is not an
   actuator but the *feasibility set* `F_t` itself.
6. **Does `Phi_t` appear smooth and convex enough for a battery-only DDP
   method?** Where it is defined, yes, on both systems: first derivatives stable
   across two decades of step size (worst relative spread `7.3e-06` on IEEE123,
   `9.3e-05` on IEEE2522) and curvature positive in every usable probe (137/137
   and 132/132). IEEE2522 shows many more active-set changes (40/132 vs 10/137)
   without losing that stability. The real caveat is not smoothness but
   **domain**: on IEEE2522 `Phi_t` is not defined on the whole box, and the
   centralized optimum *rides* the boundary of its domain at 5 of the 24 hours
   (the voltage floor is active at t = 1, 2, 17, 18, 19). Any outer method must
   be able to sit on that boundary, not merely avoid it.
7. **What must the inner solve return to reproduce FilterDDP's backward
   recursion?** At minimum `Phi_t` and `dPhi_t/dP_B` (the battery-bus balance
   duals, free from the solve). A second-order outer method also needs
   `d2Phi_t/dP_B2`; that was obtained here only by finite differences, and
   getting it analytically would require a sensitivity solve against the inner
   KKT system -- exactly the cost structure the reduced space is meant to avoid,
   and the open question for phase 2.
8. **Exact reformulation or approximation?** The *elimination* is exact on both
   systems -- no constraint is dropped, the eliminated variables are recovered to
   `1e-9`-`1e-7`, and the gradient is analytic. But the reformulation is only
   equivalent to the original problem if the outer layer also represents `F_t`.
   Stated as originally proposed -- outer keeps only battery energy, dispatch,
   dynamics and bounds -- it is **not** an exact reformulation on IEEE2522,
   because it admits dispatches the network cannot serve.

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
