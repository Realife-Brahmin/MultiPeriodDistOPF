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

1. **The boundary is undervoltage under charging.** Restoration depths reach
   `7.4e-2` in `v` units at `t=6` `all_max_charge`, i.e. about 0.91 pu voltage
   against a 0.95 floor. An earlier version of this note added "and it tightens
   with demand"; **that is wrong** -- see the horizon-independence section below,
   where the low-*net*-load hours turn out to be harder than the medium ones
   because they are the high-PV hours, and PV output consumes the DER reactive
   headroom that would otherwise hold voltage up.
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

## Soft voltage limits: fixed vs adaptive penalty (IEEE2522, T = 3/6/12)

Replacing the hard voltage box with an L1 penalty on violation slacks makes
`Phi_t` defined on the whole battery-power box instead of only on `F_t`.

**A fixed coefficient does not work, and fails in the dangerous direction.**
Over 324 paired solves: at `rho = 1e2`, 10 of 90 feasible cases are not
reproduced -- the penalty buys cost reductions with real voltage violations (up
to `2.4e-4`, `Phi` off by 0.127 USD) and still reports a converged solve. At
`rho = 1e4` all 90 are exact. At `rho = 1e6` the answer is physically right but
numerically degraded by conditioning. The usable window is about two decades
wide and its lower edge is roughly `max|voltage dual|`, which moves with price
level, network and operating point -- so a hard-coded `rho` is per-system tuning
whose failure mode is a *silent* voltage violation.

**Convergence is universal, not cleaner.** Every penalised solve converges (0
failures vs 18 infeasible under hard limits), but at the exact coefficient each
solve costs more: mean iterations 35.8 -> 43.4, worst 52 -> 94, mean time
0.90 -> 1.64 s.

**The adaptive scheme removes the constant and reproduces the hard verdict
exactly.** `inner_opf_adaptive` raises `rho` until the violation either vanishes
(feasible) or stops moving between rounds (genuinely infeasible), using the
measured fact that the minimum violation is `rho`-independent for a truly
infeasible dispatch. Validated against the hard solve on 108 decisions:

| metric | result |
|---|---|
| verdict agreement | **108 / 108** |
| false "feasible" (silent voltage violation) | **0** |
| false "infeasible" (legal dispatch discarded) | **0** |
| inconclusive | **0** |
| max `Phi` gap where both feasible | `2.0e-07` USD |
| inner solves per decision | 1 .. 2, mean **1.26** |
| wall time per decision | 2.13 s (vs 0.90 s for one hard solve) |

The mechanism is visible in the round counts: 80 of 90 feasible cases settle in
one round at `rho = 1e2`; the 10 that need a second round at `rho = 1e3` are
*exactly* the 10 that a fixed `rho = 1e2` would have silently relaxed. All 18
infeasible cases take two rounds and are identified by violation stabilisation,
with settled violations spanning `0.55` to `64.6` -- a usable measure of how far
outside `F_t` a dispatch sits.

**What this still does not do.** It is a reliable *interface*, not yet an
algorithm: it tells an outer layer whether a dispatch is servable and by how
much it misses, but not what to do about it. And the L1 penalty is nonsmooth
exactly at the constraint boundary -- where the centralized optimum sits at 5 of
24 hours -- which cannot be tested until something iterates on top.

## `F_t` is horizon-independent; `T` only selects which snapshots get tested

Worth stating plainly because it governs how far these results generalise. The
inner problem is a single snapshot: network equations, load and PV at instant
`t`, `P_B^t` fixed. Nothing in it references `t-1`, `t+1`, `B^{t-1}` or `T`. So
`F_t` is fixed by (network, load at `t`, PV at `t`) alone -- it is not even
price-dependent, since price moves where in `F_t` the optimum sits, not where
the boundary is. `T` enters only by choosing which snapshots the tADMM profiles
resample to.

Verified rather than assumed. Across `T = 3/6/12` the nine sampled slots are
only **eight distinct snapshots**: `T=3` is degenerate (its "medium" and "high"
are the same 2.7541 pu, and price is a constant 0.1400 -- the sampling artifact
`CLAUDE.md` warns about at `T=3`). Three slots at 2.7541 pu -- `T=3` `t=1`,
`T=3` `t=3`, `T=6` `t=1` -- return the identical infeasible pattern and the
identical violation `v0 = 7.81`. Same snapshot, three horizons, same answer.

**Difficulty is not monotone in demand, by either measure.** Ordered by net load:

| net | gross load | PV | sum reactive headroom | `all_max_charge` violation |
|---|---|---|---|---|
| 1.9011 | 3.2312 | 1.3301 | 0.8826 | 19.8 |
| 1.9567 | 3.2868 | 1.3301 | 0.8826 | 27.2 |
| 1.9917 | 3.3218 | 1.3301 | 0.8826 | 32.2 |
| 2.6667 | 2.6667 | 0 | 1.5963 | 2.28 |
| 2.7541 | 2.7541 | 0 | 1.5963 | 7.81 |
| 3.1043 | 3.1043 | 0 | 1.5963 | 58.6 |
| 3.1393 | 3.1393 | 0 | 1.5963 | 64.5 |

The low-net-load hours are the **highest gross load in the set** -- PV is masking
1.33 pu of it -- and since `qmax = sqrt(S_D_R^2 - p_D^2)`, that same PV output
consumes **45% of the DER reactive headroom** precisely when the feeder needs it
for voltage support. Two effects fight (local PV injection helps voltage, lost
reactive capability hurts it), and the result is monotone in neither net nor
gross load.

Consequence, and it is the same shape as the aggregate-battery-power result:
**no scalar load index is a sufficient statistic for `F_t`.** An outer layer
cannot pre-screen hours as "easy" or "hard" from a demand number; difficulty
depends on the full spatial snapshot including inverter reactive headroom.

## Feasibility cuts: representing `F_t` without writing it down

The adaptive penalty is a diagnostic -- it says whether a dispatch is servable
and by how much it misses, not which way to move. A Benders feasibility cut
supplies the direction. At an infeasible `P_B^0`, solve the *pure feasibility*
inner problem (`pure_feasibility=true`: minimise violation alone, no priced
substation term, so the balance duals mean `dv/dP_B` rather than a
cost-plus-violation mixture) and linearise:

```
v0 + g' (P_B - P_B^0) <= 0,    g = d v / d P_B = real-power balance duals
```

The inner problem is the BFM SOCP, so `v` is convex in `P_B` and its
linearisation is a global under-estimator: the cut cannot remove a dispatch the
network can serve. That is the claim, and it is what makes the route worth
having -- the cut is built from duals alone, with no tuned constant and no
per-system bus list.

Tested on IEEE2522 at `T = 3/6/12`: 18 cuts generated at the 18 hard-infeasible
dispatches, each evaluated against all 12 patterns at its own `(T, t)` -- 216
evaluations.

| check | result |
|---|---|
| **validity** -- cuts excluding a hard-feasible dispatch | **0 of 168** |
| margin at the closest feasible dispatch | `-4.54` (median `-71.7`) |
| self-exclusion -- cut removes the point that generated it | 18 / 18 |
| **strength** -- other infeasible dispatches also removed | 13 of 30 (**43%**) |
| gradient vs central differences | max rel err `6.9e-07` |
| cut cost | 1.82 s, one extra inner solve |

Validity holds with room, not marginally: no feasible dispatch comes within
`4.5` of being cut off. Strength is the meaningful number -- **each cut already
excludes 43% of the other infeasible dispatches it never saw**, so the cuts are
learning the shape of `F_t` rather than memorising points one at a time.

`g` has 249 nonzero entries out of 250 batteries. The single zero is the
root-bus battery, which is the known transcription defect recorded below, not a
property of the cut.

**This also sidesteps the L1 kink.** The nonsmoothness concern applies to a
penalty sitting in the objective at the boundary where the optimum lives. Cuts
are *linear constraints* in the outer problem, kept out of the objective, so an
outer method sits on them the way it sits on any active linear constraint. If
the outer loop is built on cuts rather than on the penalty, the kink question
does not arise.

Still not established: nothing here iterates. Cut *accumulation* -- whether a
growing bundle converges, and how many cuts a full horizon needs -- requires the
outer loop, which does not exist yet.

## Does the decomposition actually pay? Not yet -- it is currently dominated

The reduced-space run is *correct* (see above: it reproduces the full-space
FilterDDP optimum to `6.0e-09`). It is not *faster*. Measured, `T = 3`:

| system | full-space `nu` | full-space wall | iters | per-iter | reduced `nu` | reduced wall |
|---|---|---|---|---|---|---|
| ieee123C_1ph | 791 | **18.4 s** | 48 | 0.38 s | 102 | **42.6 s** (2.3x slower) |
| ieee2522C_1ph | 13358 | **99.7 s** | 67 | 1.49 s | 500 | **926.8 s** (9.3x slower) |
| large10kC_1ph | 54665 | **2569 s** | 126 | 20.4 s | 2040 | projected ~7 h |

**Correcting an overstatement made earlier in this work.** The claim that the
reduced space does "465x less work per iteration" came from comparing `nu^3`,
i.e. assuming a dense backward pass. The measured per-iteration cost scales like
`nu^1.3` (123 -> 2522) and `nu^1.85` (2522 -> 10k), nowhere near cubic: a dense
factorisation at large10k would be 22.3 GB per stage and roughly 540 s per
factorisation, against 2.25 s measured, so the production backward pass is
plainly sparse (consistent with the UMFPACK multi-RHS artefacts in
`ddp/results/network_filterddp/`). The theoretical advantage being chased was
therefore much smaller than quoted.

**Where the cost actually goes.** Two separate problems, both fatal on their own:

1. *The finite-difference Hessian is O(nB) Ipopt solves per stage.* Even frozen
   (computed once per stage rather than per backward pass) it costs
   `T * nB * c_inner` = 10 s at ieee123, **1335 s** at ieee2522, **~24500 s** at
   large10k -- that is 13x and 9.5x the entire full-space solve, just for
   curvature.
2. *Even with a free Hessian, the value/gradient solves lose.* The reduced
   method needs `T` inner solves per outer iteration = 5.3 s at ieee2522,
   against 1.49 s for a full-space iteration covering all stages. A single cold
   inner Ipopt solve (~35 interior-point iterations) already costs more than one
   full-space FilterDDP iteration.

The ieee2522 run confirms the projection and is otherwise a success on
correctness: converged in **46** outer iterations (fewer than full-space), gap
`1.29e-08` relative, max `|dP_B|` `0.0279` kW, with **95%** of the 926.8 s spent
inside Ipopt (750 Hessian solves + 187 value/gradient at ~0.935 s each).

**Zero infeasible trial dispatches, on both systems.** This is worth recording
because it was not the expectation: ieee2522 rejects 8 of 111 *arbitrary*
dispatches, yet FilterDDP never once proposed a dispatch outside `F_t` on the
path from `P_B = 0` to the optimum. The adaptive-penalty fallback and the
feasibility cuts are therefore insurance rather than load-bearing machinery on
this trajectory -- they have not yet been exercised by a real run, and should
not be described as validated in situ.

**This does not invalidate the decomposition; it identifies what has to change.**
Two standard fixes, both untried here:

* **Warm-start the inner solves.** Every inner solve is currently cold, from a
  flat-voltage start. Between consecutive outer iterations `P_B` moves very
  little, so a warm start should cut ~35 interior-point iterations to a handful.
* **Replace finite differences with NLP sensitivity.** One factorisation of the
  inner KKT system at the inner solution yields all `nB` Hessian columns by
  back-substitution, turning `O(nB)` *solves* into one solve plus `nB` cheap
  back-solves. This is the standard sIPOPT construction.

The arithmetic says both are needed and neither suffices alone. At ieee2522,
with a free Hessian the method still costs 175 s of value/gradient solves plus
51 s overhead against 99.7 s full-space -- 2.3x slower. With a free Hessian
*and* a 10x warm-start gain it lands near 68 s, which finally beats full-space.
Warm-starting attacks both at once, and should help the Hessian most of all:
those perturbed solves differ from the base solve by only `h = 1e-6`, so a warm
start ought to converge them in one or two interior-point iterations instead of
roughly 35.

### Both fixes for the inner solve are now in; the gap narrowed but did not close

Warm start (bound multipliers included -- primals alone leave the solve at 22-28
interior-point iterations instead of ~6) plus a persistent model updated with
`set_normalized_rhs` instead of rebuilt per call:

| `T = 3` | full-space | reduced before | reduced now | still |
|---|---|---|---|---|
| ieee123C_1ph | 18.4 s | 42.6 s | **25.4 s** | 1.4x slower |
| ieee2522C_1ph | 99.7 s | 926.8 s | **369.6 s** | 3.7x slower |

Answers are unchanged: ieee2522 returns the identical objective
(`8515.8804065391`, gap `1.29e-08` relative, max `|dP_B|` `0.0279` kW).

**The warm start partly pays for itself, which is worth recording.** Warm-started
solves return slightly looser duals (dual residual `8.9e-08` against `4.2e-09`
cold), and the outer method then works harder: 50 outer iterations instead of 46
and **442** value/gradient solves instead of 187. The per-solve saving still wins
by 2.5x overall, but the naive "cheaper solve = proportionally cheaper run"
accounting is wrong here.

**Even a free Hessian would not close it at `T = 3`.** Of the 319 s now spent
inside Ipopt, roughly 201 s is Hessian and 118 s value/gradient; zeroing the
Hessian entirely leaves ~168 s against 99.7 s. The floor is that the reduced
method needs about three inner solves per stage per outer iteration (2.4 s at
`T = 3`) versus 1.49 s for one full-space iteration covering all stages.

**So the case for the decomposition is now specifically a large-system case.**
Full-space per-iteration cost grows steeply with `nu` (0.38 -> 1.49 -> 20.4 s
across the three systems) while the reduced method's cost grows only with the
inner solve. Extrapolating the measured inner cost to large10k gives roughly
10.8 s per outer iteration against 20.4 s full-space -- the first size at which
the reduced method would be ahead. That projection is untested and the
one-time Hessian would still dominate unless the low-rank route below works.

**Next, and the reason not to give up on the Hessian.** `H*d` for an arbitrary
direction `d` costs exactly ONE inner solve (perturb along `d`, difference the
gradients), so randomized low-rank SVD applies directly. The structure probe
already measured rank-25 reproducing `H` to 2.6%, comfortably inside the ~6%
that `frozen` is known to tolerate -- so ~35 probes per stage should replace
`nB`. At ieee2522 that is 28 s of Hessian instead of 201 s.

### The low-rank Hessian: right idea, and it does not work as stated

`H*d` for any direction costs one inner solve, so a randomised Nystrom sketch is
available, and the required rank genuinely does NOT grow with the system:

| | `nB` | effective rank @1% | ideal rank-25 error | Nystrom 25+10 | Nystrom 50+15 |
|---|---|---|---|---|---|
| ieee123 | 51 | 36 | 2.6% | 3.58% (35 solves) | -- |
| ieee2522 | 250 | **59** | 2.3% | 5.36% (35 solves) | 2.74% (65 solves) |

`nB` grows 5x while the effective rank goes 36 -> 59. So the sketch cost is
roughly size-independent, which is exactly the property the large systems need.

**But it does not converge.** rank-50+15 on ieee2522 `T = 3` (2.74% Hessian
error) stalls at the 200-iteration cap, status 8, dual residual `1.0e-02` -- the
same failure as `battery_only`, not a slowdown.

**An inference made earlier in this work was wrong and is withdrawn.** From
"`frozen` works while the Hessian drifts 1.1-6.3% between points" it was
concluded that roughly 6% Hessian error is tolerable. That does not follow.
Smooth drift perturbs the whole matrix consistently; low-rank truncation zeroes
*specific directions*, and a Newton model with no curvature in a direction takes
an unbounded step along it.

**Why this system is unusually unforgiving, measured.** The network cases carry
`C_B = 1.4e-07` and `dt = 8` (NOT the `C_B = 0.05` recorded in `CLAUDE.md`,
which is the copper-plate value), so the battery term adds just **2.24** to the
diagonal while `d2Phi` eigenvalues run **1.19 .. 2482**. This is the tADMM
regime. The consequence is structural: **`d2Phi` is the entire curvature model**,
there is no well-conditioned term underneath it to regularise the step, and any
approximation that drops directions is therefore fatal rather than merely
inaccurate. It also explains `battery_only` -- `2.24*I` against a matrix reaching
2482 is not an approximation of anything.

The implied fix -- an eigenvalue FLOOR replacing the truncated (~0) eigenvalues
with the smallest captured one -- was implemented (`REDUCED_HESS_FLOOR`) and
**tested: it helps but does not fix it.** Same rank-50+15 configuration on
ieee2522 `T = 3`:

| | iterations | dual residual | wall | verdict |
|---|---|---|---|---|
| Nystrom, no floor | 200 (cap) | `1.04e-02` | 447 s | stalled |
| Nystrom, floored | 200 (cap) | `6.34e-04` | 640 s | stalled |

The floor improves the dual residual 16x -- the diagnosis was right -- but does
not reach the `1e-7` tolerance, and it costs more inner solves (2178 vs 810)
because the over-stated curvature shortens every step. Confirming the diagnosis
is not the same as fixing the method.

### Curvature summary: only the exact per-stage Hessian converges

| curvature model | cost | result |
|---|---|---|
| exact, recomputed per backward pass | `nB` solves per stage per iteration | converges (409 s, ieee123) |
| **exact, frozen per stage** | `nB` solves per stage, once | **converges** (25.4 s / 369.6 s) |
| Nystrom low-rank, floored | ~65 solves per stage, once | **stalls** at cap |
| Nystrom low-rank, unfloored | ~65 solves per stage, once | **stalls** at cap |
| battery term only | free | **stalls** at cap |

Every approximation tried fails, and the reason is the `C_B = 1.4e-07` regime
above: `d2Phi` is the whole curvature model, so there is nothing to fall back on
when it is approximated. The remaining idea is therefore not a better
approximation but a cheaper route to the EXACT Hessian: NLP sensitivity, where
one factorisation of the inner KKT system at the inner solution yields all `nB`
columns by back-substitution. The repository already contains multi-RHS
UMFPACK machinery of exactly this shape
(`ddp/results/network_filterddp/large10k_umfpack_parallel_rhs.csv`).

**Recommendation on large10k: do not run it yet.** With the only convergent
curvature model, its Hessian alone costs `1020 * 3 * 1.316 s = 4028 s` against
**2569 s** for the entire full-space solve -- roughly 80 minutes to produce a
result that loses. It becomes worth the machine time only once the Hessian is
cheap and exact.

Until a curvature approximation is found that actually converges, the honest
statement is that the two-stage workflow is exact and convergent only with a
per-stage exact Hessian, and on that basis it is still slower than solving the
problem whole at both system sizes tested.

## Sizing `C_B`, and why the charge/discharge split is not needed here

`C_B` enters the stage Hessian only as a perfectly-conditioned diagonal
`d = 2*C_B*S^2*dt` added to every eigenvalue of `d2Phi`, so it can be sized from
measurement rather than tuned:

* **Floor (mandatory):** `d > -lambda_min(d2Phi)`, or the stage Hessian is not
  positive definite. Measured on ieee2522, `d2Phi` spans `-0.148 .. 21194`, so
  `C_B > 9.3e-9`. The exported `C_B = 1.4e-07` clears this only barely.
* **Conditioning target:** `d ~ lambda_max/kappa`. `kappa = 100` gives
  `C_B ~ 1.3e-5`; `kappa = 10` gives `C_B ~ 1.3e-4`.

This is usable because `lambda_max` is cheap: `H*d` costs one inner solve, so
about ten power iterations price it. No per-system constant is involved.

**A flat direction exists even with losses.** `lambda_min(d2Phi) = -0.148`
against `lambda_max = 21194` on a 2522-bus network: network losses do NOT fully
break the degeneracy that makes battery dispatch non-unique on a copper plate.
So the regularisation has a real job at network scale, not merely a
justification by appeal to operating cost. (The value is within about an order
of magnitude of the finite-difference noise floor, so read it as "flat", not as
established nonconvexity.)

**Charge/discharge split: not needed in this formulation.** The model carries no
round-trip efficiency -- there are no `eta`-like keys in the exported network
data and the dynamics are `B^t = B^{t-1} - dt*P_B^t` -- so `P_C`/`P_D` would be
pure redundancy. Two facts decide this generally:

* `Phi` depends only on the net `P_B = P_D - P_C`, so in split coordinates
  `d2Phi = [[H, -H], [-H, H]]`, rank `nB` out of `2nB`. The null direction
  `(delta, delta)` is exactly the simultaneous-charge-and-discharge direction:
  the network cannot see it at all.
* A LINEAR throughput penalty `c*(P_C + P_D)` has zero Hessian. It suppresses
  SCD at first order but contributes nothing to curvature, so on its own it
  would leave the stage Hessian singular in `nB` directions.
* A QUADRATIC on the SPLIT variables, `C_B*(P_C^2 + P_D^2)`, does both jobs:
  minimising it subject to `P_D - P_C = p` with both nonnegative lands on
  `P_C = 0` or `P_D = 0` automatically. The same quadratic written on the NET
  `P_B^2` cannot, since SCD leaves `P_B` unchanged.

Decision rule: **no efficiency -> net `P_B` with a quadratic `C_B`, no split.
With `eta != 1` the split becomes mandatory** (the SOC dynamics are not
expressible in net power), and the quadratic must then be written on
`P_C^2 + P_D^2` rather than on the net.

## Measured crossover (2026-09-15), and the large10k result verified

All at `T = 3` with BOTH formulations re-run at the same `C_B = 1e-3`, since a
reference at a different `C_B` is a different problem:

| | full-space | reduced (low-rank) | outcome |
|---|---|---|---|
| ieee123 (`nu` 791 -> 102) | 14.4 s / 46 it | 21.7 s / 26 it | 1.5x slower |
| ieee2522 (`nu` 13358 -> 500) | 107.3 s / 56 it | 107.1 s / 33 it | parity |
| large10k (`nu` 54665 -> 2040) | 1689.98 s / 100 it | **802.49 s / 13 it** | **2.07x faster** |

**large10k verified, not just timed.** The full-space reference was regenerated
at `C_B = 1e-3` (filename tagged so the `C_B = 1.4e-07` reference survives) and
compared directly:

| | value |
|---|---|
| full-space objective | `2998133.704778` USD |
| reduced objective | `2998133.703480` USD |
| objective gap | `4.330e-10` relative |
| max `|dP_B|` | `6.0e-11` pu (0.0000 kW) |
| max `|dB|` | `4.7e-10` pu |

So the speedup is on the same solution, not a cheaper wrong one.

Two trends worth carrying forward. Outer iteration counts move in opposite
directions with size -- full-space 46 / 56 / 100, reduced 26 / 33 / **13** --
because the reduced problem's dimension is set by battery count, not network
size. And the bottleneck has moved: at ieee123 inner solves were ~70% of wall,
at large10k they are **25%**, so further gains now come from the outer solve.

## The method does not survive longer horizons (2026-09-15)

Sweep over `T` and `C_B` on ieee2522, both formulations at matching `C_B`
(`ddp/results/reduced_space/overnight/sweep.csv`, raw logs beside it):

| `T` | `C_B` | full-space | reduced |
|---|---|---|---|
| 12 | 1e-5 | 545 s, 78 it, **status 0** | 3529 s, 200 it, **status 8** |
| 12 | 1e-4 | 542 s, 78 it, **status 0** | 3857 s, 200 it, **status 8**, NaN objective |
| 12 | 1e-3 | 544 s, 78 it, **status 0** | 2692 s, 200 it, **status 8** |
| 24 | 1e-5 | 1112 s, 84 it, **status 0** | 3989 s, 200 it, **status 8**, NaN objective |
| 24 | 1e-4 | 1107 s, 83 it, **status 0** | (operator-killed, no result) |

Dual residuals `7.3e+02`, `1.0e+01`, `1.6e+04` -- diverging, not nearly
converged. Full space is untouched by this and converges in 78-84 iterations
throughout.

**Cause: the feasibility machinery finally engages, and it is not sound.**

```
INFEASIBLE trial dispatches (penalty fallback used):   80   (T=12)
                                                      108   (T=12)
                                                      157   (T=24)
```

At `T = 3` this count was **zero on every system**, which is why the earlier
runs looked healthy. This note previously recorded the adaptive penalty and the
feasibility cuts as "insurance rather than load-bearing machinery ... not
validated in situ". At `T >= 12` they become load-bearing and the insurance
fails.

The defect is not accuracy, it is **consistency**: where a dispatch is servable
the outer method is handed `Phi`, and where it is not it is handed
`Phi + rho*violation` with **`rho` chosen per call by the adaptive scheme**. The
objective therefore changes definition between evaluations, and no Newton-type
method can converge against that. The enormous dual residuals are the signature.

**`C_B` is not the culprit here.** All three values fail alike, which confirms
the `Delta t`-invariance algebra rather than contradicting it. The error was
inferring that `T`-invariant *conditioning* implied `T`-invariant *behaviour*: a
separate mechanism scales with horizon, since longer horizons force deeper
battery cycling and push trial dispatches outside `F_t`.

**Consequences for the two candidate fixes.** A single FIXED `rho` at least
makes `Phi_rho` one well-defined function everywhere, which is the minimum bar
for a Newton method; it does not remove the L1 kink at the boundary. The
feasibility cuts already prototyped (0/168 validity violations, 43%
transferable) are the principled route precisely because they keep the inner
problem hard-constrained and exclude bad dispatches with linear constraints,
leaving the objective a single smooth `Phi` on `F_t`. Neither is yet wired into
the outer loop.

Until one of them is, the `T = 3` results below stand and nothing above `T = 3`
does.

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
