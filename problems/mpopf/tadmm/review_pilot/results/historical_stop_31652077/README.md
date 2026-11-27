# Historical RMS stopping: large10k T24

The user-selected rule stops at **48 outer iterations**: **50.5359 s core** and
**110.4559 s loop time including diagnostics**. These are exploratory single-run
observations after three discarded compile steps and fresh model/state reset;
they exclude startup/warmup and are not warmed benchmark medians.

Settings: MA57, exact locals, July C_B=8.013987384856764e-8, current modified
network and periodic profiles, fixed rho=48, free terminal SOC, proved bound
screening, corrected branch-flow losses. Source/status checks pass. Initial
start hashes and the full 48-step numerical prefix match the previous exact-local
assessment. This exit changes the termination point, not the trajectory.

The rule is owned-period primal RMS <=1e-4 AND (unique consensus-change dual
RMS <=1e-2 OR objective range/max <1e-3 over five iterations). At exit, primal
RMS is 7.33274594e-5 and dual RMS is 1.20867148e-2: the objective stagnation
branch supplies the second condition. This is the selected historical variant,
not the later relaxed July defaults.

Endpoint objective is 2605228.033277688 versus centralized 2605186.166625776:
relative gap 1.60705e-5 (0.00160705%). Maximum original-model violation is
**8.53850406e-3**. Strict convergence and the separate physical NO criterion
are false. OpenDSS converges with voltage agreement 1.48e-10 pu; that agreement
does not establish temporal storage feasibility.

A sub-100 exit is therefore available under the selected rule. Penalty tuning
should assess both this stopping count and the solution quality reached there.
Local cutoff tuning is useful primarily on medium: prior large counts are
median 11/max 24 Ipopt steps, versus medium median 26/max 37 with persistent
harder periods 18-20. See the [cutoff assessment](../july_weight_comparison/README.md).

- [Stopping/feasibility plot](historical-stop.png) and [PDF](historical-stop.pdf).
- [Validation summary](summary.json) and [settings/endpoint](assessment.toml).
- historical-stop-31652077.tgz preserves source, traces, statuses, initial
  starts, OpenDSS logs, hardware/affinity, gate outputs and scheduler states.
  Serialized dispatches/checkpoints remain on Kamiak in
  tadmm_review_nov08_historical_stop.
- Frozen source cad4367702a45bb879c2fea3553e28ece751d4bd; gate 31652076 and
  assessment 31652077 both completed with exit 0.
- Archive SHA256:
  2d2d6ab83353264b2d8de510d1ea7a70ccf0707d2cef9d067ec0e904728b28b4.
