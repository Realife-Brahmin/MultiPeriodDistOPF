# Large10k T24/T48/T96 rho sweep

Frozen campaign source: 492451e3a3dfdf2cdc48273db75e71656374a8a8.
Input/model gate 31652339 passed. Jobs: T24 31652340 (completed), T48 31652341
and T96 31652342 (running at this checkpoint). All candidate outcomes will be
retained; a completed collection validator does not imply physical feasibility.
The [predeclared campaign plan](../../HORIZON_SWEEP_NOV08.md) states settings,
resource limits, stopping rules and the extrapolated T96 anchor.

## Completed T24

Both measured trajectories reach owned AND every-copy RMS <=1e-5 and unique
dual RMS <=1e-2. Whole warmup/measured numerical trajectories match after fresh
primal/dual resets. Physical validation is outside core solve timers.

- rho0=15000: final RMS goal at 103 iterations / 97.5177 s core; relative
  objective gap 6.38549e-7; original maximum violation 5.50819e-4.
- rho0=60000: 101 iterations / 92.3924 s core; objective gap 4.07301e-7;
  original maximum violation 5.51150e-4.
- Historical paper-style owned RMS <=1e-3 AND 0.5% objective gap occurs at
  38 / 32.9811 s and 40 / 32.7552 s respectively. These points have original
  violations 8.84612e-3 and 1.37849e-2; they are not strict physical NO points.
- With owned AND every-copy RMS <=1e-4 AND 0.1% cost gap, first crossings are
  62 / 55.5470 s and 64 / 54.9562 s. Dual residual and physical error are
  separately recorded, rather than silently appended to that criterion.
- Neither trajectory reaches the separate physical NO target (0.1% objective
  gap AND maximum original violation <=1e-5) before its RMS termination.

The complete central Ipopt/MA57 solve is 33.8596 s (same tight tolerances and
full warmup/reset). Its original-printed-infeasibility posthoc candidate for
0.1% cost and 1e-5 primal appears at iteration 11 / 16.7376 s. This candidate
is not an independently checked early central dispatch; full centralized
convergence and its validated final dispatch remain the primary reference.

All 19 endpoint/unique-checkpoint OpenDSS log records report converged flows,
zero voltage violations, and maximum optimizer/DSS voltage difference <=1.48e-10
pu. This does not eliminate battery SOC inconsistency. At the two final RMS
points the maximum cone equality gaps are about 7.86e-6 in the model's squared
per-unit quantities, and fictitious line loss is about 0.0002534 kWh. These
are measured relaxation diagnostics, not a general proof of exactness.

## Battery error in physical units

A passive audit of the actual saved arrays confirms that the final largest
bound violation is reconstructed battery energy, at bus 88006, period 21:
0.550819 / 0.551150 kWh for the two starts, or 0.041316 / 0.041340 percentage
points of SOC on a 1333.2-kWh battery. Consensus-energy box bounds and battery
power bounds have zero violation. The per-step SOC dynamics mismatch is
0.224601 / 0.249176 kWh; the largest integrated/consensus energy disagreement
is 0.572461 / 0.572795 kWh. These quantities describe different checks.
The worst overshoot is 95.0413% versus the 95% operating upper bound, under
the lossless battery model. The [SOC trajectory and error plot](T24-SOC.png)
([PDF](T24-SOC.pdf)) shows both the full trajectory and the deviation scale.

The absolute engineering error is small, but these points do not satisfy the
pre-existing strict 1e-5 physical threshold. Any engineering tolerance or repair
step must be declared separately; this report does not relabel them as exact
feasible solutions. [Audit results](T24-energy-audit.toml) and compressed actual
endpoint arrays are included so the interpretation is independently reviewable.
The passive audit ran locally with Julia 1.12.6; numerical solves used Kamiak
Julia 1.12.4. Reconstructed errors exactly reproduce the saved trace metrics.

## Resource and tuning evidence

T24 central measured-phase peak aggregate PSS is 3.04 GiB. The two tADMM
measured phases peak at 22.44 and 23.07 GiB, with 25 processes. Summed RSS
peaks are 28.09/28.76 GiB and double-count shared pages; PSS is reported
separately. The current process implementation does not show an aggregate
memory advantage at T24. Both methods have the same 64-CPU/128GB allocation;
central MA57 uses one thread, while tADMM has 24 pinned workers and a coordinator.
Sampling is pinned to a spare CPU outside the 25-CPU compute mask.

The complete T24 workflow takes 1072.0 s including central reference,
construction, full warmups, both candidates and saved-point validation.
That is not substituted for core solve time, nor is tuning treated as free.
Only one measured trajectory per candidate is reported; no final repeat-median
benchmark or broadly optimal rho heuristic is claimed.

## Files

- [T24 convergence plot](analysis-T24/horizon-convergence.png),
  [PDF](analysis-T24/horizon-convergence.pdf), [summary](analysis-T24/summary.json),
  [all first crossings](analysis-T24/milestones.csv).
- horizon-T24-31652340.tgz: source/input/hardware, complete logs/traces,
  validation, checkpoint hashes and OpenDSS records; large serialized files
  remain on Kamiak.
- T24-energy-audit-inputs.tgz: completed input instance and two saved endpoint
  arrays for the passive physical-unit audit.

SHA256:

- horizon-T24-31652340.tgz:
  307ffcf2d8a384d54a1aba1e784dbb7d94555cf2986064feeaa7d5810a44d5c4
- T24-energy-audit-inputs.tgz:
  7577c4a82ca2afd8ce03f3e5096300259e95cd2e20f8eeb0725c03a1c696533b
