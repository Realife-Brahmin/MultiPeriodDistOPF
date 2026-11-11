# Med2522 cap30 and residual-gated comparison

Cap30 (31681328, audit 31681329) hit 1000 rounds in 243.15 s core
and 279.88 s loop versus 63.42 s fully solved centralized MA57 (63 Ipopt
iterations). The cap activated twice in 96,000 local solves. The endpoint
reproduced the exact fixed-48 result, with a 38.97% objective gap, 96.61 kWh
maximum requested SOC-bound error, and 32.69% maximum chronological SOC
disagreement as a share of battery capacity. Maximum absolute cone gap
was 8.20e-7 p.u. squared versus 7.79e-7 centralized. The p95 normalized
gaps were 0.3810% and 0.2394%.

On identical inputs, the exact residual-gated run (31681398, audit
31681399) also reached 1000 rounds. It took 244.52 s core and 284.85 s
loop versus 63.57 s centralized (63 Ipopt iterations). It reduced the
objective gap to 0.4203%, the largest requested SOC-bound error to 2.771
kWh and maximum chronological SOC disagreement to 17.89%. Maximum
actual dispatch disagreement was still 6.385 kW. The maximum cone gap
was 8.75e-7 p.u. squared and p95 normalized gap was 0.3496%. The
controller improved agreement, but no audited candidate reached the
engineering target. The first saved historical residual-plus-cost
milestone was round 670 at 164.12 s core with an original-model
residual of 0.003224 in the model convention.

Mean per-round median/slowest local solve times were 0.0953/0.1374 s
for cap30 and 0.0870/0.1404 s for residual-gated. Measured tADMM
sampled RSS/PSS/private peaks were 88.10/64.98/64.73 GiB and
87.54/64.77/64.52 GiB, respectively. Centralized PSS peaked at 2.91
and 2.90 GiB. Both used 96 distinct workers on cn180. The source and
input manifests, Slurm accounting, worker traces and chronological
audits are preserved in the two archives.

## July-input initialization failure and retry

Job 31681415 stopped during centralized model construction because the
idle-battery backward/forward load-flow sweep returned NaN at period 1.
Its dependent audit 31681416 was cancelled. The original logs, source
manifest and initializer are archived in july-init-failure-nov11.tgz.
This failure occurred before Ipopt, so it does not establish OPF
infeasibility.

Retry 31681491 (audit 31681492) uses a shared flat network initialization
at 1.05 p.u. for centralized and the tADMM local models. This changes
only the start in the isolated July source copy. The release voltage,
input hashes, OPF equations, residual-gated controller and solver
tolerances are retained. Any comparison of iteration counts or timing
against the current-input run must account for the initialization change.
