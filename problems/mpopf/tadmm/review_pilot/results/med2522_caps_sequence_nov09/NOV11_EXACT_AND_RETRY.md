# Med2522 T96 single-node follow-up, November 11

Both exact-local comparisons used cn180, one fully allocated 128-core Xeon
6760P node with 96 distinct single-thread workers. The source was frozen at
87a0c749. The serialized optimization inputs had identical SHA256 hashes.

- MA57 (31660388) centralized solved in 63.59 s and 63 Ipopt iterations.
  tADMM reached 1000 rounds in 241.06 s core and 278.94 s loop.
- MUMPS (31660390) centralized solved in 106.60 s and 63 Ipopt iterations.
  tADMM reached 1000 rounds in 378.29 s core and 422.76 s loop.

Both tADMM endpoints missed the residual criterion and had 38.97% relative
objective disagreement and a maximum 0.09661 per-unit original-model
residual. Chronological OpenDSS replays 31681324 and 31681325 found maximum
SOC disagreement of 32.691% of affected battery capacity for both endpoints.
Peak dispatch disagreement was 6.385 kW. Independent maximum absolute cone
equality gaps were 7.792e-7 p.u. squared for centralized and 8.197e-7 for
the endpoints. Relative gap p95 was 0.2394% for centralized and 0.3810% for
the endpoints. The replays and the preserved original audits are outside
the benchmark clocks. These endpoints establish no speedup.

Average per-round median/slowest local solve times were 0.0950/0.1367 s
for MA57 and 0.2021/0.2695 s for MUMPS. Sampled summed RSS/PSS/private
peaks during the measured tADMM phases were 87.41/64.28/64.03 GiB
(MA57) and 89.98/66.28/66.05 GiB (MUMPS). Centralized-phase PSS peaked
at 2.72 and 2.86 GiB. RSS counts shared pages once per process.

Original cap10 job 31667880 exited after recording its allocation because
the inner wrapper required CPU 6438Y+ on the requested 6760P node. Its
dependent cap30 job 31667882 was cancelled. Their original logs and frozen
launcher are archived in cap-launch-failure-nov11.tgz. The retry wrapper
checks the original source manifest and changes only the CPU check in a
temporary launcher copy. The input, model, local caps, and solver settings
are unchanged. Cap10 was resubmitted as 31681326 with audit 31681327 and
cap30 as 31681328 with audit 31681329.

The two exact-run archives preserve original logs, manifests, traces,
expanded chronological and SOCP audits, and Slurm accounting. Full dispatch
files and replay tables remain in their original Kamiak result directories.
