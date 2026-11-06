# Preserved diagnostic-budget failure (v3)

Job31612252 on cn179 ended naturally at the immutable50-minute process cap:
Slurm FAILED, exit124:0, elapsed50m09s. It was not cancelled or modified.
No complete optimizer solution, scientific completion gate, independent physical/
OpenDSS result, or clean timing pair was produced. Do not use these timings to
judge the optimizer or cache speed.

The terminal traceback is inside the coordinator fresh-shadow reference,
OpenBLAS parallel Cholesky (dpotrf_U_parallel). The coordinator was pinned to one
physical CPU while this diagnostic requested two BLAS threads. Four completed
warmup passes show reference-including round walls near500s versus actual
worker merge maxima summed across rounds of2.5-4.7s after first compilation.
This diagnoses a validation execution problem rather than a measured algorithm
regression. The completed source proof and frozen pilot remain separate.

Before the scratch warmup log was lost at process termination, read-only short
Slurm steps preserved partial progress logs. Latest preserved trace contains
four completed warmup passes (k0-3), all168 actual stage negative controls,
and6252 observed debug comparisons, all with zero absolute/scaled/relative error.
There was no completed eligible hit or full solve in this attempt. Partial
checks cannot replace missing full scientific validation.

Raw source/harness checks before and after pass. Full raw timeout traceback,
Slurm accounting/hardware/worker masks, code, manifests, unit logs and partial
warmup snapshots are retained in raw/ and evidence.tgz; runtime-source.tgz
inside raw preserves the actual private source. Binary inputs remain at pinned
cluster locations. Archive SHA256:
1e400a4e81ca27af5f4ff9505d1b23fd1ddcd96f4eaa04e42601f7b99ddbd064.

The new immutable v4 moves only the diagnostic reference to an existing worker
with two assigned physical CPUs, on independent copies of the synchronized
leaves. It also preserves diagnostic warmup logs, resets negative flags between
solves, and warms/retires cache payloads through the measured exit method.
No cache mathematics, numerical gates or clean baseline configuration changed.
It was submitted only after this job ended and evidence was archived.
