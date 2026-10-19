# Kamiak large10k T=96 three-way value-pass diagnostic

Date: 2026-10-15. Frozen snapshot:
`/home/aryan.jha/kamiak_jobs/matched_benchmark_oct15_large96_threeway_v3`.

## Question

The measured one-node pairwise run at T=96 was only 1.53 times faster than
centralized Ipopt. We tested whether this was caused by assigning 96 stages to
only 63 workers, or by the value-function pass itself.

All three large10k runs used Gold 6438Y+ nodes, the same physical starting
point and stopping rules, and 16 FilterDDP iterations. The two-node jobs used
96 workers plus the coordinator. An IEEE123 T=6 smoke test first verified that
the ordinary and explicit-transfer pairwise paths reproduced the same answer.

## Result

| arrangement | nodes | workers | solve (s) | value pass per backward pass (s) | result |
|---|---:|---:|---:|---:|---|
| ordinary value recursion | 1 | 63 | 220.376 | 8.31 | OpenDSS pass |
| ordinary value recursion | 2 | 96 | 219.030 | 8.45 | OpenDSS pass |
| pairwise explicit-transfer recursion | 2 | 96 | 76.062 | 2.19 | OpenDSS pass |

Giving every stage its own worker reduces independent stage work, but it does
not shorten the ordinary serial value recursion. Consequently, the two-node
ordinary run is only 0.61% faster than the one-node run. Pairwise recursion
removes that chain and is 2.88 times faster than the two-node ordinary run.

The existing one-node SharedArray pairwise result (66.482 s median) remains
14.4% faster than the 76.062 s two-node explicit-transfer result. The extra
node therefore does not improve T=96: cross-node transfer and synchronization
cost more than the last stage batch they eliminate. The implementation
bottleneck is the value-function dependency and its data movement, not merely
the number of workers available for stage-local Hessian/Jacobian/KKT work.

All three large10k runs reached objective 2704807.066934211645 and primal
residual 6.559214901225e-06. The tiny difference in the pairwise dual residual
(6.384104506434e-10 versus 6.384101446222e-10) is floating-point rounding.
Every result passed OpenDSS replay without voltage violations.

## Evidence

Raw solve logs, metadata, summaries, OpenDSS replays, completion markers, and
snapshot checksums are archived under
`ddp/results/kamiak/large10k_T96_threeway_oct15/`. The reported worker RSS sums
are sums of individually sampled worker maxima, not simultaneous node peaks.
