# Kamiak parallel correctness gate

Job 31643439 completed successfully on November 7, 2026. This is correctness
evidence from a shared-node allocation, not benchmark timing evidence.

The archive contains the frozen source/input identities, allocation details,
solver traces, captured source files, resource samples, validation reports and
independent OpenDSS logs for IEEE123 T=6: serial MA57, two-worker MA57 and
two-worker MA97. All cases completed full warm-up and measured restarts.
Large Julia serialized dispatch/instance files are omitted from this backup;
they remain in the original Kamiak campaign directory.

All three cases converged in 73 iterations. Parallel MA57's final objective was
identical to the serial reference; MA97's relative difference was
6.151977859692316e-13. Both parallel cases passed worker-layout, source snapshot,
warm-up/restart, optimization and independent OpenDSS validation.

Archive SHA256:
`a474da1de270af513da8f4c5559e5dc94e06ff6305183e8e6fc7f2d010c2ed41`

The archived orchestration scripts were captured after the gate and identify
the subsequent system-verification and benchmark scheduling protocol. They do
not imply that those later runs have completed. The authoritative gate launcher
is the captured numerical source's `kamiak_parallel_gate.sh`.
