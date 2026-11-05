# Dense merge kernel and worker/thread study, November 5, 2026

The right-sided kernel and the 31-worker/two-thread configuration both improve this frozen Phase B replay. Their combination lowers the median from1.217836s to0.875896s, a28.08% reduction (1.390x). All first-pass and saved late-state accuracy checks pass.

| Kernel | Workers x BLAS threads | Phase B median, s | Range, s | Mean complete merge, ms |
|---|---:|---:|---:|---:|
| current | 63 x 1 | 1.217836 | 1.211392-1.222948 | 134.52 |
| right | 63 x 1 | 1.126903 | 1.122226-1.136242 | 124.33 |
| current | 31 x 1 | 1.590151 | 1.589714-1.590938 | 132.15 |
| current | 31 x 2 | 0.968472 | 0.962320-0.977411 | 81.10 |
| current | 15 x 1 | 2.938739 | 2.922396-2.955443 | 128.95 |
| current | 15 x 4 | 1.088239 | 1.083693-1.089119 | 48.51 |
| current | 31 x 2 | 0.968723 | 0.961563-0.992746 | 81.11 |
| right | 31 x 2 | 0.875896 | 0.873898-0.885279 | 73.73 |

The second31x2 current row is the repeated baseline within the combined-kernel comparison. All values use five clean timed passes, after a warm pass, with alternating order within each same-process pair. Mean merge time sums worker kernel elapsed time and divides by256 merges; concurrent merge latencies do not add to Phase B wall time.

The right kernel wins all five63x1 pairs. The current31x2 configuration beats31x1 in all five pairs, as does15x4 versus15x1. The combined right31x2 kernel beats its own current31x2 baseline in all five pairs. The configuration results therefore include their single-thread controls; a cross-worker-count comparison alone would not isolate BLAS threading.

## What changed

The current left matrix solve is ldiv!(U', Y1) after explicitly transposing Mab into Y1. For this strided Float64 matrix, Julia1.12.7 dispatches through LAPACK.trtrs!; previous native stack samples identify OpenBLAS dtrsm underneath it. The new kernel copies Mab into Y1 and explicitly calls BLAS.trsm!('R','U','N','N',1,U.data,Y1), so Y1=Mab/U. The left Hessian update becomes syrk with trans=N, the cross block becomes -Y1*Y2, and the left gradient becomes ma-Y1*g. The right Hessian and right gradient retain their original Y2 orientation. Only upper triangles are used and mirrored, including asymmetric raw inputs. Both failed-Cholesky branches call the original pint_combine and copy back the same fields.

The speedup covers the complete merge and all surrounding Phase B dispatch/barriers, including matrix assembly and copies. It cannot be assigned exclusively to removing the transpose because the TRSM side and subsequent BLAS orientations also change.

## Accuracy

All96 Hessians and gradients are validated after every warm, timed, and allocation-diagnostic pass. Across64 replay passes the largest scaled infinity-norm error is below3.1e-14, versus the unchanged1e-11 gate. The saved lateT24 trajectory is fully reconstructed with all24 leaves and the true future suffix. At BLAS widths1 and2, value errors are at most5.75e-11 (gate1e-7); recovered control/equality/bound errors are at most8.75e-13 scaled and1.48e-6 relative (gates1e-7 and1e-5). Saved balanced control/equality probes reproduce within2.1e-16 scaled. All24 battery state/power directions pass, with at most2.06e-14 scaled and1.53e-10 relative error. The synthetic regression includes12 cases spanning dimensions1/3/17, both boundary types, raw asymmetry, and six original-fallback cases.

## Allocation and placement

Each worker is pinned to a disjoint group of physical cores within one socket; one core is reserved for the coordinator. The allocations use64 cores for63x1,63 cores for31x2, and61 cores for15x4, including the coordinator. CPU groups, physical topology, actual BLAS thread counts, process IDs and OS thread placement are archived for every mode/pass. The backend is ILP64 OpenBLAS. Node cn154 had no other Slurm jobs at startup. NUMA page placement and first touch were not independently controlled.

The same2,398,648,320-byte shared-slot storage and three-matrix/one-vector PairWork scratch layout are retained. Scratch contains24,977,760 coefficient bytes per worker; the supplementary setup probe measures24,977,976 bytes including object metadata. The initial SCRATCH_OBJECT_BYTES0 log was emitted before worker setup and is invalid. Its raw record is preserved. No scratch or resident-memory saving follows from the kernel change.

Allocation measurement is outside clean timing passes, with its own warm pass. Each configuration/kernel allocates10,752 Julia bytes in all256 merge calls combined and records zero GC and zero Julia compilation/recompilation during that diagnostic. This measures Julia allocation inside the merge, not native BLAS-internal allocation, serializer allocation, or simultaneous resident RAM. GC is explicitly disabled within replay passes and runs between passes, matching the pinned replay protocol. Lifetime worker maxRSS includes shared mappings. Slurm batch maxRSS is69,009,060KiB for the overall study, not a per-configuration peak.

## Decision and scope

The strongest next full-solve candidate is the right kernel with31 workers and2 BLAS threads. A paired optimizer test must include Phase A/C, since fewer workers and wider BLAS teams may change those phases and erase part of the Phase B benefit. A63x1 right-kernel control would isolate the kernel change if the combined candidate regresses elsewhere. No full optimizer run or production integration is included in this batch.

This is the original frozen large10k T96 first backward pass, on one socket-pair node. The late trajectory checks support numerical accuracy at a different iterate/horizon; they do not establish trajectory-wide optimization robustness. The study measures component/replay improvements, not full-solve or near-optimality speedups. Earlier grouped-versus-balanced benefits are separate and should not be multiplied into these ratios without a matched comparison.

## Provenance and reproduction

Base commit dcc3e3096a2539ca1e10609dc8704aa608d36682 on ddp/nov04-grouped-phase-b. No newer remote commit was present when checked. Job31574957 completed0:0 in9m59s. Original T96 capture under stage_b_replay_oct15/capture; all96 manifest hashes passed. Timed numerical source comes from stage_b_replay_oct15/repo with Julia1.12.7. Late validation uses the pinned compact_campaign_capture_oct23_v3 source and Julia1.12.4, as in the earlier late check. Source, input, runtime, native-dispatch, hardware, affinity, Slurm accounting and raw replay logs are included in the evidence archive.

The evidence archive retains the tested source exactly. The readable replay harness fixes only the placement of the scratch-size metadata print, after setup; the original misplaced print and supplementary correction are archived. The collector was corrected to accept Julia's spaces around equals signs in metadata. Scientific gates and raw measurements were unchanged.

Extract the evidence archive and scratch supplement into the same folder, then run analyze.py EXTRACTED/dense_merge_study_nov05 OUTPUT. summary.csv provides the short comparison table; rounds.csv includes round wall times and full task distributions; diagnostics.csv contains allocation/GC/compilation readings; errors.csv contains every first-pass reference check. All experiment code is opt-in.
