# Fresh tree-elimination ablation (20 October 2026)

This is an **experimental diagnostic**, not a replacement for the validated
October 16 paper campaign. Never copy these results into the paper automatically.

## Questions and controls

1. **Topology-specific gain:** radial elimination versus generic sparse UMFPACK
   on the network block. Both retain the exact same battery Schur interface,
   second `2nB -> nB` elimination, exact curvature, policies, temporal suffix
   scan, initialization, screening, and tolerances. The UMFPACK backend forms
   the selected Schur complement through batches of 16 RHS columns, not an
   explicit inverse. It is a reproducible generic baseline, not a claim of
   best-possible tuning across every sparse solver or ordering.
2. **Broader selected-computation gain:** `replay.jl` additionally compares a
   full sparse KKT factorization with blocked full network sensitivity solves,
   on identical frozen matrices. This is a local linear-algebra diagnostic;
   its ratio is not a complete-optimizer speedup or the tree-only benefit.

The production default remains `FILTERDDP_NETWORK_BACKEND=tree`. The opt-in
`umfpack` route rejects unsupported diagonal/non-halved paths rather than
silently running another algorithm. Both process-parallel arms explicitly
use `FILTERDDP_PROCS_BALANCED=1`: combinations are assigned evenly across
workers, independently of time-step ownership. The original scheduling
remains the default when that experimental flag is absent.

The new generic route still uses the existing interface ordering metadata.
It does not factor network blocks along the feeder tree, assume independent
feeders, truncate battery couplings, approximate Hessians, or remove batteries.
It supports batteries at both a parent and a child.

## Local correctness smoke

Use Julia **1.12.4** and `envs/spar_ddp`, per the current DDP4OPF compatibility.
Instantiate in a private checkout/depot; never modify an active MSOPF environment.

```bash
julia --startup-file=no --project=envs/spar_ddp \
  spar_ddp/examples/power_system/tree_ablation_smoke.jl
```

The branching synthetic system has batteries at adjacent parent/child buses.
The test adds coupled future battery curvature and checks the complete Newton
correction, selected responses, policy feedback and value-curvature update
against the tree path and an independent full UMFPACK solve. This checks algebra,
not real-network performance. A separate IEEE123 process-parallel pair and
OpenDSS replay must pass before larger timing runs.

## Kamiak campaign

From this committed branch on Kamiak, with Julia loaded:

```bash
bash spar_ddp/experiments/tree_ablation/submit.sh \
  /absolute/path/to/current/periodic/exports \
  /absolute/path/to/a/new/tree_ablation_oct20_campaign
```

The first directory must contain the three current `T96_periodic` exports
and the IEEE123 `T3_periodic` export. They are copied and hashed, not regenerated.
The launcher creates a private source/results snapshot and Julia depot and
submits **one exclusive node** (64 allocated CPUs, 63 workers plus coordinator).
It never cancels, edits, or shares a compute node with Claude's MSOPF jobs.
It uses the site's default account/partition; normal `SBATCH_*` environment
settings may be set by the user before submission.

The job runs the synthetic test, a two-worker IEEE123 T3 complete-solve pair,
and then three T96 systems x three interleaved repeats x two backends. Arm
order reverses on repeat 2. Both arms run to full convergence, use full warm-up,
and require finite solutions, final residuals, solver success, close physical
trajectories, and independent OpenDSS convergence with zero voltage violations.
All pairs run sequentially in the same allocation; compilation, setup and
OpenDSS are outside the optimizer timer. Slurm timeout/OOM/failure must not be
reported as a completed timing row.

The campaign currently reports **full-convergence optimizer time**, phase A/B/C
wall-time components and iteration counts. It does not invent posthoc time to
near-optimality from final-only residuals. A paper-ready near-optimality ablation
must reuse/extend the frozen campaign's complete per-iterate original-feasibility
logging. Comparing backend ratios at full convergence is useful bookkeeping
but remains distinct from the paper's time-to-near-optimality headline.

Outputs: raw logs, per-arm saved dispatches, OpenDSS reports, pair validation,
input/source hashes, Julia manifests, CPU/allocation details, `summary.csv`,
`results/per_run.csv`, and an explicit `campaign.complete` marker. Summed worker
peak RSS is not simultaneous node memory. CSV replay `allocated_bytes` is
cumulative allocation, not peak RAM. No memory-saving claim is made from it.

## Frozen matrix replay

A separate diagnostic process can set `TREE_ABLATION_CAPTURE_DIR` to a fresh
private directory (one Julia thread/process). `TREE_ABLATION_CAPTURE_CALL=1`
selects that process's first network factorization. The captured matrix is
**before future-curvature injection**; it is not a complete optimized trajectory.
Capture I/O contaminates timings, so never enable it for a timed arm.

```bash
julia --startup-file=no --project=envs/spar_ddp \
  spar_ddp/experiments/tree_ablation/replay.jl capture.jls fresh-output.csv 5
```

All three routes warm up before alternating repetitions and must agree on the
full correction and retained sensitivity rows. Layout construction is shared
and outside this kernel timer; extraction, factorization, Schur construction,
battery solves and recovery are inside it. Preserve the capture's SHA256 and
metadata. Keep local smoke measurements separate from Kamiak timing results.
