# Bounded coupling-retaining replay (no optimization)

replay.jl consumes one hash-verified selected snapshot, original model input,
new output directory and window length (0 means full captured horizon).
replay_job.sh gates on all capture coverage validators, then tests all five
IEEE123T3 snapshots before mediumT24 and original1020 four-period windows.
Every process is capped at8minutes; batch capped at45minutes. No new capture,
optimizer solve, production hook or benchmark rerun is involved.

The reference is ordinary dense suffix recursion over all24 original leaves.
Window closure uses that reference at period4's exit, NOT zero. The first
gate reconstructs saved full control/equality policies from deep-copied local
factors at identical zero/sine perturbations (<1e-7 unit-floor Frobenius error).
Only then run the actual suffix-scan association pattern, preserving every
stage's starting interval. Compare all three matrices and both vectors at
each merge and each final suffix, plus battery state/power directions and
selected full control/equality/bound-multiplier directions. Dummy lambda in
policy reconstruction is never reported as an original dual residual.

Exact assessment precedes approximate budgets: zero-truncation hierarchical
storage preserves all numerical entries (dense fallback if factors cost more).
The saved KKT layout records local parent interfaces and battery/energy rows.
A small per-period interface does NOT prove bounded temporal rank: composing
multiple periods can concatenate independent separator responses, and Schur
elimination can introduce fill. Without a closed exact separator algebra and
rank bound, this experiment establishes only exact reconstruction via fallback,
not a new exact compact composition algorithm. Full cross-block spectra in
zero-truncation mode support that assessment; numerical rank alone is not proof
of algebraic exactness. General battery-weighted graph separators descend below
the root; no feeder-number-specific grouping is used.

Adaptive candidates truncate only hierarchy off-diagonal blocks of the joint
two-boundary Hessian at global relative Frobenius budgets1e-8/1e-5 per element.
All local blocks and both vectors remain full; cross terms are retained when
rank factors are useful, otherwise dense fallback. A separate dropped mode
zeros inter-group blocks at every element merge for comparison. No physical
battery/state/network/location/limit is changed.

Outputs: errors.csv (absolute/reference-relative/unit-floor errors), ranks.csv,
spectra.csv, inertia.csv (raw eig extrema, machine-scaled zero threshold,
negative/near-zero counts, absolute condition), costs.csv (combination plus
construction/materialization/allocation/GC), directions.csv, probes.csv,
separators.csv and metadata with input/leaf hashes and actual sampling labels.
The global Hessian budget does NOT automatically bound Newton-direction error,
preserve inertia, or imply convergence; those are explicit candidate decisions.

Important limits: SVD construction is diagnostic and can dominate; this code
holds BOTH the payload and dense expansion. Payload coefficients exclude
indices/tree/object overhead; no process-RAM saving is claimed. Dense expansion
plus production dense combination is NOT a compact-kernel speedup. Timings
include cold compilation and single observations, are serial diagnostic clocks,
and are not Phase-B parallel wall times or reportable solve speedups. Closed
suffix values are currently dense; no storage saving is claimed for them.
No approximate full-solve campaign until real replay accuracy/inertia/storage/
cost evidence justifies one. Near-optimal and full convergence remain separate.
