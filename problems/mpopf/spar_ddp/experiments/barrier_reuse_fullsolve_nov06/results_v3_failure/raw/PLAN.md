# Opt-in process-parallel barrier reuse: full-solve study

Continuation explicitly authorized by Aryan through the DDP Paper and Theory
thread after the completed PR185 pilot: obtain actual full-optimizer results,
including a negative result. No production default or paper/result-row changes.
Start e6775b2c; remote fetched, no newer work, no active Kamiak jobs. Preserve the
two unrelated battery reference files. This is separate from the reported replay.

Use the validated grouped b2/right kernel with a fresh toggle-only control.
T24: 24 workers x2 BLAS; only after correctness passes, T96: 48x2. One otherwise
unoccupied 128-physical-core Xeon6760P allocation, 600GiB, three hours, no requeue.
Two alternating paired repeats (base,cache / cache,base), equivalent full
warmups. Third pair only for unresolved overlap. No central reruns or sweeps.
Reuse frozen input/start/Cb2.4e-5/gamma420/screening/passGC/original gates.

Private overlay only: the actual solve loop emits solve_start, accepted and pure
barrier events. Every event except barrier advances a source-owned epoch. A hit
requires that explicit event, identical epoch/k/context, consecutive j, reduced
mu, and every worker's complete cache. Unknown state/config/model/terminal/
regularization/path updates must invalidate via the API. The supported private
fixed driver owns its model/data after construction and does not expose mutable
callbacks during solve. Direct external mutation is unsupported. Debug stage
fingerprints audit the source proof; clean timing never hashes/serializes S.

The supported _tp_pre route applies no shift: exact Hessian, identity linear
state dynamics, tree/network factors; unsupported derivatives stop the solve.
Do not infer applied regularization from d.reg_last, which is reset downstream.
Phase A/current network factors remain fresh. Only Phase B matrix algebra is
reused. The right-sided kernel retains U, right-whitened BL, and Y2 or final V.
No gradients or leaf curvature are retained. Current leaf vectors, each current
merge h/solve/two vector updates, all current policies and the full forward path
remain. Re-population reuses storage but rewrites every factor/block. A missing
preflight entry forces whole-pass fresh population; unexpected mid-scan loss
fails explicitly rather than fresh-combining unmaterialized intermediates.
Failed Cholesky follows unchanged right-kernel fallback and invalidates the next
hit. Warmup validity and payloads are dropped before warmup collection. Cache
retirement/full collection is charged inside measured candidate solve exit.

Separate debug T24 before clean runs: fingerprints, actual fresh right-kernel
shadow with every intermediate, population/hit suffixes, current recovery,
full control/equality/lower/upper actions at multiple scales and actual full
forward trajectory. Strict scaled1e-7 and relative1e-5/floor1e-10, bit exact
matrix/vector algebra wherever expected. Explicit invalidation and stale-vector
negative controls. Debug instrumentation is excluded from clean timings.

Validate every full solve (numerical/physical/bounds1e-5, independent OpenDSS,
thread masks/backend/start/source/input/runtime hashes) and compare same-repeat
solutions/iterations at unchanged1e-10. Disclose any trajectory differences;
score existing shared objective0.1% plus original-primal1e-5 target posthoc using
completed same128node central reference. Full optimizer, target, A/B/C, cache
counts, allocations/GC/RSS limitations, setup/build/warmup/process wall remain
separate. Archive raw logs/CSV/manifests/negative controls, commit/push/report.
Never modify an active/queued immutable snapshot; preserve failed attempts.
