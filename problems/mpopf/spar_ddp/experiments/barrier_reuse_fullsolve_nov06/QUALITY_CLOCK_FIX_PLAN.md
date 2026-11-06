# Correct target-time accounting without repeating cache validation

The original extractor sets remainder=full optimizer minus sum ITER_TIMING
wall times and adds ALL of it before the first quality point. It includes final
solve-exit/retirementGC. This is a conservative prefix estimate, not an actual
wall-clock target timestamp. Cache retirement is larger, so comparing those
upper estimates introduced a directional bias. Original raw measurements stay
untouched; audit_quality_clock.py computes logged-prefix work and honest bounds.
Do not call T24+6.13% or T96-3.16% directly measured target changes.

Full-convergence13.909->12.4145 and27.0085->24.3585 are unchanged. Logged prefix
work improves4.62%/9.22%, but that alone is not total elapsed target wall time.
T24's exact wall-target direction is unresolved; T96 bounds still imply earlier
attainment but do not establish an exact3.16% gain. Before T24 target, three
hits save~1.21s scan work, fresh cache population adds~0.14s scan overhead and
pass-start/GC adds~0.74s. Net logged prefix is~0.33s lower. The legacy estimator
then adds~0.78s extra unallocated solve time upfront, reversing the apparent sign.

Current T144 must finish untouched. No new cache diagnostic, central solve,
worker sweep or repetition of completed benchmarks. Report its legacy estimate,
logged prefix and bounded wall-time interpretation if no direct timestamp exists.

Before submitting not-yet-queued T288, prepare a NEW immutable snapshot from its
prepared source with ONLY cheap passive timestamp instrumentation. Record an
elapsed clock after each backward pass, before convergence/barrier/forward
branches. Start its timer at the same outer optimizer-call boundary as the full
clock (including trajectory initialization). Set/warm the marker for full
warmup too, preventing first measured print/method compilation. Both arms use
IDENTICAL clock-instrumented source. Do not change algorithm, gates, inputs,
settings, worker policy, warmup or timer boundaries. Preserve the old prepared
T288 source/manifest and exact small diff; prove all files except the two clock
hook sites unchanged, plus unchanged model/curvature source. Pin new manifests
before queueing, parse/static sanity check only; no new diagnostic solve.

Suggested implementation in private copied solve.jl: benchmark-start Ref{UInt64}
and one FILTERDDP_ELAPSED_CLOCK marker with k/j and elapsed_s computed at end of
backward pass. The copied outer ieee123c_filterddp driver sets the reference
immediately at both warm and measured solve-call starts. Align marker keys to
original primal/objective trace in POSTHOC collector. Preserve legacy raw CSV;
use actual marker for elapsed target only when complete/monotonic/key-matched,
otherwise report prefix/bounds. Do not assume missing final cleanup happens
before the target. Include the same instrumentation in baseline and candidate,
with all clock overhead charged. No new warmup/validation campaign.

Continue per-case evidence commits and original scientific checks. Correct PR
and saved follow-up prompt so legacy target claims are not repeated. No main
benchmark row or production-default extractor edits are authorized/needed.
