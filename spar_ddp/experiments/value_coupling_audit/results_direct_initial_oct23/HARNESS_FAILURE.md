# Preserved first attempt

Job31369046/overlay143a5681 FAILED1:0 in1m16 before warmed timings.
Full source/logs/partial CSV/hardware/accounting are preserved in
failed_harness_evidence.tgz, SHA256
912199f8117151563624dcc7f3e8889ece0d5c390b07935d99a9efc4d5a01b59.
No performance result is attributed to this attempt.

Numerical gates had passed: existing references/independent consumer,
direct leaf generation, poisoning dense Sh0 independence, chained/prefix
values and full directions/negative control. However a nested value-logging
closure assigned d=norm(A-B), capturing the outer d used for the dense chain.
The subsequent warmup di_dense_prefix therefore received a Float64 instead
of a dense-reference NamedTuple (FieldError: Float64 has no field value).
This is a logging/warmup harness bug, not compression/optimizer failure.

Replacement moves that diagnostic computation into a top-level helper with
an explicit local result in the closure. Unit regression checks dense-reference
identity and logging result. No generator arithmetic, factor budgets, consumer
residual/direction gates, cases, hardware or scope changes. New immutable
snapshot only; failed job/source/partial files stay unchanged.
