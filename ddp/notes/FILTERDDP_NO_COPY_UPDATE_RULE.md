# Removing the duplicate FilterDDP update-rule copy

After factoring the bound sensitivities, the backward pass already owns the
arrays used to construct each stage's `UpdateRule`. The original converting
constructor nevertheless copied every array again with `Vector{T}(...)` and
`Matrix{T}(...)`. For large10k, the copied `beta` and `omega` alone total
754.603 MiB per stage.

The no-copy patch adds a more-specific constructor for already-owned
`Vector{T}` and `Matrix{T}` inputs. It stores those arrays directly. The old
abstract-array constructor remains as a compatibility fallback and still
converts external views or differently typed arrays. The ownership audit found
no mutation of these arrays after assignment; forward rollout only reads them.

On a bounded large10kC `T = 3` iteration, warm-stage update allocation fell
from 4282.176 MiB to 3525.163 MiB, a reduction of 757.013 MiB (17.68%) per
stage. The reduction equals one complete factored update rule to measurement
precision. Persistent storage is unchanged: this removes temporary copying,
not the retained `beta` or `omega` maps.

The full IEEE123C `T = 3` regression remained identical: 48 iterations,
objective `2808.924770488853`, primal residual `6.227371345322e-8`, dual
residual `2.303870314635e-8`, and complementarity `1.624842969030e-9`.
The one-iteration large10k wall time is not used as a speed claim because a
single iteration is noisy; a future full run is needed before attributing a
runtime change to this allocation reduction.

## Full runtime regressions

Strict IEEE2522C cold-start regressions isolate the no-copy constructor from
the preceding factored-bound implementation:

| Horizon | Factored only | No-copy | Change | Iterations | Practical `1e-6` iteration |
|---:|---:|---:|---:|---:|---:|
| 3 | 245.175 s | 237.108 s | -3.29% | 56 | 52 |
| 12 | 1183.039 s | 1131.022 s | -4.40% | 79 | 76 |
| 24 | 2330.122 s | 2249.267 s | -3.47% | 84 | 80 |

Every trace is numerically identical to its factored-only predecessor at the
stored precision, including objectives, all infeasibilities, step sizes, and
line-search decisions. The consistent 3.3--4.4% reduction across three full
runs supports a modest runtime benefit from lower allocation/garbage-collection
pressure. It does not justify calling copying the dominant bottleneck.

Another full large10k `T = 3` run would likely resolve a few-percent effect but
cost roughly two hours. The IEEE2522 series already establishes the direction
and scale of the benefit, so large10k should be rerun only if a large-system
runtime number is specifically needed for publication.

Apply `ddp/patches/no_copy_update_rule.patch` after both
`dynamic_network_scaling.patch` and `factor_bound_sensitivities.patch`.
Raw allocation rows are in `no_copy_update_rule_memory.csv`.
Full runtime rows are in `no_copy_update_rule_runtime.csv`, with complete
traces stored as `no_copy_updates_*_trace.csv`.
