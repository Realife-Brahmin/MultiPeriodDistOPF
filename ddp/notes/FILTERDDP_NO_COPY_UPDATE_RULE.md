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

Apply `ddp/patches/no_copy_update_rule.patch` after both
`dynamic_network_scaling.patch` and `factor_bound_sensitivities.patch`.
Raw allocation rows are in `no_copy_update_rule_memory.csv`.
