# Blocked FilterDDP value-update right-hand sides

## Purpose

The factor-backed policy mode removes the dense `beta` and `omega` maps from
persistent horizon storage, but the backward pass still used one dense
`(nu + nc) x (nx + 1)` right-hand side.  On large10k this workspace occupied
`755.343 MiB`.

The optional blocked mode factorizes the same KKT matrix once, solves the
affine column, and then solves the state-sensitivity columns in blocks.  Each
solved block immediately contributes its rows to

```text
Vxx = C + beta' * B + omega' * cx
Vx  = lx + beta' * Qu + omega' * c + fx' * Vx_next + cx' * phi.
```

The block is then overwritten.  No sensitivity columns are omitted or
approximated.  The mode requires the exact factor-backed active-row network
path and is enabled with:

```text
FILTERDDP_FACTOR_BACKED_POLICY=1
FILTERDDP_BLOCKED_VALUE_RHS=1
FILTERDDP_VALUE_BLOCK_WIDTH=128
```

## Exactness gates

- IEEE123 `T = 3`, width 16: the complete 49-row trace is byte-for-byte
  identical to the dense-map baseline and converges in 48 iterations.
- IEEE2522 `T = 3`, width 128: the complete 57-row trace is byte-for-byte
  identical to the stored baseline and converges in 56 iterations.  The solve
  took `108.459 s`; this is recorded as functional evidence rather than a
  controlled timing comparison.

## Large10k memory profile

For one bounded large10k `T = 3` iteration with width 128:

- dense RHS workspace: `755.343 MiB -> 94.695 MiB` (87.46% reduction);
- retained factor-backed policy: `34.888 MiB/stage`, unchanged;
- first-stage maximum-RSS increase: `1589.223 MiB -> 242.652 MiB` relative to
  factor-backed policy without blocking (84.73% reduction);
- KKT assembly allocation: `863.783 MiB -> 188.601 MiB` at the first stage;
- bounded solve time: `28.109 s -> 37.179 s` (32.27% slower).

Cumulative allocation during all block solves remains high (`787.943 MiB` on
the warm stages), because every sensitivity value is still computed.  The
gain is peak workspace, not arithmetic complexity.

This is therefore a minimum-memory mode for horizons or systems that would
otherwise exhaust RAM.  The faster factor-backed-only and dense-map modes
remain preferable when memory is available.  A controlled runtime comparison
can be repeated when the machine is idle, but cannot change the exact storage
result.

## Idle-machine runtime comparison

A sequential IEEE2522 `T = 12` comparison on an otherwise idle machine used
the same strict `1e-7` settings for all modes.  Dense policies took
`523.357 s`, factor-backed policies took `474.091 s`, and factor-backed plus
width-128 blocking took `516.285 s`.  All three complete traces have the same
SHA-256 hash and therefore the same 79-iteration trajectory.  Blocking was
8.90% slower than factor-backed alone, but 1.35% faster than dense in this run.
It remains the minimum-memory mode; the test shows its memory reduction need
not imply a large full-solve penalty.
