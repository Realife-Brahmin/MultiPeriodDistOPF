# `ddp2026` — the single environment for all DDP work

Created 2026-08-07. This is the one environment for everything DDP-related:
the FilterDDP evaluation, the centralized JuMP/Ipopt reference, and any further
formulation added later.

```powershell
julia --startup-file=no --project=envs/ddp2026 <script>
```

## What it holds

| | |
|---|---|
| `FilterDDP` v0.6.0 | dev'd from the pinned clone at `ddp/external/FilterDDP.jl` (`513a104`) |
| `JuMP` v1.31.1 + `Ipopt` v1.15.0 | the centralized reference formulation |
| `StaticArrays` v1.9.18 | required by the FilterDDP stage models |

FilterDDP and JuMP coexist without conflict — verified by a clean resolve, which
settles `OrderedCollections` on 1.8.2.

## What it supersedes, and why

**`envs/ddp/Project.toml`** (2025-11) — the environment for the user's own
first-order Differential Dynamic Programming scheme. Unused for ~9 months and not trusted to
resolve as it once did. Its `Plots`, `Crayons`, `LaTeXStrings`, `Parameters` and
`Revise` dependencies existed because verification then meant a human reading
formatted output; the checks here are numerical and self-asserting, so none of
them are carried over. The old directory is left in place as a **read-only
reference** for what that algorithm did — it is not expected to be run again.

**`ddp/env`** (2026-08-03) — FilterDDP-only. **Removed** 2026-08-07. Stages 5 and
6 were re-run here first and reproduce its recorded logs exactly, and its
`Manifest.toml` was never tracked, so it preserved no version information this
environment lacks. Recoverable from git history.

`Gurobi` is deliberately omitted: it needs a licensed local install and nothing
in this stack requires it. Add it back only if the user's own DDP code is run
from here.

## Note on editing `Project.toml`

Pkg rewrites this file on every `add`/`resolve` and **strips all comments**, so
rationale lives in this README instead.
