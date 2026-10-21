# SPAR-DDP Julia environment

Use Julia 1.12 for the active DDP4OPF solver and centralized JuMP/Ipopt
references. Run from the repository root:

```bash
julia --project=envs/spar_ddp -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=envs/spar_ddp <script>
```

[Project.toml](Project.toml) declares dependencies and points DDP4OPF to the
local source at `../../spar_ddp/DDP4OPF.jl`. This folder was previously
`envs/ddp2026`; the rename changes paths, not package identity or algorithms.
No Manifest is committed here, so this README does not promise fixed resolved
package versions. For published timing reproduction use the recorded campaign
metadata and environment at its evidence commit.

OpenDSS instance export uses `envs/tadmm`, whose tracked Manifest is retained.
HSL libraries and licenses are external to this environment; see the recorded
configuration for MA57/MA97 runs. Gurobi is optional and not declared here.

The retired 2025 first-order implementation is in
[archives/first_order_ddp](../../archives/first_order_ddp/README.md); it is
not a second supported environment. Pkg strips Project.toml comments during
resolve/add operations, so keep explanatory notes in this README.
