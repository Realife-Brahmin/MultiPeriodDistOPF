# DDP4OPF.jl

Filter line-search **differential dynamic programming** for multi-period optimal
power flow on distribution networks.

Derived from [FilterDDP.jl](https://github.com/mingu6/FilterDDP.jl) by Mingda Xu
(MIT). See [NOTICE.md](NOTICE.md) for provenance, attribution and the list of
changes, and [LICENSE](LICENSE) for terms. Please cite the original method papers
(arXiv 2504.08278, 2606.01487) when using it.

## Use

It is developed as a local package inside the `envs/ddp2026` environment:

```bash
julia --startup-file=no --project=envs/ddp2026 <script>
```

```julia
using DDP4OPF
```

The public API is upstream's: `Objective`, `Dynamics`, `EqualityConstraints`,
`ControlLimits`, `build_ocp`, `Solver`, `Options`, `solve!`, `get_trajectory`.
Stage structs can be built directly from hand-written callables, bypassing the
`Symbolics`-based convenience constructors -- which is what the network and
reduced-space drivers under `ddp/examples/power_system/` do.
