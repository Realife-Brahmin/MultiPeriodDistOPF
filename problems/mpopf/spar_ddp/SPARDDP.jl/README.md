# SPARDDP.jl

Filter line-search **differential dynamic programming** for multi-period optimal
power flow on distribution networks.

Derived from [FilterDDP.jl](https://github.com/mingu6/FilterDDP.jl) by Mingda Xu
(MIT). See [NOTICE.md](NOTICE.md) for provenance, attribution and the list of
changes, and [LICENSE](LICENSE) for terms. Please cite the original method papers
(arXiv 2504.08278, 2606.01487) when using it.

## Use

It is developed as a local package inside the `problems/mpopf/spar_ddp` environment:

```bash
julia --startup-file=no --project=problems/mpopf/spar_ddp <script>
```

```julia
using SPARDDP
```

The public API is upstream's: `Objective`, `Dynamics`, `EqualityConstraints`,
`ControlLimits`, `build_ocp`, `Solver`, `Options`, `solve!`, `get_trajectory`.
Stage structs can be built directly from hand-written callables, bypassing the
`Symbolics`-based convenience constructors -- which is what the network and
reduced-space drivers under `problems/mpopf/spar_ddp/examples/power_system/` do.
