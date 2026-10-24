# Historical tADMM plotting tools

These scripts plot saved tADMM convergence, trajectories and scaling results.
They are not SPAR-DDP entry points. They previously lived in root `plots/`.
Run from the repository root with the tADMM environment, for example:

```bash
julia --project=problems/mpopf/tadmm problems/mpopf/tadmm/plots/plot_convergence_from_csv.jl <csv> <output>
```

Several scripts retain study-specific inputs and sibling-paper output paths.
Inspect them before execution. Historical results are preserved for reference;
retention does not establish that they meet the current matched-comparison
protocol. Fresh tADMM testing/research is deferred.
