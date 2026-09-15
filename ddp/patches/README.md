# Development history of the DDP4OPF fork -- NOT a build recipe

These patches record, in the order they were developed, the changes made to
FilterDDP.jl (upstream commit `513a104`) while adapting it to multi-period OPF.
Each has a companion note in `ddp/notes/`.

**They do not rebuild the solver.** Checked 2026-09-15: applied to a clean
`513a104` in the order the old README documented, three of the nine listed
patches failed to apply (`in_place_kkt_rhs`, `reuse_kkt_rhs_workspace`,
`active_B_rows`), and the resulting tree differed from the solver that had
actually produced the results in four source files. Several patches here were
also never part of that recipe at all.

The solver is now committed as source at [`ddp/DDP4OPF.jl`](../DDP4OPF.jl),
verified bit-identical to the working clone the results came from. Read these
patches to understand *why* a given change was made; read `ddp/DDP4OPF.jl/src`
to see *what the solver is*.
