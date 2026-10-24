# Distribution-system optimization problems

Each project keeps its Julia environment beside its implementation, drivers,
studies and results. Choose a problem first, then an algorithm.

- **Multi-period OPF (MPOPF)**
  - [SPAR-DDP](mpopf/spar_ddp/README.md): sparse parallel differential dynamic
    programming; instantiate `--project=problems/mpopf/spar_ddp`.
  - [tADMM](mpopf/tadmm/README.md): temporal ADMM and a monolithic baseline;
    instantiate `--project=problems/mpopf/tadmm`.
- **[Multi-source OPF (MSOPF)](msopf/README.md)**: explicit source angles and
  multi-period scheduling; instantiate `--project=problems/msopf`.
  Its optional SDP environment is under `environments/sdp/`.

The shared [OpenDSS networks and parser](../test_networks/README.md) are outside
the algorithms. SPAR-DDP includes its own OpenDSS dependency: no tADMM environment
is needed for its export or validation workflow. Models retain their separate
assumptions, including different terminal-SOC formulations.

See the [layout migration guide](../docs/problem-layout-migration.md) for old
paths and the distinction between runnable sources and frozen evidence.
