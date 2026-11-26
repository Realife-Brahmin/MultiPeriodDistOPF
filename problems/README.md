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

[Browse the test networks](../test_networks/README.md) or follow the
[setup guide](../docs/first-time-installation-guide.md) to get started.
