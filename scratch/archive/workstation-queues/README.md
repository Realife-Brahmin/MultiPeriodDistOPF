# Superseded workstation queues

Seven queue/monitor scripts were moved here unchanged from `scratch/` during
the October 23 merge cleanup. They describe old timing-matrix/Hessian queues,
hardcoded checkouts, private Julia locations and unattended run behavior.
They are not portable entry points or jobs that should now be resumed.

Their bytes are preserved for history, but moving them intentionally makes
their old `$PSScriptRoot` assumptions unsuitable for execution here. For exact
historical behavior use the recorded pre-cleanup commit/path, with explicit
review of side effects and inputs. Do not fix their paths and run them blindly.
Reusable extractors, individual case/regression helpers and drawing tools
remain in `scratch/`; supported SPAR-DDP drivers are documented in
[the workstream guide](../../../spar_ddp/README.md).
