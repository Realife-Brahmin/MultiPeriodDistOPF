# Keeping master a usable research baseline

The repository hosts SPAR-DDP, tADMM and MSOPF. A project may retain useful
experimental work without making it a supported feature. The directory layout
must distinguish implementation, evidence and unfinished local exploration.

## Placement and retention

- Maintained solver/model code belongs to its workstream. Environments declare
  dependencies; do not create an unrelated root Julia package.
- A new bounded study belongs in `spar_ddp/experiments/<study>/` (or its owning
  project's equivalent), with a question, scope, input assumptions, source,
  checks and findings/status. Keep failures that explain scientific decisions.
- Keep benchmark evidence immutable and indexed. Do not move hash-pinned
  harnesses or erase negative results merely to reduce file count. Historical
  paths are interpreted at their evidence commits, not rewritten as current.
- Shared decks/provenance stay in `rawData/` and `test_systems/`; figures and
  maintained general tools have their own homes. Do not delete another
  workstream's results or source models as part of a DDP merge.
- Temporary scripts go in `scratch/` until promoted, archived or removed after
  a reference check. Superseded machine-specific queues are archived separately
  and are not supported launchers. Private libraries, package depots, binary
  captures, editor files and Git transfer bundles remain local/ignored.

## Before merging

1. Bring in current master; preserve concurrent workstream changes.
2. Document what is maintained, opt-in, historical, rejected or unresolved.
3. Check navigation, syntax and relevant small tests. Do not launch expensive
   simulations simply to tidy documentation or folders.
4. Preserve source/input manifests and failed attempts. Check a rename's
   references; do not relocate shared model helpers without updating consumers.
5. Merge through the workstream PR; do not change scientific defaults, published
   rows or the separate paper as an incidental cleanup.

## October 23 cleanup scope

The battery-count/value-function branch was reconciled with current master,
including MSOPF. Added workstream/study/evidence navigation and the final
compact-campaign status; corrected the root guide to use `envs/msopf`.
Seven superseded workstation queues moved unchanged below
`scratch/archive/workstation-queues/`; individual reusable helpers remain.
Git transfer bundles are ignored, not deleted from local checkouts.
No frozen DDP evidence, numerical source, feeder model, benchmark row or paper
was removed or changed by this organization pass.
Obsolete root-package CI/CompatHelper/TagBot templates were replaced by
navigation and parse-only checks. Removed templates remain recoverable in Git
history; these checks deliberately make no scientific-validation claim.
