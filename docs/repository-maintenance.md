# Keeping master a usable research baseline

## Public documentation and coding-agent guidance

The root README is for visitors: what the repository offers, how to get an
OpenDSS network or try a method, and how to cite it. Campaign status belongs
with the study; workflow safeguards belong here or in developer instructions,
not in the main introduction.

Full shared agent instructions live in `.agents/AGENTS.md`. The root
`AGENTS.md` is only a short pointer because Codex discovers guidance along the
root-to-working-directory path, not by recursively loading arbitrary hidden
folders. Claude's project entry lives in `.claude/CLAUDE.md` and imports the
shared instructions; local sessions/worktrees remain ignored. Workstream
Claude notes remain scoped to their own folders. No user-global settings changed.
See [Codex discovery](https://learn.chatgpt.com/docs/agent-configuration/agents-md)
and [Claude instruction locations](https://code.claude.com/docs/en/memory).

The network catalogue and physical adapted decks live together in `test_networks/`;
published source decks are preserved in `test_networks/sources/`.
Standalone network ZIPs use `test_networks/<system>/` and include shared
profiles, license and citation. The problem-first migration renamed rawData
without changing electrical inputs; historical harnesses use their recorded commits.
The repository citation is a request for scholarly credit, not a license change.

The repository hosts SPAR-DDP, tADMM and MSOPF. A project may retain useful
experimental work without making it a supported feature. The directory layout
must distinguish implementation, evidence and unfinished local exploration.

## Placement and retention

- Maintained solver/model code belongs to its workstream. Environments declare
  dependencies; do not create an unrelated root Julia package.
- A new bounded study belongs in `problems/mpopf/spar_ddp/experiments/<study>/` (or its owning
  project's equivalent), with a question, scope, input assumptions, source,
  checks and findings/status. Keep failures that explain scientific decisions.
- Keep benchmark evidence immutable and indexed. Do not move hash-pinned
  harnesses or erase negative results merely to reduce file count. Historical
  paths are interpreted at their evidence commits, not rewritten as current.
- Shared decks/provenance stay in `test_networks/` and `test_networks/sources/`; figures and
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
compact-campaign status; corrected the root guide to use `problems/msopf`.
Seven superseded workstation queues moved unchanged below
`scratch/archive/workstation-queues/`; individual reusable helpers remain.
Git transfer bundles are ignored, not deleted from local checkouts.
No frozen DDP evidence, numerical source, feeder model, benchmark row or paper
was removed or changed by this organization pass.
Obsolete root-package CI/CompatHelper/TagBot templates were replaced by
navigation and parse-only checks. Removed templates remain recoverable in Git
history; these checks deliberately make no scientific-validation claim.
