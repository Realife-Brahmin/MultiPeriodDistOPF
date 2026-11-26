# Problem-first layout migration

This migration is stacked on README PR176 and includes the merged MSOPF PR177.
The baseline immediately before moving files is commit
`0efa8015b8e7d4a59dee23c7173fd1494a2a0963`. No paper, solver equations,
benchmark rows or OpenDSS electrical inputs are changed.

## New homes

- `spar_ddp/` and `envs/spar_ddp/` → `problems/mpopf/spar_ddp/`.
  The environment's Project.toml is beside SPARDDP, drivers and evidence.
- `envs/tadmm/` and root `results/` → `problems/mpopf/tadmm/`.
- `envs/multi_poi/` and the current `envs/msopf/` environment →
  `problems/msopf/`. The superseded multi_poi Project.toml is retained in
  `archive/legacy-environment/`, not used as the main environment.
- `envs/msopf_sdp/` → `problems/msopf/environments/sdp/`.
- `rawData/` → `test_networks/` (the public catalogue already lives there).
- `test_systems/` → `test_networks/sources/`.

The duplicate rawData landing README was removed; the catalogue is the single
network landing page. Public source decks and adapted deck files remain intact.
No symlinks or machine-specific junctions are needed by a fresh clone.

## Running code

From the repository root, choose one project:

```bash
julia --project=problems/mpopf/spar_ddp -e 'using Pkg; Pkg.instantiate()'
julia --project=problems/mpopf/tadmm -e 'using Pkg; Pkg.instantiate()'
julia --project=problems/msopf -e 'using Pkg; Pkg.instantiate()'
```

These are alternatives, not a requirement to install every environment.
SPAR-DDP's exporter and OpenDSS validator both use its own environment.
The numerical parser is shared at `test_networks/parse_opendss.jl`; tADMM's
old include filename delegates to it without changing its API or parsing.
Environment consolidation is not a new set of matched benchmark versions:
use pinned campaign snapshots and their manifests for exact historical runs.

Local generated inputs/solutions have moved with their algorithm where present.
Re-export fresh inputs under the new tree as needed; do not overwrite frozen
inputs to impersonate an older campaign. Existing worktrees/jobs continue at
their original commits and are not migrated in place.

## Evidence and validation boundaries

Logs, CSVs, saved manifests, archives, figures, certificates and raw source decks
move without rewriting their recorded contents. Notes and frozen harnesses may
therefore mention old paths. Reproduce them at their recorded commit; those
commands are not promised to run against the reorganized current checkout.
Current public drivers, environment references and maintained navigation are
updated separately. Retired first-order code and archived workstation queues
remain historical, not supported entry points.

`python tools/check_layout.py --strict-migration` checks preserved Git blob
identities and unchanged maintained solver/parser source against the merged
master commit5742a5d1, and validates environment/include/root paths. The normal
CI invocation omits the strict flag so future deliberate solver development
is not prohibited; immutable evidence/data remain checked. Both require Git
history (not only a shallow snapshot). Navigation checks,
Python tests and Julia parsing provide additional structural coverage. These
checks are not scientific revalidation of the stored benchmark/certificate rows.
No cluster jobs or benchmark campaigns are required by this migration.

## Package name

The local solver package/module is now `SPARDDP` in `SPARDDP.jl/`, formerly
`DDP4OPF`. Its UUID, authors, API functions, upstream MIT notice and numerical
method are retained. Current drivers import `SPARDDP`; frozen study scripts
and serialized solver captures still use their recorded package name/commit.
The `FILTERDDP_*` environment variables and historical log formats are not
renamed. If upgrading a local checkout with an old generated Manifest, preserve
it and regenerate that local environment; fresh clones have the new source map.

The strict audit accounts for the entry-file rename by allowing exactly
`module DDP4OPF` → `module SPARDDP`, not arbitrary source differences.

## Validation checkpoint

The relocation/namespace audit checks4,924 unchanged evidence/data/figure/solver
blobs plus the single exact module-declaration change,
19 literal repository roots and86 literal includes. Maintained navigation and
all16 adapted-network packaging cases pass. The Julia parse check covers215
files across the projects, shared networks and tools.

A bounded runtime smoke using the SPAR-DDP environment checks parser equivalence
against the baseline on IEEE123 T3 and exports a periodic instance with OpenDSS.
Its isolated exports do not overwrite campaign inputs. This is import/export
coverage, not a re-run of either algorithm's optimization benchmarks.
Fresh SPAR-DDP resolution is explicitly distinct from recorded campaign pins.
The runtime smoke passed with Julia1.12.7, OpenDSSDirect0.9.9 and Parameters0.12.3,
including SPARDDP, JuMP and Ipopt imports. The extended smoke also checks a
one-variable Ipopt/MUMPS solve without HSL/Gurobi; no MPOPF benchmark is run.

The old Inkscape and evidence bytes are unchanged. tADMM plotting defaults now
write to its own results/figures directory; explicitly set TADMM_FIGURE_DIR to
export elsewhere, rather than implicitly writing into a sibling paper checkout.
