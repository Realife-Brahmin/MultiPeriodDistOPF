# Maintained repository utilities

`check_layout.py` verifies relocated evidence/data/solver identities and key
entry-point paths; it needs the baseline Git history (see the
[migration guide](../docs/problem-layout-migration.md)).
`smoke_project_layout.jl` is an optional small runtime import/parser/export
check in the SPAR-DDP environment; it does not run an optimizer and is not a
benchmark. It writes its export into a temporary directory.

`package_test_network.py` creates a standalone OpenDSS ZIP containing one
adapted deck, its shared profiles, license and citation metadata. It does not
change inputs or overwrite existing archives. See the
[network catalogue](../test_networks/README.md) for a command and OpenDSS usage.
Run `python -m unittest discover -s tools -p 'test_*.py' -v` for packaging checks.

`fetch_resources.sh` downloads missing papers listed in resources/MANIFEST.txt
files. Existing files are left untouched. Run `bash tools/fetch_resources.sh`
from the repository root. Historical experiment launchers live in scratch/.

`fetch_test_systems.sh` retrieves public source decks; `summarize_test_systems.jl`
checks those decks with OpenDSS. Read [their provenance](../test_networks/sources/README.md)
before running either; the latter performs solves and is not a cleanup check.

Lightweight checks, run from any checkout:

```bash
python tools/check_repository.py
julia --startup-file=no tools/check_julia_syntax.jl
```

These check curated maintained navigation and Python/Julia syntax only, without
installing packages or launching optimizations. They do not establish solver
correctness, feasibility or benchmark reproducibility. CI runs these checks
on Julia 1.12.4 instead of treating this multi-project repository as a root
Julia package. The obsolete root-package CompatHelper/TagBot templates were
retired; a future software release needs an explicitly scoped package workflow.
