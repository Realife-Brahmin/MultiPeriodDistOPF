# Maintained repository utilities

`fetch_resources.sh` downloads missing papers listed in resources/MANIFEST.txt
files. Existing files are left untouched. Run `bash tools/fetch_resources.sh`
from the repository root. Historical experiment launchers live in scratch/.

`fetch_test_systems.sh` retrieves public source decks; `summarize_test_systems.jl`
checks those decks with OpenDSS. Read [their provenance](../test_systems/README.md)
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
