---
name: Known Issues
description: Ongoing bugs and technical debt in the project
type: project
---

- `git reset --hard` / `git revert` fails due to long filenames in `archives/` directory
- GR backend hangs on `savefig` in headless Windows — fix: set `ENV["GKSwstype"] = "png"` before `using Plots`, or pass `GKSwstype=png` as env var when running Julia
- OpenDSSDirect precompilation warnings (harmless, Julia 1.12 compat)
- `.jls` files (Serialization) can't be deserialized outside the full tadmm environment (MathOptInterface dependency)
- Serialization-based save/load (commit `fb732c69`) is overengineered — user prefers simple txt/csv
