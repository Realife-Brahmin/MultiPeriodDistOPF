# First-Time Installation Guide for MultiPeriodDistOPF

## Prerequisites
1. Install Julia from the official website: https://julialang.org/downloads/
   - Download the current stable release (e.g., Julia 1.11.x) for your operating system.
   - Run the installer and finish setup.

2. Install Julia Extension in VS Code
   - Open VS Code.
   - Go to the Extensions Marketplace (left sidebar, square icon).
   - Search for "Julia" and install the official extension by julialang.

---

## 1. Open the Project in VS Code
1. Clone the repository from GitHub (or unzip if shared as a folder).
   ```bash
   git clone <repo-url>
   ```
   or just copy the folder.
2. In VS Code, click File → Open Folder… and select the MultiPeriodDistOPF folder.

---

## 2. Activate and Install Dependencies
From the terminal in VS Code, start Julia in project mode:

```bash
julia --project=.
```

Then in the Julia REPL (prompt looks like `julia>`), run:

```julia
using Pkg
Pkg.instantiate()
```

This downloads and installs all packages listed in Project.toml and Manifest.toml.
Only needs to be done once (unless dependencies change later).

---

## 3. Automatic Environment Activation Inside Scripts
To make running scripts easier, environment activation can be done inside the script itself.  
For example, at the top of `main.jl` (which lives in the repo root), add:

```julia
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
```

- `Pkg.activate(@__DIR__)` ensures the environment in the repo root is used.  
- `Pkg.instantiate()` ensures dependencies are installed if missing.  
- This allows running the script directly with:
  ```bash
  julia main.jl
  ```
  without needing to manually activate the environment each time.

If the script is located in a subfolder (e.g., `src/`), use:
```julia
Pkg.activate(joinpath(@__DIR__, ".."))
```

---

## 4. Running the Project
Run any of the scripts in the repo, for example:

```julia
include("main.jl")
```
or
```julia
include("admm_temporal_copper_plate.jl")
```

---

## 5. Restarting Later
Every time Julia restarts, either:
- Start Julia with the project environment:
  ```bash
  julia --project=.
  ```
- Or just run a script that has the activation snippet included at the top.

---

## 6. Tips
- If VS Code asks to choose a Julia environment, select the MultiPeriodDistOPF folder.
- If issues occur, delete Manifest.toml and run `Pkg.instantiate()` again.
- To keep all dependencies inside this repo instead of the global Julia folder, set `JULIA_DEPOT_PATH=$PWD/.julia_depot` before running Julia.
