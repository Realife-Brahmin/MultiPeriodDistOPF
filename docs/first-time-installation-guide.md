# Installing MultiPeriodDistOPF

Install Julia 1.12 from https://julialang.org/downloads/. VS Code and its
Julia extension are optional; open the repository folder directly. Personal
`.code-workspace` files are local and ignored.

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

Use the environment for your workstream; there is no root Project.toml
or universal main.jl. For SPAR-DDP, instantiate its colocated environment,
which includes OpenDSS export and validation:

```bash
julia --project=problems/mpopf/spar_ddp -e 'using Pkg; Pkg.instantiate()'
```

Use Bash for these commands on Windows. The tracked `problems/mpopf/tadmm/Manifest.toml`
preserves its dependency versions; do not update or remove it during setup.
The SPAR-DDP environment declares its local SPARDDP source in Project.toml.

See [SPAR-DDP](../problems/mpopf/spar_ddp/README.md), [tADMM](../problems/mpopf/tadmm/README.md) or
[MSOPF](../problems/msopf/README.md) for actual drivers and required inputs.
SPAR-DDP requires no proprietary solver license. The standard Ipopt.jl reference
works with bundled MUMPS. Optional HSL MA57/MA97 and Gurobi were used for more
competitive centralized comparisons; they require extra setup/access but are
not prerequisites. A default Ipopt/MUMPS install is not the recorded Kamiak
MA57/MA97 timing configuration.

Optional reference downloads: `bash tools/fetch_resources.sh`. This reads
resources/MANIFEST.txt files and leaves existing files untouched.

Windows paths can be long. Use a short checkout location if necessary.
No optimization jobs need to be launched merely to install or inspect code.
