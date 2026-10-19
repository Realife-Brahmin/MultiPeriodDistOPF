# Installing MultiPeriodDistOPF

Install Julia 1.12 from https://julialang.org/downloads/. VS Code and its
Julia extension are optional; open the repository folder directly. Personal
`.code-workspace` files are local and ignored.

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

Use the environment for your workstream; there is no root Project.toml
or universal main.jl. For SPAR-DDP, instantiate both the solver environment
and the OpenDSS export environment:

```bash
julia --project=envs/spar_ddp -e 'using Pkg; Pkg.instantiate()'
julia --project=envs/tadmm -e 'using Pkg; Pkg.instantiate()'
```

Use Bash for these commands on Windows. The tracked `envs/tadmm/Manifest.toml`
preserves its dependency versions; do not update or remove it during setup.
The SPAR-DDP environment declares its local DDP4OPF source in Project.toml.

See [SPAR-DDP](../spar_ddp/README.md), [tADMM](../envs/tadmm/README.md) or
[MSOPF](../envs/multi_poi/README.md) for actual drivers and required inputs.
Gurobi and HSL are optional separately licensed backends; access/configuration
is not supplied by cloning this repository. A default Ipopt install does not
by itself recreate the recorded Kamiak MA57/MA97 benchmark.

Optional reference downloads: `bash tools/fetch_resources.sh`. This reads
resources/MANIFEST.txt files and leaves existing files untouched.

Windows paths can be long. Use a short checkout location if necessary.
No optimization jobs need to be launched merely to install or inspect code.
