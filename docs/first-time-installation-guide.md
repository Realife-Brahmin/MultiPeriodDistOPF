# Getting started

**Just exploring the networks?** Start with [test networks](../test_networks/README.md)
and open a deck in [OpenDSS](https://sourceforge.net/projects/electricdss/).
Julia and the optimization environments are unnecessary for that route.

**Running optimization code?** Install Julia and the packages for your chosen
project below. SPAR-DDP requires no proprietary solver license. Its environment
includes OpenDSS export and validation dependencies; the standard
[Ipopt.jl installation](https://jump.dev/JuMP.jl/stable/packages/Ipopt/) provides
MUMPS for centralized reference solves. HSL and Gurobi are optional additions
for making centralized comparisons more competitive.

## Install Julia

Use [Juliaup, the official Julia installer](https://julialang.org/downloads/).
On Windows, run this in PowerShell:

```powershell
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

On Linux or macOS, the official installer is:

```bash
curl -fsSL https://install.julialang.org | sh
```

Review the installer and follow your organization's software policy. Open a
new terminal after installation. The helpers support Julia 1.12.4 or newer
within the 1.12 series; 1.12.7 has been checked for local setup. Install that
version without changing your global default:

```text
juliaup add 1.12.7
```

Git is convenient for obtaining the repository; downloading and extracting
the [repository ZIP](https://github.com/Realife-Brahmin/MultiPeriodDistOPF/archive/refs/heads/master.zip)
also works. VS Code and its Julia extension are optional.

```bash
git clone --depth 1 https://github.com/Realife-Brahmin/MultiPeriodDistOPF.git
cd MultiPeriodDistOPF
```

## Install project packages

From the repository directory, this command works in PowerShell, Bash, and
other terminals with Juliaup available:

```text
julia +1.12.7 --startup-file=no tools/setup.jl spar-ddp
```

Or use the wrapper for your terminal:

```powershell
.\tools\setup.ps1 -Project spar-ddp -JuliaChannel 1.12.7
```

```bash
MPOPF_JULIA_CHANNEL=1.12.7 bash tools/setup.sh spar-ddp
```

Replace `spar-ddp` with `tadmm` or `msopf` to install that project's packages.
Each environment is beside its implementation under `problems/`; there is no
root Julia environment. The helpers instantiate the selected environment,
preserve existing dependency pins, and do not launch an optimization run.
Internet access and disk space are needed for packages and binary artifacts.
First use may spend several minutes compiling packages; add `--precompile`
to the Julia/Bash command or `-Precompile` to the PowerShell wrapper to do
that during setup.

To check Julia's version and the project path before installing packages,
add `--check-only` (PowerShell: `-CheckOnly`). With a standalone Julia binary,
omit `+1.12.7`; wrapper overrides are `-JuliaExecutable` and `MPOPF_JULIA`.
If an institutional PowerShell policy blocks scripts, use the Julia command
directly or consult your administrator.

Next choose a driver and its inputs in [SPAR-DDP](../problems/mpopf/spar_ddp/README.md),
[tADMM](../problems/mpopf/tadmm/README.md), or [MSOPF](../problems/msopf/README.md).
The small IEEE123 examples and network exports are the best starting points
on a workstation. Large networks and long horizons can require substantial
RAM, worker/thread configuration, and cluster scheduling; the largest recorded
experiments used Kamiak. Installation alone does not reproduce those timings.

## Optional centralized solvers

- **HSL MA57 / MA97:** obtain the appropriate software and license through
  [STFC's Coin-HSL page](https://licences.stfc.ac.uk/product/coin-hsl).
  Follow the [Ipopt.jl HSL instructions](https://jump.dev/JuMP.jl/stable/packages/Ipopt/#HSL)
  for Julia integration, including `hsllib` and `linear_solver`. Choose a library
  matching your operating system and architecture. Building from source may
  require a Fortran compiler and supporting numerical libraries; follow the
  provider's build instructions.
- **Gurobi:** use the [download center](https://www.gurobi.com/product/download-center)
  and [license instructions](https://support.gurobi.com/hc/en-us/articles/12684663118993-How-do-I-obtain-a-Gurobi-license).
  [Academic options](https://www.gurobi.com/academics) depend on eligibility.
  Follow [Gurobi.jl's installation guide](https://jump.dev/JuMP.jl/stable/packages/Gurobi/)
  for the Julia interface; it can also obtain the solver binary. A usable license
  still needs activation.

MA57 and MA97 are linear solvers used inside Ipopt; Gurobi is a separate
optimizer supported by relevant centralized formulations. None is required
for the default SPAR-DDP path or an Ipopt/MUMPS reference installation.
The setup helpers do not acquire or activate vendor licenses.

## Prefer an AI assistant?

An assistant with terminal access can usually help install Julia, obtain the
repository, run these helpers, configure library paths, and diagnose package
errors. Give it this guide and the README for your chosen project. Ask it to
start with a small example appropriate to your machine.

You still need to handle vendor accounts, license requests, institutional
verification, and any required installation permissions. Supply approved
library locations after installation; keep license keys and credentials out
of chat messages and Git. An assistant cannot grant software access or make
large experiments fit hardware that lacks the required resources.

Optional reference-paper downloads are available through
`bash tools/fetch_resources.sh`; these are unrelated to package installation.
On Windows, a short checkout path can help avoid path-length limitations.
