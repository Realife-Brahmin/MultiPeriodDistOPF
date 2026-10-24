# Install only the selected project's declared dependencies. No solves or licenses.
using Pkg
const SETUP_PROJECTS = Dict(
    "spar-ddp" => "problems/mpopf/spar_ddp",
    "tadmm" => "problems/mpopf/tadmm",
    "msopf" => "problems/msopf",
)

function setup_main(args)
    if any(x -> x in ("--help", "-h"), args)
        println("Usage: julia tools/setup.jl [spar-ddp|tadmm|msopf] [--check-only] [--precompile]")
        println("Default: spar-ddp. Does not install Julia or acquire/activate vendor licenses.")
        return
    end
    flags = filter(x -> startswith(x, "--"), args)
    all(x -> x in ("--check-only", "--precompile"), flags) || error("Unknown setup flag; use --help")
    names = filter(x -> !startswith(x, "--"), args)
    length(names) <= 1 || error("Choose one project; use --help")
    name = isempty(names) ? "spar-ddp" : only(names)
    haskey(SETUP_PROJECTS, name) || error("Unknown project '$name'; choose spar-ddp, tadmm or msopf")
    v"1.12.4" <= VERSION < v"1.13" || error("Use Julia 1.12.4+ in the 1.12 series; see the setup guide")
    root = normpath(joinpath(@__DIR__, ".."))
    project = joinpath(root, SETUP_PROJECTS[name])
    isfile(joinpath(project, "Project.toml")) || error("Project.toml missing: extract or clone the complete repository")
    println("Julia version: ", VERSION)
    println("Selected project: ", name)
    println("Environment: ", project)
    if "--check-only" in flags
        println("SETUP_CHECK_PASS (paths/version only; no packages installed or licenses checked)")
        return
    end
    Pkg.activate(project)
    Pkg.instantiate(; allow_autoprecomp=false)
    "--precompile" in flags && Pkg.precompile()
    println("SETUP_COMPLETE: declared packages installed; no optimization or benchmark run.")
    println("First use may compile packages. Optional HSL/Gurobi need separate activation/configuration.")
end
setup_main(ARGS)
