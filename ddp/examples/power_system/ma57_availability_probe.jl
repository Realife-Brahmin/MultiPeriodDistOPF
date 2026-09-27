# Is MA57 (or any HSL solver) usable from Julia or Ipopt on this machine?
#
# Records, without installing anything:
#   * versions of every solver-related package in envs/ddp2026
#   * for each linear_solver value Ipopt knows, what Ipopt does when asked to use
#     it on a 2-variable QP (its own journal is kept, so the exact error survives)
#   * whether any HSL shared library is loadable by name
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/ma57_availability_probe.jl [out.txt]

using JuMP, Ipopt, Pkg, Printf
using Base.Libc.Libdl: find_library

# Ipopt holds its journal open until the problem is freed, hence the function
# scope before the parent process reads the file.
function try_solver(solver, journal)
    m = Model(Ipopt.Optimizer)
    set_attribute(m, "linear_solver", solver)
    set_attribute(m, "output_file", journal)
    set_attribute(m, "file_print_level", 5)
    set_attribute(m, "print_level", 0)
    @variable(m, x); @variable(m, y)
    @constraint(m, x + y == 1)
    @objective(m, Min, x^2 + y^2)
    optimize!(m)
    return string(termination_status(m)), raw_status(m)
end

if length(ARGS) >= 3 && ARGS[1] == "--one"
    status, raw = try_solver(ARGS[2], ARGS[3])
    println("PROBE_RESULT status=$status raw=$raw")
    exit(0)
end

out = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "..", "results", "kkt_ordering", "ma57_availability_julia.txt")
mkpath(dirname(out))
io = open(out, "w")
say(args...) = (println(io, args...); println(args...))

say("MA57 / HSL availability probe, ", gethostname(), ", Julia ", VERSION)
say("\n== packages (envs/ddp2026) ==")
for (_, d) in sort(collect(Pkg.dependencies()); by = p -> p[2].name)
    d.name in ("Ipopt", "Ipopt_jll", "MUMPS", "MUMPS_seq_jll", "MUMPS_jll", "METIS_jll",
               "PARMETIS_jll", "SCOTCH_jll", "SPRAL_jll", "SuiteSparse_jll", "OpenBLAS32_jll",
               "JuMP", "DDP4OPF", "HSL", "HSL_jll") || continue
    say(rpad(d.name, 16), d.version)
end
for name in ("HSL", "HSL_jll")
    say(rpad(name, 16), Base.find_package(name) === nothing ? "not installed" : "installed")
end

say("\n== HSL shared libraries loadable by name ==")
for lib in ("libhsl", "libcoinhsl", "libma57", "libma27", "libhsl_ma57", "libhsl_ma97", "libmwma57")
    found = find_library([lib])
    say(rpad(lib, 14), isempty(found) ? "not found" : found)
end

say("\n== Ipopt ", Ipopt.Ipopt_jll.libipopt, " ==")

# Each solver runs in its own process: in one process, the fifth consecutive
# failed HSL load crashed Julia (access violation) instead of returning an error.
# SPRAL is retried with the OpenMP settings its documentation requires.
cases = [(s, Pair{String,String}[]) for s in
         ("mumps", "ma27", "ma57", "ma77", "ma86", "ma97", "pardiso", "pardisomkl", "spral", "wsmp")]
push!(cases, ("spral", ["OMP_CANCELLATION" => "TRUE", "OMP_PROC_BIND" => "TRUE"]))
for (solver, envs) in cases
    journal = tempname() * ".txt"
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $(@__FILE__) --one $solver $journal`
    buf = IOBuffer()
    proc = run(pipeline(ignorestatus(addenv(cmd, envs...)); stdout=buf, stderr=buf))
    text = String(take!(buf))
    m = match(r"PROBE_RESULT status=(\S+) raw=(\S+)", text)
    status, raw = m === nothing ? ("process exit $(proc.exitcode)", "") : (m[1], m[2])
    lines = isfile(journal) ? readlines(journal) : String[]
    msg = filter(l -> occursin(r"(?i)hsl|ma\d\d|pardiso|wsmp|spral|linear solver|not available|obtain|error|exception", l), lines)
    label = isempty(envs) ? solver : solver * " (" * join(first.(envs), ",") * ")"
    say(@sprintf("%-40s status=%-20s raw=%s", label, status, raw))
    if m === nothing
        err = filter(l -> occursin(r"(?i)error|exception|invalid|violation", l), split(text, '\n'))
        for l in first(err, 3)
            say("            | ", strip(l))
        end
    end
    for l in first(unique(msg), 4)
        say("            | ", strip(l))
    end
end
close(io)
println("wrote ", out)
