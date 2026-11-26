# The driver with every period in a process of its own (procs_pass.jl).
# Usage: as ieee123c_filterddp.jl. FILTERDDP_PROCS_WORKERS: the number of worker
# processes (default: one per period); each runs one Julia thread and one BLAS
# thread. The coordinator is one more process.
using Distributed
const PS = @__DIR__
const PROCS_DRIVER = joinpath(PS, "ieee123c_filterddp.jl")
const PROCS_PASS = joinpath(PS, "procs_pass.jl")
let T = parse(Int, ARGS[2])
    W = min(T, parse(Int, get(ENV, "FILTERDDP_PROCS_WORKERS", string(T))))
    t0 = time()
    addprocs(W; exeflags=`--startup-file=no --project=$(Base.active_project()) --threads=1`,
             env=["OPENBLAS_NUM_THREADS" => "1", "JULIA_NUM_THREADS" => "1"])
    println("PROCS workers started: ", W, " in ", round(time() - t0; digits=1), " s")
end
@everywhere begin
    include($PROCS_DRIVER)                # defines main(); does not run it
    include($PROCS_PASS)
end
let system = ARGS[1], T = parse(Int, ARGS[2])
    ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
    data = deserialize(joinpath(REPO, "problems", "mpopf", "spar_ddp", "results", "network_filterddp", "network_data_$(system)_T$(T)$(ptag).jls"))
    procs_install(system, T, length(data[:Bset]))
end
main()
rmprocs(workers())
