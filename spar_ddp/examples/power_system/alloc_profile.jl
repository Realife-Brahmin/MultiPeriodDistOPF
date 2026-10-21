# Where does a FilterDDP solve allocate? (2026-10-11: the collector takes 17-22%
# of the timed solve on med2522 and large10k, 33% on ieee123 T=96.)
#
# Runs the driver's main() under Julia's allocation profiler and attributes
# every sampled allocation to the innermost frame in this repository's code.
# Untimed; set the environment as run_parsim.sh does (FILTERDDP_WARMUP=0).
#
#   ALLOC_SAMPLE_RATE=0.01 julia --project=envs/spar_ddp alloc_profile.jl <system> <T> solve
#
# Prints the sites with the most allocations and the most bytes. The model
# build and compilation are in the profile too; they show up under
# build_model and under Base/Core frames and are listed for what they are.

using Profile

include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))      # defines main(); does not run it

# With FILTERDDP_PAIRED_SWEEP set, the pass of paired_pass.jl is the one profiled.
if get(ENV, "FILTERDDP_PAIRED_SWEEP", "0") != "0"
    include(joinpath(@__DIR__, "paired_pass.jl"))
    let system = ARGS[1], T = parse(Int, ARGS[2])
        ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
        data = deserialize(joinpath(REPO, "spar_ddp", "results", "network_filterddp", "network_data_$(system)_T$(T)$(ptag).jls"))
        idx, _ = control_layout(data)
        install_paired_pass(data, idx)
    end
end

function alloc_site(st)
    for f in st
        file = String(f.file)
        (occursin("DDP4OPF", file) || occursin("power_system", file)) || continue
        occursin("alloc_profile", file) && continue
        return string(f.func, "  ", basename(file), ":", f.line)
    end
    # no frame of ours: name the innermost frames instead
    return "(outside) " * join((string(f.func, " ", basename(String(f.file)), ":", f.line) for f in st[1:min(3, length(st))]), " < ")
end

function alloc_report(rate)
    res = Profile.Allocs.fetch()
    count = Dict{String, Int}(); bytes = Dict{String, Int}()
    for a in res.allocs
        k = alloc_site(a.stacktrace) * "  [" * string(a.type) * "]"
        count[k] = get(count, k, 0) + 1
        bytes[k] = get(bytes, k, 0) + a.size
    end
    n = length(res.allocs); b = sum(a.size for a in res.allocs; init=0)
    println("ALLOC_PROFILE samples=", n, " sample_rate=", rate, " estimated allocations=", round(Int, n / rate),
            " estimated GiB=", round(b / rate / 2^30, digits=2))
    println("-- by number of allocations")
    for (k, v) in sort!(collect(count); by=last, rev=true)[1:min(30, length(count))]
        println(rpad(string(round(100v / n, digits=1), "%"), 8), rpad(string(round(100bytes[k] / b, digits=1), "% bytes"), 14), k)
    end
    println("-- by bytes")
    for (k, v) in sort!(collect(bytes); by=last, rev=true)[1:min(20, length(bytes))]
        println(rpad(string(round(100v / b, digits=1), "%"), 8), rpad(string(round(100count[k] / n, digits=1), "% count"), 14), k)
    end
end

let rate = parse(Float64, get(ENV, "ALLOC_SAMPLE_RATE", "0.01"))
    # Once without the profiler, so that the profiled run compiles nothing.
    get(ENV, "ALLOC_SKIP_FIRST", "0") == "0" && redirect_stdout(devnull) do
        main()
    end
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=rate main()
    alloc_report(rate)
end
