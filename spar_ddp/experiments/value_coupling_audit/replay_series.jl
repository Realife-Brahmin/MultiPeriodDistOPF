# Multiple independent frozen points in one process to amortize compilation.
# Point order does NOT evolve an optimizer; every point reloads its own state.
length(ARGS)==5 || error("replay_series.jl CAPTURE_ROOT SYSTEM T WINDOW OUTPUT_CASE")
include(joinpath(@__DIR__,"replay.jl"))
function series(args)
    cap,sys,ts,ws,out=args;T=parse(Int,ts);window=parse(Int,ws)
    base=joinpath(cap,"results","$(sys)_T$(T)_dense_capture")
    input=joinpath(cap,"source","spar_ddp/results/network_filterddp","network_data_$(sys)_T$(T)_periodic.jls")
    records=[split(line,',') for line in readlines(joinpath(base,"captures/coverage.csv"))[2:end]]
    @assert Set(first.(records))==Set(("early","middle","late","barrier_before","barrier_after"))
    sys=="ieee2522C_1ph" && sort!(records;by=r->r[1]=="late" ? 0 : 1)
    mkpath(out)
    for r in records
        tag,snapshot=r[1:2]
        println("SERIES_POINT_BEGIN system=$sys tag=$tag snapshot=$snapshot");flush(stdout)
        run_replay(joinpath(base,"captures",snapshot),input,joinpath(out,tag);window=window)
        GC.gc()
    end
    println("FROZEN_SERIES_COMPLETE system=$sys T=$T window=$window diagnostic_only=true")
end
series(ARGS)
