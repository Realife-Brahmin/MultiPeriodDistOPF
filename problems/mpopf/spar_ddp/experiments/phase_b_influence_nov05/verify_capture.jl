using SHA,TOML
@assert length(ARGS)==1
root=ARGS[1];caps=joinpath(root,"captures");log=read(joinpath(root,"solve.log"),String)
states=Dict{String,Any}()
for dir in readdir(caps;join=true)
    isdir(dir)||continue
    m=TOML.parsefile(joinpath(dir,"metadata.toml"))
    @assert m["stages"]==24&&m["full_stage_records"]
    @assert m["schedule"]=="original_balanced"&&m["kernel"]=="original_inplace"
    @assert length(m["tags"])==1
    tag=only(m["tags"]);@assert !haskey(states,tag)
    @assert occursin("INFLUENCE_CAPTURE tag=$tag ",log)
    @assert isfinite(m["regularization_before_pass"])&&isfinite(m["regularization_reported_after_pass"])
    @assert bytes2hex(sha256(read(joinpath(dir,"leaves.jls"))))==m["leaves_sha256"]
    files=m["stage_files"];@assert length(files)==24
    @assert sort([parse(Int,match(r"stage_([0-9]+).jls",r["file"])[1]) for r in files])==collect(1:24)
    for r in files
        path=joinpath(dir,r["file"])
        @assert filesize(path)==r["bytes"]&&bytes2hex(sha256(read(path)))==r["sha256"]
    end
    states[tag]=m
end
@assert length(states)==8
for name in ("early","middle","late","barrier")
    a=states[name*"_before"];b=states[name*"_after"]
    @assert b["pass"]==a["pass"]+1
    if name=="barrier"
        @assert b["iteration"]==a["iteration"]&&b["barrier_iteration"]==a["barrier_iteration"]+1&&b["mu"]<a["mu"]
    else
        @assert b["iteration"]==a["iteration"]+1&&b["barrier_iteration"]==a["barrier_iteration"]&&b["mu"]==a["mu"]
        line=only(filter(x->startswith(x,"FILTERDDP_ITER_TIMING iteration=$(a["iteration"]) barrier_iteration=$(a["barrier_iteration"]) "),split(log,'\n')))
        @assert occursin("outcome=accepted",line)
    end
end
@assert states["early_after"]["pass"]<states["middle_before"]["pass"]<states["late_before"]["pass"]
open(joinpath(root,"capture_validation.toml"),"w") do io
    TOML.print(io,Dict("coverage_pass"=>true,"states"=>8,"full_stage_records"=>192,
        "ordinary_consecutive_pairs"=>3,"barrier_consecutive_pairs"=>1,"all_hashes_pass"=>true))
end
println("INFLUENCE_CAPTURE_COVERAGE_PASS states=8 full_stage_records=192 consecutive_ordinary=3 consecutive_barrier=1")

