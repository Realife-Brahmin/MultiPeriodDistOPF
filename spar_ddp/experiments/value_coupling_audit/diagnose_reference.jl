# One frozen point, no optimizer, no gate relaxation. Complete all comparisons
# before requiring production-path reproduction of the saved direction probes.
using Serialization,LinearAlgebra,SHA,TOML
length(ARGS)==3 || error("diagnose_reference.jl ORIGINAL_SOURCE CAPTURE_DIR OUTPUT")
source,capture,out=abspath.(ARGS);@assert !ispath(out);mkpath(out)
include(joinpath(source,"spar_ddp/experiments/value_coupling_audit/replay.jl"))
include(joinpath(@__DIR__,"production_reference.jl"))
meta=TOML.parsefile(joinpath(capture,"metadata.toml"))
leafpath=joinpath(capture,"leaves.jls")
@assert bytes2hex(sha256(read(leafpath)))==meta["leaves_sha256"]
payload=deserialize(leafpath);elems=payload.elems;n=length(first(elems).ma)
terminal=PintVal(zeros(n,n),zeros(n))
ord=ordinary(elems,terminal)
alloc=Any[deepcopy(e) for e in elems];push!(alloc,deepcopy(terminal))
_pp_rounds!(alloc,false,false)
production=production_values(elems,terminal)
maximum_prod_error=0.0
open(joinpath(out,"directions.csv"),"w") do io
    emit(io,"path,stage,probe,field,absolute_fro,reference_fro,unit_floor_scaled,relative_fro")
    open(joinpath(out,"values.csv"),"w") do vv
        emit(vv,"path,stage,field,absolute_fro,production_fro,unit_floor_scaled")
        for (mode,values) in (("ordinary",ord),("allocating_scan",alloc),("production_inplace",production))
            for t in 1:length(elems),f in (:V,:v)
                A=getproperty(values[t],f);B=getproperty(production[t],f)
                emit(vv,mode,t,f,norm(A-B),norm(B),scaled_error(A,B))
            end
            for r in meta["stage_files"]
                file=joinpath(capture,r["file"])
                @assert bytes2hex(sha256(read(file)))==r["sha256"]
                packet=deserialize(file);t=packet.stage
                rule=policy(packet,values[t+1],payload.mu)
                for (i,pr) in enumerate(packet.probes)
                    d=directions(rule,pr.dx)
                    for f in (:du,:dphi)
                        A=getproperty(d,f);B=getproperty(pr,f);err=scaled_error(A,B)
                        emit(io,mode,t,i,f,norm(A-B),norm(B),err,norm(A-B)/max(norm(B),floatmin(Float64)))
                        global maximum_prod_error=mode=="production_inplace" ? max(maximum_prod_error,err) : maximum_prod_error
                    end
                end
            end
        end
    end
end
println("REFERENCE_DIAGNOSIS production_max_saved_probe_error=$maximum_prod_error same_gate=1e-7")
@assert maximum_prod_error<1e-7 "production-path dense reproduction still fails; do not relax the gate"
println("REFERENCE_DIAGNOSIS_PASS")
