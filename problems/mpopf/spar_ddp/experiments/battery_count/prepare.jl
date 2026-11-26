using Serialization, SHA, Printf
const C=abspath(ARGS[1])
const original=joinpath(C,"inputs")
const horizons=(3,24,96,144,288)
const fields=(:B0_pu,:B_R,:B_R_pu,:P_B_R,:P_B_R_pu,:S_B_R_pu,:soc_min,:soc_max)
base=deserialize(joinpath(original,"network_data_large10kC_1ph_T96_periodic.jls"))
@assert length(base[:Bset])==1020
root=base[:substationBus]
function feeder(b)
    while base[:parent][b] != root; b=base[:parent][b]; end
    return b
end
groups=Dict{Int,Vector{Int}}()
for b in sort(base[:Bset]); push!(get!(groups,feeder(b),Int[]),b); end
ordered=Int[]
for rank in 1:maximum(length,values(groups)), g in sort(collect(keys(groups)))
    rank<=length(groups[g]) && push!(ordered,groups[g][rank])
end
@assert length(unique(ordered))==1020
open(joinpath(C,"battery_selection.csv"),"w") do io
    println(io,"rank,bus,feeder,power_kw,energy_kwh,in_X,in_Y")
    for (k,b) in enumerate(ordered)
        println(io,join((k,b,feeder(b),base[:P_B_R][b],base[:B_R][b],k<=250,k<=100),','))
    end
end
open(joinpath(C,"data_dictionary.txt"),"w") do io
    for k in sort(collect(keys(base));by=string); println(io,k," ",typeof(base[k])); end
end
open(joinpath(C,"variants.csv"),"w") do io
    println(io,"variant,T,batteries,feeder_coverage,total_kw,total_kwh,shared_entries_per_time_step,shared_bytes_all_time_steps,input_sha256")
    for (label,count) in (("10k_X",250),("10k_Y",100))
        chosen=Set(ordered[1:count]); src=joinpath(C,"variants",label,"source")
        @assert !ispath(src); mkpath(dirname(src)); cp(joinpath(C,"source"),src)
        datadir=joinpath(src,"spar_ddp/results/network_filterddp");mkpath(datadir)
        for T in horizons
            filename="network_data_large10kC_1ph_T$(T)_periodic.jls"
            a=deserialize(joinpath(original,filename))
            @assert a[:Bset]==base[:Bset]
            for f in fields; @assert isequal(a[f],base[f]); end
            # Refuse silently missing or unexpected battery metadata.
            @assert all(haskey(a,f) && a[f] isa AbstractDict for f in fields)
            b=deepcopy(a); b[:Bset]=[x for x in a[:Bset] if x in chosen];b[:n_B]=count
            for f in fields; filter!(p->first(p) in chosen,b[f]); end
            for k in keys(a)
                if k in fields
                    @assert Set(keys(b[k]))==chosen
                    @assert all(isequal(b[k][j],a[k][j]) for j in chosen)
                elseif k ∉ (:Bset,:n_B)
                    @assert isequal(a[k],b[k]) "non-battery field changed: $k"
                end
            end
            @assert length(b[:Bset])==b[:n_B]==count && issubset(chosen,Set(a[:Bset]))
            dest=joinpath(datadir,filename);serialize(dest,b)
            entries=3count^2+2count
            println(io,join((label,T,count,length(unique(feeder.(b[:Bset]))),sum(values(b[:P_B_R])),sum(values(b[:B_R])),entries,8T*entries,bytes2hex(sha256(read(dest)))),','))
        end
        storage=joinpath(src,"rawData/large10kC_1ph/Storage.dss")
        lines=readlines(storage); kept=String[];seen=Int[]
        for line in lines
            if occursin(r"(?i)^\s*new\s+storage\.",line)
                m=match(r"(?i)\bbus1\s*=\s*(\d+)",line);@assert !isnothing(m)
                bus=parse(Int,m[1]);bus in chosen || continue;push!(seen,bus)
            end
            push!(kept,line)
        end
        @assert length(seen)==count && Set(seen)==chosen
        write(storage,join(kept,'\n')*"\n")
        @assert read(joinpath(src,"rawData/large10kC_1ph/PVSystem.dss"))==read(joinpath(C,"source/rawData/large10kC_1ph/PVSystem.dss"))
        println(label," constructed and validated; feeder coverage=",length(unique(feeder.(collect(chosen)))))
    end
end
@assert issubset(Set(ordered[1:100]),Set(ordered[1:250]))
write(joinpath(C,"construction.complete"),"validated nested subsets and unchanged non-battery fields\n")
