# Late frozen correctness replay using the previously validated policy probes.
include(joinpath(ENV["VALUE_FOLLOWUP_BASE"],"replay.jl"))
include(joinpath(@__DIR__,"schedule.jl"))
function check_late(cap)
    meta=TOML.parsefile(joinpath(cap,"metadata.toml"))
    leaf=joinpath(cap,"leaves.jls")
    @assert bytes2hex(open(sha256,leaf)) == meta["leaves_sha256"] == "cd7922cf75df0743cf5e3fd38e2e87e4e65a861346e89bc66512cf1ffb90f416"
    saved=deserialize(leaf)
    elems=saved.elems
    N=length(elems)
    n=length(first(elems).ma)
    @assert N==24 && n==1020 && "late" in meta["tags"]
    terminal=PintVal(zeros(n,n),zeros(n))
    balanced=production_values(elems,terminal)
    serial=ordinary(elems,terminal)
    function grouped_values(elems,terminal)
        xs=Any[deepcopy(e) for e in elems]
        push!(xs,deepcopy(terminal))
        work=PairWork(length(first(elems).ma))
        for round in replay_rounds(length(elems),2)
            for (i,j) in round
                xs[i]=_pp_combine!(xs[i],xs[j],work)
            end
        end
        @assert all(x->x isa PintVal,xs)
        xs
    end
    grouped=grouped_values(elems,terminal)
    accurate=true
    for (label,values,refs) in (("balanced_serial",balanced,serial),("grouped_serial",grouped,serial),("grouped_balanced",grouped,balanced))
        for t in 1:N, f in (:v,:V)
            err=scaled_error(getproperty(values[t],f),getproperty(refs[t],f))
            println("LATE_VALUE mode=$label stage=$t field=$f error=$err")
            accurate &= isfinite(err) && err<1e-7
        end
    end
    firstpacket=nothing
    for row in meta["stage_files"]
        path=joinpath(cap,row["file"])
        @assert filesize(path)==row["bytes"] && bytes2hex(open(sha256,path))==row["sha256"]
        packet=deserialize(path)
        t=packet.stage
        t==1 && (firstpacket=packet)
        rules=[policy(packet,values[t+1],saved.mu) for values in (balanced,grouped)]
        for (k,probe) in enumerate(packet.probes)
            ref,got=[directions(rule,probe.dx) for rule in rules]
            for f in (:du,:dphi,:dzl,:dzu)
                a,b=getproperty(got,f),getproperty(ref,f)
                e=scaled_error(a,b)
                rel=norm(a-b)/max(1e-10,norm(b))
                println("LATE_DIRECTION stage=$t probe=$k field=$f error=$e relative=$rel")
                accurate &= isfinite(e) && e<1e-7 && rel<1e-5
                if f in (:du,:dphi)
                    e_saved=scaled_error(b,getproperty(probe,f))
                    println("LATE_SAVED_PROBE stage=$t probe=$k field=$f error=$e_saved")
                    accurate &= isfinite(e_saved) && e_saved<1e-7
                end
            end
        end
    end
    @assert firstpacket!==nothing
    ref=battery_path(elems,balanced,firstpacket.f)
    got=battery_path(elems,grouped,firstpacket.f)
    for t in 1:N,f in (:dx,:dp)
        a,b=getproperty(got[t],f),getproperty(ref[t],f)
        e=scaled_error(a,b)
        rel=norm(a-b)/max(1e-10,norm(b))
        println("LATE_BATTERY stage=$t field=$f error=$e relative=$rel")
        accurate &= isfinite(e) && e<1e-7 && rel<1e-5
    end
    println("LATE_COMPLETE accurate=$accurate")
    @assert accurate
end
check_late(only(ARGS))
