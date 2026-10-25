# Matched IEEE9500 primary/secondary and dense SF follow-up, T=24 only.
const MATCH_ARGS=copy(ARGS)
ENV["SECONDARY_LIBRARY_ONLY"]="1"
include(joinpath(@__DIR__,"secondary_8500_study.jl"))
using Serialization
const MATCH_ROOT=normpath(get(ENV,"MATCHED_ROOT",joinpath(@__DIR__,"..","results","matched_secondary_oct25")))
const LABEL=MATCH_ARGS[1]
const PHASE=MATCH_ARGS[2]
const N=parse(Int,string(last(LABEL)))
function matched_network(label)
    n=parse(Int,string(last(label)));names=["subs$i" for i in 1:n]
    if startswith(label,"sf")
        deck="smartds_p1u_xf_1ph";extra=variant_switches(deck,"Master_4sub_typeIII.dss")
    else
        deck=startswith(label,"secondary") ? "ieee9500xf_1ph" : "ieee9500_1ph"
        extra=["Line.$t.enabled=true" for t in TIES_9500]
        if n==4
            # Preserve the historical 2705 fourth-source site and its positive delivery impedance.
            append!(extra,["New Line.poi4 phases=1 bus1=4s bus2=l2804256 r1=1e-6 x1=1e-6 c1=0 c0=0",
             "New Vsource.subs4 phases=1 bus1=4s basekv=7.1996 pu=1.05 angle=0 R1=0.104186 X1=0.617339 R0=0.104186 X0=0.617339"])
        end
    end
    net,pvs=load_sources(names;system=deck,extra)
    # Same delivery-point voltage contract in both representations, including direct-MV loads.
    bs=union(customer_buses(net),Set(l.bus for l in net.loads))
    s=setup(net,pvs;T=24,feas=parse(Float64,get(ENV,"MATCHED_FEAS","1e-7")),band_buses=bs)
    @assert Set(x.name for x in net.sources)==Set(names)
    say("$label: $(length(net.buses)) buses, $(length(bs)) delivery buses, $(length(net.batteries)) batteries; T=24")
    return s,names,deck
end
function matchdirs(n)
    n==2 && return [[-1.0],[1.0]]
    n==3 && return [[cos(2pi*i/32),sin(2pi*i/32)] for i in 0:31]
    return fib_dirs(78)
end
function checkdirect(s,case,x,angles,stage)
    if !(x.status in RAY_OK)
        appendrow(joinpath(ROOT8500,"direct_replays.csv"),"case,stage,status,r,replay_ok,dP_kW,dQ_kvar,dV_pu,dB_kW,vmin,vmax,min_import_kW,terminal_kWh,original_violation",[case,stage,x.status,x.r,false,fill(NaN,9)...]);return false
    end
    rp=replay_opendss(s.net,x,angles;load_factor=s.lam,pvs=s.P,pv_factor=s.F);dg=diagnostics8500(s,x)
    good=rp.converged && rp.pq_band && rp.dP<1e-3 && rp.dQ<1e-3 && rp.dV<1e-7 && rp.dB<1e-3 && max(dg.terminal,dg.socviol,dg.dyn,dg.bpow,dg.qviol)<.01
    appendrow(joinpath(ROOT8500,"direct_replays.csv"),"case,stage,status,r,replay_ok,dP_kW,dQ_kvar,dV_pu,dB_kW,vmin,vmax,min_import_kW,terminal_kWh,original_violation",[case,stage,x.status,x.r,good,rp.dP,rp.dQ,rp.dV,rp.dB,dg.vmin,dg.vmax,dg.minp,dg.terminal,max(dg.socviol,dg.dyn,dg.bpow,dg.qviol)])
    say(@sprintf("  ray replay %s dP %.2e dV %.2e | V %.6f %.6f minP %.4f",good ? "PASS" : "FAIL",rp.dP,rp.dV,dg.vmin,dg.vmax,dg.minp))
    good || error("Direct ray failed independent replay/original constraints")
    return dg.vmin>=.95-1e-6 && dg.vmax<=1.05+1e-6 && dg.minp>=-.01
end
function feasible_origin(s,case,start)
    x=solve_multiperiod(s.net,equal_angles(s.net);problem(s,case)...,objective=:cost,no_backflow=true,start,max_iter=500)
    good=checkdirect(s,case,x,equal_angles(s.net),"hard_origin")
    if good
        d=diagnostics8500(s,x)
        room=min(d.minp,20000*(d.vmin^2-.95^2),20000*(1.05^2-d.vmax^2))
        return merge(x,(;margin_kW=room))
    end
    return merge(x,(;margin_kW=NaN))
end
if get(ENV,"MATCHED_LIBRARY_ONLY","0")!="1"
s,names,deck=matched_network(LABEL)
cache=joinpath(MATCH_ROOT,LABEL,"state.jls")
if PHASE=="origin_test"
    packet=deserialize(cache)
    for case in Symbol.(split(get(ENV,"MATCHED_CASES","cyclic"),','))
        x=feasible_origin(s,case,packet.sens.x0)
        say("hard origin $case: status=$(x.status) feasible=$(feasible8500(x)) headroom=$(x.margin_kW)")
        feasible8500(x) || error("Origin pilot failed")
    end
elseif PHASE=="repair"
    packet=deserialize(cache);x0s=packet.x0s
    for case in (:nobatt,:cyclic)
        feasible8500(x0s[case]) && continue
        x=feasible_origin(s,case,packet.sens.x0);x0s[case]=x
        dg=x.status in RAY_OK ? diagnostics8500(s,x) : nothing
        appendrow(BASE8500,"configuration,case,T,status,margin_kW,vmin,vmax,min_import_kW,origin_feasible",[LABEL,case,24,x.status,x.margin_kW,dg===nothing ? NaN : dg.vmin,dg===nothing ? NaN : dg.vmax,dg===nothing ? NaN : dg.minp,feasible8500(x)])
    end
    backup=joinpath(MATCH_ROOT,LABEL,"state_initial.jls")
    isfile(backup) || cp(cache,backup)
    serialize(cache,merge(packet,(;x0s)))
elseif PHASE=="prep"
    sens=sensitivity8500(s,LABEL,names);sens===nothing && error("Sensitivity failure")
    serialize(joinpath(MATCH_ROOT,LABEL,"sensitivity_state.jls"),sens)
    x0s=Dict{Symbol,Any}()
    for case in (:nobatt,:cyclic)
        # The campaign needs an admissible origin, not a maximizer of spare headroom.
        # This uses hard network, no-backflow, voltage and cyclic-SOC constraints.
        x=feasible_origin(s,case,sens.x0)
        x0s[case]=x
        dg=x.status in RAY_OK ? diagnostics8500(s,x) : nothing
        appendrow(BASE8500,"configuration,case,T,status,margin_kW,vmin,vmax,min_import_kW,origin_feasible",[LABEL,case,24,x.status,x.margin_kW,dg===nothing ? NaN : dg.vmin,dg===nothing ? NaN : dg.vmax,dg===nothing ? NaN : dg.minp,feasible8500(x)])
    end
    open(joinpath(ROOT8500,"profiles.csv"),"w") do io
        println(io,"period,load_factor,pv_factor,price");for t in 1:24;println(io,join([t,s.lam[t],s.F[t],s.c[t]],','));end
    end
    serialize(cache,(;sens,x0s,names,deck))
else
    packet=deserialize(cache);@assert packet.names==names && packet.deck==deck
    sens,x0s=packet.sens,packet.x0s
    ids=parse.(Int,split(PHASE,','));dirs=matchdirs(N)
    cases=Symbol.(split(get(ENV,"MATCHED_CASES","nobatt,cyclic"),','))
    for case in cases
        x0=x0s[case]
        if case===:cyclic && get(ENV,"MATCHED_INITIAL","cost")=="idle"
            x0=x0s[:nobatt]
            @assert maximum((maximum(abs,v) for v in values(x0.P_B_kW));init=0.0)<1e-6
            say("  cyclic retry starts from the independently admissible idle trajectory")
        end
        feasible8500(x0) || (say("$case origin unresolved/negative local margin; no rays");continue)
        for i in ids
            u=dirs[i];dir=Dict(names[j] => (j==1 ? 0.0 : u[j-1]) for j in eachindex(names));at(r)=Dict(k=>r*v for (k,v) in dir)
            pred,lim=predicted_reach(sens,Dict(zip(names[2:end],u));day=case===:cyclic)
            x=solve_multiperiod(s.net,equal_angles(s.net);problem(s,case)...,objective=:ray,ray_dir=dir,ray_limit=RMAX8500,start=x0,max_iter=500)
            if get(ENV,"MATCHED_SAVE_ATTEMPTS","0")=="1"
                serialize(joinpath(ROOT8500,"$(case)_$(i)_direct.jls"),x)
            end
            r=x.r;method="unresolved";upper=NaN;best=x0
            if x.status in RAY_OK && isfinite(r) && checkdirect(s,case,x,at(r),"ray_$i")
                best=x
                # Hard-constrained direct maximum, then independently checked local sign bracket.
                inside=checked8500(s,LABEL,case,at(max(0.0,r-.002));start=x)
                if feasible8500(inside)
                    best=inside;lower=max(0.0,r-.002)
                    if r>=RMAX8500-1e-4
                        method="at_limit";r=lower;upper=Inf
                    else
                        outside=checked8500(s,LABEL,case,at(r+.01);start=x)
                        if outside.status in RAY_OK && outside.margin_kW<0
                            method="direct_checked";upper=r+.01;r=lower
                        elseif feasible8500(outside)
                            # A direct local maximum can be suboptimal. Fall back to fixed-angle search.
                            b=edge8500(s,LABEL,case,names,u,r+.01,x0);r=b.r;upper=b.upper;best=b.x;method=b.method
                        else
                            r=lower
                        end
                    end
                end
            elseif !(x.status in RAY_OK)
                checkdirect(s,case,x,at(0.0),"ray_$i")
                b=edge8500(s,LABEL,case,names,u,pred,x0);r=b.r;upper=b.upper;best=b.x;method=b.method
            end
            best === x0 && (r=0.0)
            dg=diagnostics8500(s,best);imports=Dict(k=>minimum(best.P_subs_kW[k]) for k in names);bind=dg.minp<=1 ? argmin(imports) : (dg.vmin<.9501 || dg.vmax>1.0499 ? "voltage" : "unknown")
            appendrow(EDGES8500,"configuration,case,T,direction,u1,u2,u3,predicted_deg,measured_deg,upper_deg,method,binding,min_import_kW,vmin,vmax,imports",[LABEL,case,24,i,[j<=length(u) ? u[j] : 0.0 for j in 1:3]...,pred,r,upper,method,bind,dg.minp,dg.vmin,dg.vmax,join(["$k=$(imports[k])" for k in names],';')])
            say("EDGE $LABEL $case $i/$(length(dirs)): $r ($method), predicted $pred, bound $bind")
        end
    end
end
say("COMPLETED $LABEL $PHASE")
end # executable driver; shared model and direction definitions remain reusable
