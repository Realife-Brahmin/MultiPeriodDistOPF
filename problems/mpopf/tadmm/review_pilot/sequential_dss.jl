module SequentialDSS
using OpenDSSDirect, Serialization, TOML, SHA, Printf
const D=OpenDSSDirect
const REPO=normpath(joinpath(@__DIR__,"..","..","..",".."))
include(joinpath(@__DIR__,"..","..","spar_ddp","examples","power_system","control_layout.jl"))
busnum(s)=parse(Int,first(split(strip(s),'.')))
prop(s)=D.Properties.Value(String(s))
number(s)=parse(Float64,prop(s))
csv(io,x)=println(io,join(x,','))
energy_delta(p,dt,ec,ed)=dt*(ec*max(-p,0.0)-max(p,0.0)/ed)
bound_error(e,lo,hi)=max(0.0,lo-e,e-hi)

"""One explicitly clocked step. Never overwrites stored energy."""
function storage_step!(names,command,dt,ec,ed)
    before=Float64[]
    for (name,p) in zip(names,command)
        D.Circuit.SetActiveElement("Storage."*name)
        push!(before,number("kWhStored"))
        D.Text.Command("Edit Storage.$name kW=$p kvar=0")
    end
    h0=D.Solution.DblHour()
    D.Solution.SolveSnap()
    actual=Float64[]
    for (i,name) in enumerate(names)
        D.Circuit.SetActiveElement("Storage."*name)
        @assert isapprox(number("kWhStored"),before[i];atol=1e-7,rtol=1e-12) "SolveSnap unexpectedly advanced energy"
        push!(actual,-real(first(D.CktElement.Powers())))
    end
    return (;before,actual,h0,converged=D.Solution.Converged())
end

function finish_storage_step!(names,step,dt,ec,ed)
    D.Solution.FinishTimeStep() # exactly one energy update and clock increment
    @assert isapprox(D.Solution.DblHour(),step.h0+dt;atol=1e-8,rtol=0) "Wrong time advancement"
    after=Float64[];balance=Float64[]
    for (i,name) in enumerate(names)
        D.Circuit.SetActiveElement("Storage."*name)
        push!(after,number("kWhStored"))
        push!(balance,after[end]-step.before[i]-energy_delta(step.actual[i],dt,ec[i],ed[i]))
    end
    return (;after,balance,h1=D.Solution.DblHour())
end

"""Chronological Storage replay with full device/period discrepancy tables."""
function replay(instance,solution,out;efficiency="model",repo=REPO)
    efficiency in ("model","dss") || error("efficiency must be model or dss")
    ispath(out) && error("Refusing to overwrite replay output")
    mkpath(out);out=abspath(out)
    data=deserialize(instance);sol=deserialize(solution);ix,nu=control_layout(data)
    T=data[:T];dt=data[:delta_t_h];base=data[:kVA_B];root=data[:substationBus]
    @assert length(sol[:u])==T && all(length(u)==nu for u in sol[:u])
    buses=data[:Nset];lines=data[:Lset];bs=data[:Bset];ds=data[:Dset]
    bp=Dict(j=>n for (n,j) in enumerate(buses));bat=Dict(j=>n for (n,j) in enumerate(bs));der=Dict(j=>n for (n,j) in enumerate(ds))
    D.Text.Command("clear")
    D.Text.Command("Redirect \"$(joinpath(repo,"test_networks",sol[:system],"Master.dss"))\"")
    D.Text.Command("Set DataPath=\"$out\"")
    D.Text.Command("Edit Vsource.source pu=1.05")
    D.Text.Command("BatchEdit Line..* c1=0 c0=0")
    # Explicit series are applied below; remove any daily/yearly multiplier.
    loads=Tuple{String,Float64,Float64}[]
    for name in D.Loads.AllNames()
        D.Loads.Name(name);push!(loads,(name,D.Loads.kW(),D.Loads.kvar()))
        D.Text.Command("Edit Load.$name daily=none yearly=none duty=none Vminpu=0.3 Vmaxpu=1.7")
    end
    pvs=Tuple{Int,Float64,String}[]
    for name in D.PVsystems.AllNames()
        name=="NONE" && continue
        D.Circuit.SetActiveElement("PVSystem."*name)
        push!(pvs,(busnum(prop("bus1")),number("Pmpp"),prop("kV")))
        D.Text.Command("Edit PVSystem.$name enabled=no")
    end
    filter!(p->p[1]!=root,pvs)
    @assert sort(first.(pvs))==sort(filter(!=(root),ds))
    for (j,_,kv) in pvs
        D.Text.Command("New Load.seq_pv$j bus1=$j.1 phases=1 model=1 kV=$kv kW=0 kvar=0 daily=none Vminpu=0.3 Vmaxpu=1.7")
    end
    names=String[];buses_b=Int[];capacity=Float64[];ec=Float64[];ed=Float64[];original_ec=Float64[];original_ed=Float64[]
    for name in D.Storages.AllNames()
        name=="NONE" && continue
        D.Circuit.SetActiveElement("Storage."*name);j=busnum(prop("bus1"))
        if j==root
            D.Text.Command("Edit Storage.$name enabled=no");continue
        end
        @assert haskey(bat,j) && D.CktElement.NumPhases()==1
        cap=number("kWhRated");@assert isapprox(cap,data[:B_R_pu][j]*base;atol=1e-6,rtol=1e-9)
        @assert isapprox(number("kWrated"),data[:P_B_R_pu][j]*base;atol=1e-6,rtol=1e-9)
        @assert number("%IdlingkW")==0 "Nonzero idle losses need an explicit energy-accounting model"
        @assert lowercase(strip(prop("EffCurve"))) in ("","none") "Inverter efficiency curve not supported silently"
        push!(original_ec,number("%EffCharge")/100);push!(original_ed,number("%EffDischarge")/100)
        if efficiency=="model"
            D.Text.Command("Edit Storage.$name %EffCharge=100 %EffDischarge=100")
        end
        push!(ec,number("%EffCharge")/100);push!(ed,number("%EffDischarge")/100)
        # Initialize once only. Operating upper SOC is checked, not replaced by
        # a fake smaller kWhRated value. Native reserve clamping is also audited.
        e0=data[:B0_pu][j]*base
        D.Text.Command("Edit Storage.$name DispMode=External model=1 kWhStored=$e0 %reserve=$(100data[:soc_min][j]) %cutin=0 %cutout=0 Vminpu=0.3 Vmaxpu=1.7 daily=none yearly=none duty=none")
        @assert lowercase(prop("DispMode"))=="external"
        push!(names,name);push!(buses_b,j);push!(capacity,cap)
    end
    @assert sort(buses_b)==sort(filter(!=(root),bs))
    @assert length(unique(buses_b))==length(buses_b)
    linemap=Dict{Tuple{Int,Int},Tuple{String,Int}}();line_set=Set(lines)
    for name in D.Lines.AllNames()
        D.Lines.Name(name);ends=busnum.(D.CktElement.BusNames());nc=D.CktElement.NumConductors()
        @assert D.CktElement.NumPhases()==1
        edge=(ends[1],ends[2])
        if edge in line_set
            linemap[edge]=(name,1)
        else
            edge=(ends[2],ends[1]);@assert edge in line_set;linemap[edge]=(name,nc+1)
        end
    end
    @assert length(linemap)==length(lines)
    D.Text.Command("Set mode=daily stepsize=$(3600dt)s number=1 controlmode=off hour=0 sec=0 tolerance=1e-9 maxiterations=200")
    requested=[data[:B0_pu][j]*base for j in buses_b];model_requested=copy(requested)
    expected_before=copy(requested)
    stats=Dict{String,Any}("efficiency_mode"=>efficiency,"T"=>T,"dt_h"=>dt,"batteries"=>length(names),
        "engine_version"=>D.Basic.Version(),"instance_sha256"=>bytes2hex(sha256(read(instance))),
        "solution_sha256"=>bytes2hex(sha256(read(solution))),"source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "original_charge_efficiency_min"=>minimum(original_ec),"original_discharge_efficiency_min"=>minimum(original_ed),
        "used_charge_efficiency_min"=>minimum(ec),"used_discharge_efficiency_min"=>minimum(ed),
        "converged"=>true,"max_voltage_difference_pu"=>0.0,"voltage_violations"=>0,
        "max_power_tracking_kw"=>0.0,"max_engine_energy_balance_kwh"=>0.0,
        "max_engine_SOC_bound_kwh"=>0.0,"max_requested_SOC_bound_kwh"=>0.0,
        "max_engine_vs_opf_energy_kwh"=>0.0,"max_substation_difference_kw"=>0.0,
        "max_line_P_difference_kw"=>0.0,"max_line_Q_difference_kvar"=>0.0,"max_ell_difference_pu"=>0.0,
        "reserve_or_capacity_clamp_records"=>0,"power_tracking_records"=>0,"opf_objective"=>get(sol,:objective,NaN),
        "note"=>"Live Storage chronology; energy initialized once. Snapshot checks are separate. No thermal upper compliance is claimed. PV injections, zero line capacitance, source1.05pu and constant-power voltage range match the optimization. Efficiency=dss retains battery efficiencies, not every original network setting.")
    objective=0.0;command_objective=0.0;files=Dict{String,String}()
    ios=Dict(name=>open(joinpath(out,name*".csv"),"w") for name in ("battery","bus","branch","period"))
    csv(ios["battery"],("t","bus","name","command_kw","actual_kw","power_difference_kw","energy_before_kwh","energy_after_kwh","opf_energy_kwh","engine_minus_opf_kwh","requested_efficiency_energy_kwh","requested_model_energy_kwh","engine_energy_balance_error_kwh","SOC_engine_percent","SOC_opf_percent","SOC_low_kwh","SOC_high_kwh","engine_bound_error_kwh","requested_bound_error_kwh","engine_state"))
    csv(ios["bus"],("t","bus","opf_voltage_pu","dss_voltage_pu","difference_pu","lower_pu","upper_pu"))
    csv(ios["branch"],("t","from_bus","to_bus","opf_P_kw","dss_P_kw","P_difference_kw","opf_Q_kvar","dss_Q_kvar","Q_difference_kvar","opf_ell_pu","dss_ell_pu","ell_difference_pu","dss_line_loss_kw"))
    csv(ios["period"],("t","hour_start","hour_end","converged","opf_substation_kw","dss_substation_kw","substation_difference_kw","opf_substation_kvar","dss_substation_kvar","opf_line_loss_kw","dss_line_loss_kw","command_battery_kw","actual_battery_kw","max_voltage_difference_pu","max_dispatch_difference_kw","max_engine_energy_balance_kwh"))
    try
        for t in 1:T
            @assert isapprox(D.Solution.DblHour(),(t-1)*dt;atol=1e-8)
            u=sol[:u][t]
            for (name,p,q) in loads
                D.Text.Command("Edit Load.$name kW=$(p*data[:LoadShapeLoad][t]) kvar=$(q*data[:LoadShapeLoad][t])")
            end
            for (j,pmpp,_) in pvs
                qmax=sqrt(max(0.0,data[:S_D_R][j]^2-data[:p_D_pu][j,t]^2))
                D.Text.Command("Edit Load.seq_pv$j kW=$(-pmpp*data[:LoadShapePV][t]) kvar=$(-qmax*u[ix.qnorm[der[j]]]*base)")
            end
            command=[u[ix.pb[bat[j]]]*base for j in buses_b]
            step=storage_step!(names,command,dt,ec,ed)
            @assert all(isapprox.(step.before,expected_before;atol=1e-7,rtol=1e-12)) "Energy was reset between periods"
            stats["converged"] &= step.converged
            actual_sub=-D.Circuit.TotalPower();ps=real(actual_sub);qs=imag(actual_sub)
            nodes=D.Circuit.AllNodeNames();vraw=D.Circuit.AllBusVMag()./(data[:kV_B]*1000)
            vd=Dict(busnum(n)=>vraw[i] for (i,n) in enumerate(nodes));vdelta=0.0
            for j in buses
                vo=sqrt(max(0.0,u[ix.v[bp[j]]]))
                diff=vd[j]-vo;vdelta=max(vdelta,abs(diff))
                stats["voltage_violations"]+=j!=root && (vd[j]<data[:Vminpu][j]-1e-4 || vd[j]>data[:Vmaxpu][j]+1e-4)
                csv(ios["bus"],(t,j,vo,vd[j],diff,data[:Vminpu][j],data[:Vmaxpu][j]))
            end
            loss=0.0;opf_loss=0.0
            for (n,edge) in enumerate(lines)
                name,pos=linemap[edge];D.Circuit.SetActiveElement("Line."*name)
                s=D.CktElement.Powers()[pos];amp=abs(D.CktElement.Currents()[pos])
                ell=(amp/(base/data[:kV_B]))^2;lp=real(first(D.CktElement.Losses()))/1000
                po=u[ix.P[n]]*base;qo=u[ix.Q[n]]*base;lo=u[ix.ell[n]]
                stats["max_line_P_difference_kw"]=max(stats["max_line_P_difference_kw"],abs(real(s)-po))
                stats["max_line_Q_difference_kvar"]=max(stats["max_line_Q_difference_kvar"],abs(imag(s)-qo))
                stats["max_ell_difference_pu"]=max(stats["max_ell_difference_pu"],abs(ell-lo))
                loss+=lp;opf_loss+=data[:rdict_pu][edge]*lo*base
                csv(ios["branch"],(t,edge[1],edge[2],po,real(s),real(s)-po,qo,imag(s),imag(s)-qo,lo,ell,ell-lo,lp))
            end
            done=finish_storage_step!(names,step,dt,ec,ed);expected_before=copy(done.after)
            track=maximum(abs.(step.actual.-command);init=0.0);bal=maximum(abs.(done.balance);init=0.0)
            for (n,j) in enumerate(buses_b)
                requested[n]+=energy_delta(command[n],dt,ec[n],ed[n]);model_requested[n]-=dt*command[n]
                lo=data[:soc_min][j]*capacity[n];hi=data[:soc_max][j]*capacity[n]
                be=bound_error(done.after[n],lo,hi);re=bound_error(requested[n],lo,hi);eo=sol[:B][bat[j],t]*base
                stats["max_engine_SOC_bound_kwh"]=max(stats["max_engine_SOC_bound_kwh"],be)
                stats["max_requested_SOC_bound_kwh"]=max(stats["max_requested_SOC_bound_kwh"],re)
                stats["max_engine_vs_opf_energy_kwh"]=max(stats["max_engine_vs_opf_energy_kwh"],abs(done.after[n]-eo))
                stats["reserve_or_capacity_clamp_records"]+=abs(done.balance[n])>base*1e-5
                stats["power_tracking_records"]+=abs(step.actual[n]-command[n])>base*1e-5
                D.Circuit.SetActiveElement("Storage."*names[n])
                csv(ios["battery"],(t,j,names[n],command[n],step.actual[n],step.actual[n]-command[n],step.before[n],done.after[n],eo,done.after[n]-eo,requested[n],model_requested[n],done.balance[n],100done.after[n]/capacity[n],100eo/capacity[n],lo,hi,be,re,prop("State")))
            end
            stats["max_voltage_difference_pu"]=max(stats["max_voltage_difference_pu"],vdelta)
            stats["max_substation_difference_kw"]=max(stats["max_substation_difference_kw"],abs(ps-u[ix.ps]*base))
            stats["max_power_tracking_kw"]=max(stats["max_power_tracking_kw"],track)
            stats["max_engine_energy_balance_kwh"]=max(stats["max_engine_energy_balance_kwh"],bal)
            battery_cost=data[:C_B]*dt*sum(abs2,step.actual)
            objective+=data[:LoadShapeCost][t]*dt*ps+battery_cost
            command_objective+=data[:LoadShapeCost][t]*dt*u[ix.ps]*base+data[:C_B]*dt*sum(abs2,command)
            csv(ios["period"],(t,step.h0,done.h1,step.converged,u[ix.ps]*base,ps,ps-u[ix.ps]*base,u[ix.qs]*base,qs,opf_loss,loss,sum(command),sum(step.actual),vdelta,track,bal))
        end
    finally
        foreach(close,values(ios))
    end
    stats["replayed_objective"]=objective;stats["command_objective"]=command_objective
    stats["terminal_hour"]=D.Solution.DblHour();stats["physical_tolerance_pu"]=1e-5
    stats["network_voltage_agreement_pass"]=stats["converged"] && stats["voltage_violations"]==0 && stats["max_voltage_difference_pu"]<=1e-4
    stats["dispatch_energy_execution_pass"]=stats["converged"] && stats["voltage_violations"]==0 &&
        maximum(stats[k]/base for k in ("max_power_tracking_kw","max_engine_energy_balance_kwh","max_engine_SOC_bound_kwh","max_requested_SOC_bound_kwh"))<=1e-5
    stats["opf_state_agreement_pass"]=stats["max_engine_vs_opf_energy_kwh"]/base<=1e-5 &&
        maximum(stats[k]/base for k in ("max_substation_difference_kw","max_line_P_difference_kw","max_line_Q_difference_kvar"))<=1e-5 && stats["network_voltage_agreement_pass"]
    stats["sequential_schedule_pass"]=stats["dispatch_energy_execution_pass"] && stats["opf_state_agreement_pass"]
    for name in keys(ios)
        path=joinpath(out,name*".csv");files[name*".csv"]=bytes2hex(sha256(read(path)))
    end
    stats["table_sha256"]=files
    open(io->TOML.print(io,stats),joinpath(out,"summary.toml"),"w")
    println("SEQUENTIAL_DSS efficiency=$efficiency pass=",stats["sequential_schedule_pass"]," max_tracking_kw=",stats["max_power_tracking_kw"]," max_energy_balance_kwh=",stats["max_engine_energy_balance_kwh"]," max_SOC_bound_kwh=",stats["max_engine_SOC_bound_kwh"])
    return stats
end
end
