# Passive saved-array audit; no model construction or optimization.
using Serialization,TOML,SHA
const AUDIT_LAYOUT=get(ENV,"TADMM_AUDIT_LAYOUT",joinpath(@__DIR__,"..","..","spar_ddp/examples/power_system/control_layout.jl"))
include(AUDIT_LAYOUT)
function energy_audit(folder,out)
    d=deserialize(joinpath(folder,"instance.jls"));ix,_=control_layout(d)
    report=Dict{String,Any}("julia"=>string(VERSION),"T"=>d[:T],"batteries"=>length(d[:Bset]),
        "base_kva"=>d[:kVA_B],"instance_sha256"=>bytes2hex(sha256(read(joinpath(folder,"instance.jls")))),
        "audit_source_sha256"=>bytes2hex(sha256(read(@__FILE__))),"layout_sha256"=>bytes2hex(sha256(read(AUDIT_LAYOUT))))
    bp=Dict(j=>n for (n,j) in enumerate(d[:Nset]))
    parent=[bp[e[1]] for e in d[:Lset]]
    for (root,_,files) in walkdir(folder),file in files
        endswith(file,".jls") && file!="instance.jls" || continue
        path=joinpath(root,file)
        occursin("warmup",relpath(path,folder)) && continue
        sol=deserialize(path)
        haskey(sol,:u) && haskey(sol,:B) || continue
        energy=Float64[d[:B0_pu][j] for j in d[:Bset]]
        bound=0.0;socpp=0.0;dyn=0.0;drift=0.0;power=0.0;consensusbound=0.0
        bus=0;period=0;socbus=0;soctime=0
        lower_binding=0;upper_binding=0;power_binding=0;energy_binding=0;normalized_cone=0.0
        for t in 1:d[:T]
            pb=sol[:u][t][ix.pb]
            u=sol[:u][t];v=u[ix.v];p=u[ix.P];q=u[ix.Q];ell=u[ix.ell]
            for (n,j) in enumerate(d[:Nset])
                j==d[:substationBus] && continue
                lower_binding+=abs(v[n]-d[:Vminpu][j]^2)<=1e-6
                upper_binding+=abs(v[n]-d[:Vmaxpu][j]^2)<=1e-6
            end
            for n in eachindex(parent)
                ve=v[parent[n]]*ell[n];pq=p[n]^2+q[n]^2
                normalized_cone=max(normalized_cone,abs(ve-pq)/max(1.0,abs(ve),pq))
            end
            for (n,j) in enumerate(d[:Bset])
                cap=d[:B_R_pu][j];lo=d[:soc_min][j]*cap;hi=d[:soc_max][j]*cap
                energy[n]-=d[:delta_t_h]*pb[n]
                error=max(lo-energy[n],energy[n]-hi,0.0)
                if error>bound;bound=error;bus=j;period=t;end
                if 100error/cap>socpp;socpp=100error/cap;socbus=j;soctime=t;end
                prev=t==1 ? d[:B0_pu][j] : sol[:B][n,t-1]
                dyn=max(dyn,abs(sol[:B][n,t]-prev+d[:delta_t_h]*pb[n]))
                drift=max(drift,abs(sol[:B][n,t]-energy[n]))
                power=max(power,abs(pb[n])-d[:P_B_R_pu][j])
                power_binding+=abs(abs(pb[n])-d[:P_B_R_pu][j])<=1e-6
                energy_binding+=min(abs(sol[:B][n,t]-lo),abs(sol[:B][n,t]-hi))/cap<=1e-6
                consensusbound=max(consensusbound,lo-sol[:B][n,t],sol[:B][n,t]-hi)
            end
        end
        key=replace(relpath(path,folder),'\\'=>'/')
        report[key]=Dict("integrated_energy_bound_pu_h"=>bound,"integrated_energy_bound_kwh"=>bound*d[:kVA_B],
            "max_SOC_bound_percentage_points"=>socpp,"max_energy_error_bus"=>bus,"max_energy_error_period"=>period,
            "max_SOC_error_bus"=>socbus,"max_SOC_error_period"=>soctime,
            "max_dynamics_error_kwh"=>dyn*d[:kVA_B],"max_integrated_consensus_difference_kwh"=>drift*d[:kVA_B],
            "consensus_energy_bound_pu_h"=>consensusbound,"max_power_bound_kw"=>power*d[:kVA_B],
            "voltage_lower_binding_nonroot"=>lower_binding,"voltage_upper_binding_nonroot"=>upper_binding,
            "battery_power_binding"=>power_binding,"consensus_SOC_binding"=>energy_binding,
            "voltage_power_binding_tolerance_pu"=>1e-6,"SOC_binding_tolerance_fraction"=>1e-6,
            "normalized_cone_gap"=>normalized_cone,
            "normalized_cone_gap_definition"=>"abs(v*ell-P^2-Q^2)/max(1,abs(v*ell),P^2+Q^2)",
            "thermal_upper_constraints"=>"not modeled in this pilot; no thermal binding count claimed",
            "snapshot_sha256"=>bytes2hex(sha256(read(path))))
        if socbus!=0
            n=findfirst(==(socbus),d[:Bset]);cap=d[:B_R_pu][socbus]
            e=d[:B0_pu][socbus];actual=[100e/cap];consensus=[100e/cap]
            for t in 1:d[:T]
                e-=d[:delta_t_h]*sol[:u][t][ix.pb[n]]
                push!(actual,100e/cap);push!(consensus,100sol[:B][n,t]/cap)
            end
            report[key]["hours"]=collect(0:d[:T]).*d[:delta_t_h]
            report[key]["integrated_SOC_percent"]=actual
            report[key]["consensus_SOC_percent"]=consensus
            report[key]["minimum_SOC_percent"]=100d[:soc_min][socbus]
            report[key]["maximum_SOC_percent"]=100d[:soc_max][socbus]
        end
        println(key," bound_kwh=",bound*d[:kVA_B]," SOC_bound_percentage_points=",socpp)
    end
    open(io->TOML.print(io,report),out,"w")
end
energy_audit(ARGS...)
