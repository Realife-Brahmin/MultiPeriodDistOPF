# Definitions fixed before the sweep; raw crossings are not convergence claims.
function sweep_flags(owned,allcopies,dual,gap,physical)
    flags=Dict{String,Bool}()
    for (label,eps) in (("1e-3",1e-3),("1e-4",1e-4),("1e-5",1e-5))
        flags["owned_raw_"*label]=owned<=eps
        flags["all_raw_"*label]=allcopies<=eps
        flags["owned_dual_"*label]=owned<=eps && dual<=1e-2
        flags["all_dual_"*label]=owned<=eps && allcopies<=eps && dual<=1e-2
        flags["paper_cost_"*label]=owned<=eps && gap<=.005
        flags["current_cost_"*label]=owned<=eps && allcopies<=eps && gap<=.001
    end
    flags["physical_NO"]=physical<=1e-5 && gap<=.001
    flags["final_rms"]=owned<=1e-5 && allcopies<=1e-5 && dual<=1e-2
    return flags
end

function sweep_observer(folder,d,system,reference;save_points=true)
    milestones=Dict{String,Any}()
    function observe(x)
        gap=abs(x.metrics.objective-reference)/max(1.0,abs(reference))
        flags=sweep_flags(x.owned_r,x.rr,x.owned_s,gap,x.metrics.original_primal)
        fresh=sort!([key for (key,yes) in flags if yes && !haskey(milestones,key)])
        point="point_"*lpad(string(x.k),5,'0')*".jls"
        if !isempty(fresh) && save_points
            serialize(joinpath(folder,point),Dict(:system=>system,:T=>d[:T],:objective=>x.metrics.objective,
                :u=>x.us,:B=>x.bb,:status=>"SAVED_SWEEP_MILESTONE",:iteration=>x.k,:labels=>fresh))
        end
        digest=!isempty(fresh) && save_points ? bytes2hex(sha256(read(joinpath(folder,point)))) : ""
        for key in fresh
            milestones[key]=Dict("iteration"=>x.k,"core_s"=>x.elapsed,"snapshot"=>point,"snapshot_sha256"=>digest,
                "owned_r_rms"=>x.owned_r,"all_r_rms"=>x.rr,"owned_s_rms"=>x.owned_s,
                "objective"=>x.metrics.objective,"relative_objective_gap"=>gap,
                "original_primal"=>x.metrics.original_primal,"maxeq"=>x.metrics.maxeq,
                "maxbound"=>x.metrics.maxbound,"maxcone"=>x.metrics.maxcone,
                "cone_equality_gap"=>x.metrics.soc_equality_gap,"fictitious_loss_kwh"=>x.metrics.fictitious_loss_kwh)
            @printf("SWEEP_MILESTONE label=%s k=%d core_s=%.6f owned=%.3e all=%.3e dual=%.3e gap=%.3e original=%.3e\n",
                key,x.k,x.elapsed,x.owned_r,x.rr,x.owned_s,gap,x.metrics.original_primal)
        end
        if !isempty(fresh)
            open(io->TOML.print(io,milestones),joinpath(folder,"milestones.toml"),"w")
            flush(stdout)
        end
        return flags["final_rms"]
    end
    return observe,milestones
end
