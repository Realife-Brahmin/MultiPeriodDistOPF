include("model.jl")
using TOML, Statistics
reference,out=ARGS
mkpath(out); report=Dict{String,Any}()
for (system,job) in (("ieee2522C_1ph","31650588"),("large10kC_1ph","31650589"))
    case=joinpath(reference,"results/$(system)_T24_july_pair_$(job)")
    d=deserialize(joinpath(case,"instance.jls"))
    sol=deserialize(joinpath(case,"exact/final_dispatch.jls"))
    us=sol[:u]; bb=sol[:B]; ix,_=control_layout(d)
    c=typed_physical_layout(d)
    probes=("endpoint","negative_ell","negative_import","voltage_lower","voltage_upper",
        "battery_power","battery_energy","branch_flow","cone","reactive_bound","integrated_soc")
    for probe in probes
        u=deepcopy(us); b=copy(bb)
        if probe=="negative_ell"; u[1][first(ix.ell)]=-0.01
        elseif probe=="negative_import"; u[1][ix.ps]=-0.01
        elseif probe=="voltage_lower"; u[1][last(ix.v)]=0.1
        elseif probe=="voltage_upper"; u[1][last(ix.v)]=2.0
        elseif probe=="battery_power"; u[1][first(ix.pb)]=d[:P_B_R_pu][first(d[:Bset])]+1.0
        elseif probe=="battery_energy"; b[1,1]+=1.0
        elseif probe=="branch_flow"; u[1][first(ix.P)]+=1.0; u[1][first(ix.Q)]-=0.5
        elseif probe=="cone"; u[1][first(ix.ell)]=0.0
        elseif probe=="reactive_bound"; u[1][first(ix.qnorm)]=2.0
        elseif probe=="integrated_soc"; for t in 1:d[:T]; u[t][first(ix.pb)]=d[:P_B_R_pu][first(d[:Bset])]; end
        end
        a=physical_metrics_cached(d,u,b,physical_layout(d)); z=typed_physical_metrics(u,b,c)
        for key in keys(a)
            @assert isapprox(a[key],z[key];atol=1e-10,rtol=1e-13) "$system $probe $key: $(a[key]) != $(z[key])"
        end
    end
    # Alternate implementations after both have compiled; no solver calls.
    oldtimes=Float64[]; newtimes=Float64[]
    for _ in 1:7
        push!(oldtimes,@elapsed physical_metrics_cached(d,us,bb,physical_layout(d)))
        push!(newtimes,@elapsed typed_physical_metrics(us,bb,c))
    end
    # Weight changes must select a fresh numeric cache.
    d[:C_B]*=2
    c2=typed_physical_layout(d)
    @assert c2!==c && c2.cb==2c.cb
    d[:C_B]/=2
    @assert typed_physical_layout(d)===c
    report[system]=Dict("probes"=>collect(probes),"all_metrics_match"=>true,
        "cached_median_s"=>median(oldtimes),"typed_median_s"=>median(newtimes),
        "isolated_speedup"=>median(oldtimes)/median(newtimes),"original_primal"=>typed_physical_metrics(us,bb,c).original_primal)
    println(system," ",report[system]); flush(stdout)
end
open(io->TOML.print(io,report),joinpath(out,"diagnostic_profile.toml"),"w")
println("TYPED_DIAGNOSTICS_VALIDATED")
