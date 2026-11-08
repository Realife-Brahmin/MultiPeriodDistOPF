# Validate exact DDP screening and cached diagnostics against frozen C4 code.
using Serialization, TOML, SHA, JuMP
newsource,reference,out=ARGS
mkpath(out)
old=Module(:FrozenReference); fresh=Module(:DDPTransferred)
Base.include(old,joinpath(reference,"source/problems/mpopf/tadmm/review_pilot/model.jl"))
Base.include(fresh,joinpath(newsource,"problems/mpopf/tadmm/review_pilot/model.jl"))
results=Dict{String,Any}()
for (system,job) in (("ieee2522C_1ph","31645746"),("large10kC_1ph","31645790"))
    case=joinpath(reference,"results/$(system)_T24_ma57_firstpass_$(job)_a0")
    isfile(joinpath(case,"COMPLETE")) || error("Reference not fully validated")
    d=deserialize(joinpath(case,"instance.jls"))
    sol=deserialize(joinpath(case,"tadmm_dispatch.jls")); us=sol[:u]; B=sol[:B]
    ix,_=fresh.control_layout(d)
    checks=String[]
    for probe in ("saved_dispatch","negative_ell","negative_import","voltage_upper","battery_power","battery_energy")
        u=deepcopy(us); b=copy(B)
        if probe=="negative_ell"; u[1][first(ix.ell)]=-0.001
        elseif probe=="negative_import"; u[1][ix.ps]=-0.01
        elseif probe=="voltage_upper"; u[1][last(ix.v)]=d[:Vmaxpu][last(d[:Nset])]^2+0.01
        elseif probe=="battery_power"; u[1][first(ix.pb)]=d[:P_B_R_pu][first(d[:Bset])]+0.01
        elseif probe=="battery_energy"; b[1,1]=d[:soc_max][first(d[:Bset])]*d[:B_R_pu][first(d[:Bset])]+0.01
        end
        a=old.physical_metrics(d,u,b); z=fresh.physical_metrics(d,u,b)
        for key in keys(a)
            isapprox(a[key],z[key];atol=1e-10,rtol=1e-13) || error("$system $probe changed $key: $(a[key]) vs $(z[key])")
        end
        probe=="saved_dispatch" || z.maxbound>=0.000999 || error("Original-bound violation was missed")
        push!(checks,probe)
    end
    fresh.physical_layout(d)===fresh.physical_layout(d) || error("Topology not cached")
    # All optional rules remain auditable; actual defaults are selected using
    # the October 20 ablation (medium keeps nonroot voltage upper bounds).
    ENV["TADMM_SCREEN"]="substation,vupper,ell,psubs"
    allscreen=fresh.tadmm_screen(d)
    shared=fresh.voltage_screen(d,d[:T])
    allscreen.keep_upper==[any(@view shared.keep_up[n,:]) for n in eachindex(d[:Nset])] || error("Shared screen drift")
    ENV["TADMM_SCREEN"]=system=="ieee2522C_1ph" ? "substation,ell" : "substation,vupper,ell"
    rho=system=="ieee2522C_1ph" ? 3800.0 : 48.0
    period=system=="ieee2522C_1ph" ? 19 : 13
    s=fresh.make_model(d,[period];local_window=true,rho=rho)
    all(!has_lower_bound(s.ell[e,period]) for e in d[:Lset]) || error("ell lower bound retained")
    all(has_lower_bound(s.v[j,period]) for j in d[:Nset] if j!=d[:substationBus]) || error("Positive voltage proof lost")
    z=[d[:B0_pu][j] for j in d[:Bset],t in 1:d[:T]]
    fresh.stage_step!(s,d,z,zeros(length(d[:Bset]),length(s.times)),rho,1)
    local_feas=maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0)
    local_feas<=1e-8 || error("Transferred local model inaccurate")
    minell=minimum(value(s.ell[e,period]) for e in d[:Lset])
    minell>=-1e-8 || error("SOC did not enforce implied current nonnegativity")
    upper_violation=maximum(value(s.v[j,period])-d[:Vmaxpu][j]^2 for j in d[:Nset])
    upper_violation<=1e-8 || error("Screened original voltage bound violated")
    results[system]=Dict("checks"=>checks,"original_metrics_match"=>true,
        "positive_voltage_preserved"=>true,"ell_lower_bounds_removed"=>length(d[:Lset]),
        "local_primal"=>local_feas,"minimum_ell"=>minell,"original_upper_violation"=>upper_violation,
        "model_sha256"=>bytes2hex(sha256(read(joinpath(newsource,"problems/mpopf/tadmm/review_pilot/model.jl")))))
    s=nothing; GC.gc()
end
open(io->TOML.print(io,results),joinpath(out,"validation.toml"),"w")
println("DDP_TRANSFER_VALIDATED ",out)
