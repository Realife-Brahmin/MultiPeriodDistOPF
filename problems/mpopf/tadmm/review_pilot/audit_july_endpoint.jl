# Saved-point reconstruction plus two controlled next-step local solves.
using Serialization, SHA, TOML
source,case,out=ARGS
include(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))
mkpath(out)
d=deserialize(joinpath(case,"instance.jls")); state=deserialize(joinpath(case,"exact/trace_checkpoint.jls"))
@assert bytes2hex(sha256(read(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))))==state[:model_sha256]
meta=TOML.parsefile(joinpath(case,"comparison.toml")); ENV["TADMM_SCREEN"]=meta["screening"]
ENV["TADMM_LOCAL_SOLVE_POLICY"]="exact"
ix,nu=control_layout(d); us=Vector{Vector{Float64}}(undef,d[:T]); worst=(error=-Inf,period=0,battery=0,time=0)
for t in 1:d[:T]
    s=make_model(d,[t];local_window=true,rho=state[:rho]); vals=state[:stage_primals][t]
    @assert length(vals)==length(s.vars) && all(index(v).value==n for (n,v) in enumerate(s.vars))
    at(v)=vals[index(v).value]
    u=zeros(nu); u[ix.ps]=at(s.ps[t]); u[ix.qs]=at(s.qs[t])
    for (n,e) in enumerate(d[:Lset]); u[ix.P[n]]=at(s.p[e,t]);u[ix.Q[n]]=at(s.q[e,t]);u[ix.ell[n]]=at(s.ell[e,t]);end
    for (n,j) in enumerate(d[:Nset]);u[ix.v[n]]=at(s.v[j,t]);end
    for (n,j) in enumerate(d[:Bset]);u[ix.pb[n]]=at(s.pb[j,t]);end
    for (n,j) in enumerate(d[:Dset])
        qm=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2));u[ix.qnorm[n]]=qm>0 ? at(s.qd[j,t])/qm : 0.0
    end
    for (n,j) in enumerate(d[:Bset]),tt in s.times
        err=abs(at(s.b[j,tt])-state[:z][n,tt])
        if err>worst.error;global worst=(error=err,period=t,battery=j,time=tt);end
    end
    us[t]=u;s=nothing;GC.gc()
end
m=physical_metrics(d,us,state[:z]); lines=readlines(joinpath(case,"exact/trace.csv"))
header=Dict(k=>n for (n,k) in enumerate(split(first(lines),','))); lastrow=split(last(lines),',')
@assert isapprox(m.objective,parse(Float64,lastrow[header["objective"]]);rtol=1e-13,atol=1e-9)
@assert isapprox(m.original_primal,parse(Float64,lastrow[header["original_primal"]]);atol=1e-10,rtol=1e-12)
dyn=0.0;driftbound=0.0;energy=[d[:B0_pu][j] for j in d[:Bset]]
for t in 1:d[:T],(n,j) in enumerate(d[:Bset])
    prev=t==1 ? d[:B0_pu][j] : state[:z][n,t-1]
    global dyn=max(dyn,abs(state[:z][n,t]-prev+d[:delta_t_h]*us[t][ix.pb[n]]))
    energy[n]-=d[:delta_t_h]*us[t][ix.pb[n]]
    global driftbound=max(driftbound,d[:soc_min][j]*d[:B_R_pu][j]-energy[n],energy[n]-d[:soc_max][j]*d[:B_R_pu][j])
end
result=Dict{String,Any}("saved_iteration"=>state[:iteration],"objective_trace_matches"=>true,"primal_trace_matches"=>true,
    "additional_optimization_solves"=>2,"probe_next_iteration"=>state[:iteration]+1,
    "original_primal"=>m.original_primal,"max_equation"=>m.maxeq,"max_bound"=>m.maxbound,"max_cone"=>m.maxcone,
    "battery_dynamic_equation"=>dyn,"integrated_soc_bound"=>driftbound,"worst_copy_period"=>worst.period,
    "worst_copy_battery"=>worst.battery,"worst_copy_time"=>worst.time,"worst_copy_error"=>worst.error)
copies=Dict{String,Matrix{Float64}}()
for variant in ("default","absolute_dual")
    s=make_model(d,[worst.period];local_window=true,rho=state[:rho])
    set_start_value.(s.vars,state[:stage_primals][worst.period])
    for (n,j) in enumerate(d[:Bset]),(a,t) in enumerate(s.times)
        set_objective_coefficient(s.m,s.b[j,t],state[:rho]*(state[:dual][worst.period][n,a]-state[:z][n,t]))
        set_objective_coefficient(s.m,s.b[j,t],s.b[j,t],state[:rho]/2)
    end
    set_optimizer_attribute(s.m,"max_iter",150)
    set_optimizer_attribute(s.m,"output_file",joinpath(out,variant*".log"))
    set_optimizer_attribute(s.m,"file_print_level",5)
    variant=="absolute_dual" && set_optimizer_attribute(s.m,"dual_inf_tol",1e-8)
    optimize!(s.m)
    status=string(termination_status(s.m)); @assert has_values(s.m)
    copies[variant]=[value(s.b[j,t]) for j in d[:Bset],t in s.times]
    primal=maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0)
    maxdual=has_duals(s.m) ? maximum(abs(JuMP.dual(c)) for c in all_constraints(s.m;include_variable_in_set_constraints=true);init=0.0) : NaN
    result[variant]=Dict("status"=>status,"ipopt_iterations"=>MOI.get(s.m,MOI.BarrierIterations()),
        "local_primal"=>primal,"max_raw_multiplier"=>maxdual,"nlp_objective"=>objective_value(s.m))
    s=nothing;GC.gc()
end
result["max_energy_difference_default_vs_absolute_dual"]=maximum(abs.(copies["default"].-copies["absolute_dual"]))
open(io->TOML.print(io,result),joinpath(out,"audit.toml"),"w")
println("JULY_ENDPOINT_AUDITED ",out)
