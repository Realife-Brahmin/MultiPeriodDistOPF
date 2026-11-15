# Why does the repair solve stop at Ipopt's acceptable level in some periods,
# and which remedy gives a strict solve? Runs the storage-coordinated iteration
# serially to a given count, repairs every period, and for each period that is
# not strictly solved tries a few variants.
include("storage_coordinator.jl")
system=ARGS[1]; T=parse(Int,ARGS[2]); stop_at=parse(Int,ARGS[3])
d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
apply_review_cb!(d,system)
bs=d[:Bset]; nb=length(bs)
prated=sort([d[:P_B_R_pu][j] for j in bs]); unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*d[:kVA_B]*d[:delta_t_h]/prated[cld(nb,2)]
stages=[make_power_stage(d,t;rho=2unit) for t in 1:T]; coord=make_coordinator(d)
z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T); rho=2unit
for k in 1:stop_at
    global rho,z
    next=2unit*(0.0625/2)^(min(k-1,20)/20); u.*=rho/next; rho=next
    for t in 1:T; x[:,t],_,_=power_stage_step!(stages[t],d,view(z,:,t),view(u,:,t),rho,k); end
    z,_=coordinator_step!(coord,x.+u,k); u.+=x.-z
end
@printf("after %d iterations: max |x-z| = %.3e pu\n",stop_at,maximum(abs,x.-z))

function attempt(t,label;from_stage=false,options=Pair{String,Any}[],band=0.0)
    s=make_power_stage(d,t;rho=0.0)
    for (k,v) in options; set_optimizer_attribute(s.m,k,v); end
    from_stage && set_start_value.(s.vars,value.(stages[t].vars))
    for (n,j) in enumerate(bs)
        if band>0
            set_lower_bound(s.pb[j],max(z[n,t]-band,-d[:P_B_R_pu][j])); set_upper_bound(s.pb[j],min(z[n,t]+band,d[:P_B_R_pu][j]))
        else
            fix(s.pb[j],z[n,t];force=true)
        end
    end
    optimize!(s.m)
    viol=has_values(s.m) ? maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0) : NaN
    cone=has_values(s.m) ? maximum(value(s.v[e[1]])*value(s.ell[e])-value(s.p[e])^2-value(s.q[e])^2 for e in d[:Lset]) : NaN
    @printf("  period %2d %-34s status=%-28s iterations=%4d violation=%.2e largest cone slack=%.2e\n",t,label,raw_status(s.m),MOI.get(s.m,MOI.BarrierIterations()),viol,cone)
    return termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
end

for t in 1:T
    attempt(t,"as now (fixed, cold start)") && continue
    attempt(t,"start from the period's own solve";from_stage=true)
    attempt(t,"adaptive barrier";options=["mu_strategy"=>"adaptive"])
    attempt(t,"own-solve start + adaptive";from_stage=true,options=["mu_strategy"=>"adaptive"])
    attempt(t,"batteries within 1e-7 of schedule";band=1e-7)
    attempt(t,"band 1e-7 + own-solve start";band=1e-7,from_stage=true)
end
println("DIAGNOSE_DONE")
