# Probe for reviewer items AE-4 and R2-3: a battery with separate charging and
# discharging powers and efficiencies. In the storage-coordinated scheme the
# battery model appears only in Update II, so the period OPFs are untouched.
#
#   julia --project=<conic_env> battery_efficiency_probe.jl system T eta [iterations]
#
# Model: P_B = P_dis - P_ch, both nonnegative, and
#   B_t = B_{t-1} + dt*(eta*P_ch - P_dis/eta).
# Simultaneous charging and discharging is not forbidden (that would be
# non-convex); it is checked on the result. Serial run, iteration counts only.
# The centralized reference is the explicit cone program solved by Clarabel.
include("storage_coordinator.jl")
using Clarabel

function central_with_efficiency(d,eta)
    T=d[:T]; ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    loads=Set(d[:NLset]); batteries=Set(bs); ders=Set(ds); root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]
    m=Model(optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,"max_threads"=>1))
    @variable(m,ps[1:T]>=0); @variable(m,p[es,1:T]); @variable(m,q[es,1:T]); @variable(m,v[ns,1:T]); @variable(m,ell[es,1:T]>=0)
    @variable(m,qd[ds,1:T]); @variable(m,pc[bs,1:T]>=0); @variable(m,pdis[bs,1:T]>=0); @variable(m,b[bs,1:T])
    pb=Dict((j,t)=>pdis[j,t]-pc[j,t] for j in bs,t in 1:T)
    @objective(m,Min,sum(d[:LoadShapeCost][t]*base*dt*ps[t] for t in 1:T)+sum(d[:C_B]*base^2*dt*(pdis[j,t]-pc[j,t])^2 for j in bs,t in 1:T))
    for t in 1:T
        @constraint(m,ps[t]==sum(p[e,t] for e in d[:L1set]))
        for j in d[:Nm1set]
            e=(d[:parent][j],j)
            pl=j in loads ? d[:p_L_pu][j,t] : 0.0; ql=j in loads ? d[:q_L_pu][j,t] : 0.0; pd=j in ders ? d[:p_D_pu][j,t] : 0.0
            @constraint(m,sum(p[(j,k),t] for k in d[:children][j];init=0.0)-p[e,t]+d[:rdict_pu][e]*ell[e,t]==(j in batteries ? pb[(j,t)] : 0.0)+pd-pl)
            @constraint(m,sum(q[(j,k),t] for k in d[:children][j];init=0.0)-q[e,t]+d[:xdict_pu][e]*ell[e,t]==(j in ders ? qd[j,t] : 0.0)-ql)
        end
        for e in es
            i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
            @constraint(m,v[j,t]==v[i,t]-2*(r*p[e,t]+x*q[e,t])+(r^2+x^2)*ell[e,t])
            @constraint(m,[0.5*v[i,t],ell[e,t],p[e,t],q[e,t]] in RotatedSecondOrderCone())
        end
        @constraint(m,v[root,t]==1.05^2)
        for j in ns
            j==root && continue
            set_lower_bound(v[j,t],d[:Vminpu][j]^2); set_upper_bound(v[j,t],d[:Vmaxpu][j]^2)
        end
        for j in ds
            qmax=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2)); set_lower_bound(qd[j,t],-qmax); set_upper_bound(qd[j,t],qmax)
        end
    end
    for j in bs, t in 1:T
        set_upper_bound(pc[j,t],d[:P_B_R_pu][j]); set_upper_bound(pdis[j,t],d[:P_B_R_pu][j])
        set_lower_bound(b[j,t],d[:soc_min][j]*d[:B_R_pu][j]); set_upper_bound(b[j,t],d[:soc_max][j]*d[:B_R_pu][j])
        @constraint(m,b[j,t]==(t==1 ? d[:B0_pu][j] : b[j,t-1])+dt*(eta*pc[j,t]-pdis[j,t]/eta))
    end
    optimize!(m)
    both=maximum(min(value(pc[j,t]),value(pdis[j,t])) for j in bs,t in 1:T)
    return (;objective=objective_value(m),status=string(termination_status(m)),simultaneous=both)
end

# Update II with efficiencies: the schedule closest to the target that the
# batteries can follow. The small charge on throughput picks, among equally close
# schedules, one that does not charge and discharge at once.
function make_coordinator_eff(d,eta;throughput=1e-6)
    T=d[:T]; bs=d[:Bset]; dt=d[:delta_t_h]; nb=length(bs)
    m=Model(Ipopt.Optimizer); set_silent(m)
    set_optimizer_attribute(m,"tol",1e-10); set_optimizer_attribute(m,"bound_relax_factor",0.0)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","ma57"))
    for kv in split(get(ENV,"TADMM_IPOPT_OPTIONS",""),';';keepempty=false)
        k,v=strip.(split(kv,'=';limit=2)); set_optimizer_attribute(m,String(k),something(tryparse(Int,v),tryparse(Float64,v),String(v)))
    end
    @variable(m,z[1:nb,1:T]); @variable(m,pc[1:nb,1:T]>=0); @variable(m,pdis[1:nb,1:T]>=0); @variable(m,b[1:nb,1:T])
    for (n,j) in enumerate(bs),t in 1:T
        set_upper_bound(pc[n,t],d[:P_B_R_pu][j]); set_upper_bound(pdis[n,t],d[:P_B_R_pu][j])
        set_lower_bound(b[n,t],d[:soc_min][j]*d[:B_R_pu][j]); set_upper_bound(b[n,t],d[:soc_max][j]*d[:B_R_pu][j])
        set_start_value(b[n,t],d[:B0_pu][j]); set_start_value(pc[n,t],0.0); set_start_value(pdis[n,t],0.0); set_start_value(z[n,t],0.0)
        @constraint(m,z[n,t]==pdis[n,t]-pc[n,t])
        @constraint(m,b[n,t]==(t==1 ? d[:B0_pu][j] : b[n,t-1])+dt*(eta*pc[n,t]-pdis[n,t]/eta))
    end
    @objective(m,Min,sum(z[n,t]^2 for n in 1:nb,t in 1:T)+throughput*sum(pc[n,t]+pdis[n,t] for n in 1:nb,t in 1:T))
    return (;m,z,pc,pdis,b,vars=all_variables(m))
end

function main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); eta=parse(Float64,ARGS[3]); maxiter=length(ARGS)>3 ? parse(Int,ARGS[4]) : 60
    d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    apply_review_cb!(d,system); bs=d[:Bset]; nb=length(bs); dt=d[:delta_t_h]; base=d[:kVA_B]
    ref=central_with_efficiency(d,eta)
    @printf("EFF_CENTRAL system=%s T=%d eta=%.3f status=%s objective=%.9f largest simultaneous charge and discharge=%.3e\n",system,T,eta,ref.status,ref.objective,ref.simultaneous)
    prated=sort([d[:P_B_R_pu][j] for j in bs]); unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*base*dt/prated[cld(nb,2)]
    rho=2unit; stages=[make_power_stage(d,t;rho) for t in 1:T]; repair=[make_power_stage(d,t;rho=0.0) for t in 1:T]; coord=make_coordinator_eff(d,eta)
    z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T); first_ok=0
    for k in 1:maxiter
        next=2unit*(0.0625/2)^(min(k-1,20)/20); u.*=rho/next; rho=next
        for t in 1:T; x[:,t],_,_=power_stage_step!(stages[t],d,z[:,t],u[:,t],rho,k); end
        if k>1
            vals=value.(coord.vars); set_start_value.(coord.vars,vals)
        end
        target=x.+u
        for i in eachindex(target); set_objective_coefficient(coord.m,coord.z[i],-2*target[i]); end
        optimize!(coord.m); termination_status(coord.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL) || error("battery step failed")
        zold=z; z=value.(coord.z); u.+=x.-z
        both=maximum(min.(value.(coord.pc),value.(coord.pdis)))
        r_rms=sqrt(sum(abs2,x.-z)/length(z))
        # Repair with the batteries at the battery-side schedule; cost from the period solves.
        cost=0.0; strict=true
        for t in 1:T
            s=repair[t]
            for (n,j) in enumerate(bs); fix(s.pb[j],z[n,t];force=true); end
            optimize!(s.m); strict&=termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
            strict || break
            cost+=d[:LoadShapeCost][t]*base*dt*value(s.ps)+d[:C_B]*base^2*dt*sum(abs2,z[:,t])
        end
        gap=strict ? abs(cost-ref.objective)/abs(ref.objective) : NaN
        strict && first_ok==0 && gap<=1e-3 && (first_ok=k)
        (k<=3 || k%5==0 || k==first_ok) && (@printf("EFF_ITER k=%d r_rms=%.3e repaired_gap=%.3e repair_strict=%s simultaneous=%.3e\n",k,r_rms,gap,strict,both); flush(stdout))
    end
    @printf("EFF_DONE system=%s T=%d eta=%.3f first iteration with a repaired dispatch within 0.1%% of centralized=%d\n",system,T,eta,first_ok)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
