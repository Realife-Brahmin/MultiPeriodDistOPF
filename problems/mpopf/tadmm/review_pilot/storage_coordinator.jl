# Prototype: temporal ADMM with a whole-horizon storage coordinator.
#
# Same parallel step as the overlapping-window tADMM (one branch-flow OPF per
# period), but the consensus variable is each battery's POWER at each period,
# and the coordination step projects those powers onto the set of trajectories
# that satisfy the battery dynamics and energy limits over the whole horizon.
# Every round therefore carries information across all T periods, whereas
# averaging neighbouring energy copies moves it one period per round.
#
#   x-step (parallel in t):  min cost_t + rho/2 * |pb_t - z_t + u_t|^2
#                            s.t. the branch-flow model of period t
#   z-step (one small QP):   min |z - (x + u)|^2
#                            s.t. B_t = B_{t-1} - dt z_t, energy and power limits
#   u-step:                  u += x - z
#
# Verification only: this file measures round counts and solution quality, not
# publishable wall-clock times.
include("model.jl")

function make_power_stage(d,t;rho)
    ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    loads=Set(d[:NLset]); batteries=Set(bs); ders=Set(ds)
    root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]
    screen=tadmm_screen(d); layout=physical_layout(d)
    m=Model(Ipopt.Optimizer); set_silent(m)
    set_optimizer_attribute(m,"tol",1e-9)
    set_optimizer_attribute(m,"constr_viol_tol",1e-9)
    set_optimizer_attribute(m,"compl_inf_tol",1e-8)
    set_optimizer_attribute(m,"acceptable_iter",0)
    set_optimizer_attribute(m,"max_iter",3000)
    set_optimizer_attribute(m,"bound_relax_factor",0.0)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","ma57"))
    for kv in split(get(ENV,"TADMM_IPOPT_OPTIONS",""),';';keepempty=false)
        k,v=strip.(split(kv,'=';limit=2))
        set_optimizer_attribute(m,String(k),something(tryparse(Int,v),tryparse(Float64,v),String(v)))
    end
    @variable(m,ps>=0); @variable(m,qs)
    @variable(m,p[es]); @variable(m,q[es]); @variable(m,v[ns]); @variable(m,ell[es]>=0)
    @variable(m,qd[ds]); @variable(m,pb[bs])
    "ell" in screen.flags && foreach(e->delete_lower_bound(ell[e]),es)
    screen.drop_ps && delete_lower_bound(ps)
    cb=d[:C_B]*base^2*dt
    @objective(m,Min,d[:LoadShapeCost][t]*base*dt*ps+sum((cb+rho/2)*pb[j]^2 for j in bs))
    @constraint(m,ps==sum(p[e] for e in d[:L1set]))
    @constraint(m,qs==sum(q[e] for e in d[:L1set]))
    for j in d[:Nm1set]
        e=(d[:parent][j],j)
        pl=j in loads ? d[:p_L_pu][j,t] : 0.0
        ql=j in loads ? d[:q_L_pu][j,t] : 0.0
        pd=j in ders ? d[:p_D_pu][j,t] : 0.0
        pbat=j in batteries ? pb[j] : 0.0
        qder=j in ders ? qd[j] : 0.0
        @constraint(m,sum(p[(j,k)] for k in d[:children][j];init=0.0)-p[e]+d[:rdict_pu][e]*ell[e]==pbat+pd-pl)
        qcap=haskey(get(d,:qC_pu,Dict{Int,Float64}()),j) ? d[:qC_pu][j]*v[j] : 0.0   # fixed capacitor, q = qC * v
        @constraint(m,sum(q[(j,k)] for k in d[:children][j];init=0.0)-q[e]+d[:xdict_pu][e]*ell[e]==qder-ql+qcap)
    end
    for e in es
        i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
        @constraint(m,v[j]==v[i]-2*(r*p[e]+x*q[e])+(r^2+x^2)*ell[e])
        get(ENV,"SC_AC","0")=="1" ? @constraint(m,p[e]^2+q[e]^2==v[i]*ell[e]) : @constraint(m,p[e]^2+q[e]^2<=v[i]*ell[e])   # SC_AC=1: nonconvex equality
    end
    @constraint(m,v[root]==1.05^2)
    for (n,j) in enumerate(ns)
        j==root && continue
        set_lower_bound(v[j],d[:Vminpu][j]^2)
        (!("vupper" in screen.flags) || screen.keep_upper[n]) && set_upper_bound(v[j],d[:Vmaxpu][j]^2)
    end
    for j in ds
        qmax=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
        set_lower_bound(qd[j],-qmax); set_upper_bound(qd[j],qmax); set_start_value(qd[j],0.0)
    end
    vv,pp,qq,ll,_=loadflow_initial_values(d,t)
    for (n,j) in enumerate(ns); set_start_value(v[j],vv[n]); end
    for (n,e) in enumerate(es)
        set_start_value(p[e],pp[n]); set_start_value(q[e],qq[n]); set_start_value(ell[e],ll[n])
    end
    lpos=Dict(e=>n for (n,e) in enumerate(es))
    set_start_value(ps,sum(pp[lpos[e]] for e in d[:L1set]))
    set_start_value(qs,sum(qq[lpos[e]] for e in d[:L1set]))
    for j in bs
        set_lower_bound(pb[j],-d[:P_B_R_pu][j]); set_upper_bound(pb[j],d[:P_B_R_pu][j])
        set_start_value(pb[j],0.0)
    end
    return (;m,t,ps,qs,p,q,v,ell,qd,pb,cb,vars=all_variables(m),layout)
end

function power_stage_step!(s,d,zt,ut,rho,k;start=nothing)
    if start!==nothing
        # The model was last used for a repair solve: resume from the saved iterate.
        set_start_value.(s.vars,start)
    elseif k>1
        vals=value.(s.vars); set_start_value.(s.vars,vals)
    end
    for (n,j) in enumerate(d[:Bset])
        set_objective_coefficient(s.m,s.pb[j],rho*(ut[n]-zt[n]))
        set_objective_coefficient(s.m,s.pb[j],s.pb[j],s.cb+rho/2)
    end
    optimize!(s.m)
    termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL) || error("stage $(s.t) failed: $(termination_status(s.m))")
    return [value(s.pb[j]) for j in d[:Bset]],solve_time(s.m),MOI.get(s.m,MOI.BarrierIterations())
end

function stage_dispatch(s,d)
    ix=s.layout.ix; u=zeros(s.layout.nu); t=s.t
    u[ix.ps]=value(s.ps); u[ix.qs]=value(s.qs)
    for (n,e) in enumerate(d[:Lset])
        u[ix.P[n]]=value(s.p[e]); u[ix.Q[n]]=value(s.q[e]); u[ix.ell[n]]=value(s.ell[e])
    end
    for (n,j) in enumerate(d[:Nset]); u[ix.v[n]]=value(s.v[j]); end
    for (n,j) in enumerate(d[:Bset]); u[ix.pb[n]]=value(s.pb[j]); end
    for (n,j) in enumerate(d[:Dset])
        qm=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
        u[ix.qnorm[n]]=qm>0 ? value(s.qd[j])/qm : 0.0
    end
    return u
end

# Whole-horizon battery step. Separable by battery; solved here as one sparse QP.
function make_coordinator(d;idx=eachindex(d[:Bset]))
    T=d[:T]; bs=d[:Bset][idx]; dt=d[:delta_t_h]; nb=length(bs)
    m=Model(Ipopt.Optimizer); set_silent(m)
    set_optimizer_attribute(m,"tol",1e-10)
    set_optimizer_attribute(m,"bound_relax_factor",0.0)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","ma57"))
    for kv in split(get(ENV,"TADMM_IPOPT_OPTIONS",""),';';keepempty=false)
        k,v=strip.(split(kv,'=';limit=2))
        set_optimizer_attribute(m,String(k),something(tryparse(Int,v),tryparse(Float64,v),String(v)))
    end
    @variable(m,z[1:nb,1:T]); @variable(m,b[1:nb,1:T]); scd=AffExpr[]
    for (n,j) in enumerate(bs),t in 1:T
        set_lower_bound(z[n,t],-d[:P_B_R_pu][j]); set_upper_bound(z[n,t],d[:P_B_R_pu][j])
        set_lower_bound(b[n,t],d[:soc_min][j]*d[:B_R_pu][j]); set_upper_bound(b[n,t],d[:soc_max][j]*d[:B_R_pu][j])
        set_start_value(z[n,t],0.0); set_start_value(b[n,t],d[:B0_pu][j])
        eta=sc_eta()
        if eta<1   # charging and discharging powers; a small throughput term rules out simultaneous use
            pc=@variable(m,lower_bound=0,upper_bound=d[:P_B_R_pu][j]); pd=@variable(m,lower_bound=0,upper_bound=d[:P_B_R_pu][j])
            @constraint(m,z[n,t]==pd-pc)
            @constraint(m,b[n,t]==(t==1 ? d[:B0_pu][j] : b[n,t-1])+dt*(eta*pc-pd/eta))
            push!(scd,pc+pd)
        else
            @constraint(m,b[n,t]==(t==1 ? d[:B0_pu][j] : b[n,t-1])-dt*z[n,t])
        end
    end
    @objective(m,Min,sum(z[n,t]^2 for n in 1:nb,t in 1:T)+1e-6*sum(scd;init=0.0))
    return (;m,z,b,vars=all_variables(m))
end

function coordinator_step!(c,target,k;start=nothing)
    if start!==nothing
        set_start_value.(c.vars,start)
    elseif k>1
        vals=value.(c.vars); set_start_value.(c.vars,vals)
    end
    for i in eachindex(target); set_objective_coefficient(c.m,c.z[i],-2*target[i]); end
    optimize!(c.m)
    termination_status(c.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL) || error("coordinator failed: $(termination_status(c.m))")
    return value.(c.z),solve_time(c.m)
end

function run_storage_coordinator(d;rho,maxiter=300,reference=NaN,no_tol=1e-5,alpha=1.0,adapt=false,label="sc",io=nothing,rho_end=rho,ramp=0,spec=string(rho),check_repair=false)
    T=d[:T]; bs=d[:Bset]; nb=length(bs); dt=d[:delta_t_h]
    stages=[make_power_stage(d,t;rho) for t in 1:T]; coord=make_coordinator(d)
    z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T)
    no_iteration=0; strict=0; ideal=0.0; serial=0.0; k=0; metrics=nothing; rmax=Inf; smax=Inf
    us=Vector{Vector{Float64}}(undef,T)
    repair=check_repair ? [make_power_stage(d,t;rho=0.0) for t in 1:T] : nothing; repair_no=0
    rho_start=rho
    for iteration in 1:maxiter
        k=iteration; slowest=0.0; iters=0
        if ramp>0
            # Geometric schedule from rho_start to rho_end over `ramp` rounds, then held.
            next=rho_start*(rho_end/rho_start)^(min(k-1,ramp)/ramp)
            u.*=rho/next; rho=next
        end
        for t in 1:T
            x[:,t],st,it=power_stage_step!(stages[t],d,view(z,:,t),view(u,:,t),rho,k)
            slowest=max(slowest,st); serial+=st; iters=max(iters,it)
        end
        zold=copy(z)
        xh=alpha.*x.+(1-alpha).*zold
        z,ct=coordinator_step!(coord,xh.+u,k)
        u.+=xh.-z
        ideal+=slowest+ct; serial+=ct
        r=x.-z; s=rho.*(z.-zold)
        rmax=maximum(abs,r); smax=maximum(abs,s)
        rr=sqrt(sum(abs2,r)/length(r)); ss=sqrt(sum(abs2,s)/length(s))
        for t in 1:T; us[t]=stage_dispatch(stages[t],d); end
        bb=zeros(nb,T)
        for n in 1:nb
            e=d[:B0_pu][bs[n]]
            for t in 1:T; e-=dt*x[n,t]; bb[n,t]=e; end
        end
        metrics=physical_metrics(d,us,bb)
        gap=abs(metrics.objective-reference)/max(1.0,abs(reference))
        if no_iteration==0 && gap<=1e-3 && metrics.original_primal<=no_tol
            no_iteration=k
            @printf("SC_NO label=%s spec=%s rho=%g k=%d gap=%.3e original=%.3e ideal_s=%.2f\n",label,spec,rho,k,gap,metrics.original_primal,ideal)
        end
        # Repair check: fix every battery at the coordinator's schedule (which
        # satisfies the energy dynamics and limits) and re-solve each period's
        # network. If all periods solve, the result violates nothing.
        repair_gap=NaN; repair_primal=NaN; repair_s=NaN
        if repair!==nothing
            repair_s=0.0; solved=true
            for t in 1:T
                ps=repair[t]
                for (n,j) in enumerate(bs); fix(ps.pb[j],z[n,t];force=true); end
                optimize!(ps.m); repair_s=max(repair_s,solve_time(ps.m))
                solved &= termination_status(ps.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
                solved || break
            end
            if solved
                zb=zeros(nb,T)
                for n in 1:nb
                    e=d[:B0_pu][bs[n]]
                    for t in 1:T; e-=dt*z[n,t]; zb[n,t]=e; end
                end
                pm=physical_metrics(d,[stage_dispatch(repair[t],d) for t in 1:T],zb)
                repair_gap=abs(pm.objective-reference)/max(1.0,abs(reference)); repair_primal=pm.original_primal
                if repair_no==0 && repair_gap<=1e-3 && repair_primal<=no_tol
                    repair_no=k
                    @printf("SC_REPAIR_NO label=%s spec=%s rho=%g k=%d gap=%.3e original=%.3e ideal_s=%.2f repair_s=%.2f\n",label,spec,rho,k,repair_gap,repair_primal,ideal,repair_s)
                end
            end
        end
        io===nothing || println(io,join((label,spec,k,metrics.objective,gap,metrics.original_primal,rr,rmax,ss,smax,ideal,serial,iters,rho,slowest,ct,repair_gap,repair_primal,repair_s),','))
        if k<=5 || k%10==0
            @printf("SC_ITER label=%s rho=%g k=%d gap=%.3e original=%.3e rmax=%.3e smax=%.3e\n",label,rho,k,gap,metrics.original_primal,rmax,smax); flush(stdout)
        end
        if rmax<=2.5e-7 && smax<=1e-4 && metrics.original_primal<=5e-7
            strict=k; break
        end
        if adapt && k%5==0
            old=rho
            rr>10ss && (rho*=2); ss>10rr && (rho/=2)
            old==rho || (u.*=old/rho)
        end
    end
    gap=abs(metrics.objective-reference)/max(1.0,abs(reference))
    @printf("SC_DONE label=%s spec=%s rho_final=%g rounds=%d strict=%d no=%d repair_no=%d gap=%.3e original=%.3e rmax=%.3e ideal_s=%.2f serial_s=%.2f\n",label,spec,rho,k,strict,no_iteration,repair_no,gap,metrics.original_primal,rmax,ideal,serial)
    flush(stdout)
    return (;k,strict,no_iteration,gap,metrics,ideal,serial)
end

function main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    profile=apply_review_cb!(d,system)
    no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    reffile=joinpath(out,"reference_$(profile).txt")
    if isfile(reffile)
        reference=parse(Float64,read(reffile,String))
    else
        central=make_model(d,collect(1:T)); tick=time_ns(); solve_checked!(central)
        cu,cbb=dispatch([central],d); cm=physical_metrics(d,cu,cbb); reference=cm.objective
        @printf("CENTRAL objective=%.12f primal=%.3e iterations=%d wall_s=%.2f\n",reference,cm.original_primal,MOI.get(central.m,MOI.BarrierIterations()),(time_ns()-tick)/1e9)
        write(reffile,string(reference))
    end
    # Penalty at which the day's price spread moves a median battery by its rating in one round.
    prated=sort([d[:P_B_R_pu][j] for j in d[:Bset]]); pmed=prated[cld(length(prated),2)]
    rho_unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*d[:kVA_B]*d[:delta_t_h]/pmed
    @printf("CASE system=%s T=%d C_B=%.6g profile=%s batteries=%d reference=%.12f rho_unit=%.6g median_P_rated_pu=%.6g\n",system,T,d[:C_B],profile,length(d[:Bset]),reference,rho_unit,pmed); flush(stdout)
    maxiter=parse(Int,get(ENV,"SC_MAX_ITER","300"))
    alpha=parse(Float64,get(ENV,"SC_ALPHA","1.0")); adapt=get(ENV,"SC_ADAPT","0")=="1"
    tag=get(ENV,"SC_TAG","sweep")
    open(joinpath(out,"trace_$(profile)_$(tag).csv"),"w") do io
        println(io,"label,spec,iteration,objective,gap,original_primal,r_rms,r_max,s_rms,s_max,ideal_parallel_estimate_s,serial_s,max_local_ipopt_iterations,rho,slowest_stage_s,coordinator_s,repair_gap,repair_primal,repair_slowest_stage_s")
        # SC_RHO entries: a number, or "a:b:n" for a geometric ramp from a to b over n rounds.
        # With SC_RHO_UNITS=1 the numbers multiply rho_unit.
        scale=get(ENV,"SC_RHO_UNITS","0")=="1" ? rho_unit : 1.0
        for spec in split(get(ENV,"SC_RHO","100"),',')
            parts=split(spec,':')
            rho=scale*parse(Float64,parts[1])
            rho_end=length(parts)==3 ? scale*parse(Float64,parts[2]) : rho
            ramp=length(parts)==3 ? parse(Int,parts[3]) : 0
            run_storage_coordinator(d;rho,maxiter,reference,no_tol,alpha,adapt,label=tag,io,rho_end,ramp,spec=String(spec),check_repair=get(ENV,"SC_REPAIR","0")=="1")
        end
    end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
