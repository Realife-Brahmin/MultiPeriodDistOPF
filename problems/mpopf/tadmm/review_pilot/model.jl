# A common physical model for centralized and overlapping-window tADMM solves.
# Each stage owns its network and cost. Copies of neighboring battery energies
# and powers enforce local dynamics; consensus identifies every energy copy.
using JuMP, Ipopt, LinearAlgebra, Printf, Serialization
const MOI = JuMP.MOI
const SHARED = normpath(joinpath(@__DIR__, "..", "..", "spar_ddp", "examples", "power_system"))
include(joinpath(SHARED, "terminal_soc_penalty.jl"))
include(joinpath(SHARED, "loadflow_start.jl"))
include(joinpath(SHARED, "control_layout.jl"))
include(joinpath(SHARED, "voltage_screening.jl"))
include("experiment_settings.jl")

mutable struct StageReuse
    rho::Float64
    warm_start::String
    objective_update::String
    energy_vars::Vector{VariableRef}
    coefficients::Vector{Float64}
    constraints::Vector{ConstraintRef}
    quadratic_updates::Int
end

# Reuse SPAR-DDP's proved, data-only screen; never infer an active set from a
# previous optimum. Cache topology and screening once per exported instance.
function tadmm_screen(d)
    spec=get(ENV,"TADMM_SCREEN","substation,ell")
    cache=get!(d,:_tadmm_screen_cache,Dict{String,Any}())
    return get!(cache,spec) do
        flags=Set(String.(strip.(split(spec,',';keepempty=false))))
        isempty(setdiff(flags,Set(["substation","ell","vupper","psubs"]))) || error("Unknown TADMM_SCREEN item")
        push!(flags,"substation") # fixed-root bounds are always redundant
        # Positive voltage is needed for SOC => ell>=0. Keep all nonroot
        # lower-voltage bounds, as in the validated DDP screening studies.
        all(j->d[:Vminpu][j]>0,d[:Nset]) || error("ell screening requires positive voltage")
        vs=("vupper" in flags || "psubs" in flags) ? voltage_screen(d,d[:T]) : nothing
        keep_upper=vs===nothing ? trues(length(d[:Nset])) : [any(@view vs.keep_up[n,:]) for n in eachindex(d[:Nset])]
        drop_ps="psubs" in flags && minimum(vs.psubs_min)>0
        (;flags,keep_upper,drop_ps,psubs_min=vs===nothing ? NaN : minimum(vs.psubs_min))
    end
end

function physical_layout(d)
    return get!(d,:_tadmm_physical_layout) do
        ix,nu=control_layout(d)
        (;ix,nu,bp=Dict(j=>n for (n,j) in enumerate(d[:Nset])),
            ep=Dict(e=>n for (n,e) in enumerate(d[:Lset])),
            bat=Dict(j=>n for (n,j) in enumerate(d[:Bset])),
            der=Dict(j=>n for (n,j) in enumerate(d[:Dset])),
            root=d[:substationBus],loads=Set(d[:NLset]))
    end
end

function make_model(d, owned; local_window=false, rho=0.0)
    T = d[:T]; ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    loads=Set(d[:NLset]); batteries=Set(bs); ders=Set(ds)
    root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]
    times = local_window ? collect(max(1,only(owned)-1):min(T,only(owned)+1)) : collect(1:T)
    screen=tadmm_screen(d)
    layout=physical_layout(d)
    # Only dynamics owned by t and its successor are needed in each window.
    dyn = local_window ? collect(only(owned):min(T,only(owned)+1)) : collect(1:T)
    m=Model(Ipopt.Optimizer)
    set_silent(m)
    set_optimizer_attribute(m,"tol",1e-9)
    set_optimizer_attribute(m,"constr_viol_tol",1e-9)
    # tol is scaled: huge multipliers can otherwise permit inaccurate proximal
    # solves. Enforce an explicit absolute complementarity tolerance as well.
    set_optimizer_attribute(m,"compl_inf_tol",1e-8)
    set_optimizer_attribute(m,"acceptable_iter",0)
    local_policy=local_window ? get(ENV,"TADMM_LOCAL_SOLVE_POLICY","exact") : "exact"
    local_policy in ("exact","cap30","warm_retry","early_cap30") || error("Unknown local-solve policy")
    report_setting=get(ENV,"TADMM_LOCAL_PRIMAL_REPORT","1")
    report_setting in ("0","1") || error("TADMM_LOCAL_PRIMAL_REPORT must be 0 or 1")
    local_primal_report=report_setting=="1"
    set_optimizer_attribute(m,"max_iter",local_policy=="cap30" ? 30 : 3000)
    set_optimizer_attribute(m,"bound_relax_factor",0.0)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","ma57"))
    for kv in split(get(ENV,"TADMM_IPOPT_OPTIONS",""),';';keepempty=false)
        k,v=strip.(split(kv,'=';limit=2))
        val=something(tryparse(Int,v),tryparse(Float64,v),String(v))
        set_optimizer_attribute(m,String(k),val)
    end
    @variable(m, ps[owned]>=0)
    @variable(m, qs[owned])
    @variable(m, p[es,owned])
    @variable(m, q[es,owned])
    @variable(m, v[ns,owned])
    @variable(m, ell[es,owned]>=0)
    @variable(m, qd[ds,owned])
    @variable(m, pb[bs,dyn])
    @variable(m, b[bs,times])
    if "ell" in screen.flags
        for e in es,t in owned; delete_lower_bound(ell[e,t]); end
    end
    if screen.drop_ps
        for t in owned; delete_lower_bound(ps[t]); end
    end
    @expression(m,cost,
        sum(d[:LoadShapeCost][t]*base*dt*ps[t] for t in owned)+
        sum(d[:C_B]*base^2*dt*pb[j,t]^2 for j in bs,t in owned))
    @objective(m,Min,cost+rho/2*sum(b[j,t]^2 for j in bs,t in times))
    for t in owned
        # Root-device convention matches the shared SPAR-DDP reference.
        @constraint(m,ps[t]==sum(p[e,t] for e in d[:L1set]))
        @constraint(m,qs[t]==sum(q[e,t] for e in d[:L1set]))
        for j in d[:Nm1set]
            e=(d[:parent][j],j)
            pl=j in loads ? d[:p_L_pu][j,t] : 0.0
            ql=j in loads ? d[:q_L_pu][j,t] : 0.0
            pd=j in ders ? d[:p_D_pu][j,t] : 0.0
            pbat=j in batteries ? pb[j,t] : 0.0
            qder=j in ders ? qd[j,t] : 0.0
            @constraint(m,sum(p[(j,k),t] for k in d[:children][j];init=0.0)-p[e,t]+d[:rdict_pu][e]*ell[e,t]==pbat+pd-pl)
            @constraint(m,sum(q[(j,k),t] for k in d[:children][j];init=0.0)-q[e,t]+d[:xdict_pu][e]*ell[e,t]==qder-ql)
        end
        for e in es
            i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
            @constraint(m,v[j,t]==v[i,t]-2*(r*p[e,t]+x*q[e,t])+(r^2+x^2)*ell[e,t])
            @constraint(m,p[e,t]^2+q[e,t]^2<=v[i,t]*ell[e,t])
        end
        @constraint(m,v[root,t]==1.05^2)
        for (n,j) in enumerate(ns)
            if j==root
                # The equality fixes root voltage at its permitted maximum.
                # An additional upper bound is redundant and yielded opposing
                # ~1e9 multipliers, distorting Ipopt's scaled stopping test.
                d[:Vminpu][j]<=1.05<=d[:Vmaxpu][j] || error("Source voltage outside original bounds")
                continue
            end
            set_lower_bound(v[j,t],d[:Vminpu][j]^2)
            if !("vupper" in screen.flags) || screen.keep_upper[n]
                set_upper_bound(v[j,t],d[:Vmaxpu][j]^2)
            end
        end
        for j in ds
            qmax=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
            set_lower_bound(qd[j,t],-qmax); set_upper_bound(qd[j,t],qmax)
        end
        vv,pp,qq,ll,_=loadflow_initial_values(d,t)
        for (n,j) in enumerate(ns); set_start_value(v[j,t],vv[n]); end
        for (n,e) in enumerate(es)
            set_start_value(p[e,t],pp[n]); set_start_value(q[e,t],qq[n]); set_start_value(ell[e,t],ll[n])
        end
        lpos=Dict(e=>n for (n,e) in enumerate(es))
        set_start_value(ps[t],sum(pp[lpos[e]] for e in d[:L1set]))
        set_start_value(qs[t],sum(qq[lpos[e]] for e in d[:L1set]))
        for j in ds; set_start_value(qd[j,t],0.0); end
    end
    for j in bs
        for t in times
            set_lower_bound(b[j,t],d[:soc_min][j]*d[:B_R_pu][j])
            set_upper_bound(b[j,t],d[:soc_max][j]*d[:B_R_pu][j])
            set_start_value(b[j,t],d[:B0_pu][j])
        end
        for t in dyn
            set_lower_bound(pb[j,t],-d[:P_B_R_pu][j]); set_upper_bound(pb[j,t],d[:P_B_R_pu][j])
            set_start_value(pb[j,t],0.0)
            prev=t==1 ? d[:B0_pu][j] : b[j,t-1]
            @constraint(m,b[j,t]==prev-dt*pb[j,t])
        end
    end
    vars=all_variables(m); initial=start_value.(vars)
    warm_start=local_window ? get(ENV,"TADMM_WARM_START","primal") : "primal"
    objective_update=get(ENV,"TADMM_OBJECTIVE_UPDATE","legacy")
    warm_start in ("primal","primal_dual") || error("Unknown warm-start policy")
    objective_update in ("legacy","batched") || error("Unknown objective update policy")
    energy_vars=local_window ? VariableRef[b[j,t] for t in times for j in bs] : VariableRef[]
    constraints=warm_start=="primal_dual" ? all_constraints(m;include_variable_in_set_constraints=true) : ConstraintRef[]
    reuse=StageReuse(rho,warm_start,objective_update,energy_vars,zeros(length(energy_vars)),constraints,0)
    return (;m,ps,qs,p,q,v,ell,qd,pb,b,cost,times,dyn,owned,vars,initial,layout,local_policy,local_primal_report,reuse)
end

function solve_checked!(s; reuse=false,allow_iteration_limit=false)
    if reuse
        vals=value.(s.vars)
        set_start_value.(s.vars,vals)
    end
    optimize!(s.m)
    status=termination_status(s.m)
    solved=status in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
    truncated=allow_iteration_limit && status==MOI.ITERATION_LIMIT && has_values(s.m)
    solved || truncated || error("subproblem failed: $status")
    truncated && !all(isfinite,value.(s.vars)) && error("Nonfinite truncated local primal values")
    return solve_time(s.m)
end

function dispatch_period(s,d,t)
        ix=s.layout.ix; u=zeros(s.layout.nu); bb=zeros(length(d[:Bset]))
        u[ix.ps]=value(s.ps[t]); u[ix.qs]=value(s.qs[t])
        for (n,e) in enumerate(d[:Lset])
            u[ix.P[n]]=value(s.p[e,t]); u[ix.Q[n]]=value(s.q[e,t]); u[ix.ell[n]]=value(s.ell[e,t])
        end
        for (n,j) in enumerate(d[:Nset]); u[ix.v[n]]=value(s.v[j,t]); end
        for (n,j) in enumerate(d[:Bset])
            u[ix.pb[n]]=value(s.pb[j,t]); bb[n]=value(s.b[j,t])
        end
        for (n,j) in enumerate(d[:Dset])
            qm=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
            u[ix.qnorm[n]]=qm>0 ? value(s.qd[j,t])/qm : 0.0
        end
    return u,bb
end

function dispatch(stages,d; z=nothing)
    T=d[:T]; us=Vector{Vector{Float64}}(undef,T); bb=zeros(length(d[:Bset]),T)
    for t in 1:T
        s=length(stages)==1 ? only(stages) : stages[t]
        us[t],own_b=dispatch_period(s,d,t)
        bb[:,t]=z===nothing ? own_b : z[:,t]
    end
    return us,bb
end

function stage_step!(s,d,z,u,rho,k)
    tick=time_ns()
    warm_tick=time_ns()
    if k>1
        # Cache every result before setting any start or modifying the model.
        vals=value.(s.vars)
        multipliers=s.reuse.warm_start=="primal_dual" ? dual.(s.reuse.constraints) : Float64[]
        all(isfinite,vals) && all(isfinite,multipliers) || error("Nonfinite warm start")
        set_start_value.(s.vars,vals)
        if s.reuse.warm_start=="primal_dual"
            set_dual_start_value.(s.reuse.constraints,multipliers)
            set_optimizer_attribute(s.m,"warm_start_init_point","yes")
            set_optimizer_attribute(s.m,"warm_start_bound_push",1e-8)
            set_optimizer_attribute(s.m,"warm_start_bound_frac",1e-8)
            set_optimizer_attribute(s.m,"warm_start_slack_bound_push",1e-8)
            set_optimizer_attribute(s.m,"warm_start_slack_bound_frac",1e-8)
            set_optimizer_attribute(s.m,"warm_start_mult_bound_push",1e-8)
            set_optimizer_attribute(s.m,"mu_init",1e-6)
        end
    end
    warm_start_s=(time_ns()-warm_tick)/1e9
    update_tick=time_ns()
    if s.reuse.objective_update=="batched"
        nbat=length(d[:Bset])
        for (a,t) in enumerate(s.times),n in 1:nbat
            s.reuse.coefficients[(a-1)*nbat+n]=rho*(u[n,a]-z[n,t])
        end
        set_objective_coefficient(s.m,s.reuse.energy_vars,s.reuse.coefficients)
        if rho!=s.reuse.rho
            set_objective_coefficient(s.m,s.reuse.energy_vars,s.reuse.energy_vars,fill(rho/2,length(s.reuse.energy_vars)))
            s.reuse.quadratic_updates+=1
        end
    else
        for (n,j) in enumerate(d[:Bset]),(a,t) in enumerate(s.times)
            set_objective_coefficient(s.m,s.b[j,t],rho*(u[n,a]-z[n,t]))
            set_objective_coefficient(s.m,s.b[j,t],s.b[j,t],rho/2)
        end
        s.reuse.quadratic_updates+=1
    end
    s.reuse.rho=rho
    objective_update_s=(time_ns()-update_tick)/1e9
    cold_retry=false;warm_attempt_iterations=0
    inexact_accepted=false;cap_recovery=false;cap_attempt_iterations=0;local_primal_gate_s=0.0
    gated_primal=NaN
    cap_end=parse(Int,get(ENV,"TADMM_CAP_END_ITERATION","30"))
    cap_limit=parse(Int,get(ENV,"TADMM_CAP_ITERATIONS","30"))
    cap_tol=parse(Float64,get(ENV,"TADMM_CAP_PRIMAL_TOL","1e-6"))
    cap_screen=get(ENV,"TADMM_CAP_SCREEN","feasible")
    if s.local_policy=="early_cap30" && 1<k<=cap_end
        cap_limit>0 && cap_end>=2 && isfinite(cap_tol) && cap_tol>0 || error("Invalid early cap settings")
        cap_screen in ("feasible","finite") || error("Invalid early cap screen")
        set_optimizer_attribute(s.m,"max_iter",cap_limit)
        solver_s=solve_checked!(s;allow_iteration_limit=true)
        iterations=MOI.get(s.m,MOI.BarrierIterations())
        if termination_status(s.m)==MOI.ITERATION_LIMIT
            cap_attempt_iterations=iterations
            gate_tick=time_ns()
            gated_primal=maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0)
            multipliers=s.reuse.warm_start=="primal_dual" ? dual.(s.reuse.constraints) : Float64[]
            reusable=all(isfinite,multipliers)
            inexact_accepted=isfinite(gated_primal) && reusable && (cap_screen=="finite" || gated_primal<=cap_tol)
            if !inexact_accepted
                # Continue from the interrupted primal (and finite multipliers).
                # Rejected attempt, feasibility gate and exact recovery all count.
                cap_recovery=true
                nextvals=value.(s.vars)
                set_start_value.(s.vars,nextvals)
                if s.reuse.warm_start=="primal_dual" && reusable
                    set_dual_start_value.(s.reuse.constraints,multipliers)
                else
                    set_optimizer_attribute(s.m,"warm_start_init_point","no")
                end
            end
            local_primal_gate_s=(time_ns()-gate_tick)/1e9
            if cap_recovery
                set_optimizer_attribute(s.m,"max_iter",3000)
                solver_s+=solve_checked!(s)
                iterations+=MOI.get(s.m,MOI.BarrierIterations())
            end
        end
    elseif s.local_policy=="warm_retry" && k>1 && s.reuse.warm_start=="primal_dual"
        cutoff=parse(Int,get(ENV,"TADMM_WARM_RETRY_LIMIT","20"));@assert cutoff>0
        interrupted=Ref(false);enabled=Ref(true)
        function callback(alg::Cint,it::Cint,obj::Float64,pr::Float64,du::Float64,mu::Float64,
                dn::Float64,reg::Float64,ad::Float64,ap::Float64,ls::Cint)
            stop=enabled[] && it>=cutoff
            interrupted[] |= stop
            return !stop
        end
        MOI.set(s.m,Ipopt.CallbackFunction(),callback)
        optimize!(s.m)
        solver_s=solve_time(s.m);iterations=MOI.get(s.m,MOI.BarrierIterations())
        if interrupted[] && !(termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL))
            # Restart from the last fully solved primal, not the interrupted
            # iterate. The recovery solve must meet all original tolerances.
            cold_retry=true;warm_attempt_iterations=iterations;enabled[]=false
            set_start_value.(s.vars,vals)
            set_optimizer_attribute(s.m,"warm_start_init_point","no")
            solver_s+=solve_checked!(s)
            iterations+=MOI.get(s.m,MOI.BarrierIterations())
        else
            termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL) || error("Warm solve failed")
        end
    else
        s.local_policy=="early_cap30" && set_optimizer_attribute(s.m,"max_iter",3000)
        solver_s=solve_checked!(s;allow_iteration_limit=s.local_policy=="cap30")
        iterations=MOI.get(s.m,MOI.BarrierIterations())
    end
    iteration_limit_hit=termination_status(s.m)==MOI.ITERATION_LIMIT
    local_primal_evaluated=inexact_accepted || (iteration_limit_hit && s.local_primal_report)
    local_primal=inexact_accepted ? gated_primal : (local_primal_evaluated ? maximum(values(primal_feasibility_report(s.m;atol=0.0));init=0.0) : NaN)
    copy_b=[value(s.b[j,t]) for j in d[:Bset],t in s.times]
    return (;period=only(s.owned),copy_b,solver_s,iterations,status=string(termination_status(s.m)),
        iteration_limit_hit,local_primal,local_primal_evaluated,warm_start_s,objective_update_s,
        quadratic_updates=s.reuse.quadratic_updates,cold_retry,warm_attempt_iterations,
        inexact_accepted,cap_recovery,cap_attempt_iterations,cap_attempt_primal=gated_primal,
        local_primal_gate_s,elapsed_s=(time_ns()-tick)/1e9)
end

# Evaluate the assembled dispatch against original physical equations and bounds,
# independently of JuMP's per-subproblem convergence/status and ADMM residuals.
include("typed_diagnostics.jl")
function physical_metrics(d,us,bb)
    mode=get(ENV,"TADMM_DIAGNOSTICS","cached")
    mode in ("cached","typed") || error("Unknown diagnostic mode")
    mode=="typed" && return typed_physical_metrics(us,bb,typed_physical_layout(d))
    return physical_metrics_cached(d,us,bb,physical_layout(d))
end

function physical_metrics_cached(d,us,bb,layout)
    (;ix,bp,ep,bat,der,root,loads)=layout
    maxeq=0.0; maxbound=0.0; maxcone=0.0; gap=0.0; loss_gap=0.0; obj=0.0
    recovered_energy=[d[:B0_pu][j] for j in d[:Bset]]
    for t in 1:d[:T]
        u=us[t]; p=u[ix.P]; q=u[ix.Q]; ell=u[ix.ell]; v=u[ix.v]; pb=u[ix.pb]
        maxeq=max(maxeq,abs(u[ix.ps]-sum(p[ep[e]] for e in d[:L1set])),
            abs(u[ix.qs]-sum(q[ep[e]] for e in d[:L1set])),abs(v[bp[root]]-1.05^2))
        maxbound=max(maxbound,-u[ix.ps])
        for j in d[:Nm1set]
            e=(d[:parent][j],j); n=ep[e]
            pl=j in loads ? d[:p_L_pu][j,t] : 0.0
            ql=j in loads ? d[:q_L_pu][j,t] : 0.0
            pd=haskey(der,j) ? d[:p_D_pu][j,t] : 0.0
            pbat=haskey(bat,j) ? pb[bat[j]] : 0.0
            qd=haskey(der,j) ? u[ix.qnorm[der[j]]]*sqrt(max(0.0,d[:S_D_R][j]^2-pd^2)) : 0.0
            maxeq=max(maxeq,abs(sum(p[ep[(j,k)]] for k in d[:children][j];init=0.0)-p[n]+d[:rdict_pu][e]*ell[n]-pbat-pd+pl),
                abs(sum(q[ep[(j,k)]] for k in d[:children][j];init=0.0)-q[n]+d[:xdict_pu][e]*ell[n]-qd+ql))
        end
        for (n,e) in enumerate(d[:Lset])
            i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
            maxeq=max(maxeq,abs(v[bp[j]]-v[bp[i]]+2*(r*p[n]+x*q[n])-(r^2+x^2)*ell[n]))
            cg=v[bp[i]]*ell[n]-p[n]^2-q[n]^2
            maxcone=max(maxcone,-cg); gap=max(gap,abs(cg))
            loss_gap+=r*cg/v[bp[i]]*d[:kVA_B]*d[:delta_t_h]
            maxbound=max(maxbound,-ell[n])
        end
        for (n,j) in enumerate(d[:Nset])
            maxbound=max(maxbound,d[:Vminpu][j]^2-v[n],v[n]-d[:Vmaxpu][j]^2)
        end
        for (n,j) in enumerate(d[:Bset])
            # A dispatch handed to OpenDSS contains battery POWERS. Check the
            # energy obtained by integrating those powers as well as consensus
            # energies; small per-period errors must not hide accumulated drift.
            recovered_energy[n]-=d[:delta_t_h]*pb[n]
            prev=t==1 ? d[:B0_pu][j] : bb[n,t-1]
            maxeq=max(maxeq,abs(bb[n,t]-prev+d[:delta_t_h]*pb[n]))
            maxbound=max(maxbound,abs(pb[n])-d[:P_B_R_pu][j],
                d[:soc_min][j]*d[:B_R_pu][j]-bb[n,t],bb[n,t]-d[:soc_max][j]*d[:B_R_pu][j],
                d[:soc_min][j]*d[:B_R_pu][j]-recovered_energy[n],recovered_energy[n]-d[:soc_max][j]*d[:B_R_pu][j])
        end
        maxbound=max(maxbound,maximum(abs.(u[ix.qnorm]);init=0.0)-1)
        obj+=d[:LoadShapeCost][t]*d[:kVA_B]*d[:delta_t_h]*u[ix.ps]+
            d[:C_B]*d[:kVA_B]^2*d[:delta_t_h]*sum(abs2,pb)
    end
    return (;objective=obj,original_primal=max(maxeq,maxbound,maxcone),maxeq,maxbound,maxcone,soc_equality_gap=gap,fictitious_loss_kwh=loss_gap)
end
