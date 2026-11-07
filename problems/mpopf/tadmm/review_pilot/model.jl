# A common physical model for centralized and overlapping-window tADMM solves.
# Each stage owns its network and cost. Copies of neighboring battery energies
# and powers enforce local dynamics; consensus identifies every energy copy.
using JuMP, Ipopt, LinearAlgebra, Printf, Serialization
const MOI = JuMP.MOI
const SHARED = normpath(joinpath(@__DIR__, "..", "..", "spar_ddp", "examples", "power_system"))
include(joinpath(SHARED, "terminal_soc_penalty.jl"))
include(joinpath(SHARED, "loadflow_start.jl"))
include(joinpath(SHARED, "control_layout.jl"))

function make_model(d, owned; local_window=false, rho=0.0)
    T = d[:T]; ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]
    times = local_window ? collect(max(1,only(owned)-1):min(T,only(owned)+1)) : collect(1:T)
    # Only dynamics owned by t and its successor are needed in each window.
    dyn = local_window ? collect(only(owned):min(T,only(owned)+1)) : collect(1:T)
    m=Model(Ipopt.Optimizer)
    set_silent(m)
    set_optimizer_attribute(m,"tol",1e-9)
    set_optimizer_attribute(m,"constr_viol_tol",1e-9)
    set_optimizer_attribute(m,"acceptable_iter",0)
    set_optimizer_attribute(m,"max_iter",3000)
    set_optimizer_attribute(m,"bound_relax_factor",0.0)
    set_optimizer_attribute(m,"linear_solver",get(ENV,"TADMM_LINEAR_SOLVER","mumps"))
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
            pl=j in d[:NLset] ? d[:p_L_pu][j,t] : 0.0
            ql=j in d[:NLset] ? d[:q_L_pu][j,t] : 0.0
            pd=j in ds ? d[:p_D_pu][j,t] : 0.0
            pbat=j in bs ? pb[j,t] : 0.0
            qder=j in ds ? qd[j,t] : 0.0
            @constraint(m,sum(p[(j,k),t] for k in d[:children][j];init=0.0)-p[e,t]+d[:rdict_pu][e]*ell[e,t]==pbat+pd-pl)
            @constraint(m,sum(q[(j,k),t] for k in d[:children][j];init=0.0)-q[e,t]+d[:xdict_pu][e]*ell[e,t]==qder-ql)
        end
        for e in es
            i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
            @constraint(m,v[j,t]==v[i,t]-2*(r*p[e,t]+x*q[e,t])+(r^2+x^2)*ell[e,t])
            @constraint(m,p[e,t]^2+q[e,t]^2<=v[i,t]*ell[e,t])
        end
        @constraint(m,v[root,t]==1.05^2)
        for j in ns
            set_lower_bound(v[j,t],d[:Vminpu][j]^2)
            set_upper_bound(v[j,t],d[:Vmaxpu][j]^2)
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
    return (;m,ps,qs,p,q,v,ell,qd,pb,b,cost,times,dyn,owned,vars,initial)
end

function solve_checked!(s; reuse=false)
    if reuse
        vals=value.(s.vars)
        set_start_value.(s.vars,vals)
    end
    optimize!(s.m)
    status=termination_status(s.m)
    status in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL) || error("subproblem failed: $status")
    return solve_time(s.m)
end

function dispatch(stages,d; z=nothing)
    T=d[:T]; ix,nu=control_layout(d); us=[zeros(nu) for _ in 1:T]
    bb=zeros(length(d[:Bset]),T)
    for t in 1:T
        s=length(stages)==1 ? only(stages) : stages[t]; u=us[t]
        u[ix.ps]=value(s.ps[t]); u[ix.qs]=value(s.qs[t])
        for (n,e) in enumerate(d[:Lset])
            u[ix.P[n]]=value(s.p[e,t]); u[ix.Q[n]]=value(s.q[e,t]); u[ix.ell[n]]=value(s.ell[e,t])
        end
        for (n,j) in enumerate(d[:Nset]); u[ix.v[n]]=value(s.v[j,t]); end
        for (n,j) in enumerate(d[:Bset])
            u[ix.pb[n]]=value(s.pb[j,t]); bb[n,t]=z===nothing ? value(s.b[j,t]) : z[n,t]
        end
        for (n,j) in enumerate(d[:Dset])
            qm=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
            u[ix.qnorm[n]]=qm>0 ? value(s.qd[j,t])/qm : 0.0
        end
    end
    return us,bb
end

# Evaluate the assembled dispatch against original physical equations and bounds,
# independently of JuMP's per-subproblem convergence/status and ADMM residuals.
function physical_metrics(d,us,bb)
    ix,_=control_layout(d); bp=Dict(j=>n for (n,j) in enumerate(d[:Nset]))
    ep=Dict(e=>n for (n,e) in enumerate(d[:Lset])); bat=Dict(j=>n for (n,j) in enumerate(d[:Bset]))
    der=Dict(j=>n for (n,j) in enumerate(d[:Dset])); root=d[:substationBus]
    maxeq=0.0; maxbound=0.0; maxcone=0.0; gap=0.0; loss_gap=0.0; obj=0.0
    recovered_energy=[d[:B0_pu][j] for j in d[:Bset]]
    for t in 1:d[:T]
        u=us[t]; p=u[ix.P]; q=u[ix.Q]; ell=u[ix.ell]; v=u[ix.v]; pb=u[ix.pb]
        maxeq=max(maxeq,abs(u[ix.ps]-sum(p[ep[e]] for e in d[:L1set])),
            abs(u[ix.qs]-sum(q[ep[e]] for e in d[:L1set])),abs(v[bp[root]]-1.05^2))
        maxbound=max(maxbound,-u[ix.ps])
        for j in d[:Nm1set]
            e=(d[:parent][j],j); n=ep[e]
            pl=j in d[:NLset] ? d[:p_L_pu][j,t] : 0.0
            ql=j in d[:NLset] ? d[:q_L_pu][j,t] : 0.0
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
