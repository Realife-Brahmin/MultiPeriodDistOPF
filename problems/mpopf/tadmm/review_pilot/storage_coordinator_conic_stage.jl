# tADMM's own solves with a conic solver instead of Ipopt
# (SC_STAGE_SOLVER=mosek|gurobi|clarabel). Included by storage_coordinator_actor.jl.
#
# The period OPF states the cone explicitly (rotated second-order cone) and keeps
# every original bound; the battery step is the same small quadratic program.
# A conic interior-point solve cannot start from the previous iterate, so each
# iteration solves every period from scratch.
const SC_STAGE_SOLVER=get(ENV,"SC_STAGE_SOLVER","ipopt")
SC_STAGE_SOLVER=="mosek" && @eval using MosekTools
SC_STAGE_SOLVER=="gurobi" && @eval using Gurobi
SC_STAGE_SOLVER=="clarabel" && @eval using Clarabel

# Same tolerances as the centralized conic baseline (conic_baseline.jl).
function stage_optimizer(;threads=1)
    tol=parse(Float64,get(ENV,"CONIC_TOL","1e-8"))
    if SC_STAGE_SOLVER=="mosek"
        return optimizer_with_attributes(Main.MosekTools.Optimizer,"QUIET"=>true,"MSK_IPAR_NUM_THREADS"=>threads,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>tol,"MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>tol,"MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>tol,
            "MSK_DPAR_INTPNT_QO_TOL_REL_GAP"=>tol,"MSK_DPAR_INTPNT_QO_TOL_PFEAS"=>tol,"MSK_DPAR_INTPNT_QO_TOL_DFEAS"=>tol)
    elseif SC_STAGE_SOLVER=="gurobi"
        return optimizer_with_attributes(Main.Gurobi.Optimizer,"OutputFlag"=>0,"Threads"=>threads,"Method"=>2,"Crossover"=>0,
            "BarConvTol"=>tol,"BarQCPConvTol"=>tol)
    elseif SC_STAGE_SOLVER=="clarabel"
        return optimizer_with_attributes(Main.Clarabel.Optimizer,"verbose"=>false,"max_threads"=>threads,
            "tol_gap_abs"=>tol,"tol_gap_rel"=>tol,"tol_feas"=>tol,"max_iter"=>500)
    end
    error("no conic stage solver selected")
end
conic_ok(m)=termination_status(m) in (MOI.OPTIMAL,MOI.ALMOST_OPTIMAL,MOI.LOCALLY_SOLVED)
conic_strict(m)=termination_status(m) in (MOI.OPTIMAL,MOI.LOCALLY_SOLVED)

function make_power_stage_conic(d,t;rho)
    ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    loads=Set(d[:NLset]); batteries=Set(bs); ders=Set(ds)
    root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]
    m=Model(stage_optimizer())
    @variable(m,ps>=0); @variable(m,qs)
    @variable(m,p[es]); @variable(m,q[es]); @variable(m,v[ns]); @variable(m,ell[es]>=0)
    @variable(m,qd[ds]); @variable(m,pb[bs])
    @constraint(m,ps==sum(p[e] for e in d[:L1set]))
    @constraint(m,qs==sum(q[e] for e in d[:L1set]))
    for j in d[:Nm1set]
        e=(d[:parent][j],j)
        pl=j in loads ? d[:p_L_pu][j,t] : 0.0; ql=j in loads ? d[:q_L_pu][j,t] : 0.0
        pd=j in ders ? d[:p_D_pu][j,t] : 0.0
        @constraint(m,sum(p[(j,k)] for k in d[:children][j];init=0.0)-p[e]+d[:rdict_pu][e]*ell[e]==(j in batteries ? pb[j] : 0.0)+pd-pl)
        @constraint(m,sum(q[(j,k)] for k in d[:children][j];init=0.0)-q[e]+d[:xdict_pu][e]*ell[e]==(j in ders ? qd[j] : 0.0)-ql)
    end
    for e in es
        i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
        @constraint(m,v[j]==v[i]-2*(r*p[e]+x*q[e])+(r^2+x^2)*ell[e])
        @constraint(m,[0.5*v[i],ell[e],p[e],q[e]] in RotatedSecondOrderCone())
    end
    @constraint(m,v[root]==1.05^2)
    for j in ns
        j==root && continue
        set_lower_bound(v[j],d[:Vminpu][j]^2); set_upper_bound(v[j],d[:Vmaxpu][j]^2)
    end
    for j in ds
        qmax=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
        set_lower_bound(qd[j],-qmax); set_upper_bound(qd[j],qmax)
    end
    for j in bs
        set_lower_bound(pb[j],-d[:P_B_R_pu][j]); set_upper_bound(pb[j],d[:P_B_R_pu][j])
    end
    cost=@expression(m,d[:LoadShapeCost][t]*base*dt*ps)
    return (;m,t,ps,qs,p,q,v,ell,qd,pb,cb=d[:C_B]*base^2*dt,cost,vars=all_variables(m),layout=physical_layout(d),conic=true)
end

function power_stage_step_conic!(s,d,zt,ut,rho)
    bs=d[:Bset]
    @objective(s.m,Min,s.cost+sum((s.cb+rho/2)*s.pb[j]^2+rho*(ut[n]-zt[n])*s.pb[j] for (n,j) in enumerate(bs)))
    optimize!(s.m)
    conic_ok(s.m) || error("stage $(s.t) failed: $(termination_status(s.m))")
    return [value(s.pb[j]) for j in bs],solve_time(s.m),0
end

function repair_stage_conic!(s,d,zt)
    bs=d[:Bset]
    for (n,j) in enumerate(bs); fix(s.pb[j],zt[n];force=true); end
    @objective(s.m,Min,s.cost)
    optimize!(s.m)
    strict=conic_strict(s.m); ok=has_values(s.m) && all(isfinite,value.(s.vars))
    u=ok ? stage_dispatch(s,d) : Float64[]; status=string(termination_status(s.m)); solver_s=solve_time(s.m)
    for j in bs
        unfix(s.pb[j]); set_lower_bound(s.pb[j],-d[:P_B_R_pu][j]); set_upper_bound(s.pb[j],d[:P_B_R_pu][j])
    end
    return (;ok,strict,status,u,solver_s)
end

function make_coordinator_conic(d;idx=eachindex(d[:Bset]))
    T=d[:T]; bs=d[:Bset][idx]; dt=d[:delta_t_h]; nb=length(bs)
    m=Model(stage_optimizer())
    @variable(m,z[1:nb,1:T]); @variable(m,b[1:nb,1:T])
    for (n,j) in enumerate(bs),t in 1:T
        set_lower_bound(z[n,t],-d[:P_B_R_pu][j]); set_upper_bound(z[n,t],d[:P_B_R_pu][j])
        set_lower_bound(b[n,t],d[:soc_min][j]*d[:B_R_pu][j]); set_upper_bound(b[n,t],d[:soc_max][j]*d[:B_R_pu][j])
        @constraint(m,b[n,t]==(t==1 ? d[:B0_pu][j] : b[n,t-1])-dt*z[n,t])
    end
    return (;m,z,b,vars=all_variables(m),conic=true)
end

function coordinator_step_conic!(c,d,idx,target;pinned=Dict{Int,Vector{Float64}}())
    limits=[d[:P_B_R_pu][j] for j in d[:Bset][idx]]
    for (t,p) in pinned, n in eachindex(p); fix(c.z[n,t],p[n];force=true); end
    @objective(c.m,Min,sum(c.z[i]^2-2*target[i]*c.z[i] for i in eachindex(target)))
    optimize!(c.m)
    ok=conic_ok(c.m); z=ok ? value.(c.z) : zeros(0,0); solver_s=solve_time(c.m)
    for (t,p) in pinned, n in eachindex(p)
        unfix(c.z[n,t]); set_lower_bound(c.z[n,t],-limits[n]); set_upper_bound(c.z[n,t],limits[n])
    end
    return (;ok,z,solver_s)
end
