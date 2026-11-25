# Centralized MPOPF as an explicit second-order cone program, solved by a
# dedicated conic solver, next to the Ipopt/MA57 solve of the same instance.
#
#   julia --project=<conic_env> conic_baseline.jl system T output-directory solver[,solver...]
#
# solvers: clarabel, gurobi, mosek, ipopt. Same data, objective and limits as
# model.jl's centralized model, without its screening: every original bound is
# kept and the cone is passed as a rotated second-order cone
#   2*(v_i/2)*ell_ij >= P_ij^2 + Q_ij^2.
# Each solver gets a warm-up solve; the second solve is the timed one.
include("model.jl")
using Dates
# Load the requested conic solvers before any function that uses them is compiled.
const CONIC_SOLVERS=length(ARGS)>=4 ? split(ARGS[4],',') : String[]
"clarabel" in CONIC_SOLVERS && @eval using Clarabel
"gurobi" in CONIC_SOLVERS && @eval using Gurobi
"mosek" in CONIC_SOLVERS && @eval using MosekTools

function make_conic_model(d,optimizer)
    T=d[:T]; ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]
    loads=Set(d[:NLset]); batteries=Set(bs); ders=Set(ds)
    root=d[:substationBus]; dt=d[:delta_t_h]; base=d[:kVA_B]; owned=1:T
    m=Model(optimizer)
    # CONIC_SCREEN=1: drop the same provably inactive bounds as the Ipopt model (model.jl):
    # upper voltage limits the data-only screen rules out, and P_Subs >= 0 when it is implied.
    screened=get(ENV,"CONIC_SCREEN","0")=="1"; screen=screened ? tadmm_screen(d) : nothing
    @variable(m,ps[owned]); @variable(m,qs[owned])
    (screened && screen.drop_ps) || for t in owned; set_lower_bound(ps[t],0.0); end
    @variable(m,p[es,owned]); @variable(m,q[es,owned]); @variable(m,v[ns,owned]); @variable(m,ell[es,owned]>=0)
    @variable(m,qd[ds,owned]); @variable(m,pb[bs,owned]); @variable(m,b[bs,owned])
    @objective(m,Min,sum(d[:LoadShapeCost][t]*base*dt*ps[t] for t in owned)+sum(d[:C_B]*base^2*dt*pb[j,t]^2 for j in bs,t in owned))
    for t in owned
        @constraint(m,ps[t]==sum(p[e,t] for e in d[:L1set]))
        @constraint(m,qs[t]==sum(q[e,t] for e in d[:L1set]))
        for j in d[:Nm1set]
            e=(d[:parent][j],j)
            pl=j in loads ? d[:p_L_pu][j,t] : 0.0; ql=j in loads ? d[:q_L_pu][j,t] : 0.0
            pd=j in ders ? d[:p_D_pu][j,t] : 0.0
            pbat=j in batteries ? pb[j,t] : 0.0; qder=j in ders ? qd[j,t] : 0.0
            qcap=haskey(get(d,:qC_pu,Dict{Int,Float64}()),j) ? d[:qC_pu][j]*v[j,t] : 0.0   # fixed capacitor, q = qC * v
            @constraint(m,sum(p[(j,k),t] for k in d[:children][j];init=0.0)-p[e,t]+d[:rdict_pu][e]*ell[e,t]==pbat+pd-pl)
            @constraint(m,sum(q[(j,k),t] for k in d[:children][j];init=0.0)-q[e,t]+d[:xdict_pu][e]*ell[e,t]==qder-ql+qcap)
        end
        for e in es
            i,j=e; r=d[:rdict_pu][e]; x=d[:xdict_pu][e]
            @constraint(m,v[j,t]==v[i,t]-2*(r*p[e,t]+x*q[e,t])+(r^2+x^2)*ell[e,t])
            if get(ENV,"SC_CONIC_CONE","rsoc")=="quad"   # as in the tADMM period model
                @constraint(m,p[e,t]^2+q[e,t]^2<=v[i,t]*ell[e,t])
            else
                @constraint(m,[0.5*v[i,t],ell[e,t],p[e,t],q[e,t]] in RotatedSecondOrderCone())
            end
        end
        @constraint(m,v[root,t]==1.05^2)
        for (n,j) in enumerate(ns)
            j==root && continue
            set_lower_bound(v[j,t],d[:Vminpu][j]^2)
            (screened && "vupper" in screen.flags && !screen.keep_upper[n]) || set_upper_bound(v[j,t],d[:Vmaxpu][j]^2)
        end
        for j in ds
            qmax=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
            set_lower_bound(qd[j,t],-qmax); set_upper_bound(qd[j,t],qmax)
        end
    end
    for j in bs, t in owned
        set_lower_bound(b[j,t],d[:soc_min][j]*d[:B_R_pu][j]); set_upper_bound(b[j,t],d[:soc_max][j]*d[:B_R_pu][j])
        set_lower_bound(pb[j,t],-d[:P_B_R_pu][j]); set_upper_bound(pb[j,t],d[:P_B_R_pu][j])
        prev=t==1 ? d[:B0_pu][j] : b[j,t-1]; eta=sc_eta()
        if eta<1   # charging and discharging powers with efficiency eta each way, as in model.jl
            pc=@variable(m,lower_bound=0,upper_bound=d[:P_B_R_pu][j]); pd=@variable(m,lower_bound=0,upper_bound=d[:P_B_R_pu][j])
            @constraint(m,pb[j,t]==pd-pc); @constraint(m,b[j,t]==prev+dt*(eta*pc-pd/eta))
        else
            @constraint(m,b[j,t]==prev-dt*pb[j,t])
        end
    end
    return (;m,ps,qs,p,q,v,ell,qd,pb,b,layout=physical_layout(d))
end

function optimizer_for(name)
    if name=="clarabel"
        return optimizer_with_attributes(Main.Clarabel.Optimizer,"verbose"=>false,
            "max_threads"=>1,
            "tol_gap_abs"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),"tol_gap_rel"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),
            "tol_feas"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),"max_iter"=>500)
    elseif name=="gurobi"
        return optimizer_with_attributes(Main.Gurobi.Optimizer,"OutputFlag"=>0,"Threads"=>parse(Int,get(ENV,"CONIC_THREADS","1")),
            "Method"=>2,"Crossover"=>0,"BarConvTol"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),"BarQCPConvTol"=>parse(Float64,get(ENV,"CONIC_QCPTOL",get(ENV,"CONIC_TOL","1e-8"))))
    elseif name=="mosek"
        return optimizer_with_attributes(Main.MosekTools.Optimizer,"QUIET"=>true,"MSK_IPAR_NUM_THREADS"=>parse(Int,get(ENV,"CONIC_THREADS","1")),
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),"MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")),
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>parse(Float64,get(ENV,"CONIC_TOL","1e-8")))
    end
    error("unknown solver $name")
end

function conic_main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out); solvers=split(ARGS[4],',')
    d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    apply_review_cb!(d,system)
    println("CONIC_HOST host=",gethostname()," time=",now()," julia=",VERSION," system=",system," T=",T," eta=",sc_eta(),
        " screen=",get(ENV,"CONIC_SCREEN","0")," tol=",get(ENV,"CONIC_TOL","1e-8")," qcptol=",get(ENV,"CONIC_QCPTOL",get(ENV,"CONIC_TOL","1e-8"))," cone=",get(ENV,"SC_CONIC_CONE","rsoc"))
    "gurobi" in solvers && println("CONIC_GUROBI version=",join(let a=Ref{Cint}(),b=Ref{Cint}(),c=Ref{Cint}(); Main.Gurobi.GRBversion(a,b,c); (a[],b[],c[]) end,'.'))
    rows=String[]
    for name in solvers
        if name=="ipopt"
            s=make_model(d,collect(1:T))
            solve_checked!(s); set_start_value.(s.vars,s.initial)
            tick=time_ns(); solve_checked!(s); wall=(time_ns()-tick)/1e9
            us,bb=dispatch([s],d); pm=physical_metrics(d,us,bb)
            status="LOCALLY_SOLVED"; iterations=MOI.get(s.m,MOI.BarrierIterations())
        else
            # Warm-up on one model, timed solve on a second, freshly built one: a
            # solver may return at once when asked to re-solve an unchanged model.
            warm=make_conic_model(d,optimizer_for(name)); optimize!(warm.m); warm=nothing; GC.gc()
            build0=time_ns(); s=make_conic_model(d,optimizer_for(name)); build_s=(time_ns()-build0)/1e9
            tick=time_ns(); optimize!(s.m); wall=(time_ns()-tick)/1e9
            status=string(termination_status(s.m)); iterations=try MOI.get(s.m,MOI.BarrierIterations()) catch; -1 end
            if !has_values(s.m)
                @printf("CONIC solver=%s status=%s wall_s=%.3f no solution\n",name,status,wall); continue
            end
            us,bb=dispatch([s],d); pm=physical_metrics(d,us,bb)
            wall=max(wall,0.0); solver_s=try solve_time(s.m) catch; NaN end
            @printf("CONIC_DETAIL solver=%s build_s=%.2f reported_solve_s=%.3f\n",name,build_s,solver_s)
        end
        @printf("CONIC solver=%s status=%s iterations=%d wall_s=%.3f objective=%.9f violation=%.3e cone_gap=%.3e maxrss_gb=%.2f\n",
            name,status,iterations,wall,pm.objective,pm.original_primal,pm.soc_equality_gap,Sys.maxrss()/2^30); flush(stdout)
        push!(rows,join((system,T,name,status,iterations,wall,pm.objective,pm.original_primal,pm.soc_equality_gap),','))
        serialize(joinpath(out,"$(name)_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>pm.objective,:u=>us,:B=>bb,:status=>status))
        s=nothing; GC.gc()
    end
    open(joinpath(out,"conic_summary.csv"),"w") do io
        println(io,"system,T,solver,status,iterations,wall_s,objective,violation,cone_gap"); foreach(r->println(io,r),rows)
    end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && conic_main()
