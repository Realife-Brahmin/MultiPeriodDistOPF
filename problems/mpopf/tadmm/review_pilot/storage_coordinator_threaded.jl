# Storage-coordinated tADMM on the threads of ONE Julia process, to avoid the
# fixed memory of a Julia process per worker. Same algorithm and stopping rule
# as storage_coordinator_lean.jl: one model per period, repair only when the
# residual test is met, pinned repair for periods without slack.
#
#   julia --threads=N storage_coordinator_threaded.jl system T output-directory
#
# Ipopt with MA57 from an external HSL library is not documented as safe to call
# from several threads; results must be checked against the process version.
using Dates
include("storage_coordinator_actor.jl")
BLAS.set_num_threads(1)
stamp(name)=(@printf("SC_PHASE %s %.3f\n",name,time()); flush(stdout))

function repair_stage!(s,d,zt,iterate;max_iter)
    set_start_value.(s.vars,iterate)
    for (n,j) in enumerate(d[:Bset])
        set_objective_coefficient(s.m,s.pb[j],0.0)
        set_objective_coefficient(s.m,s.pb[j],s.pb[j],s.cb)
        fix(s.pb[j],zt[n];force=true)
    end
    set_optimizer_attribute(s.m,"max_iter",max_iter)
    optimize!(s.m)
    strict=termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
    ok=has_values(s.m) && all(isfinite,value.(s.vars))
    u=ok ? stage_dispatch(s,d) : Float64[]
    set_optimizer_attribute(s.m,"max_iter",3000)
    for j in d[:Bset]
        unfix(s.pb[j]); set_lower_bound(s.pb[j],-d[:P_B_R_pu][j]); set_upper_bound(s.pb[j],d[:P_B_R_pu][j])
    end
    return (;ok,strict,u)
end

function pinned_step!(c,d,idx,target,pinned,iterate)
    limits=[d[:P_B_R_pu][j] for j in d[:Bset][idx]]
    set_start_value.(c.vars,iterate)
    for i in eachindex(target); set_objective_coefficient(c.m,c.z[i],-2*target[i]); end
    for (t,p) in pinned, n in eachindex(p); fix(c.z[n,t],p[n];force=true); end
    optimize!(c.m)
    ok=termination_status(c.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
    z=ok ? value.(c.z) : zeros(0,0)
    for (t,p) in pinned, n in eachindex(p)
        unfix(c.z[n,t]); set_lower_bound(c.z[n,t],-limits[n]); set_upper_bound(c.z[n,t],limits[n])
    end
    return (;ok,z)
end

function run_threaded(d,system;rho,rho_end,ramp,maxiter,reference,no_tol,phase,out,eps_pri,eps_dual,balance,balance_step,repair_max_iter,repair_backoff)
    T=d[:T]; bs=d[:Bset]; nb=length(bs); dt=d[:delta_t_h]; N=Threads.nthreads()
    C=min(N,nb); chunks=[collect(c:C:nb) for c in 1:C]
    build0=time_ns()
    stages=Vector{Any}(undef,T); coords=Vector{Any}(undef,C)
    Threads.@threads :static for t in 1:T; stages[t]=make_power_stage(d,t;rho); end
    Threads.@threads :static for c in 1:C; coords[c]=make_coordinator(d;idx=chunks[c]); end
    build_s=(time_ns()-build0)/1e9
    resume=Vector{Union{Nothing,Vector{Float64}}}(nothing,T); coord_resume=Vector{Union{Nothing,Vector{Float64}}}(nothing,C)
    z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T); zrep=z
    algorithm=0.0; k=0; rho_start=rho; rho_next=rho; attempts=0; next_attempt=0; pinned_passes=0; pinned_accepted=0
    converged=0; gap=NaN; violation=NaN
    io=open(joinpath(out,"threaded_trace_$(phase).csv"),"w")
    # A collection requested by one thread stalls every thread until all of them
    # reach a safe point, and a thread inside a long Ipopt call does not. With
    # SC_GC_EVERY=n the collector is off during the threaded loops and runs once
    # every n iterations, between them and inside the clock.
    gc_every=parse(Int,get(ENV,"SC_GC_EVERY","0"))
    gc_every>0 && GC.enable(false)
    println(io,"iteration,r_rms,s_relative,rho,algorithm_s,iteration_s,repair_attempted,repair_strict")
    for iteration in 1:maxiter
        k=iteration; tick=time_ns()
        if k<=ramp+1 || balance<=1
            next=rho_start*(rho_end/rho_start)^(min(k-1,ramp)/max(ramp,1)); u.*=rho/next; rho=next
        elseif rho_next!=rho
            u.*=rho/rho_next; rho=rho_next
        end
        Threads.@threads :static for t in 1:T
            x[:,t],_,_=power_stage_step!(stages[t],d,z[:,t],u[:,t],rho,k;start=resume[t]); resume[t]=nothing
        end
        zold=copy(z); target=x.+u
        Threads.@threads :static for c in 1:C
            z[chunks[c],:],_=coordinator_step!(coords[c],target[chunks[c],:],k;start=coord_resume[c]); coord_resume[c]=nothing
        end
        u.+=x.-z
        r_rms=sqrt(sum(abs2,x.-z)/length(z)); s_rms=rho*sqrt(sum(abs2,z.-zold)/length(z))
        s_rel=s_rms/max(rho*sqrt(sum(abs2,u)/length(u)),1e-12)
        rho_next=rho
        if balance>1 && k>ramp
            a=r_rms/eps_pri; b=s_rel/eps_dual
            a>balance*b && (rho_next=rho*balance_step)
            b>balance*a && (rho_next=rho/balance_step)
        end
        attempted=false; strict=false; rus=Vector{Vector{Float64}}(undef,T)
        if k>ramp && k>=next_attempt && r_rms<=eps_pri && s_rel<=eps_dual
            attempted=true; attempts+=1; zrep=z
            flags=falses(T)
            Threads.@threads :static for t in 1:T
                resume[t]=value.(stages[t].vars)
                row=repair_stage!(stages[t],d,z[:,t],resume[t];max_iter=repair_max_iter)
                row.ok && (rus[t]=row.u); flags[t]=row.ok && row.strict
            end
            strict=all(flags)
            if !strict
                tight=findall(!,flags); pinned_passes+=1
                if length(tight)<=max(1,T÷4)
                    zp=copy(z); feasible=fill(true,C)
                    Threads.@threads :static for c in 1:C
                        coord_resume[c]=value.(coords[c].vars)
                        row=pinned_step!(coords[c],d,chunks[c],z[chunks[c],:],Dict(t=>x[chunks[c],t] for t in tight),coord_resume[c])
                        feasible[c]=row.ok; row.ok && (zp[chunks[c],:]=row.z)
                    end
                    if all(feasible)
                        Threads.@threads :static for t in 1:T
                            row=repair_stage!(stages[t],d,zp[:,t],resume[t];max_iter=repair_max_iter)
                            row.ok && (rus[t]=row.u); flags[t]=row.ok && row.strict
                        end
                        strict=all(flags); strict && (zrep=zp; pinned_accepted=length(tight))
                    end
                end
            end
            strict || (next_attempt=k+repair_backoff)
        end
        if gc_every>0 && k%gc_every==0
            GC.enable(true); GC.gc(false); GC.enable(false)
        end
        iteration_s=(time_ns()-tick)/1e9; algorithm+=iteration_s
        println(io,join((k,r_rms,s_rel,rho,algorithm,iteration_s,attempted,strict),',')); flush(io)
        if attempted && strict
            zb=zeros(nb,T)
            for n in 1:nb
                e=d[:B0_pu][bs[n]]
                for t in 1:T; e-=dt*zrep[n,t]; zb[n,t]=e; end
            end
            pm=physical_metrics(d,rus,zb); violation=pm.original_primal
            gap=abs(pm.objective-reference)/max(1.0,abs(reference))
            if violation<=no_tol
                converged=k
                phase=="timed" && serialize(joinpath(out,"converged_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>pm.objective,:u=>rus,:B=>zb,:status=>"SC_THREADED_CONVERGED",:iteration=>k))
                break
            end
        end
    end
    close(io)
    gc_every>0 && GC.enable(true)
    @printf("SC_THREADED_DONE phase=%s threads=%d iterations=%d converged=%d repair_attempts=%d pinned_passes=%d pinned_periods=%d algorithm_s=%.3f gap=%.3e violation=%.3e build_s=%.3f\n",
        phase,N,k,converged,attempts,pinned_passes,pinned_accepted,algorithm,gap,violation,build_s); flush(stdout)
    return (;k,converged,attempts,algorithm,gap,violation,build_s)
end

function threaded_main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    datafile=abspath(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    d=deserialize(datafile); profile=apply_review_cb!(d,system)
    serialize(joinpath(out,"instance.jls"),d)
    no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    println("SC_HOST host=",gethostname()," time=",now()," threads=",Threads.nthreads()," julia=",VERSION," pid=",getpid())
    stamp("central_build_begin")
    central=make_model(d,collect(1:T))
    set_optimizer_attribute(central.m,"inf_pr_output","original")
    trace=Tuple{Int,Float64,Float64,Float64}[]; clock=Ref(time_ns())
    function callback(alg_mod::Cint,iter::Cint,obj::Float64,inf_pr::Float64,inf_du::Float64,mu::Float64,
                      dnorm::Float64,reg::Float64,alpha_du::Float64,alpha_pr::Float64,ls::Cint)
        push!(trace,(Int(iter),obj,inf_pr,(time_ns()-clock[])/1e9)); return true
    end
    MOI.set(central.m,Ipopt.CallbackFunction(),callback)
    stamp("central_warmup_begin"); solve_checked!(central); stamp("central_warmup_end")
    set_start_value.(central.vars,central.initial); empty!(trace)
    stamp("central_timed_begin"); clock[]=time_ns(); tick=time_ns(); solve_checked!(central); central_wall=(time_ns()-tick)/1e9; stamp("central_timed_end")
    cu,cb=dispatch([central],d); cm=physical_metrics(d,cu,cb); reference=cm.objective
    cno=findfirst(r->abs(r[2]-trace[end][2])<=1e-3*max(1.0,abs(trace[end][2])) && r[3]<=no_tol,trace)
    @printf("CENTRAL objective=%.12f violation=%.3e iterations=%d wall_s=%.3f within_0.1%%_s=%.3f\n",reference,cm.original_primal,trace[end][1],central_wall,trace[cno][4])
    serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>cu,:B=>cb,:status=>"LOCALLY_SOLVED"))
    central=nothing; cu=nothing; cb=nothing; GC.gc(); GC.gc()
    stamp("central_released")
    tadmm_screen(d); physical_layout(d)   # fill the caches before any thread reads d
    prated=sort([d[:P_B_R_pu][j] for j in d[:Bset]]); rho_unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*d[:kVA_B]*d[:delta_t_h]/prated[cld(length(prated),2)]
    parts=split(get(ENV,"SC_RHO","2:0.0625:20"),':'); rho=rho_unit*parse(Float64,parts[1])
    rho_end=length(parts)==3 ? rho_unit*parse(Float64,parts[2]) : rho; ramp=length(parts)==3 ? parse(Int,parts[3]) : 0
    common=(;rho,rho_end,ramp,maxiter=parse(Int,get(ENV,"SC_MAX_ITER","150")),reference,no_tol,out,
        eps_pri=parse(Float64,get(ENV,"SC_EPS_PRI","1e-5")),eps_dual=parse(Float64,get(ENV,"SC_EPS_DUAL","1e-3")),
        balance=parse(Float64,get(ENV,"SC_BALANCE","2")),balance_step=parse(Float64,get(ENV,"SC_BALANCE_STEP","1.5")),
        repair_max_iter=parse(Int,get(ENV,"SC_REPAIR_MAX_ITER","60")),repair_backoff=parse(Int,get(ENV,"SC_REPAIR_BACKOFF","3")))
    @printf("CASE system=%s T=%d profile=%s batteries=%d threads=%d rho_unit=%.6g\n",system,T,profile,length(d[:Bset]),Threads.nthreads(),rho_unit); flush(stdout)
    stamp("tadmm_warmup_begin"); warm=run_threaded(d,system;common...,phase="warmup"); stamp("tadmm_warmup_end")
    GC.gc(); GC.gc()
    stamp("tadmm_timed_begin"); timed=run_threaded(d,system;common...,phase="timed"); stamp("tadmm_timed_end")
    warm.k==timed.k || println("SC_WARNING warm-up and timed passes differ: ",warm.k," vs ",timed.k)
    @printf("SC_THREADED_SUMMARY system=%s T=%d threads=%d | centralized converged_s=%.3f within_0.1%%_s=%.3f | tadmm converged_s=%.3f iterations=%d repair_attempts=%d gap=%.2e violation=%.2e | ratio=%.2f | process_maxrss_gb=%.2f\n",
        system,T,Threads.nthreads(),central_wall,trace[cno][4],timed.algorithm,timed.k,timed.attempts,timed.gap,timed.violation,central_wall/timed.algorithm,Sys.maxrss()/2^30)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && threaded_main()
