# Lean run of the storage-coordinated tADMM, as it would run in service:
# no per-iteration diagnostics, each worker holds only its own periods' data and
# one model per period, and the repair is solved only when the residual test is
# met. Centralized (Ipopt/MA57) is timed first on the same node.
#
#   julia storage_coordinator_lean.jl system T output-directory
#
# Settings as in storage_coordinator_parallel.jl (SC_NWORKERS, SC_RHO, SC_BALANCE,
# SC_EPS_PRI, SC_EPS_DUAL). SC_PHASE lines carry wall-clock stamps so that an
# external memory sampler can attribute its samples to a phase.
#
# tADMM clock: every iteration with its messages and residual test, plus every
# repair attempt, successful or not.
using Distributed, Dates
include("storage_coordinator_parallel.jl")

stamp(name)=(@printf("SC_PHASE %s %.3f\n",name,time()); flush(stdout))

function run_lean(d,system;rho,rho_end,ramp,maxiter,reference,no_tol,phase,out,eps_pri,eps_dual,balance,balance_step,repair_max_iter,repair_backoff)
    T=d[:T]; bs=d[:Bset]; nb=length(bs); dt=d[:delta_t_h]; W=length(SC_WIDS)
    groups=[collect(w:W:T) for w in 1:W]; chunks=[collect(w:W:nb) for w in 1:W]
    build0=time_ns()
    layout=fetch.([remotecall(sc_prepare_lean!,SC_WIDS[w],slim_instance(d,groups[w]),groups[w],chunks[w],rho) for w in 1:W])
    build_s=(time_ns()-build0)/1e9
    phase=="warmup" && length(unique((i.host,i.affinity) for i in layout))!=W && Sys.islinux() && get(ENV,"SC_PIN","1")=="1" && error("Workers share CPUs")
    z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T)
    tight_periods=Int[]; pinned_set=Set{Int}(); zrep=z; pinned_passes=0; pinned_accepted=0
    algorithm=0.0; k=0; rho_start=rho; rho_next=rho; attempts=0; next_attempt=0; rho_floor=0.0; refuse_rho=parse(Float64,get(ENV,"SC_REFUSE_RHO","1")); refuse_cap=parse(Float64,get(ENV,"SC_REFUSE_RHO_MAX","16")); rho_refused=0.0; converged=0;
    # SC_PERIOD_RHO=f: each period t carries its own penalty rho*prm[t]. After the ramp, a period whose
    # mismatch exceeds f times the RMS mismatch (and eps_pri) doubles its multiplier, up to
    # SC_PERIOD_RHO_MAX. SC_REFUSE_PERIOD=1: a refused repair raises only the refused periods (by
    # SC_REFUSE_RHO) instead of the common penalty. Multipliers never decrease, so they settle.
    prm=ones(T); period_f=parse(Float64,get(ENV,"SC_PERIOD_RHO","0")); period_cap=parse(Float64,get(ENV,"SC_PERIOD_RHO_MAX","16"))
    refuse_period=get(ENV,"SC_REFUSE_PERIOD","0")=="1"
    raise!(ts,f)=for t in ts; nv=min(prm[t]*f,period_cap); u[:,t].*=prm[t]/nv; prm[t]=nv; end
    gap=NaN; violation=NaN; objective=NaN
    io=open(joinpath(out,"lean_trace_$(phase).csv"),"w")
    println(io,"iteration,r_rms,s_relative,rho,algorithm_s,iteration_s,slowest_stage_s,slowest_battery_step_s,repair_attempted,repair_strict,r_max,r_max_period,periods_raised,max_multiplier")
    for iteration in 1:maxiter
        k=iteration; tick=time_ns()
        if k<=ramp+1 || balance<=1
            # SC_RHO_STEP=s changes rho only every s iterations (same start and end), e.g. s=4 with
            # 2:0.0625:20 halves rho at iterations 5, 9, 13, 17 and 21.
            step=parse(Int,get(ENV,"SC_RHO_STEP","1"))
            next=rho_start*(rho_end/rho_start)^((step*div(min(k-1,ramp),step))/max(ramp,1)); u.*=rho/next; rho=next
        elseif rho_next!=rho
            u.*=rho/rho_next; rho=rho_next
        end
        jobs=[remotecall(sc_xstep_lean!,SC_WIDS[w],Dict(t=>z[:,t] for t in groups[w]),Dict(t=>u[:,t] for t in groups[w]),Dict(t=>rho*prm[t] for t in groups[w]),k) for w in 1:W]
        slowest=0.0
        for answer in fetch.(jobs), row in answer; x[:,row.t]=row.x; slowest=max(slowest,row.elapsed_s); end
        zold=copy(z); target=x.+u; battery_s=0.0
        zjobs=[(w,remotecall(sc_zstep!,SC_WIDS[w],target[chunks[w],:],k;weights=all(==(1.0),prm) ? nothing : prm)) for w in 1:W if !isempty(chunks[w])]
        for (w,job) in zjobs
            answer=fetch(job); z[chunks[w],:]=answer.z; battery_s=max(battery_s,answer.elapsed_s)
        end
        u.+=x.-z
        rmax,imax=findmax(abs.(x.-z)); r_max_period=imax[2]   # largest single battery-period mismatch
        r_rms=sqrt(sum(abs2,x.-z)/length(z)); s_rms=rho*sqrt(sum(abs2,(z.-zold).*prm')/length(z))
        s_rel=s_rms/max(rho*sqrt(sum(abs2,u.*prm')/length(u)),1e-12)
        rho_next=rho
        if balance>1 && k>ramp
            a=r_rms/eps_pri; b=s_rel/eps_dual
            a>balance*b && (rho_next=rho*balance_step)
            b>balance*a && (rho_next=rho/balance_step)
        end
        rho_next=max(rho_next,rho_floor)
        if period_f>0 && k>ramp
            rt=vec(sqrt.(sum(abs2,x.-z;dims=1)./nb))
            raise!([t for t in 1:T if rt[t]>period_f*r_rms && rt[t]>eps_pri && prm[t]<period_cap],2.0)
        end
        attempted=false; strict=false; rus=Vector{Vector{Float64}}(undef,T)
        # SC_EPS_MAX=e: also require every battery-period mismatch to be below e before a repair.
        eps_max=parse(Float64,get(ENV,"SC_EPS_MAX","Inf"))
        if k>ramp && k>=next_attempt && r_rms<=eps_pri && s_rel<=eps_dual && rmax<=eps_max
            attempted=true; attempts+=1; strict=true
            snap=(;k,z=copy(z),x=copy(x),u=copy(u),rho)   # the state the repair starts from (SC_SAVE_REPAIR)
            empty!(tight_periods); zrep=z
            for answer in fetch.([remotecall(sc_repair_shared!,SC_WIDS[w],Dict(t=>z[:,t] for t in groups[w]);max_iter=repair_max_iter) for w in 1:W]), row in answer
                strict&=row.ok && row.strict; row.ok && (rus[row.t]=row.u); (row.ok && row.strict) || push!(tight_periods,row.t)
                (row.ok && row.strict) || get(ENV,"SC_REPAIR_LOG","0")=="1" && phase=="timed" &&
                    println("SC_REPAIR_REFUSED k=$k t=$(row.t) ok=$(row.ok) status=$(row.status) solver_s=$(round(row.solver_s,digits=3))")
            end
            get(ENV,"SC_REPAIR_LOG","0")=="1" && phase=="timed" &&
                println("SC_REPAIR_ATTEMPT k=$k refused=$(length(tight_periods)) of $T"); flush(stdout)
            # SC_REPAIR_BAND=b: the refused periods try again with the batteries free to move by
            # at most b p.u. (and a loss penalty); the schedule then takes their repaired powers.
            band=parse(Float64,get(ENV,"SC_REPAIR_BAND","0"))
            if !strict && band>0
                zrep=copy(z); ipb=physical_layout(d).ix.pb; ok2=true; refused=sort(tight_periods)
                for answer in fetch.([remotecall(sc_repair_shared!,SC_WIDS[w],Dict(t=>z[:,t] for t in groups[w] if t in refused);max_iter=repair_max_iter,band) for w in 1:W]), row in answer
                    ok2&=row.ok && row.strict
                    row.ok && (rus[row.t]=row.u; zrep[:,row.t]=row.u[ipb])
                end
                phase=="timed" && @printf("SC_REPAIR_BAND k=%d periods=%d strict=%s max_shift=%.3e\n",k,length(refused),ok2,maximum(abs.(zrep.-z)))
                ok2 ? (strict=true; empty!(tight_periods)) : (zrep=z)
            end
            # Periods whose network cannot absorb the battery schedule (no slack in
            # the voltage band) keep the powers their own solve asked for. The
            # battery step is re-solved with those periods pinned, and the repair
            # is run once more with the resulting schedule.
            if !strict
                # SC_PIN_UNION=1: a period that refused once stays pinned in every later attempt.
                get(ENV,"SC_PIN_UNION","0")=="1" ? union!(pinned_set,tight_periods) : (empty!(pinned_set); union!(pinned_set,tight_periods))
                tight=sort!(collect(pinned_set)); pinned_passes+=1
                if !isempty(tight) && length(tight)<=max(1,T÷4)
                    pjobs=[(w,remotecall(sc_zstep_pinned!,SC_WIDS[w],z[chunks[w],:],Dict(t=>x[chunks[w],t] for t in tight))) for w in 1:W if !isempty(chunks[w])]
                    zp=copy(z); feasible=true
                    for (w,job) in pjobs
                        answer=fetch(job); feasible&=answer.ok; answer.ok && (zp[chunks[w],:]=answer.z)
                    end
                    if feasible
                        strict=true
                        for answer in fetch.([remotecall(sc_repair_shared!,SC_WIDS[w],Dict(t=>zp[:,t] for t in groups[w]);max_iter=repair_max_iter) for w in 1:W]), row in answer
                            strict&=row.ok && row.strict; row.ok && (rus[row.t]=row.u)
                            (row.ok && row.strict) || get(ENV,"SC_PIN_UNION","0")=="1" && push!(pinned_set,row.t)
                        end
                        get(ENV,"SC_REPAIR_LOG","0")=="1" && phase=="timed" && println("SC_PINNED_ATTEMPT k=$k pinned=$(length(tight)) strict=$strict")
                        strict && (zrep=zp; pinned_accepted=length(tight))
                    end
                end
            end
            # After a refused repair, iterate a few times before trying again.
            strict || (next_attempt=k+repair_backoff)
            # SC_REFUSE_RHO=f: a refused repair multiplies the penalty by f, and balancing may not take it below that again.
            if !strict && refuse_rho>1
                refuse_period ? raise!(sort(collect(union(Set(tight_periods),pinned_set))),refuse_rho) :
                    (rho_refused==0 && (rho_refused=rho); rho_floor=min(rho*refuse_rho,refuse_cap*rho_refused); rho_next=max(rho_next,rho_floor))
            end
        end
        iteration_s=(time_ns()-tick)/1e9; algorithm+=iteration_s
        # Outside the clock: save the pre-repair state, so that repair variants can be tried offline
        # from it (agreed schedule z, period powers x, scaled duals u, penalty, and each period
        # model's primal iterate before the repair). SC_SAVE_REPAIR=0 turns it off.
        if attempted && phase=="timed" && get(ENV,"SC_SAVE_REPAIR","1")!="0"
            try
                dir=joinpath(out,"repair_snapshots"); mkpath(dir)
                iterates=merge(fetch.([remotecall(()->Dict(t=>SC[:resume][t] for t in SC[:periods] if haskey(SC[:resume],t)),SC_WIDS[w]) for w in 1:W])...)
                serialize(joinpath(dir,"k$(k).jls"),Dict(:system=>system,:T=>T,:k=>snap.k,:z=>snap.z,:x=>snap.x,:u=>snap.u,:rho=>snap.rho,
                    :eta=>sc_eta(),:refused=>sort(collect(tight_periods)),:strict=>strict,:period_iterates=>iterates,:env=>Dict(k=>v for (k,v) in ENV if startswith(k,"SC_") || startswith(k,"TADMM_"))))
            catch err
                println("SC_WARNING repair snapshot not saved: ",sprint(showerror,err))
            end
        end
        println(io,join((k,r_rms,s_rel,rho,algorithm,iteration_s,slowest,battery_s,attempted,strict,rmax,r_max_period,count(>(1.0),prm),maximum(prm)),',')); flush(io)
        if attempted && strict
            # Outside the clock: score the reported dispatch against the original model.
            zb=zeros(nb,T)
            for n in 1:nb
                e=d[:B0_pu][bs[n]]
                for t in 1:T; e-=dt*battery_drain(zrep[n,t]); zb[n,t]=e; end
            end
            pm=physical_metrics(d,rus,zb); objective=pm.objective; violation=pm.original_primal
            gap=abs(objective-reference)/max(1.0,abs(reference))
            violation>no_tol && phase=="timed" && @printf("SC_CHECK_FAILED k=%d eq=%.3e bound=%.3e cone=%.3e pinned=%d\n",k,pm.maxeq,pm.maxbound,pm.maxcone,length(pinned_set))
            if violation>no_tol   # the repaired dispatch fails the original model: same as a refusal
                next_attempt=k+repair_backoff
                refuse_rho>1 && (rho_refused==0 && (rho_refused=rho); rho_floor=min(rho*refuse_rho,refuse_cap*rho_refused); rho_next=max(rho_next,rho_floor))
            end
            if violation<=no_tol
                converged=k
                phase=="timed" && serialize(joinpath(out,"converged_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>objective,:u=>rus,:B=>zb,:eta=>sc_eta(),:status=>"SC_LEAN_CONVERGED",:iteration=>k))
                break
            end
        end
    end
    close(io)
    @printf("SC_LEAN_DONE phase=%s workers=%d iterations=%d converged=%d repair_attempts=%d pinned_passes=%d pinned_periods=%d algorithm_s=%.3f gap=%.3e violation=%.3e build_s=%.3f\n",
        phase,W,k,converged,attempts,pinned_passes,pinned_accepted,algorithm,gap,violation,build_s); flush(stdout)
    return (;k,converged,attempts,algorithm,gap,violation,build_s)
end

function lean_main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    datafile=abspath(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    d=deserialize(datafile); profile=apply_review_cb!(d,system)
    serialize(joinpath(out,"instance.jls"),d)
    no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    nworkers=min(T,parse(Int,get(ENV,"SC_NWORKERS",string(T))))
    println("SC_HOST host=",gethostname()," time=",now()," affinity=",sc_affinity()," julia=",VERSION," pid=",getpid())
    # SC_REFERENCE=<objective> skips the centralized solve, so that the coordinator
    # process holds nothing left over from it and its own memory can be measured.
    skip_central=haskey(ENV,"SC_REFERENCE")
    reference=NaN; central_wall=NaN; central_no_s=NaN
    if skip_central
        reference=parse(Float64,ENV["SC_REFERENCE"])
        @printf("CENTRAL skipped, reference objective %.12f taken from an earlier run\n",reference)
    else
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
        serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>cu,:B=>cb,:eta=>sc_eta(),:status=>"LOCALLY_SOLVED"))
        central=nothing; cu=nothing; cb=nothing; GC.gc(); GC.gc()
        stamp("central_released")
        central_no_s=trace[cno][4]
    end
    tadmm_screen(d); physical_layout(d)   # caches the workers receive with their data
    prated=sort([d[:P_B_R_pu][j] for j in d[:Bset]]); rho_unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*d[:kVA_B]*d[:delta_t_h]/prated[cld(length(prated),2)]
    parts=split(get(ENV,"SC_RHO","2:0.0625:20"),':'); rho=rho_unit*parse(Float64,parts[1])
    rho_end=length(parts)==3 ? rho_unit*parse(Float64,parts[2]) : rho; ramp=length(parts)==3 ? parse(Int,parts[3]) : 0
    eps_pri=parse(Float64,get(ENV,"SC_EPS_PRI","1e-5")); eps_dual=parse(Float64,get(ENV,"SC_EPS_DUAL","1e-3"))
    balance=parse(Float64,get(ENV,"SC_BALANCE","2")); balance_step=parse(Float64,get(ENV,"SC_BALANCE_STEP","1.5"))
    @printf("CASE system=%s T=%d profile=%s batteries=%d workers=%d rho_unit=%.6g rho=%s balance=%g eps_pri=%g eps_dual=%g\n",
        system,T,profile,length(d[:Bset]),nworkers,rho_unit,get(ENV,"SC_RHO","2:0.0625:20"),balance,eps_pri,eps_dual); flush(stdout)
    stamp("pool_begin"); start_pool!(nworkers); stamp("pool_end")
    println("SC_PIDS ",join(fetch.([remotecall(getpid,w) for w in SC_WIDS]),' '))
    repair_max_iter=parse(Int,get(ENV,"SC_REPAIR_MAX_ITER","60")); repair_backoff=parse(Int,get(ENV,"SC_REPAIR_BACKOFF","3"))
    @printf("SC_SETTINGS repair_max_iter=%d repair_backoff=%d stage_solver=%s eta=%g\n",repair_max_iter,repair_backoff,SC_STAGE_SOLVER,sc_eta())
    common=(;rho,rho_end,ramp,maxiter=parse(Int,get(ENV,"SC_MAX_ITER","150")),reference,no_tol,out,eps_pri,eps_dual,balance,balance_step,repair_max_iter,repair_backoff)
    stamp("tadmm_warmup_begin"); warm=run_lean(d,system;common...,phase="warmup"); stamp("tadmm_warmup_end")
    stamp("tadmm_timed_begin"); timed=run_lean(d,system;common...,phase="timed"); stamp("tadmm_timed_end")
    warm.k==timed.k || println("SC_WARNING warm-up and timed passes differ: ",warm.k," vs ",timed.k)
    @printf("SC_LEAN_SUMMARY system=%s T=%d workers=%d | centralized converged_s=%.3f within_0.1%%_s=%.3f | tadmm converged_s=%.3f iterations=%d repair_attempts=%d gap=%.2e violation=%.2e | ratio=%.2f\n",
        system,T,nworkers,central_wall,central_no_s,timed.algorithm,timed.k,timed.attempts,timed.gap,timed.violation,central_wall/timed.algorithm)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    try
        lean_main()
    finally
        isempty(SC_WIDS) || rmprocs(copy(SC_WIDS))
    end
end
