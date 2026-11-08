include("model.jl")
include("parallel_workers.jl")
include("historical_stopping.jl")
include("legacy_rho.jl")
using SHA, TOML, Dates, Sockets
BLAS.set_num_threads(1)
Threads.nthreads()==1 || error("Use one Julia thread per process")

function phase_event(out,label)
    path=joinpath(out,"phase_events.csv")
    fresh=!isfile(path)
    open(path,"a") do io
        fresh && println(io,"utc,event")
        println(io,string(now(UTC))*"Z,"*label)
    end
end

function save_admm_checkpoint(path,stages,z,dual,rho,k,phase,elapsed,ideal,solver_total,diagnostic,primal_seen,high_primal_count)
    tick=time_ns()
    primals=isempty(POOL_WIDS) ? Dict(t=>value.(s.vars) for (t,s) in enumerate(stages)) : pool_checkpoint()
    state=Dict(:iteration=>k,:phase=>phase,:z=>z,:dual=>dual,:rho=>rho,:stage_primals=>primals,
        :algorithm_s=>elapsed,:ideal_s=>ideal,:solver_s=>solver_total,:diagnostic_s=>diagnostic,
        :primal_seen=>primal_seen,:high_primal_count=>high_primal_count,
        :model_sha256=>bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl")))),
        :driver_sha256=>bytes2hex(sha256(read(@__FILE__))),
        :note=>"Post-iteration state for recovery/diagnosis. A resumed timed trajectory is not a fresh contiguous benchmark.")
    temporary=path*".tmp"; serialize(temporary,state); mv(temporary,path;force=true)
    return (time_ns()-tick)/1e9
end

function run_admm(d; maxiter=5000, rho0=parse(Float64,get(ENV,"TADMM_RHO_INIT",string(3000.0*sqrt(d[:T]/24)))), out=nothing, phase="timed",reference_objective=NaN,no_primal_tolerance=1e-5,system_name="",stopping_rule="strict",observer=nothing)
    stopping_rule in ("strict","historical_owned_rms","observer") || error("Unknown stopping rule")
    stopping_rule=="observer" && observer===nothing && error("Observer stop requires a callback")
    T=d[:T]; bs=d[:Bset]; nb=length(bs)
    build0=time_ns()
    parallel=parse(Int,get(ENV,"TADMM_NWORKERS","0"))>0
    parallel && prepare_pool!(d,rho0)
    stages=parallel ? nothing : [make_model(d,[t];local_window=true,rho=rho0) for t in 1:T]
    windows=[collect(max(1,t-1):min(T,t+1)) for t in 1:T]
    build_s=(time_ns()-build0)/1e9
    z=[d[:B0_pu][j] for j in bs,t in 1:T]
    dual=[zeros(nb,length(window)) for window in windows]
    copies=[zeros(nb,length(window)) for window in windows]
    multiplicity=[count(window->t in window,windows) for t in 1:T]
    rho=rho0; elapsed=0.0; ideal=0.0; solver_total=0.0; diagnostic=0.0
    dispatch_s=0.0; physical_check_s=0.0
    rho_policy=get(ENV,"TADMM_RHO_POLICY","phase")
    rho_policy in ("phase","residual","fixed","relative","legacy_phase") || error("unknown rho policy")
    balance_ratio=parse(Float64,get(ENV,"TADMM_RHO_BALANCE_RATIO","10"))
    isfinite(balance_ratio) && balance_ratio>1 || error("Rho balance ratio must exceed one")
    primal_seen=false; high_primal_count=0
    legacy_tau_decrease=parse(Float64,get(ENV,"TADMM_RHO_TAU_DECREASE","5"))
    legacy_watchdog=get(ENV,"TADMM_RHO_WATCHDOG","0")=="1"
    legacy_primal_gate=parse(Float64,get(ENV,"TADMM_RHO_PRIMAL_GATE","1e-4"))
    checkpoint_s=0.0; checkpoint_interval=parse(Int,get(ENV,"TADMM_CHECKPOINT_INTERVAL","25"))
    rmax=Inf; smax=Inf; ok=false
    no_iteration=0; no_wall_s=NaN; no_save_s=0.0
    objectives=Float64[]; historical_iteration=0; historical_met=false
    owned_r=Inf; owned_s=Inf; historical_stagnant=false
    tol_r=parse(Float64,get(ENV,"TADMM_PRIMAL_TOL","2.5e-7"))
    tol_s=parse(Float64,get(ENV,"TADMM_DUAL_TOL","1e-4"))
    io=out===nothing ? nothing : open(out,"w")
    details=out===nothing ? nothing : open(replace(out,".csv"=>"_subproblems.csv"),"w")
    details===nothing || println(details,"iteration,period,status,ipopt_iterations,solver_s,stage_elapsed_s,iteration_limit_hit,local_primal,local_primal_evaluated,warm_start_s,objective_update_s,quadratic_updates")
    io===nothing || println(io,"iteration,objective,original_primal,r_rms,r_max,s_rms,s_max,rho,algorithm_s,ideal_parallel_estimate_s,solver_sum_s,diagnostic_s,r_relative,s_relative,loop_wall_s,owned_r_rms,owned_s_rms,historical_stagnant,historical_stop_met,dispatch_s,physical_check_s,cone_equality_gap,fictitious_loss_kwh")
    start=time_ns(); metrics=nothing; us=nothing; bb=nothing; k=0
    for iteration in 1:maxiter
        k=iteration; iter0=time_ns(); times=Float64[]; solve_times=Float64[]; local_iterations=Int[]
        rows=parallel ? pool_step!(z,dual,rho,k,T) : [stage_step!(stages[t],d,z,dual[t],rho,k) for t in 1:T]
        for t0 in 1:T
            row=rows[t0]
            copies[t0]=row.copy_b
            push!(solve_times,row.solver_s); push!(local_iterations,row.iterations); push!(times,row.elapsed_s)
        end
        coord0=time_ns(); zold=copy(z); fill!(z,0.0)
        for t0 in 1:T,(a,t) in enumerate(windows[t0]),n in 1:nb
            z[n,t]+=copies[t0][n,a]+dual[t0][n,a]
        end
        for t in 1:T,n in 1:nb
            j=bs[n]
            z[n,t]=clamp(z[n,t]/multiplicity[t],d[:soc_min][j]*d[:B_R_pu][j],d[:soc_max][j]*d[:B_R_pu][j])
        end
        r2=0.0; s2=0.0; rmax=0.0; smax=0.0; ncopy=0; x2=0.0; z2=0.0; y2=0.0
        for t0 in 1:T,(a,t) in enumerate(windows[t0]),n in 1:nb
            r=copies[t0][n,a]-z[n,t]; s=rho*(z[n,t]-zold[n,t])
            dual[t0][n,a]+=r
            r2+=r*r; s2+=s*s; rmax=max(rmax,abs(r)); smax=max(smax,abs(s)); ncopy+=1
            x2+=copies[t0][n,a]^2; z2+=z[n,t]^2; y2+=(rho*dual[t0][n,a])^2
        end
        rr=sqrt(r2/ncopy); ss=sqrt(s2/ncopy)
        owned_r,owned_s=historical_rms(copies,z,zold,windows,rho)
        # Wohlberg (2017), Eqs. 49, 51: A selects local energy copies,
        # B=-E replicates consensus, c=0, and A' y equals the copy multipliers.
        rrel=sqrt(r2)/max(sqrt(x2),sqrt(z2),1e-12)
        srel=sqrt(s2)/max(sqrt(y2),1e-12)
        coord_s=(time_ns()-coord0)/1e9
        elapsed+=(time_ns()-iter0)/1e9
        critical=parallel ? maximum(sum(times[t] for t in group) for group in POOL_GROUPS) : maximum(times)
        ideal+=critical+coord_s
        solver_total+=sum(solve_times)
        # Physical validation is deliberately outside the algorithm timer;
        # end-to-end wall time below includes it, warm-up and setup are separate.
        diag0=time_ns(); us,bb=parallel ? pool_dispatch(d,z) : dispatch(stages,d;z=z)
        dispatch_s+=(time_ns()-diag0)/1e9
        physical0=time_ns(); metrics=physical_metrics(d,us,bb)
        physical_check_s+=(time_ns()-physical0)/1e9
        diagnostic+=(time_ns()-diag0)/1e9
        push!(objectives,metrics.objective)
        hs=historical_stop(owned_r,owned_s,objectives)
        historical_met=hs.met; historical_stagnant=hs.stagnant
        observer_stop=observer===nothing ? false : observer((;k,metrics,owned_r,owned_s,rr,rmax,smax,rho,elapsed,us,bb))
        if historical_iteration==0 && historical_met
            historical_iteration=k
            if out!==nothing
                serialize(replace(out,".csv"=>"_historical_dispatch.jls"),Dict(:system=>system_name,:T=>T,
                    :objective=>metrics.objective,:u=>us,:B=>bb,:status=>"HISTORICAL_RMS_STOP_CANDIDATE"))
            end
        end
        if no_iteration==0 && isfinite(reference_objective) && metrics.original_primal<=no_primal_tolerance &&
                abs(metrics.objective-reference_objective)/max(1.0,abs(reference_objective))<=0.001
            no_iteration=k; no_wall_s=(time_ns()-start)/1e9
            if out!==nothing
                no_tick=time_ns()
                serialize(replace(out,".csv"=>"_no_dispatch.jls"),Dict(:system=>system_name,:T=>T,
                    :objective=>metrics.objective,:u=>us,:B=>bb,:status=>"ADMM_NO_CANDIDATE"))
                no_save_s+=(time_ns()-no_tick)/1e9
            end
            @printf("TADMM_NO phase=%s k=%d objective=%.12f original=%.3e target_under_100=%s\n",phase,k,metrics.objective,metrics.original_primal,k<100)
            flush(stdout)
        end
        if details!==nothing
            for t0 in 1:T
                println(details,join((k,t0,rows[t0].status,local_iterations[t0],solve_times[t0],times[t0],
                    rows[t0].iteration_limit_hit,rows[t0].local_primal,rows[t0].local_primal_evaluated,
                    rows[t0].warm_start_s,rows[t0].objective_update_s,rows[t0].quadratic_updates),','))
            end
            flush(details)
        end
        if io!==nothing
            println(io,join((k,metrics.objective,metrics.original_primal,rr,rmax,ss,smax,rho,elapsed,ideal,solver_total,diagnostic,rrel,srel,(time_ns()-start)/1e9,owned_r,owned_s,historical_stagnant,historical_met,dispatch_s,physical_check_s,metrics.soc_equality_gap,metrics.fictitious_loss_kwh),',')); flush(io)
        end
        if k==1 || k%25==0
            @printf("TADMM_ITER phase=%s k=%d objective=%.9f physical=%.3e rmax=%.3e smax=%.3e rho=%.5g algorithm_s=%.3f\n",phase,k,metrics.objective,metrics.original_primal,rmax,smax,rho,elapsed); flush(stdout)
        end
        ok=rmax<=tol_r && smax<=tol_s && metrics.original_primal<=2tol_r
        if (stopping_rule=="strict" && ok) || (stopping_rule=="historical_owned_rms" && historical_met) || (stopping_rule=="observer" && observer_stop)
            if out!==nothing && checkpoint_interval>0
                checkpoint_s+=save_admm_checkpoint(replace(out,".csv"=>"_checkpoint.jls"),stages,z,dual,rho,k,phase,elapsed,ideal,solver_total,diagnostic,primal_seen,high_primal_count)
            end
            break
        end
        # July's phase-aware policy, now with ALL copies monitored. Never weaken
        # consensus enforcement merely because the dimensional dual residual is
        # numerically larger while primal consensus remains infeasible.
        if rho_policy=="legacy_phase"
            rho_old=rho
            if legacy_watchdog
                rho,high_primal_count=legacy_rho_watchdog(rho,owned_r,high_primal_count;eps_pri=legacy_primal_gate)
                if rho!=rho_old
                    @printf("LEGACY_WATCHDOG k=%d rho_old=%.8g rho_new=%.8g owned_r=%.3e\n",k,rho_old,rho,owned_r)
                end
            end
            rho,primal_seen=legacy_rho_update(rho,owned_r,owned_s,k,primal_seen;tau_decrease=legacy_tau_decrease,eps_pri=legacy_primal_gate)
            if rho!=rho_old
                foreach(u->u.*=rho_old/rho,dual)
                @printf("LEGACY_RHO k=%d rho_old=%.8g rho_new=%.8g owned_r=%.3e owned_s=%.3e primal_seen=%s\n",k,rho_old,rho,owned_r,owned_s,primal_seen)
                flush(stdout)
            end
        elseif rho_policy=="phase"
            rho_old=rho
            high_primal_count=rmax>2tol_r ? high_primal_count+1 : (rmax<tol_r/2 ? 0 : high_primal_count)
            if high_primal_count>=15
                rho=min(1e6,2.5rho); high_primal_count=0
            end
            if k%5==0
                primal_seen |= rmax<=tol_r
                if !primal_seen && rmax>tol_r
                    rho=min(1e6,2rho)
                elseif rr>5ss
                    rho=min(1e6,2rho)
                elseif ss>5rr && rmax<=tol_r
                    rho=max(1.0,rho/5)
                end
            end
            rho==rho_old || foreach(u->u.*=rho_old/rho,dual)
        elseif rho_policy=="relative" && k%5==0
            rho_old=rho
            if rrel>balance_ratio*srel; rho=min(rho*2,1e6)
            elseif srel>balance_ratio*rrel; rho=max(rho/2,1e-6)
            end
            rho==rho_old || foreach(u->u.*=rho_old/rho,dual)
        elseif rho_policy=="residual" && k%10==0
            rho_old=rho
            if rr>balance_ratio*ss; rho=min(rho*2,1e8)
            elseif ss>balance_ratio*rr; rho=max(rho/2,1e-6)
            end
            rho==rho_old || foreach(u->u.*=rho_old/rho,dual)
        end
        if out!==nothing && checkpoint_interval>0 && k%checkpoint_interval==0
            checkpoint_s+=save_admm_checkpoint(replace(out,".csv"=>"_checkpoint.jls"),stages,z,dual,rho,k,phase,elapsed,ideal,solver_total,diagnostic,primal_seen,high_primal_count)
        end
    end
    io===nothing || close(io)
    details===nothing || close(details)
    wall_s=(time_ns()-start)/1e9
    @printf("TADMM_DONE phase=%s converged=%s iterations=%d objective=%.12f primal=%.3e rmax=%.3e smax=%.3e algorithm_s=%.6f ideal_parallel_estimate_s=%.6f wall_with_diagnostics_s=%.6f\n",phase,ok,k,metrics.objective,metrics.original_primal,rmax,smax,elapsed,ideal,wall_s); flush(stdout)
    return (;ok,k,metrics,us,bb,rmax,smax,elapsed,ideal,solver_total,diagnostic,wall_s,build_s,checkpoint_s,no_iteration,no_wall_s,no_save_s,
        stopping_rule,historical_iteration,historical_met,owned_r,owned_s,historical_stagnant,dispatch_s,physical_check_s)
end

function main()
    length(ARGS)>=3 || error("usage: run.jl system T output-directory")
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    datafile=joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls")
    get(ENV,"TADMM_LOCAL_SOLVE_POLICY","exact")=="exact" || error("Use compare_july.jl for explicit inexact local-solve assessments")
    d=deserialize(datafile); cb_profile=apply_review_cb!(d,system)
    serialize(joinpath(out,"instance.jls"),d)
    screen=tadmm_screen(d)
    metadata=Dict("system"=>system,"T"=>T,"timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),
        "C_B"=>d[:C_B],"cb_profile"=>cb_profile,"local_solve_policy"=>"exact",
        "terminal_soc_penalty"=>0.0,"periodic"=>true,"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),"hostname"=>gethostname(),
        "data_sha256"=>bytes2hex(sha256(read(datafile))),"linear_solver"=>get(ENV,"TADMM_LINEAR_SOLVER","ma57"),
        "rho_policy"=>get(ENV,"TADMM_RHO_POLICY","phase"),"ipopt_extra_options"=>get(ENV,"TADMM_IPOPT_OPTIONS",""),
        "rho_balance_ratio"=>parse(Float64,get(ENV,"TADMM_RHO_BALANCE_RATIO","10")),
        "rho_initial"=>parse(Float64,get(ENV,"TADMM_RHO_INIT",string(3000.0*sqrt(T/24)))),
        "screening"=>join(sort!(collect(screen.flags)),','),
        "ell_lower_bounds_removed"=>("ell" in screen.flags),
        "nonroot_voltage_lower_bounds_retained"=>true,
        "import_lower_bound_removed"=>screen.drop_ps,
        "nonroot_voltage_upper_bounds_removed"=>("vupper" in screen.flags ? count(n->d[:Nset][n]!=d[:substationBus] && !screen.keep_upper[n],eachindex(d[:Nset])) : 0),
        "no_objective_relative_tolerance"=>0.001,
        "no_original_primal_tolerance"=>(system=="ieee123C_1ph" ? 1e-6 : 1e-5),
        "no_iteration_target"=>100,
        "period_workers"=>min(T,parse(Int,get(ENV,"TADMM_NWORKERS","0"))),
        "execution_mode"=>parse(Int,get(ENV,"TADMM_NWORKERS","0"))>0 ? "distributed" : "serial",
        "central_openmp_threads"=>get(ENV,"OMP_NUM_THREADS","unknown"),
        "cpu_ids"=>get(ENV,"TADMM_CPU_IDS",""),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","none"),
        "primal_tolerance"=>parse(Float64,get(ENV,"TADMM_PRIMAL_TOL","2.5e-7")),
        "dual_tolerance"=>parse(Float64,get(ENV,"TADMM_DUAL_TOL","1e-4")),
        "timing_note"=>"algorithm_s is measured stage round plus coordination, including process messages in distributed mode. Diagnostics/checkpoints are separate; loop wall time includes them. ideal_parallel_estimate_s is always hypothetical.")
    for f in ("model.jl","parallel_workers.jl","worker_julia.sh","run.jl","Project.toml","Manifest.toml")
        metadata[f*"_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,f))))
    end
    open(io->TOML.print(io,metadata),joinpath(out,"metadata.toml"),"w")
    # Centralized reference: full warm-up, then exactly the original starts.
    build0=time_ns(); central=make_model(d,collect(1:T)); central_build=(time_ns()-build0)/1e9
    set_optimizer_attribute(central.m,"output_file",joinpath(out,"central_ipopt_warmup.log"))
    set_optimizer_attribute(central.m,"file_print_level",5)
    set_optimizer_attribute(central.m,"print_timing_statistics","yes")
    set_optimizer_attribute(central.m,"inf_pr_output","original")
    central_trace=Tuple{Int,Int,Float64,Float64,Float64,Float64,Float64}[]
    central_clock=Ref(time_ns())
    function central_callback(alg_mod::Cint,iter::Cint,obj::Float64,inf_pr::Float64,inf_du::Float64,mu::Float64,
                              dnorm::Float64,reg::Float64,alpha_du::Float64,alpha_pr::Float64,ls::Cint)
        push!(central_trace,(Int(alg_mod),Int(iter),obj,inf_pr,inf_du,mu,(time_ns()-central_clock[])/1e9))
        return true
    end
    MOI.set(central.m,Ipopt.CallbackFunction(),central_callback)
    phase_event(out,"central_warmup_begin")
    central_clock[]=time_ns(); warm0=time_ns(); solve_checked!(central); central_warm=(time_ns()-warm0)/1e9
    phase_event(out,"central_warmup_end")
    warm_initial_objective=central_trace[1][3]
    warm_iters=MOI.get(central.m,MOI.BarrierIterations())
    set_start_value.(central.vars,central.initial)
    set_optimizer_attribute(central.m,"output_file",joinpath(out,"central_ipopt_timed.log"))
    empty!(central_trace); central_clock[]=time_ns()
    phase_event(out,"central_measured_begin")
    central_clock[]=time_ns()
    tick=time_ns(); central_solver=solve_checked!(central); central_wall=(time_ns()-tick)/1e9
    phase_event(out,"central_measured_end")
    central_iters=MOI.get(central.m,MOI.BarrierIterations())
    abs(warm_initial_objective-central_trace[1][3])<=1e-10*max(1.0,abs(warm_initial_objective)) || error("central initial objective changed after warm-up")
    open(joinpath(out,"central_trace.csv"),"w") do io
        println(io,"algorithm_mode,iteration,objective,internal_primal,dual,barrier,elapsed_wall_s")
        for row in central_trace; println(io,join(row,',')); end
    end
    cu,cb=dispatch([central],d); cm=physical_metrics(d,cu,cb)
    cm.original_primal<=1e-6 || error("centralized reference fails original model")
    serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>cu,:B=>cb,:status=>"LOCALLY_SOLVED"))
    @printf("CENTRAL_DONE objective=%.12f iterations=%d solver_s=%.6f wall_s=%.6f primal=%.3e\n",cm.objective,central_iters,central_solver,central_wall,cm.original_primal); flush(stdout)
    central=nothing; cu=nothing; cb=nothing; GC.gc()
    maxiter=parse(Int,get(ENV,"TADMM_MAX_ITER","5000"))
    phase_event(out,"admm_warmup_begin")
    no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    warm0=time_ns(); warm=run_admm(d;maxiter=maxiter,phase="warmup",out=joinpath(out,"warmup.csv"),reference_objective=cm.objective,no_primal_tolerance=no_tol,system_name=system); warm_s=(time_ns()-warm0)/1e9
    phase_event(out,"admm_warmup_end")
    warm.ok || error("ADMM warm-up did not converge; preserve trace and diagnose before timing")
    serialize(joinpath(out,"warmup_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>warm.metrics.objective,:u=>warm.us,:B=>warm.bb,:status=>"ADMM_CONVERGED"))
    metadata["worker_layout"]=copy(POOL_LAYOUT)
    metadata["coordinator_affinity"]=affinity_string()
    open(io->TOML.print(io,metadata),joinpath(out,"metadata.toml"),"w")
    GC.gc()
    phase_event(out,"admm_measured_begin")
    result=run_admm(d;maxiter=maxiter,out=joinpath(out,"trace.csv"),reference_objective=cm.objective,no_primal_tolerance=no_tol,system_name=system)
    phase_event(out,"admm_measured_end")
    result.ok || error("ADMM timed solve did not converge")
    result.k==warm.k || error("warm-up and timed ADMM iterations differ")
    result.no_iteration==warm.no_iteration || error("warm-up and timed NO crossings differ")
    gap=abs(result.metrics.objective-cm.objective)/max(1.0,abs(cm.objective))
    serialize(joinpath(out,"tadmm_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,:u=>result.us,:B=>result.bb,:status=>"ADMM_CONVERGED"))
    summary=Dict("system"=>system,"T"=>T,"converged"=>true,"iterations"=>result.k,
        "central_objective"=>cm.objective,"tadmm_objective"=>result.metrics.objective,"relative_objective_gap"=>gap,
        "central_solver_s"=>central_solver,"central_wall_s"=>central_wall,"central_warmup_s"=>central_warm,"central_build_s"=>central_build,
        "central_warmup_iterations"=>warm_iters,"central_timed_iterations"=>central_iters,
        "execution_mode"=>metadata["execution_mode"],"period_workers"=>metadata["period_workers"],
        "tadmm_algorithm_s"=>result.elapsed,"ideal_parallel_estimate_s"=>result.ideal,
        "tadmm_wall_with_diagnostics_s"=>result.wall_s,"tadmm_warmup_s"=>warm_s,"tadmm_build_s"=>result.build_s,
        "optimizer_loop_wall_s"=>result.wall_s,
        "checkpoint_s"=>result.checkpoint_s,
        "no_iteration"=>result.no_iteration,"no_target_under_100"=>(0<result.no_iteration<100),
        "no_arrival_wall_s"=>result.no_wall_s,"no_dispatch_save_s"=>result.no_save_s,
        "no_original_primal_tolerance"=>no_tol,"no_objective_relative_tolerance"=>0.001,
        "tadmm_original_primal"=>result.metrics.original_primal,"tadmm_consensus_max"=>result.rmax,"tadmm_dual_max"=>result.smax,
        "tadmm_soc_equality_gap"=>result.metrics.soc_equality_gap,"central_soc_equality_gap"=>cm.soc_equality_gap,
        "tadmm_fictitious_loss_kwh"=>result.metrics.fictitious_loss_kwh)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    gap<=1e-5 || error("converged ADMM objective disagrees with centralized reference: $gap")
    println("PILOT_OPTIMIZATION_COMPLETE gap=",gap," output=",out)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    try
        main()
    finally
        shutdown_pool!()
    end
end
