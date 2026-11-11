include("run.jl")
include("sweep_milestones.jl")
function horizon_sweep_main()
system,Targ,out,rhoarg=ARGS
T=parse(Int,Targ); rhos=parse.(Float64,split(rhoarg,','))
@assert ((system=="large10kC_1ph" && T in (24,48,96)) || (system=="ieee2522C_1ph" && T==96)) && length(rhos) in (1,2)
solver=ENV["TADMM_LINEAR_SOLVER"]
@assert solver in ("ma57","mumps") && ENV["TADMM_RHO_POLICY"] in ("legacy_phase","bootstrap_fixed")
@assert ENV["TADMM_RHO_WATCHDOG"]=="1" && ENV["TADMM_RHO_PRIMAL_GATE"]=="1e-4"
mkpath(out);register_process!(); total_start=time_ns()
datafile=joinpath(SHARED,"..","..","results/network_filterddp/network_data_$(system)_T$(T)_periodic.jls")
d=deserialize(datafile); apply_review_cb!(d,system;profile="july")
ENV["REDUCED_CB"]=repr(d[:C_B]); ENV["REDUCED_PROFILE"]="periodic"
ENV["TERMINAL_SOC_SOFT"]="0"; ENV["OPENDSS_ZERO_LINE_C"]="1"
serialize(joinpath(out,"instance.jls"),d)
write_toml(path,obj)=open(io->TOML.print(io,obj),path,"w")
maxiter=parse(Int,get(ENV,"TADMM_MAX_ITER","2000"))
function check_dss(path,label,log)
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(SHARED,"opendss_check.jl")) $system $T $path $label`
    open(log,"w") do io
        return success(pipeline(cmd;stdout=io,stderr=io))
    end
end
function check_sequence(path,label)
    folder=joinpath(out,label*"_sequential")
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(@__DIR__,"run_sequential_dss.jl")) $(joinpath(out,"instance.jls")) $path $folder model`
    ok=withenv("SEQUENTIAL_DSS_GZIP"=>"1") do
        open(joinpath(out,label*"_sequential.log"),"w") do io
            success(pipeline(cmd;stdout=io,stderr=io))
        end
    end
    ok || error("Chronological OpenDSS execution failed; see $(label)_sequential.log")
    return TOML.parsefile(joinpath(folder,"summary.toml")),relpath(joinpath(folder,"summary.toml"),out)
end
nw=parse(Int,get(ENV,"TADMM_NWORKERS","24"))
meta=Dict("system"=>system,"T"=>T,"C_B"=>d[:C_B],"rhos"=>rhos,"linear_solver"=>solver,
    "input_sha256"=>bytes2hex(sha256(read(datafile))),"workers"=>nw,"compute_cpu_budget"=>nw+1,
    "cluster_backend"=>get(ENV,"TADMM_CLUSTER_BACKEND","local"),
    "warm_start"=>get(ENV,"TADMM_WARM_START","primal"),
    "objective_update"=>get(ENV,"TADMM_OBJECTIVE_UPDATE","legacy"),
    "execution_policy"=>get(ENV,"TADMM_EXECUTION_POLICY","fixed24"),
    "rho_policy"=>ENV["TADMM_RHO_POLICY"],"rho_floor"=>parse(Float64,get(ENV,"TADMM_RHO_FLOOR","48")),
    "bootstrap_decrease"=>parse(Float64,get(ENV,"TADMM_BOOTSTRAP_DECREASE","5")),
    "bootstrap_gate"=>parse(Float64,get(ENV,"TADMM_BOOTSTRAP_GATE","Inf")),
    "local_solve_policy"=>get(ENV,"TADMM_LOCAL_SOLVE_POLICY","exact"),
    "warm_retry_limit"=>parse(Int,get(ENV,"TADMM_WARM_RETRY_LIMIT","20")),
    "cap_end_iteration"=>parse(Int,get(ENV,"TADMM_CAP_END_ITERATION","30")),
    "cap_iterations"=>parse(Int,get(ENV,"TADMM_CAP_ITERATIONS","30")),
    "cap_primal_tolerance"=>parse(Float64,get(ENV,"TADMM_CAP_PRIMAL_TOL","1e-6")),
    "cap_screen"=>get(ENV,"TADMM_CAP_SCREEN","feasible"),
    "nodes"=>parse(Int,get(ENV,"SLURM_JOB_NUM_NODES","1")),
    "configured_memory_gb"=>128,"maximum_iterations"=>maxiter,"controller_primal_gate"=>1e-4,
    "rho_tau_decrease"=>5.0,"rho_tau_increase"=>2.0,"rho_balance_ratio"=>5.0,
    "rho_update_interval"=>5,"rho_watchdog_window"=>15,"rho_watchdog_factor"=>2.5,
    "sequential_validation_version"=>1,
    "final_owned_rms"=>1e-5,"final_all_copy_rms"=>1e-5,"final_dual_rms"=>1e-2,
    "anchor_provenance"=>get(ENV,"TADMM_RHO_PROVENANCE",system=="large10kC_1ph" && T==96 ? "60000 is extrapolated; no July large10k T96 result" : "July paper and archived summary winner"),
    "timer_note"=>"Core stage rounds plus coordination exclude diagnostics, milestone I/O, checkpoints and validation. Actual loop time, warmup/build and tuning cost logged separately. No best-of-sweep timing reported as untuned.",
    "screening"=>join(sort!(collect(tadmm_screen(d).flags)),','),"periodic"=>true,"terminal_soc_penalty"=>0.0)
try
    phase_event(out,"central_build_begin")
    tick=time_ns(); central=make_model(d,collect(1:T)); meta["central_build_s"]=(time_ns()-tick)/1e9
    starts=bytes2hex(sha256(reinterpret(UInt8,Float64.(central.initial))))
    set_optimizer_attribute(central.m,"output_file",joinpath(out,"central_warmup.log"))
    set_optimizer_attribute(central.m,"file_print_level",5)
    set_optimizer_attribute(central.m,"inf_pr_output","original")
    central_trace=Tuple{Int,Float64,Float64,Float64,Float64,Float64}[]
    central_clock=Ref(time_ns())
    function callback(alg::Cint,k::Cint,obj::Float64,pr::Float64,du::Float64,mu::Float64,
            dn::Float64,reg::Float64,ad::Float64,ap::Float64,ls::Cint)
        push!(central_trace,(Int(k),(time_ns()-central_clock[])/1e9,obj,pr,du,mu))
        return true
    end
    MOI.set(central.m,Ipopt.CallbackFunction(),callback)
    phase_event(out,"central_warmup_begin")
    central_clock[]=time_ns()
    tick=time_ns(); warm_solver=solve_checked!(central); meta["central_warmup_wall_s"]=(time_ns()-tick)/1e9
    meta["central_warmup_solver_s"]=warm_solver
    meta["central_warmup_iterations"]=MOI.get(central.m,MOI.BarrierIterations())
    initial_objective=central_trace[1][3]
    empty!(central_trace)
    set_start_value.(central.vars,central.initial)
    @assert starts==bytes2hex(sha256(reinterpret(UInt8,Float64.(central.initial))))
    set_optimizer_attribute(central.m,"output_file",joinpath(out,"central_timed.log"))
    phase_event(out,"central_measured_begin")
    central_clock[]=time_ns()
    tick=time_ns(); meta["central_solver_s"]=solve_checked!(central); meta["central_wall_s"]=(time_ns()-tick)/1e9
    meta["central_iterations"]=MOI.get(central.m,MOI.BarrierIterations())
    phase_event(out,"central_measured_end")
    @assert isapprox(initial_objective,central_trace[1][3];atol=1e-10,rtol=1e-10)
    open(joinpath(out,"central_trace.csv"),"w") do io
        println(io,"iteration,elapsed_s,objective,internal_primal,dual,barrier")
        for row in central_trace;println(io,join(row,','));end
    end
    u,b=dispatch([central],d); cm=physical_metrics(d,u,b)
    @assert cm.original_primal<=1e-6
    for key in keys(cm);meta["central_"*string(key)]=cm[key];end
    meta["central_initial_starts_sha256"]=starts
    serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>u,:B=>b,:status=>"LOCALLY_SOLVED"))
    @assert check_dss(joinpath(out,"central_dispatch.jls"),"central",joinpath(out,"central_opendss.log"))
    cs,cs_path=check_sequence(joinpath(out,"central_dispatch.jls"),"central")
    meta["central_sequential_schedule_pass"]=cs["sequential_schedule_pass"]
    meta["central_sequential_summary"]=cs_path
    write_toml(joinpath(out,"sweep.toml"),meta)
    central=nothing;u=nothing;b=nothing;GC.gc()
    previous_starts=nothing
    for rho in rhos
        ENV["TADMM_RHO_INIT"]=repr(rho)
        arm=joinpath(out,"rho"*string(Int(rho)));mkpath(arm)
        warmdir=joinpath(arm,"warmup");mkpath(warmdir)
        warmobserve,warmmarks=sweep_observer(warmdir,d,system,cm.objective;save_points=false)
        phase_event(out,"rho$(Int(rho))_warmup_begin")
        tick=time_ns();warm=run_admm(d;maxiter=maxiter,out=joinpath(warmdir,"trace.csv"),phase="warmup_$(Int(rho))",
            system_name=system,reference_objective=cm.objective,stopping_rule="observer",observer=warmobserve)
        warmwall=(time_ns()-tick)/1e9
        warmstarts=Dict(t=>h for worker in POOL_LAYOUT for (t,h) in worker["initial_starts_sha256"])
        observe,marks=sweep_observer(arm,d,system,cm.objective)
        phase_event(out,"rho$(Int(rho))_measured_begin")
        result=run_admm(d;maxiter=maxiter,out=joinpath(arm,"trace.csv"),phase="timed_$(Int(rho))",
            system_name=system,reference_objective=cm.objective,stopping_rule="observer",observer=observe)
        phase_event(out,"rho$(Int(rho))_measured_end")
        starts_now=Dict(t=>h for worker in POOL_LAYOUT for (t,h) in worker["initial_starts_sha256"])
        @assert warmstarts==starts_now
        @assert previous_starts===nothing || previous_starts==starts_now
        previous_starts=starts_now
        endpoint=joinpath(arm,"endpoint.jls")
        serialize(endpoint,Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,:u=>result.us,:B=>result.bb,:status=>"SWEEP_ENDPOINT"))
        assessment=Dict("rho_initial"=>rho,"iterations"=>result.k,"warmup_iterations"=>warm.k,
            "warmup_total_wall_s"=>warmwall,"warmup_core_s"=>warm.elapsed,"build_s"=>result.build_s,
            "core_s"=>result.elapsed,"loop_wall_s"=>result.wall_s,"final_rms_met"=>haskey(marks,"final_rms"),
            "warmup_final_rms_met"=>haskey(warmmarks,"final_rms"),"strict_criteria_met"=>result.ok,
            "original_primal"=>result.metrics.original_primal,"objective"=>result.metrics.objective,
            "objective_gap"=>abs(result.metrics.objective-cm.objective)/abs(cm.objective),
            "owned_r_rms"=>result.owned_r,"all_r_max"=>result.rmax,"owned_s_rms"=>result.owned_s,
            "dispatch_s"=>result.dispatch_s,"physical_check_s"=>result.physical_check_s,
            "worker_layout"=>deepcopy(POOL_LAYOUT))
        write_toml(joinpath(arm,"assessment.toml"),assessment)
        # Recheck all saved first-crossing dispatches, once per unique iteration.
        logs=Dict{String,Any}()
        for point in sort!(unique([row["snapshot"] for row in values(marks)]))
            saved=deserialize(joinpath(arm,point)); audited=physical_metrics(d,saved[:u],saved[:B])
            recorded=first(row for row in values(marks) if row["snapshot"]==point)
            @assert isapprox(audited.objective,recorded["objective"];rtol=1e-12,atol=1e-10)
            @assert isapprox(audited.original_primal,recorded["original_primal"];rtol=1e-12,atol=1e-10)
            label=splitext(point)[1];log=joinpath(arm,label*"_opendss.log")
            logs[point]=check_dss(joinpath(arm,point),label,log)
        end
        finalpoints=[row["snapshot"] for row in values(marks) if row["iteration"]==result.k]
        if isempty(finalpoints)
            logs["endpoint.jls"]=check_dss(endpoint,"endpoint",joinpath(arm,"endpoint_opendss.log"))
        else
            point=first(finalpoints)
            cp(joinpath(arm,splitext(point)[1]*"_opendss.log"),joinpath(arm,"endpoint_opendss.log"))
            logs["endpoint.jls"]=logs[point]
        end
        write_toml(joinpath(arm,"opendss_status.toml"),logs)
        ss,ss_path=check_sequence(endpoint,basename(arm)*"_endpoint")
        assessment["sequential_schedule_pass"]=ss["sequential_schedule_pass"]
        assessment["sequential_summary"]=ss_path
        write_toml(joinpath(arm,"assessment.toml"),assessment)
        phase_event(out,"rho$(Int(rho))_validated")
        run(`$(get(ENV,"PYTHON","python3")) $(joinpath(@__DIR__,"validate_horizon_sweep.py")) $out $(basename(arm))`)
    end
    meta["campaign_total_wall_s"]=(time_ns()-total_start)/1e9
    write_toml(joinpath(out,"sweep.toml"),meta)
    run(`$(get(ENV,"PYTHON","python3")) $(joinpath(@__DIR__,"validate_horizon_sweep.py")) $out complete`)
finally
    shutdown_pool!()
end
end
horizon_sweep_main()
