# Controlled July-weight/cap experiment on the frozen CURRENT periodic inputs.
# It is not a replay of July's network exports, profiles or stopping rules.
include("run.jl")
system,Targ,out,maxarg=ARGS
T=parse(Int,Targ); maxiter=parse(Int,maxarg)
T==24 && maxiter>=100 || error("This first comparison requires T24 and at least 100 steps")
ENV["TADMM_LINEAR_SOLVER"]=="ma57" || error("MA57 only")
mkpath(out)
datafile=joinpath(SHARED,"..","..","results/network_filterddp/network_data_$(system)_T$(T)_periodic.jls")
d=deserialize(datafile); profile=apply_review_cb!(d,system;profile="july")
# The independent checker supports an explicit numeric weight. Do not allow
# REDUCED_CB=system to silently re-score these dispatches with DDP's weight.
ENV["REDUCED_CB"]=repr(d[:C_B])
ENV["TERMINAL_SOC_SOFT"]="0"; ENV["OPENDSS_ZERO_LINE_C"]="1"; ENV["REDUCED_PROFILE"]="periodic"
serialize(joinpath(out,"instance.jls"),d)
no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
meta=Dict("system"=>system,"T"=>T,"C_B"=>d[:C_B],"cb_profile"=>profile,
    "delta_t_h"=>d[:delta_t_h],"kVA_B"=>d[:kVA_B],
    "july_cb_source"=>JULY_CB_SOURCE,"terminal_soc_penalty"=>0.0,
    "input_sha256"=>bytes2hex(sha256(read(datafile))),"profiles"=>"current periodic, unchanged",
    "linear_solver"=>"ma57","rho_initial"=>parse(Float64,ENV["TADMM_RHO_INIT"]),
    "rho_policy"=>ENV["TADMM_RHO_POLICY"],"rho_balance_ratio"=>parse(Float64,get(ENV,"TADMM_RHO_BALANCE_RATIO","10")),
    "screening"=>join(sort!(collect(tadmm_screen(d).flags)),','),"maximum_iterations"=>maxiter,
    "local_primal_report"=>(get(ENV,"TADMM_LOCAL_PRIMAL_REPORT","1")=="1"),
    "no_objective_tolerance"=>0.001,"no_original_tolerance"=>no_tol,
    "warmup_steps_per_arm"=>3,"timing_note"=>"One assessment per arm after three discarded compile steps and fresh model/state reset. Not a full-sequence-warmed benchmark or repeat median.")
function write_toml(path,object)
    open(io->TOML.print(io,object),path,"w")
end
function save_dispatch(path,result,status)
    serialize(path,Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,
        :u=>result.us,:B=>result.bb,:status=>status))
end
function dss_check(path,label,log)
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(SHARED,"opendss_check.jl")) $system $T $path $label`
    open(log,"w") do io
        run(pipeline(cmd;stdout=io,stderr=io))
    end
end
function validate(stage)
    run(`$(get(ENV,"PYTHON","python3")) $(joinpath(@__DIR__,"validate_july_comparison.py")) $out $stage`)
end
try
    ENV["TADMM_LOCAL_SOLVE_POLICY"]="exact"
    tick=time_ns(); central=make_model(d,collect(1:T)); central_build=(time_ns()-tick)/1e9
    tick=time_ns(); central_solver=solve_checked!(central); central_wall=(time_ns()-tick)/1e9
    cu,cb=dispatch([central],d); cm=physical_metrics(d,cu,cb)
    cm.original_primal<=1e-6 || error("July-weight reference fails original physics")
    serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>cu,:B=>cb,:status=>"LOCALLY_SOLVED"))
    meta["reference_objective"]=cm.objective
    meta["reference_original_primal"]=cm.original_primal
    meta["reference_solver_s"]=central_solver; meta["reference_wall_s"]=central_wall
    meta["reference_build_s"]=central_build
    write_toml(joinpath(out,"comparison.toml"),meta)
    dss_check(joinpath(out,"central_dispatch.jls"),"central",joinpath(out,"central_opendss.log"))
    validate("reference")
    central=nothing; GC.gc()
    for policy in ("exact","cap30")
        ENV["TADMM_LOCAL_SOLVE_POLICY"]=policy
        arm=joinpath(out,policy); mkpath(arm)
        # Compilation steps are discarded. prepare_pool! rebuilds all models
        # and resets every primal start, z and dual for the assessment below.
        run_admm(d;maxiter=3,out=joinpath(arm,"warmup.csv"),phase=policy*"_compile",system_name=system)
        result=run_admm(d;maxiter=maxiter,out=joinpath(arm,"trace.csv"),phase=policy,
            reference_objective=cm.objective,no_primal_tolerance=no_tol,system_name=system)
        save_dispatch(joinpath(arm,"final_dispatch.jls"),result,"ASSESSMENT_ENDPOINT")
        row=Dict("policy"=>policy,"iterations"=>result.k,"global_strict_criteria_met"=>result.ok,
            "no_iteration"=>result.no_iteration,"target_under_100"=>(0<result.no_iteration<100),
            "objective"=>result.metrics.objective,"original_primal"=>result.metrics.original_primal,
            "objective_gap"=>abs(result.metrics.objective-cm.objective)/max(1,abs(cm.objective)),
            "r_max"=>result.rmax,"s_max"=>result.smax,"algorithm_s"=>result.elapsed,
            "loop_wall_s"=>result.wall_s,"build_s"=>result.build_s,"no_wall_s"=>result.no_wall_s,
            "worker_layout"=>deepcopy(POOL_LAYOUT))
        if result.no_iteration>0
            nop=joinpath(arm,"trace_no_dispatch.jls")
            no=deserialize(nop); nm=physical_metrics(d,no[:u],no[:B])
            row["no_objective"]=nm.objective; row["no_original_primal"]=nm.original_primal
            row["no_dispatch_sha256"]=bytes2hex(sha256(read(nop)))
            dss_check(nop,policy*"_no",joinpath(arm,"no_opendss.log"))
        end
        if result.ok
            dss_check(joinpath(arm,"final_dispatch.jls"),policy*"_final",joinpath(arm,"final_opendss.log"))
        end
        write_toml(joinpath(arm,"assessment.toml"),row)
        validate(policy)
        println("JULY_ARM_RECORDED policy=",policy," NO=",result.no_iteration," strict=",result.ok)
        flush(stdout)
    end
    validate("complete")
finally
    shutdown_pool!()
end
