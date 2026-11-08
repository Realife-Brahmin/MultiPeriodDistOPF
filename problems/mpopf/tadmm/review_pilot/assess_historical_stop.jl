# One exact-local trajectory, reusing the validated July-weight reference.
include("run.jl")
system,Targ,out,reference=ARGS
T=parse(Int,Targ)
T==24 && system=="large10kC_1ph" || error("First historical-stop assessment is large10k T24")
ENV["TADMM_LINEAR_SOLVER"]=="ma57" || error("MA57 only")
isfile(joinpath(reference,"ASSESSMENT_VALIDATED")) || error("Unvalidated reference")
prior=TOML.parsefile(joinpath(reference,"comparison.toml"))
dfile=joinpath(SHARED,"..","..","results/network_filterddp/network_data_$(system)_T24_periodic.jls")
d=deserialize(dfile); apply_review_cb!(d,system;profile="july")
@assert prior["system"]==system && prior["C_B"]==d[:C_B]
@assert prior["input_sha256"]==bytes2hex(sha256(read(dfile)))
@assert prior["rho_initial"]==parse(Float64,ENV["TADMM_RHO_INIT"])
@assert prior["rho_policy"]==ENV["TADMM_RHO_POLICY"]=="fixed"
@assert prior["screening"]==join(sort!(collect(tadmm_screen(d).flags)),',')
ENV["TADMM_LOCAL_SOLVE_POLICY"]="exact"
ENV["REDUCED_CB"]=repr(d[:C_B]); ENV["TERMINAL_SOC_SOFT"]="0"
ENV["REDUCED_PROFILE"]="periodic"; ENV["OPENDSS_ZERO_LINE_C"]="1"
mkpath(out); serialize(joinpath(out,"instance.jls"),d)
cp(joinpath(reference,"exact","trace.csv"),joinpath(out,"baseline_trace.csv"))
cp(joinpath(reference,"exact","assessment.toml"),joinpath(out,"baseline_assessment.toml"))
write_toml(path,obj)=open(io->TOML.print(io,obj),path,"w")
function check_dss(path,label,log)
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(SHARED,"opendss_check.jl")) $system $T $path $label`
    open(log,"w") do io
        return success(pipeline(cmd;stdout=io,stderr=io))
    end
end
try
    # Validate the reused reference in this execution environment as well.
    @assert check_dss(joinpath(reference,"central_dispatch.jls"),"historical_reference",joinpath(out,"central_opendss.log"))
    run_admm(d;maxiter=3,out=joinpath(out,"warmup.csv"),phase="historical_compile",system_name=system)
    result=run_admm(d;maxiter=200,out=joinpath(out,"trace.csv"),phase="historical",
        reference_objective=prior["reference_objective"],no_primal_tolerance=1e-5,system_name=system,
        stopping_rule="historical_owned_rms")
    endpoint=joinpath(out,"final_dispatch.jls")
    serialize(endpoint,Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,:u=>result.us,:B=>result.bb,
        :status=>"HISTORICAL_RMS_ASSESSMENT_ENDPOINT"))
    dss_ok=check_dss(endpoint,"historical_endpoint",joinpath(out,"endpoint_opendss.log"))
    meta=Dict("system"=>system,"T"=>T,"C_B"=>d[:C_B],"linear_solver"=>"ma57",
        "input_sha256"=>prior["input_sha256"],"reference_objective"=>prior["reference_objective"],
        "reference_metadata_sha256"=>bytes2hex(sha256(read(joinpath(reference,"comparison.toml")))),
        "rho_initial"=>prior["rho_initial"],"rho_policy"=>prior["rho_policy"],
        "stopping_rule"=>result.stopping_rule,"historical_primal_tolerance"=>1e-4,
        "historical_dual_tolerance"=>1e-2,"stagnation_window"=>5,"stagnation_range_tolerance"=>1e-3,
        "iterations"=>result.k,"historical_stopped"=>result.historical_met,
        "strict_criteria_met"=>result.ok,"physical_NO_iteration"=>result.no_iteration,
        "objective"=>result.metrics.objective,"original_primal"=>result.metrics.original_primal,
        "owned_r_rms"=>result.owned_r,"owned_s_rms"=>result.owned_s,
        "algorithm_s"=>result.elapsed,"loop_wall_s"=>result.wall_s,
        "endpoint_opendss_process_success"=>dss_ok,"worker_layout"=>deepcopy(POOL_LAYOUT),
        "timing_note"=>"Single exploratory trajectory after three compile steps and fresh reset; not a warmed benchmark median.")
    write_toml(joinpath(out,"historical.toml"),meta)
    run(`$(get(ENV,"PYTHON","python3")) $(joinpath(@__DIR__,"validate_historical_assessment.py")) $out`)
finally
    shutdown_pool!()
end
