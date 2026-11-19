# Timing-only check: same cap, with redundant per-local diagnostic work disabled.
include("run.jl")
system,Targ,out,reference=ARGS
T=parse(Int,Targ); T==24 || error("T24 only")
isfile(joinpath(reference,"ASSESSMENT_VALIDATED")) || error("Reference pair not validated")
prior=TOML.parsefile(joinpath(reference,"comparison.toml"))
ENV["TADMM_LINEAR_SOLVER"]=="ma57" || error("MA57 only")
dfile=joinpath(SHARED,"..","..","results/network_filterddp/network_data_$(system)_T24_periodic.jls")
d=deserialize(dfile); apply_review_cb!(d,system;profile="july")
@assert prior["system"]==system && prior["C_B"]==d[:C_B]
@assert prior["input_sha256"]==bytes2hex(sha256(read(dfile)))
@assert prior["rho_initial"]==parse(Float64,ENV["TADMM_RHO_INIT"]) && prior["rho_policy"]==ENV["TADMM_RHO_POLICY"]
@assert prior["rho_balance_ratio"]==parse(Float64,ENV["TADMM_RHO_BALANCE_RATIO"])
@assert prior["screening"]==join(sort!(collect(tadmm_screen(d).flags)),',')
ENV["REDUCED_CB"]=repr(d[:C_B]); ENV["TERMINAL_SOC_SOFT"]="0"
ENV["REDUCED_PROFILE"]="periodic"; ENV["OPENDSS_ZERO_LINE_C"]="1"
ENV["TADMM_LOCAL_SOLVE_POLICY"]="cap30"; ENV["TADMM_LOCAL_PRIMAL_REPORT"]="0"
mkpath(out); mkpath(joinpath(out,"comparison_reference")); mkpath(joinpath(out,"cap30"))
cp(joinpath(reference,"central_opendss.log"),joinpath(out,"central_opendss.log"))
for file in ("assessment.toml","trace.csv")
    cp(joinpath(reference,"cap30",file),joinpath(out,"comparison_reference",file))
end
cp(joinpath(reference,"comparison.toml"),joinpath(out,"comparison_reference","comparison.toml"))
prior["policies"]=["cap30"]; prior["local_primal_report"]=false
prior["reference_pair_path"]=reference
prior["reference_pair_metadata_sha256"]=bytes2hex(sha256(read(joinpath(reference,"comparison.toml"))))
write_toml(path,obj)=open(io->TOML.print(io,obj),path,"w")
write_toml(joinpath(out,"comparison.toml"),prior)
validate(stage)=run(`$(get(ENV,"PYTHON","python3")) $(joinpath(@__DIR__,"validate_july_comparison.py")) $out $stage`)
function check_dss(path,label,log)
    cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(@__DIR__) $(joinpath(SHARED,"opendss_check.jl")) $system $T $path $label`
    open(log,"w") do io
        run(pipeline(cmd;stdout=io,stderr=io))
    end
end
try
    validate("reference")
    arm=joinpath(out,"cap30")
    run_admm(d;maxiter=3,out=joinpath(arm,"warmup.csv"),phase="cap30_light_compile",system_name=system)
    result=run_admm(d;maxiter=prior["maximum_iterations"],out=joinpath(arm,"trace.csv"),phase="cap30_light",
        reference_objective=prior["reference_objective"],no_primal_tolerance=prior["no_original_tolerance"],system_name=system)
    serialize(joinpath(arm,"final_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,:u=>result.us,:B=>result.bb,:status=>"ASSESSMENT_ENDPOINT"))
    row=Dict("policy"=>"cap30","iterations"=>result.k,"global_strict_criteria_met"=>result.ok,
        "no_iteration"=>result.no_iteration,"target_under_100"=>(0<result.no_iteration<100),
        "objective"=>result.metrics.objective,"original_primal"=>result.metrics.original_primal,
        "objective_gap"=>abs(result.metrics.objective-prior["reference_objective"])/max(1,abs(prior["reference_objective"])),
        "r_max"=>result.rmax,"s_max"=>result.smax,"algorithm_s"=>result.elapsed,
        "loop_wall_s"=>result.wall_s,"build_s"=>result.build_s,"no_wall_s"=>result.no_wall_s,"worker_layout"=>deepcopy(POOL_LAYOUT))
    if result.no_iteration>0
        path=joinpath(arm,"trace_no_dispatch.jls"); no=deserialize(path); m=physical_metrics(d,no[:u],no[:B])
        row["no_objective"]=m.objective; row["no_original_primal"]=m.original_primal
        row["no_dispatch_sha256"]=bytes2hex(sha256(read(path)))
        check_dss(path,"cap30_light_no",joinpath(arm,"no_opendss.log"))
    end
    if result.ok
        check_dss(joinpath(arm,"final_dispatch.jls"),"cap30_light_final",joinpath(arm,"final_opendss.log"))
    end
    write_toml(joinpath(arm,"assessment.toml"),row)
    validate("complete")
finally
    shutdown_pool!()
end
