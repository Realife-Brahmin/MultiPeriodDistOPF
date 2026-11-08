# Bounded NO calibration, not a clean timed/full-convergence benchmark.
include("run.jl")
system,Tstr,out,reference_case=ARGS
T=parse(Int,Tstr); T==24 || error("This calibration is T=24 only")
mkpath(out)
datafile=joinpath(SHARED,"..","..","results/network_filterddp/network_data_$(system)_T$(T)_periodic.jls")
d=deserialize(datafile); d[:C_B]=battery_cb(system,"system")
refmeta=TOML.parsefile(joinpath(reference_case,"metadata.toml"))
ref=TOML.parsefile(joinpath(reference_case,"summary.toml"))
isfile(joinpath(reference_case,"COMPLETE")) || error("Reference not independently validated")
bytes2hex(sha256(read(datafile)))==refmeta["data_sha256"] || error("Calibration/reference input mismatch")
d[:C_B]==refmeta["C_B"] && refmeta["terminal_soc_penalty"]==0 || error("Calibration problem mismatch")
ENV["TADMM_LINEAR_SOLVER"]=="ma57" || error("MA57 only")
no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
try
    result=run_admm(d;maxiter=100,out=joinpath(out,"calibration.csv"),phase="NO_calibration",
        reference_objective=ref["central_objective"],no_primal_tolerance=no_tol,system_name=system)
    meta=Dict("system"=>system,"T"=>T,"linear_solver"=>"ma57",
        "rho_initial"=>parse(Float64,ENV["TADMM_RHO_INIT"]),"rho_policy"=>ENV["TADMM_RHO_POLICY"],
        "rho_balance_ratio"=>parse(Float64,get(ENV,"TADMM_RHO_BALANCE_RATIO","10")),
        "screening"=>join(sort!(collect(tadmm_screen(d).flags)),','),
        "reference_objective"=>ref["central_objective"],"reference_input_sha256"=>refmeta["data_sha256"],
        "reference_summary_sha256"=>bytes2hex(sha256(read(joinpath(reference_case,"summary.toml")))),
        "no_iteration"=>result.no_iteration,"target_under_100"=>(0<result.no_iteration<100),
        "strict_converged"=>result.ok,"iterations_executed"=>result.k,
        "no_original_tolerance"=>no_tol,"no_objective_relative_tolerance"=>0.001,
        "timing_note"=>"One capped 100-step calibration including first-execution JIT/model setup. Not warmed benchmark timing; not a full-convergence claim.")
    if result.no_iteration>0
        nodispatch=deserialize(joinpath(out,"calibration_no_dispatch.jls"))
        metrics=physical_metrics(d,nodispatch[:u],nodispatch[:B])
        meta["no_original_primal"]=metrics.original_primal
        meta["no_objective_gap"]=abs(metrics.objective-ref["central_objective"])/max(1,abs(ref["central_objective"]))
    end
    open(io->TOML.print(io,meta),joinpath(out,"calibration.toml"),"w")
    println("NO_CALIBRATION_RECORDED target_under_100=",meta["target_under_100"])
finally
    shutdown_pool!()
end
