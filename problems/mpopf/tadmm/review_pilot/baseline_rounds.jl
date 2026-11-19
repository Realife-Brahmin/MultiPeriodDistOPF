# Round counts of the overlapping-window tADMM (harness defaults) on the same
# instance and near-optimality rule as storage_coordinator.jl. Verification only.
include("run.jl")
system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
profile=apply_review_cb!(d,system)
reference=parse(Float64,read(joinpath(out,"reference_$(profile).txt"),String))
no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
r=run_admm(d;maxiter=parse(Int,get(ENV,"TADMM_MAX_ITER","1000")),out=joinpath(out,"baseline_$(profile)_trace.csv"),
    reference_objective=reference,no_primal_tolerance=no_tol,system_name=system)
@printf("BASELINE_DONE system=%s T=%d strict=%s rounds=%d no=%d ideal_s=%.2f\n",system,T,r.ok,r.k,r.no_iteration,r.ideal)
