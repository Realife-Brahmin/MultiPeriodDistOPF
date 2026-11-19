# Standalone checkpoint diagnostic; never used as a fresh benchmark trajectory.
using Serialization, SHA, TOML
source,case,out=ARGS[1:3]
profile=length(ARGS)>=4 ? ARGS[4] : "scaling"
profile in ("scaling","complementarity","robust") || error("Unknown diagnostic profile")
include(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))
BLAS.set_num_threads(1)
isdir(out) && error("Refusing existing diagnostic output: $out")
mkpath(out)
checkpoint=get(ENV,"TADMM_DIAG_CHECKPOINT",joinpath(case,"warmup_checkpoint.jls"))
cp(checkpoint,joinpath(out,"checkpoint.jls"))
state=deserialize(joinpath(out,"checkpoint.jls")); d=deserialize(joinpath(case,"instance.jls"))
bytes2hex(sha256(read(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))))==state[:model_sha256] || error("Checkpoint source mismatch")
open(joinpath(out,"comparison.csv"),"w") do io
    println(io,"backend,variant,scaling,period,checkpoint_iteration,rho,status,ipopt_iterations,objective,all_primal_difference_from_saved_max,local_constraint_violation_max")
    scalings=profile=="scaling" ? ("gradient-based","none") : ("gradient-based",)
    variants=profile=="robust" ? ("default","comp1e-8","adaptive","mc64") : (profile,)
    for backend in ("ma57","ma97"), variant in variants, scaling in scalings, t in 1:d[:T]
        backend=="ma57" && variant=="mc64" && continue
        ENV["TADMM_LINEAR_SOLVER"]=backend
        lib=joinpath(homedir(),"opt/hsl-private/install",backend=="ma57" ? "ma57-3.11.3-metis/lib/libma57.so" : "hsl_ma97-2.8.1-metis-openmp/lib/libhsl_ma97.so")
        extra=backend=="ma57" ? "ma57_pivot_order=2" : "ma97_order=amd"
        ENV["TADMM_IPOPT_OPTIONS"]="hsllib=$lib;linear_system_scaling=none;$extra;nlp_scaling_method=$scaling"
        s=make_model(d,[t];local_window=true,rho=state[:rho])
        profile=="complementarity" && set_optimizer_attribute(s.m,"compl_inf_tol",1e-10)
        if profile=="robust"
            set_optimizer_attribute(s.m,"max_iter",300)
            variant!="default" && set_optimizer_attribute(s.m,"compl_inf_tol",1e-8)
            variant=="adaptive" && set_optimizer_attribute(s.m,"mu_strategy","adaptive")
            variant=="mc64" && set_optimizer_attribute(s.m,"ma97_scaling","mc64")
        end
        saved=state[:stage_primals][t]
        set_start_value.(s.vars,saved)
        set_optimizer_attribute(s.m,"print_level",5)
        stem="$(backend)_$(variant)_$(scaling)_p$(t)"
        set_optimizer_attribute(s.m,"output_file",joinpath(out,stem*".log"))
        rho=state[:rho]
        for (n,j) in enumerate(d[:Bset]),(a,time) in enumerate(s.times)
            set_objective_coefficient(s.m,s.b[j,time],rho*(state[:dual][t][n,a]-state[:z][n,time]))
            set_objective_coefficient(s.m,s.b[j,time],s.b[j,time],rho/2)
        end
        optimize!(s.m)
        status=string(termination_status(s.m))
        iterations=MOI.get(s.m,MOI.BarrierIterations())
        copies=[value(s.b[j,time]) for j in d[:Bset],time in s.times]
        # Compare all primal variables too; the checkpoint is a previous iteration,
        # so a difference is diagnostic evidence, not an exact-replay assertion.
        report=primal_feasibility_report(s.m;atol=0.0)
        viol=maximum(values(report);init=0.0)
        delta=maximum(abs.(value.(s.vars).-saved);init=0.0)
        println(io,join((backend,variant,scaling,t,state[:iteration],rho,status,iterations,objective_value(s.m),delta,viol),',')); flush(io)
        open(joinpath(out,stem*".toml"),"w") do result
            TOML.print(result,Dict("copies"=>[collect(copies[:,i]) for i in axes(copies,2)],"times"=>s.times))
        end
        s=nothing; GC.gc()
    end
end
