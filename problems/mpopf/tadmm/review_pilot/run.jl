include("model.jl")
using SHA, TOML, Dates
BLAS.set_num_threads(1)
Threads.nthreads()==1 || error("This pilot requires one Julia thread: no parallel period solves")

function run_admm(d; maxiter=5000, rho0=3000.0*sqrt(d[:T]/24), out=nothing, phase="timed")
    T=d[:T]; bs=d[:Bset]; nb=length(bs)
    build0=time_ns()
    stages=[make_model(d,[t];local_window=true,rho=rho0) for t in 1:T]
    build_s=(time_ns()-build0)/1e9
    z=[d[:B0_pu][j] for j in bs,t in 1:T]
    dual=[zeros(nb,length(s.times)) for s in stages]
    copies=[zeros(nb,length(s.times)) for s in stages]
    multiplicity=[count(s->t in s.times,stages) for t in 1:T]
    rho=rho0; elapsed=0.0; ideal=0.0; solver_total=0.0; diagnostic=0.0
    rmax=Inf; smax=Inf; ok=false
    tol_r=parse(Float64,get(ENV,"TADMM_PRIMAL_TOL","2.5e-7"))
    tol_s=parse(Float64,get(ENV,"TADMM_DUAL_TOL","1e-4"))
    io=out===nothing ? nothing : open(out,"w")
    io===nothing || println(io,"iteration,objective,original_primal,r_rms,r_max,s_rms,s_max,rho,algorithm_s,ideal_parallel_estimate_s,solver_sum_s,diagnostic_s")
    start=time_ns(); metrics=nothing; us=nothing; bb=nothing; k=0
    for iteration in 1:maxiter
        k=iteration; iter0=time_ns(); times=Float64[]; solve_times=Float64[]
        for t0 in 1:T
            s=stages[t0]; tick=time_ns()
            # Updating an objective invalidates JuMP's cached solution: capture
            # the previous primal point BEFORE changing its coefficients.
            if iteration>1
                vals=value.(s.vars); set_start_value.(s.vars,vals)
            end
            for (n,j) in enumerate(bs),(a,t) in enumerate(s.times)
                set_objective_coefficient(s.m,s.b[j,t],rho*(dual[t0][n,a]-z[n,t]))
                set_objective_coefficient(s.m,s.b[j,t],s.b[j,t],rho/2)
            end
            push!(solve_times,solve_checked!(s))
            for (n,j) in enumerate(bs),(a,t) in enumerate(s.times)
                copies[t0][n,a]=value(s.b[j,t])
            end
            push!(times,(time_ns()-tick)/1e9)
        end
        coord0=time_ns(); zold=copy(z); fill!(z,0.0)
        for t0 in 1:T,(a,t) in enumerate(stages[t0].times),n in 1:nb
            z[n,t]+=copies[t0][n,a]+dual[t0][n,a]
        end
        for t in 1:T,n in 1:nb
            j=bs[n]
            z[n,t]=clamp(z[n,t]/multiplicity[t],d[:soc_min][j]*d[:B_R_pu][j],d[:soc_max][j]*d[:B_R_pu][j])
        end
        r2=0.0; s2=0.0; rmax=0.0; smax=0.0; ncopy=0
        for t0 in 1:T,(a,t) in enumerate(stages[t0].times),n in 1:nb
            r=copies[t0][n,a]-z[n,t]; s=rho*(z[n,t]-zold[n,t])
            dual[t0][n,a]+=r
            r2+=r*r; s2+=s*s; rmax=max(rmax,abs(r)); smax=max(smax,abs(s)); ncopy+=1
        end
        rr=sqrt(r2/ncopy); ss=sqrt(s2/ncopy)
        coord_s=(time_ns()-coord0)/1e9
        elapsed+=(time_ns()-iter0)/1e9
        ideal+=maximum(times)+coord_s
        solver_total+=sum(solve_times)
        # Physical validation is deliberately outside the algorithm timer;
        # end-to-end wall time below includes it, warm-up and setup are separate.
        diag0=time_ns(); us,bb=dispatch(stages,d;z=z); metrics=physical_metrics(d,us,bb)
        diagnostic+=(time_ns()-diag0)/1e9
        if io!==nothing
            println(io,join((k,metrics.objective,metrics.original_primal,rr,rmax,ss,smax,rho,elapsed,ideal,solver_total,diagnostic),',')); flush(io)
        end
        if k==1 || k%25==0
            @printf("TADMM_ITER phase=%s k=%d objective=%.9f physical=%.3e rmax=%.3e smax=%.3e rho=%.5g algorithm_s=%.3f\n",phase,k,metrics.objective,metrics.original_primal,rmax,smax,rho,elapsed); flush(stdout)
        end
        if rmax<=tol_r && smax<=tol_s && metrics.original_primal<=2tol_r
            ok=true; break
        end
        # Standard scaled-dual residual balancing. All copies, including edge
        # multiplicities, enter BOTH residual norms. This replaces July's
        # diagonal-only stopping check and phase-specific heuristics.
        if k%10==0
            rho_old=rho
            if rr>10ss; rho=min(rho*2,1e8)
            elseif ss>10rr; rho=max(rho/2,1e-6)
            end
            rho==rho_old || foreach(u->u.*=rho_old/rho,dual)
        end
    end
    io===nothing || close(io)
    wall_s=(time_ns()-start)/1e9
    @printf("TADMM_DONE phase=%s converged=%s iterations=%d objective=%.12f primal=%.3e rmax=%.3e smax=%.3e algorithm_s=%.6f ideal_parallel_estimate_s=%.6f wall_with_diagnostics_s=%.6f\n",phase,ok,k,metrics.objective,metrics.original_primal,rmax,smax,elapsed,ideal,wall_s); flush(stdout)
    return (;ok,k,metrics,us,bb,rmax,smax,elapsed,ideal,solver_total,diagnostic,wall_s,build_s)
end

function main()
    length(ARGS)>=3 || error("usage: run.jl system T output-directory")
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    datafile=joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls")
    d=deserialize(datafile); d[:C_B]=battery_cb(system,"system")
    serialize(joinpath(out,"instance.jls"),d)
    metadata=Dict("system"=>system,"T"=>T,"timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),
        "C_B"=>d[:C_B],"terminal_soc_penalty"=>0.0,"periodic"=>true,"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),"hostname"=>get(ENV,"COMPUTERNAME","unknown"),
        "data_sha256"=>bytes2hex(sha256(read(datafile))),"linear_solver"=>get(ENV,"TADMM_LINEAR_SOLVER","mumps"),
        "timing_note"=>"Sequential measured wall time. Ideal parallel estimate is sum of per-iteration maximum stage elapsed time plus measured coordinator; never measured parallel speedup.")
    for f in ("model.jl","run.jl","Project.toml","Manifest.toml")
        metadata[f*"_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,f))))
    end
    open(io->TOML.print(io,metadata),joinpath(out,"metadata.toml"),"w")
    # Centralized reference: full warm-up, then exactly the original starts.
    build0=time_ns(); central=make_model(d,collect(1:T)); central_build=(time_ns()-build0)/1e9
    warm0=time_ns(); solve_checked!(central); central_warm=(time_ns()-warm0)/1e9
    warm_iters=MOI.get(central.m,MOI.BarrierIterations())
    set_start_value.(central.vars,central.initial)
    tick=time_ns(); central_solver=solve_checked!(central); central_wall=(time_ns()-tick)/1e9
    central_iters=MOI.get(central.m,MOI.BarrierIterations())
    warm_iters==central_iters || error("central warm-up and timed iteration counts differ")
    cu,cb=dispatch([central],d); cm=physical_metrics(d,cu,cb)
    cm.original_primal<=1e-6 || error("centralized reference fails original model")
    serialize(joinpath(out,"central_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>cm.objective,:u=>cu,:B=>cb,:status=>"LOCALLY_SOLVED"))
    @printf("CENTRAL_DONE objective=%.12f iterations=%d solver_s=%.6f wall_s=%.6f primal=%.3e\n",cm.objective,central_iters,central_solver,central_wall,cm.original_primal); flush(stdout)
    maxiter=parse(Int,get(ENV,"TADMM_MAX_ITER","5000"))
    warm0=time_ns(); warm=run_admm(d;maxiter=maxiter,phase="warmup",out=joinpath(out,"warmup.csv")); warm_s=(time_ns()-warm0)/1e9
    warm.ok || error("ADMM warm-up did not converge; preserve trace and diagnose before timing")
    GC.gc()
    result=run_admm(d;maxiter=maxiter,out=joinpath(out,"trace.csv"))
    result.ok || error("ADMM timed solve did not converge")
    result.k==warm.k || error("warm-up and timed ADMM iterations differ")
    gap=abs(result.metrics.objective-cm.objective)/max(1.0,abs(cm.objective))
    serialize(joinpath(out,"tadmm_dispatch.jls"),Dict(:system=>system,:T=>T,:objective=>result.metrics.objective,:u=>result.us,:B=>result.bb,:status=>"ADMM_CONVERGED"))
    summary=Dict("system"=>system,"T"=>T,"converged"=>true,"iterations"=>result.k,
        "central_objective"=>cm.objective,"tadmm_objective"=>result.metrics.objective,"relative_objective_gap"=>gap,
        "central_solver_s"=>central_solver,"central_wall_s"=>central_wall,"central_warmup_s"=>central_warm,"central_build_s"=>central_build,
        "tadmm_algorithm_s"=>result.elapsed,"ideal_parallel_estimate_s"=>result.ideal,
        "tadmm_wall_with_diagnostics_s"=>result.wall_s,"tadmm_warmup_s"=>warm_s,"tadmm_build_s"=>result.build_s,
        "tadmm_original_primal"=>result.metrics.original_primal,"tadmm_consensus_max"=>result.rmax,"tadmm_dual_max"=>result.smax,
        "tadmm_soc_equality_gap"=>result.metrics.soc_equality_gap,"central_soc_equality_gap"=>cm.soc_equality_gap,
        "tadmm_fictitious_loss_kwh"=>result.metrics.fictitious_loss_kwh)
    open(io->TOML.print(io,summary),joinpath(out,"summary.toml"),"w")
    gap<=1e-5 || error("converged ADMM objective disagrees with centralized reference: $gap")
    println("PILOT_OPTIMIZATION_COMPLETE gap=",gap," output=",out)
end
main()
