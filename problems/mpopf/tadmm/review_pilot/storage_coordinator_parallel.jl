# Parallel storage-coordinated tADMM with persistent, pinned worker processes,
# and the centralized solve of the same model on the same node.
#
#   julia storage_coordinator_parallel.jl system T output-directory
#
# SC_NWORKERS (default T): period workers. With fewer workers than periods each
# worker owns several periods and solves them one after another.
# SC_RHO: one penalty specification, as in storage_coordinator.jl (default the
# decreasing schedule 2:0.0625:20 in units of rho_unit, SC_RHO_UNITS=1).
#
# Clocks. algorithm_s: Update I, II and III of every iteration, including the
# messages to and from the workers. It excludes the per-iteration quality check
# (dispatch collection, physical metrics) and the repair solve, which are timed
# separately. A full warm-up pass with the same models precedes the timed pass.
using Distributed, Dates
include("storage_coordinator_actor.jl")
BLAS.set_num_threads(1)

const SC_WIDS=Int[]

function allowed_cpus()
    Sys.islinux() || return Int[]
    cpus=Int[]
    for part in split(sc_affinity(),',')
        b=parse.(Int,split(part,'-')); append!(cpus,b[1]:b[end])
    end
    return cpus
end

function start_pool!(nworkers)
    isempty(SC_WIDS) || return
    cpus=allowed_cpus(); pin=Sys.islinux() && get(ENV,"SC_PIN","1")=="1"
    pin && length(cpus)<nworkers+1 && error("Need $(nworkers+1) CPUs, have $(length(cpus))")
    common=[String(k)=>String(v) for (k,v) in ENV if startswith(k,"TADMM_") || startswith(k,"SC_") || k in ("JULIA_DEPOT_PATH","LD_LIBRARY_PATH","LD_PRELOAD")]
    append!(common,["OPENBLAS_NUM_THREADS"=>"1","OMP_NUM_THREADS"=>"1","MKL_NUM_THREADS"=>"1","JULIA_NUM_THREADS"=>"1"])
    flags=Cmd(["--startup-file=no","--threads=1","--project=$(dirname(Base.active_project()))"])
    tasks=map(1:nworkers) do w
        @async if pin
            env=vcat(common,["TADMM_PIN_CPU"=>string(cpus[w+1]),"TADMM_JULIA_EXE"=>Base.julia_cmd().exec[1]])
            only(addprocs(1;exename=joinpath(@__DIR__,"worker_julia.sh"),exeflags=flags,env=env,dir=@__DIR__))
        else
            only(addprocs(1;exeflags=flags,env=common,dir=@__DIR__))
        end
    end
    append!(SC_WIDS,fetch.(tasks))
    fetch.([remotecall(Base.include,w,Main,joinpath(@__DIR__,"storage_coordinator_actor.jl")) for w in SC_WIDS])
    pin && run(pipeline(`taskset -apc $(cpus[1]) $(getpid())`;stdout=devnull))
end

function run_parallel(d,datafile,system;rho,rho_end=rho,ramp=0,spec,maxiter,reference,no_tol,phase,out,extra=5)
    T=d[:T]; bs=d[:Bset]; nb=length(bs); dt=d[:delta_t_h]; W=length(SC_WIDS)
    groups=[collect(w:W:T) for w in 1:W]; chunks=[collect(w:W:nb) for w in 1:W]
    build0=time_ns()
    layout=fetch.([remotecall(sc_prepare!,SC_WIDS[w],datafile,system,groups[w],chunks[w],rho) for w in 1:W])
    build_s=(time_ns()-build0)/1e9
    if phase=="warmup"
        for info in layout; println("SC_WORKER ",info); end
        length(unique((i.host,i.affinity) for i in layout))==W || get(ENV,"SC_PIN","1")!="1" || !Sys.islinux() || error("Workers share CPUs")
    end
    z=zeros(nb,T); u=zeros(nb,T); x=zeros(nb,T); us=Vector{Vector{Float64}}(undef,T)
    algorithm=0.0; no=0; no_s=NaN; repair_no=0; repair_no_s=NaN; k=0; rho_start=rho; last=0
    io=open(joinpath(out,"trace_$(phase).csv"),"w")
    println(io,"iteration,objective,gap,original_primal,r_max,s_max,rho,algorithm_s,iteration_s,slowest_stage_s,median_stage_s,slowest_battery_step_s,max_local_ipopt_iterations,repair_gap,repair_primal,repair_wall_s,check_s")
    loop0=time_ns()
    for iteration in 1:maxiter
        k=iteration; tick=time_ns()
        if ramp>0
            next=rho_start*(rho_end/rho_start)^(min(k-1,ramp)/ramp)
            u.*=rho/next; rho=next
        end
        jobs=[remotecall(sc_xstep!,SC_WIDS[w],Dict(t=>z[:,t] for t in groups[w]),Dict(t=>u[:,t] for t in groups[w]),rho,k) for w in 1:W]
        stage_s=Float64[]; iters=0
        for answer in fetch.(jobs), row in answer
            x[:,row.t]=row.x; push!(stage_s,row.elapsed_s); iters=max(iters,row.iterations)
        end
        target=x.+u; zold=copy(z)
        zjobs=[(w,remotecall(sc_zstep!,SC_WIDS[w],target[chunks[w],:],k)) for w in 1:W if !isempty(chunks[w])]
        battery_s=0.0
        for (w,job) in zjobs
            answer=fetch(job); z[chunks[w],:]=answer.z; battery_s=max(battery_s,answer.elapsed_s)
        end
        u.+=x.-z
        iteration_s=(time_ns()-tick)/1e9; algorithm+=iteration_s
        # Quality check, outside the algorithm clock.
        check0=time_ns()
        rmax=maximum(abs,x.-z); smax=rho*maximum(abs,z.-zold)
        for answer in fetch.([remotecall(sc_dispatch,w) for w in SC_WIDS]),(t,ut) in answer; us[t]=ut; end
        bb=zeros(nb,T); zb=zeros(nb,T)
        for n in 1:nb
            e=d[:B0_pu][bs[n]]; f=e
            for t in 1:T; e-=dt*x[n,t]; bb[n,t]=e; f-=dt*z[n,t]; zb[n,t]=f; end
        end
        metrics=physical_metrics(d,us,bb)
        gap=abs(metrics.objective-reference)/max(1.0,abs(reference))
        check_s=(time_ns()-check0)/1e9
        # Repair: one more parallel solve with the batteries fixed at the coordinator's schedule.
        repair0=time_ns(); rus=Vector{Vector{Float64}}(undef,T); solved=true
        for answer in fetch.([remotecall(sc_repair!,SC_WIDS[w],Dict(t=>z[:,t] for t in groups[w])) for w in 1:W]), row in answer
            solved&=row.ok; row.ok && (rus[row.t]=row.u)
        end
        repair_wall=(time_ns()-repair0)/1e9
        repair_gap=NaN; repair_primal=NaN
        if solved
            pm=physical_metrics(d,rus,zb)
            repair_gap=abs(pm.objective-reference)/max(1.0,abs(reference)); repair_primal=pm.original_primal
        end
        if repair_no==0 && solved && repair_gap<=1e-3 && repair_primal<=no_tol
            repair_no=k; repair_no_s=algorithm+repair_wall
            @printf("SC_REPAIR_NO phase=%s k=%d gap=%.3e original=%.3e algorithm_s=%.3f repair_wall_s=%.3f total_s=%.3f\n",phase,k,repair_gap,repair_primal,algorithm,repair_wall,repair_no_s)
        end
        if no==0 && gap<=1e-3 && metrics.original_primal<=no_tol
            no=k; no_s=algorithm
            @printf("SC_NO phase=%s k=%d gap=%.3e original=%.3e algorithm_s=%.3f\n",phase,k,gap,metrics.original_primal,algorithm)
        end
        sort!(stage_s)
        println(io,join((k,metrics.objective,gap,metrics.original_primal,rmax,smax,rho,algorithm,iteration_s,stage_s[end],stage_s[cld(end,2)],battery_s,iters,repair_gap,repair_primal,repair_wall,check_s),',')); flush(io)
        (k<=3 || k%10==0) && (@printf("SC_ITER phase=%s k=%d gap=%.3e original=%.3e repair_gap=%.3e rho=%g algorithm_s=%.3f\n",phase,k,gap,metrics.original_primal,repair_gap,rho,algorithm); flush(stdout))
        # Stop a few iterations after both near-optimal points are passed.
        no>0 && repair_no>0 && last==0 && (last=k+extra)
        last>0 && k>=last && break
    end
    close(io)
    @printf("SC_DONE phase=%s spec=%s workers=%d iterations=%d no=%d no_s=%.3f repair_no=%d repair_no_s=%.3f algorithm_s=%.3f loop_wall_s=%.3f build_s=%.3f\n",
        phase,spec,W,k,no,no_s,repair_no,repair_no_s,algorithm,(time_ns()-loop0)/1e9,build_s); flush(stdout)
    return (;k,no,no_s,repair_no,repair_no_s,algorithm,build_s)
end

function main()
    system=ARGS[1]; T=parse(Int,ARGS[2]); out=abspath(ARGS[3]); mkpath(out)
    datafile=abspath(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
    d=deserialize(datafile); profile=apply_review_cb!(d,system)
    no_tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    nworkers=min(T,parse(Int,get(ENV,"SC_NWORKERS",string(T))))
    println("SC_HOST host=",gethostname()," time=",now()," affinity=",sc_affinity()," julia=",VERSION)
    # Centralized reference on the same node: full warm-up, then the timed solve from the same start.
    central=make_model(d,collect(1:T))
    warm0=time_ns(); solve_checked!(central); central_warm=(time_ns()-warm0)/1e9
    set_start_value.(central.vars,central.initial)
    tick=time_ns(); central_solver=solve_checked!(central); central_wall=(time_ns()-tick)/1e9
    cu,cb=dispatch([central],d); cm=physical_metrics(d,cu,cb); reference=cm.objective
    @printf("CENTRAL objective=%.12f primal=%.3e iterations=%d solver_s=%.3f wall_s=%.3f warmup_wall_s=%.3f\n",reference,cm.original_primal,
        MOI.get(central.m,MOI.BarrierIterations()),central_solver,central_wall,central_warm); flush(stdout)
    central=nothing; cu=nothing; GC.gc()
    prated=sort([d[:P_B_R_pu][j] for j in d[:Bset]]); pmed=prated[cld(length(prated),2)]
    rho_unit=(maximum(d[:LoadShapeCost])-minimum(d[:LoadShapeCost]))*d[:kVA_B]*d[:delta_t_h]/pmed
    spec=get(ENV,"SC_RHO","2:0.0625:20"); scale=get(ENV,"SC_RHO_UNITS","1")=="1" ? rho_unit : 1.0
    parts=split(spec,':'); rho=scale*parse(Float64,parts[1])
    rho_end=length(parts)==3 ? scale*parse(Float64,parts[2]) : rho
    ramp=length(parts)==3 ? parse(Int,parts[3]) : 0
    @printf("CASE system=%s T=%d C_B=%.6g profile=%s batteries=%d workers=%d rho_unit=%.6g spec=%s rho_start=%g rho_end=%g ramp=%d\n",
        system,T,d[:C_B],profile,length(d[:Bset]),nworkers,rho_unit,spec,rho,rho_end,ramp); flush(stdout)
    pool0=time_ns(); start_pool!(nworkers); @printf("SC_POOL workers=%d start_s=%.2f coordinator_affinity=%s\n",nworkers,(time_ns()-pool0)/1e9,sc_affinity())
    maxiter=parse(Int,get(ENV,"SC_MAX_ITER","150"))
    common=(;rho,rho_end,ramp,spec,maxiter,reference,no_tol,out)
    warm=run_parallel(d,datafile,system;common...,phase="warmup")
    GC.gc()
    timed=run_parallel(d,datafile,system;common...,phase="timed")
    (timed.k,timed.no,timed.repair_no)==(warm.k,warm.no,warm.repair_no) || println("SC_WARNING warm-up and timed passes differ: ",(warm.k,warm.no,warm.repair_no)," vs ",(timed.k,timed.no,timed.repair_no))
    @printf("SC_SUMMARY system=%s T=%d workers=%d central_wall_s=%.3f tadmm_no_iteration=%d tadmm_no_s=%.3f repair_no_iteration=%d repair_no_s=%.3f speedup_no=%.2f speedup_repair_no=%.2f\n",
        system,T,nworkers,central_wall,timed.no,timed.no_s,timed.repair_no,timed.repair_no_s,central_wall/timed.no_s,central_wall/timed.repair_no_s)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    try
        main()
    finally
        isempty(SC_WIDS) || rmprocs(copy(SC_WIDS))
    end
end
