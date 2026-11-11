using Distributed, Sockets, SHA
if get(ENV,"TADMM_CLUSTER_BACKEND","local")=="slurm"
    using SlurmClusterManager
end
isdefined(Main,:make_model) || include("model.jl")
const ACTOR_DATA=Ref{Any}(nothing)
const ACTOR_MODELS=Dict{Int,Any}()
const POOL_WIDS=Int[]
const POOL_GROUPS=Vector{Int}[]
const POOL_LAYOUT=Any[]
const NODE_SAMPLER_STARTED=Ref(false)

function start_node_sampler!(folder)
    NODE_SAMPLER_STARTED[] && return
    head=get(ENV,"TADMM_COORDINATOR_HOST",gethostname())
    (myid()==1 || (gethostname()!=head && get(ENV,"SLURM_LOCALID","")=="0")) || return
    pin=get(ENV,"TADMM_MONITOR_CPU",last(split(ENV["TADMM_CPU_IDS"],',')))
    log=joinpath(dirname(folder),"sampler.$(gethostname()).log")
    open(log,"w") do io
        cmd=`taskset -c $pin python3 $(joinpath(@__DIR__,"sample_registry.py")) $folder`
        run(pipeline(detach(cmd);stdout=io,stderr=io);wait=false)
    end
    NODE_SAMPLER_STARTED[]=true
end

function register_process!()
    folder=get(ENV,"TADMM_PID_REGISTRY","")
    isempty(folder) && return
    mkpath(folder)
    write(joinpath(folder,"$(gethostname())_$(getpid()).pid"),string(getpid()))
    start_node_sampler!(folder)
end

function affinity_string()
    Sys.islinux() || return "not-linux"
    return strip(match(r"Cpus_allowed_list:\s*([^\n]+)",read("/proc/self/status",String)).captures[1])
end

function actor_prepare!(d,periods,rho,local_policy,local_primal_report,warm_start,objective_update)
    register_process!()
    BLAS.set_num_threads(1)
    ENV["TADMM_LOCAL_SOLVE_POLICY"]=local_policy
    ENV["TADMM_LOCAL_PRIMAL_REPORT"]=local_primal_report
    ENV["TADMM_WARM_START"]=warm_start
    ENV["TADMM_OBJECTIVE_UPDATE"]=objective_update
    empty!(ACTOR_MODELS); GC.gc()
    ACTOR_DATA[]=d
    for t in periods
        ACTOR_MODELS[t]=make_model(d,[t];local_window=true,rho=rho)
    end
    starts=Dict(string(t)=>bytes2hex(sha256(reinterpret(UInt8,Float64.(ACTOR_MODELS[t].initial)))) for t in periods)
    return Dict("julia_worker"=>myid(),"pid"=>getpid(),"hostname"=>gethostname(),
        "periods"=>periods,"affinity"=>affinity_string(),"blas_threads"=>BLAS.get_num_threads(),
        "julia_threads"=>Threads.nthreads(),"openmp_threads"=>get(ENV,"OMP_NUM_THREADS","unknown"),
        "local_solve_policy"=>local_policy,"local_primal_report"=>(local_primal_report=="1"),"initial_starts_sha256"=>starts,
        "expected_pin"=>get(ENV,"TADMM_PIN_CPU",""),"warm_start"=>get(ENV,"TADMM_WARM_START","primal"),
        "objective_update"=>get(ENV,"TADMM_OBJECTIVE_UPDATE","legacy"))
end

function actor_step!(z,dual_by_period,rho,k)
    return [stage_step!(ACTOR_MODELS[t],ACTOR_DATA[],z,dual_by_period[t],rho,k)
        for t in sort!(collect(keys(ACTOR_MODELS)))]
end

function actor_dispatch()
    return [(t,dispatch_period(ACTOR_MODELS[t],ACTOR_DATA[],t)[1])
        for t in sort!(collect(keys(ACTOR_MODELS)))]
end

function actor_checkpoint()
    return Dict(t=>value.(s.vars) for (t,s) in ACTOR_MODELS)
end

function prepare_pool!(d,rho)
    wanted=min(d[:T],parse(Int,get(ENV,"TADMM_NWORKERS","0")))
    wanted>0 || error("No requested period workers")
    Sys.islinux() || error("Pinned process benchmarks currently require Linux")
    cpus=parse.(Int,split(get(ENV,"TADMM_CPU_IDS",""),','))
    backend=get(ENV,"TADMM_CLUSTER_BACKEND","local")
    backend in ("local","slurm") || error("Unknown cluster backend")
    length(cpus)>=(backend=="local" ? wanted+1 : 2) || error("Insufficient reserved CPUs")
    length(unique(cpus))==length(cpus) || error("CPU IDs must be unique")
    if isempty(POOL_WIDS)
        common=[String(k)=>String(v) for (k,v) in ENV if startswith(k,"TADMM_") || k in ("JULIA_DEPOT_PATH","LD_LIBRARY_PATH")]
        append!(common,["OPENBLAS_NUM_THREADS"=>"1","OMP_NUM_THREADS"=>"1","MKL_NUM_THREADS"=>"1",
            "JULIA_NUM_THREADS"=>"1","OMP_WAIT_POLICY"=>"PASSIVE","GOMP_SPINCOUNT"=>"0"])
        if backend=="slurm"
            parse(Int,ENV["SLURM_NTASKS"])==wanted || error("Slurm task count must equal period worker count")
            ids=addprocs(SlurmManager(;launch_timeout=600.0);
                exename=joinpath(@__DIR__,"worker_slurm.sh"),
                exeflags=Cmd(["--startup-file=no","--threads=1","--project=$(dirname(Base.active_project()))"]),
                env=common,dir=@__DIR__,topology=:master_worker)
            append!(POOL_WIDS,ids)
            fetch.([remotecall(Base.include,wid,Main,joinpath(@__DIR__,"parallel_workers.jl")) for wid in ids])
            for w in 1:wanted;push!(POOL_GROUPS,collect(w:wanted:d[:T]));end
        else
        for w in 1:wanted
            env=vcat(common,["TADMM_PIN_CPU"=>string(cpus[w+1]),"TADMM_JULIA_EXE"=>Base.julia_cmd().exec[1]])
            ids=addprocs(1;exename=joinpath(@__DIR__,"worker_julia.sh"),
                exeflags=Cmd(["--startup-file=no","--threads=1","--project=$(dirname(Base.active_project()))"]),env=env,dir=@__DIR__)
            append!(POOL_WIDS,ids)
            remotecall_fetch(Base.include,only(ids),Main,joinpath(@__DIR__,"parallel_workers.jl"))
            push!(POOL_GROUPS,collect(w:wanted:d[:T]))
        end
        end
        # Pin every coordinator thread AFTER the central baseline used the full
        # explicitly selected CPU budget. Workers were independently pinned.
        run(pipeline(`taskset -apc $(cpus[1]) $(getpid())`;stdout=devnull))
    end
    length(POOL_WIDS)==wanted || error("Worker count changed between warm-up and measured solve")
    local_policy=get(ENV,"TADMM_LOCAL_SOLVE_POLICY","exact")
    local_primal_report=get(ENV,"TADMM_LOCAL_PRIMAL_REPORT","1")
    warm_start=get(ENV,"TADMM_WARM_START","primal")
    objective_update=get(ENV,"TADMM_OBJECTIVE_UPDATE","legacy")
    jobs=[remotecall(actor_prepare!,wid,d,group,rho,local_policy,local_primal_report,warm_start,objective_update) for (wid,group) in zip(POOL_WIDS,POOL_GROUPS)]
    empty!(POOL_LAYOUT); append!(POOL_LAYOUT,fetch.(jobs))
    all(info->info["julia_threads"]==1 && info["blas_threads"]==1 && info["openmp_threads"]=="1",POOL_LAYOUT) || error("Worker thread policy mismatch")
    all(info->info["warm_start"]==warm_start && info["objective_update"]==objective_update,POOL_LAYOUT) || error("Stale worker reuse policy")
    if backend=="local"
        all(i->POOL_LAYOUT[i]["affinity"]==string(cpus[i+1]),eachindex(POOL_LAYOUT)) || error("Worker CPU pin mismatch")
    else
        all(info->info["affinity"]==info["expected_pin"],POOL_LAYOUT) || error("Slurm worker CPU pin mismatch")
        length(unique((info["hostname"],info["affinity"]) for info in POOL_LAYOUT))==wanted || error("Overlapping worker CPUs")
        all(info->info["hostname"]!=gethostname() || info["affinity"]!=string(cpus[1]),POOL_LAYOUT) || error("Coordinator shares worker CPU")
    end
    @printf("PERIOD_WORKERS workers=%d coordinator_pid=%d coordinator_affinity=%s\n",wanted,getpid(),affinity_string())
    for info in POOL_LAYOUT; println("PERIOD_WORKER ",info); end
end

function pool_step!(z,duals,rho,k,T)
    jobs=[remotecall(actor_step!,wid,z,Dict(t=>duals[t] for t in group),rho,k)
        for (wid,group) in zip(POOL_WIDS,POOL_GROUPS)]
    results=Vector{Any}(undef,T); seen=falses(T)
    for answer in fetch.(jobs), row in answer
        seen[row.period] && error("Duplicate period result")
        results[row.period]=row; seen[row.period]=true
    end
    all(seen) || error("Missing period result")
    return results
end

function pool_dispatch(d,z)
    jobs=[remotecall(actor_dispatch,wid) for wid in POOL_WIDS]
    us=Vector{Vector{Float64}}(undef,d[:T])
    for answer in fetch.(jobs),(t,u) in answer; us[t]=u; end
    return us,copy(z)
end

function pool_checkpoint()
    jobs=[remotecall(actor_checkpoint,wid) for wid in POOL_WIDS]
    return merge(fetch.(jobs)...)
end

function shutdown_pool!()
    isempty(POOL_WIDS) || rmprocs(copy(POOL_WIDS))
    empty!(POOL_WIDS); empty!(POOL_GROUPS); empty!(POOL_LAYOUT)
end
