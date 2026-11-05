# Replay captured real Stage A outputs; no model construction or optimizer loop.
using Distributed, Serialization, LinearAlgebra, SharedArrays, Statistics, Printf
const ROOT = ENV["REPLAY_ROOT"]
const REPLAY_REPO = ENV["REPLAY_REPO"]
const PS = joinpath(REPLAY_REPO, "ddp", "examples", "power_system")
const N = 96
const CAP = ENV["REPLAY_CAPTURE"]
const MULTI = get(ENV, "REPLAY_MULTINODE", "0") == "1"
const MAXW = parse(Int, ENV["STUDY_WORKERS"])
if MULTI
    using SlurmClusterManager
    ENV["SLURM_NTASKS"] = string(MAXW)
    addprocs(SlurmManager(); lazy=false, exeflags=`--startup-file=no --project=$(Base.active_project()) --threads=1`)
else
    addprocs(MAXW; exeflags=`--startup-file=no --project=$(Base.active_project()) --threads=1`)
end
@everywhere begin
    using Serialization, LinearAlgebra, SharedArrays
    BLAS.set_num_threads(parse(Int, ENV["STUDY_BLAS_THREADS"]))
    include($(joinpath(PS, "ieee123c_filterddp.jl")))
    include($(joinpath(PS, "procs_pass.jl")))
    include(joinpath(ENV["REPLAY_ROOT"], "right_kernel.jl"))
    const STUDY_KERNEL = Ref("current")
    const STUDY_DIAGNOSTIC = Ref(false)
    function study_configure(kernel, diagnostic)
        parts=split(kernel,'_')
        threads=length(parts)==2 ? parse(Int,parts[2]) : parse(Int,ENV["STUDY_BLAS_THREADS"])
        STUDY_KERNEL[] = parts[1]
        STUDY_DIAGNOSTIC[] = diagnostic
        BLAS.set_num_threads(threads)
        threads
    end
    function study_pin(cpus)
        run(pipeline(`taskset -apc $cpus $(getpid())`; stdout=devnull))
        (pid=myid(),os_pid=getpid(),cpus=cpus,blas_threads=BLAS.get_num_threads(),
         backend=string(BLAS.get_config()),
         thread_status=read(`ps -L -p $(getpid()) -o pid,tid,psr,comm`,String))
    end
    const REPLAY_INPUT = Dict{Int,Any}()
    const REPLAY_E = Dict{Int,Any}()
    const REPLAY_SLOTS = Any[]
    const REPLAY_SEGMENTS = Any[]
    const REPLAY_WORK = Ref{Any}(nothing)
    const REPLAY_ZERO = Ref{Any}(nothing)
    function replay_process_info()
        mask=match(r"(?m)^Cpus_allowed_list:\s*(.+)$",read("/proc/self/status",String))
        model=match(r"(?m)^model name\s*:\s*(.+)$",read("/proc/cpuinfo",String))
        (host=gethostname(),affinity=strip(mask.captures[1]),cpu=strip(model.captures[1]),maxrss=Sys.maxrss())
    end
    function replay_setup(cap, own, slots)
        empty!(REPLAY_INPUT); empty!(REPLAY_E); empty!(REPLAY_SLOTS)
        empty!(REPLAY_SEGMENTS); append!(REPLAY_SEGMENTS,slots)
        for t in own
            REPLAY_INPUT[t] = deserialize(joinpath(cap, "stage_$(t).jls"))
        end
        n = size(deserialize(joinpath(cap, "stage_1.jls")).Maa, 1)
        for a in slots
            p=pointer(sdata(a))
            blk=i->unsafe_wrap(Array,p+(i-1)*n*n*sizeof(Float64),(n,n))
            vec=i->unsafe_wrap(Array,p+(3n*n+(i-1)*n)*sizeof(Float64),(n,))
            push!(REPLAY_SLOTS,PintElem(blk(1),blk(2),blk(3),vec(1),vec(2)))
        end
        REPLAY_WORK[] = PairWork(n)
        REPLAY_ZERO[] = PintVal(zeros(n,n), zeros(n))
        GC.enable(true); GC.gc()
        return (host=gethostname(), own=own, rss=Sys.maxrss())
    end
    function replay_reset(shared)
        GC.enable(true); GC.gc(false)
        for (t,e) in REPLAY_INPUT
            if shared
                x = REPLAY_SLOTS[t]
                for f in (:Maa,:Mab,:Mbb,:ma,:mb); copyto!(getproperty(x,f),getproperty(e,f)); end
            else
                REPLAY_E[t] = deepcopy(e)
            end
        end
        # Match once-per-pass GC: collect between passes, not within a pass.
        GC.enable(false)
        nothing
    end
    function replay_value(t, val, shared)
        e = shared ? REPLAY_SLOTS[t] : REPLAY_E[t]
        val ? (e isa PintVal ? e : PintVal(e.Maa,e.ma)) : e
    end
    function replay_combine(t, mid, val, last, owner, shared)
        ts=time_ns(); gs=Base.gc_num().total_time; comm=0.0; bytes=0
        local_read = last || owner == myid()
        y = if last
            REPLAY_ZERO[]
        elseif shared || owner == myid()
            replay_value(mid,val,shared)
        else
            tc=time_ns(); z=remotecall_fetch(replay_value,owner,mid,val,false)
            comm=(time_ns()-tc)/1e9
            bytes=8*(length(z isa PintVal ? z.V : z.Maa)*(z isa PintVal ? 1 : 3)+length(z isa PintVal ? z.v : z.ma)*(z isa PintVal ? 1 : 2))
            z
        end
        tc=time_ns()
        x=shared ? REPLAY_SLOTS[t] : REPLAY_E[t]
        if STUDY_DIAGNOSTIC[]
            measurement = @timed (STUDY_KERNEL[] == "right" ?
                right_combine!(x,y,REPLAY_WORK[]) : _pp_combine!(x,y,REPLAY_WORK[]))
            z = measurement.value
            bytes_allocated = measurement.bytes
            compile_s = measurement.compile_time + measurement.recompile_time
        else
            z = STUDY_KERNEL[] == "right" ? right_combine!(x,y,REPLAY_WORK[]) :
                _pp_combine!(x,y,REPLAY_WORK[])
            bytes_allocated = -1
            compile_s = -1.0
        end
        shared || (REPLAY_E[t]=z)
        (elapsed=(time_ns()-ts)/1e9,compute=(time_ns()-tc)/1e9,comm=comm,bytes=bytes,local_read=local_read,gc=(Base.gc_num().total_time-gs)/1e9,pid=myid(),allocated=bytes_allocated,compile_s=compile_s)
    end
    function replay_result(t,shared)
        v=replay_value(t,true,shared)
        (V=copy(v.V),v=copy(v.v))
    end
end
BLAS.set_num_threads(1)
cpulines = readlines(joinpath(ROOT,"cpu_groups.csv"))
@assert length(cpulines) == MAXW+1
run(pipeline(`taskset -apc $(cpulines[1]) $(getpid())`; stdout=devnull))
for (j,p) in enumerate(workers())
    println("AFFINITY ",remotecall_fetch(study_pin,p,cpulines[j+1]))
end
hosts=Dict(p=>remotecall_fetch(gethostname,p) for p in workers())
println("REPLAY_PLACEMENT ", hosts)
length(workers()) == MAXW || error("Incorrect worker count")
for p in procs()
    println("REPLAY_PROCESS pid=",p," info=",remotecall_fetch(replay_process_info,p))
end
# Balanced round-robin interleaving across nodes, so small W uses both nodes.
hostlist=sort(unique(collect(values(hosts))))
groups=[sort([p for p in workers() if hosts[p]==h]) for h in hostlist]
pool=Int[]
for k in 1:maximum(length.(groups)), g in groups
    k <= length(g) && push!(pool,g[k])
end
inputs=[deserialize(joinpath(CAP,"stage_$(t).jls")) for t in 1:N]
n=size(inputs[1].Maa,1)
function make_reference(inputs,n)
    w=PairWork(n); v=PintVal(zeros(n,n),zeros(n)); reference=Vector{Any}(undef,N)
    for t in N:-1:1
        v=_pp_combine!(deepcopy(inputs[t]),v,w)
        reference[t]=(V=copy(v.V),v=copy(v.v))
    end
    reference
end
reference=make_reference(inputs,n)
MULTI && error("Balanced scheduling is shared-memory one-node only")
include(joinpath(@__DIR__, "schedule.jl"))
counts=[MAXW]
modes=[true]
for shared in modes, W in counts
    # Same-node copying controls only at the larger worker counts.
    !MULTI && !shared && W < 32 && continue
    pids=pool[1:W]; owner=[pids[(t-1)%W+1] for t in 1:N]
    println("SHARED_BYTES ",8N*(3n*n+2n))
    slots=Any[]
    if shared
        for t in 1:N
            a=SharedArray{Float64}((3n*n+2n,);pids=vcat(myid(),pids))
            push!(slots,a)
        end
    end
    @sync for p in pids
        @async remotecall_fetch(replay_setup,p,CAP,findall(==(p),owner),slots)
    end
    println("SCRATCH_OBJECT_BYTES ",remotecall_fetch(p->Base.summarysize(REPLAY_WORK[]),pids[1],0))
    kernels=split(get(ENV,"STUDY_KERNELS","current,right"),',')
    totals=Dict(k=>Float64[] for k in kernels)
    for rep in vcat(collect(0:5),[-2,-1]), scheduling in (iseven(rep) ? kernels : reverse(kernels))
        threads=length(split(scheduling,'_'))==2 ? parse(Int,split(scheduling,'_')[2]) : parse(Int,ENV["STUDY_BLAS_THREADS"])
        maps=readlines(joinpath(ROOT,"cpu_groups_$threads.csv"))
        @sync for (j,p) in enumerate(pids)
            @async begin
                remotecall_fetch(study_configure,p,scheduling,rep<0)
                placement=remotecall_fetch(study_pin,p,maps[j+1])
                println("CONFIG rep=$rep mode=$scheduling info=$placement")
            end
        end
        @sync for p in pids; @async remotecall_fetch(replay_reset,p,shared); end
        valid=falses(N+1); valid[end]=true
        rounds=replay_rounds(N,2)
        @assert (sum(length,rounds),length(rounds)) == (256,8)
        roundstats=[]
        start=time_ns()
        for (round,pairs) in enumerate(rounds)
            rr=Vector{Any}(undef,length(pairs)); finished=zeros(UInt64,length(pairs)); rt=time_ns()
            assigned=[pids[(q-1)%W+1] for q in eachindex(pairs)]
            # The round's left writes are unique and no left operand is read
            # as a right operand in the same round. Any worker may write it.
            lefts=Set(first.(pairs)); rights=Set(last.(pairs))
            isempty(intersect(lefts,rights)) || error("Unsafe round dependencies")
            @sync for (q,(t,mid)) in enumerate(pairs)
                @async begin
                    rr[q]=remotecall_fetch(replay_combine,assigned[q],t,mid,valid[mid],mid==N+1,mid==N+1 ? owner[t] : owner[mid],shared)
                    finished[q]=time_ns()
                end
            end
            elapsed=(time_ns()-rt)/1e9
            push!(roundstats,(round,elapsed,assigned,rr,finished))
            for (t,mid) in pairs; valid[t]=valid[mid]; end
        end
        elapsed=(time_ns()-start)/1e9
        @assert all(valid)
        for (round,wall,assigned,rr,finished) in roundstats
            worker_time=[sum((rr[q].compute for q in eachindex(assigned) if assigned[q]==p); init=0.0) for p in pids]
            println("ROUND_DETAIL mode=$scheduling rep=$rep round=$round kernel_sum_s=$(sum(worker_time)) busiest_s=$(maximum(worker_time)) allocation_bytes=$(sum(r.allocated for r in rr)) gc_s=$(sum(r.gc for r in rr)) compile_s=$(sum(r.compile_s for r in rr)) completion_spread_s=$((maximum(finished)-minimum(finished))/1e9)")
            println("ROUND mode=$scheduling rep=$rep round=$round wall_s=$wall merges=$(length(rr)) tasks_per_worker=$(join([count(==(p),assigned) for p in pids],','))")
        end
        # Per-stage numerical check, after timing; no large gather in timing.
        grad=hess=0.0
        for t in 1:N
            got=remotecall_fetch(replay_result,owner[t],t,shared); ref=reference[t]
            hess=max(hess,norm(got.V-ref.V,Inf)/max(1,norm(ref.V,Inf)))
            grad=max(grad,norm(got.v-ref.v,Inf)/max(1,norm(ref.v,Inf)))
        end
        err=max(grad,hess)
        println("ERRORS mode=$scheduling rep=$rep grad=$grad hess=$hess")
        @printf("BALANCED_PASS scheduling=%s workers=%d rep=%d wall_s=%.6f max_relative_error=%.12e\n",scheduling,W,rep,elapsed,err)
        err < 1e-11 || error("Replay numerical mismatch: $err")
        rep>0 && push!(totals[scheduling],elapsed)
    end
    for scheduling in kernels
        times=totals[scheduling]
        @printf("BALANCED_RESULT scheduling=%s workers=%d median_s=%.6f min_s=%.6f max_s=%.6f\n",scheduling,W,median(times),minimum(times),maximum(times))
    end
    info=[remotecall_fetch(replay_process_info,p) for p in pids]
    println("REPLAY_MEMORY workers=",W," worker_maxrss_max_bytes=",maximum(x.maxrss for x in info)," worker_maxrss_sum_bytes=",sum(x.maxrss for x in info)," coordinator_maxrss_bytes=",Sys.maxrss())
end
println("BALANCED_COMPLETE")
rmprocs(workers())
