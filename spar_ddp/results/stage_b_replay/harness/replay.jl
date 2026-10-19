# Replay captured real Stage A outputs; no model construction or optimizer loop.
using Distributed, Serialization, LinearAlgebra, SharedArrays, Statistics, Printf
const ROOT = ENV["REPLAY_ROOT"]
const REPLAY_REPO = joinpath(ROOT, "repo")
const PS = joinpath(REPLAY_REPO, "ddp", "examples", "power_system")
const N = 96
const CAP = joinpath(ROOT, "capture")
const MULTI = get(ENV, "REPLAY_MULTINODE", "0") == "1"
const MAXW = MULTI ? 96 : 63
if MULTI
    using SlurmClusterManager
    ENV["SLURM_NTASKS"] = string(MAXW)
    addprocs(SlurmManager(); exeflags=`--startup-file=no --project=$(Base.active_project()) --threads=1`)
else
    addprocs(MAXW; exeflags=`--startup-file=no --project=$(Base.active_project()) --threads=1`)
end
@everywhere begin
    using Serialization, LinearAlgebra, SharedArrays
    BLAS.set_num_threads(1)
    include($(joinpath(PS, "ieee123c_filterddp.jl")))
    include($(joinpath(PS, "procs_pass.jl")))
    const REPLAY_INPUT = Dict{Int,Any}()
    const REPLAY_E = Dict{Int,Any}()
    const REPLAY_SLOTS = Any[]
    const REPLAY_SEGMENTS = Any[]
    const REPLAY_WORK = Ref{Any}(nothing)
    const REPLAY_ZERO = Ref{Any}(nothing)
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
        z=_pp_combine!(x,y,REPLAY_WORK[])
        shared || (REPLAY_E[t]=z)
        (elapsed=(time_ns()-ts)/1e9,compute=(time_ns()-tc)/1e9,comm=comm,bytes=bytes,local_read=local_read,gc=(Base.gc_num().total_time-gs)/1e9,pid=myid())
    end
    function replay_result(t,shared)
        v=replay_value(t,true,shared)
        (V=copy(v.V),v=copy(v.v))
    end
end
hosts=Dict(p=>remotecall_fetch(gethostname,p) for p in workers())
println("REPLAY_PLACEMENT ", hosts)
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
counts=MULTI ? [32,48,63,96] : [7,16,32,48,63]
modes=MULTI ? [false] : [true,false]
for shared in modes, W in counts
    # Same-node copying controls only at the larger worker counts.
    !MULTI && !shared && W < 32 && continue
    pids=pool[1:W]; owner=[pids[(t-1)%W+1] for t in 1:N]
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
    totals=Float64[]
    for rep in 0:3
        @sync for p in pids; @async remotecall_fetch(replay_reset,p,shared); end
        valid=falses(N+1); valid[end]=true
        start=time_ns()
        for (round,h) in enumerate(2 .^ (0:6))
            pairs=Tuple{Int,Int}[]
            for i in 1:2h:N+1
                mid=i+h; mid>N+1 && continue
                for t in i:mid-1; push!(pairs,(t,mid)); end
            end
            rr=Vector{Any}(undef,length(pairs)); rt=time_ns()
            @sync for (q,(t,mid)) in enumerate(pairs)
                @async rr[q]=remotecall_fetch(replay_combine,owner[t],t,mid,valid[mid],mid==N+1,mid==N+1 ? owner[t] : owner[mid],shared)
            end
            elapsed=(time_ns()-rt)/1e9
            loads=Dict(p=>sum(r.elapsed for r in rr if r.pid==p;init=0.0) for p in pids)
            cross=sum(mid<=N && hosts[owner[t]]!=hosts[owner[mid]] for (t,mid) in pairs)
            @printf("REPLAY_ROUND nodes=%d shared=%d workers=%d rep=%d round=%d tasks=%d wall_s=%.6f compute_sum_s=%.6f fetch_sum_s=%.6f gc_sum_s=%.6f busiest_worker_s=%.6f scheduling_wait_est_s=%.6f bytes=%d cross_node_tasks=%d\n",length(hostlist),shared,W,rep,round,length(pairs),elapsed,sum(r.compute for r in rr),sum(r.comm for r in rr),sum(r.gc for r in rr),maximum(values(loads)),max(0,elapsed-maximum(values(loads))),sum(r.bytes for r in rr),cross)
            for (t,mid) in pairs; valid[t]=valid[mid]; end
        end
        elapsed=(time_ns()-start)/1e9
        # Per-stage numerical check, after timing; no large gather in timing.
        err=0.0
        for t in 1:N
            got=remotecall_fetch(replay_result,owner[t],t,shared); ref=reference[t]
            err=max(err,norm(got.V-ref.V,Inf)/max(1,norm(ref.V,Inf)),norm(got.v-ref.v,Inf)/max(1,norm(ref.v,Inf)))
        end
        @printf("REPLAY_PASS nodes=%d shared=%d workers=%d rep=%d wall_s=%.6f max_relative_error=%.12e\n",length(hostlist),shared,W,rep,elapsed,err)
        err < 1e-7 || error("Replay numerical mismatch: $err")
        rep>0 && push!(totals,elapsed)
    end
    @printf("REPLAY_RESULT nodes=%d shared=%d workers=%d median_s=%.6f min_s=%.6f max_s=%.6f\n",length(hostlist),shared,W,median(totals),minimum(totals),maximum(totals))
end
println("REPLAY_COMPLETE")
rmprocs(workers())
