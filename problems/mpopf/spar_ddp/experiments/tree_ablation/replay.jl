# Same captured KKT, three implementations, interleaved repetitions.
using DDP4OPF, LinearAlgebra, SparseArrays, Serialization, Random, Printf, SHA
include(joinpath(@__DIR__, "../../examples/power_system/tree_kkt.jl"))
ENV["FILTERDDP_HALF_BLOCK"]="1";ENV["FILTERDDP_LEAN_EXACT"]="1"
BLAS.set_num_threads(1)
length(ARGS)>=2 || error("usage: replay.jl capture.jls output.csv [repeats=5]")
input,out=ARGS[1:2];repeats=length(ARGS)>2 ? parse(Int,ARGS[3]) : 5
repeats>=2 || error("at least two repetitions required")
d=deserialize(input); K,lay=d.K,d.lay
stat=tree_kkt_static(lay,K);nB=length(lay.E)÷2
rng=MersenneTwister(17);b=randn(rng,lay.n);RE=randn(rng,length(lay.E),nB)
Vc=zeros(nB,nB)
function evaluate(arm)
    if arm=="full_umfpack"
        F=lu(K);x=F\b;XE=zeros(length(lay.E),nB)
        for first in 1:16:nB
            cols=first:min(first+15,nB);R=zeros(lay.n,length(cols));R[lay.E,:]=RE[:,cols]
            XE[:,cols]=(F\R)[lay.E,:]
        end
        return x,XE
    end
    s=withenv("FILTERDDP_NETWORK_BACKEND"=>arm) do
        tree_kkt_solver(lay,stat,K;Vc=Vc)
    end
    x=copy(b);ldiv!(s,x)
    return x,tree_kkt_battery_rows(s,RE)
end
arms=["tree","umfpack","full_umfpack"]
# Warm all routes and check the same physical correction and selected response.
ref=evaluate("full_umfpack")
relative(a,b)=norm(a-b)/max(norm(b),eps())
for arm in arms
    x,Z=evaluate(arm)
    max(relative(x,ref[1]),relative(Z,ref[2]))<=1e-7 || error("$arm fails numerical agreement")
end
mkpath(dirname(abspath(out)))
ispath(out) && error("refusing to overwrite $out")
open(out,"w") do io
    println(io,"arm,repeat,elapsed_s,allocated_bytes,gc_s,n,nE,nB,nnz_K,relative_x,relative_selected_response,capture_sha256")
    digest=bytes2hex(open(sha256,input))
    for rep in 1:repeats
        for arm in (isodd(rep) ? arms : reverse(arms))
            GC.gc();r=@timed evaluate(arm);x,Z=r.value
            ex,ez=relative(x,ref[1]),relative(Z,ref[2])
            max(ex,ez)<=1e-7 || error("$arm failed during repetition $rep")
            @printf(io,"%s,%d,%.9f,%d,%.9f,%d,%d,%d,%d,%.4e,%.4e,%s\n",
                arm,rep,r.time,r.bytes,r.gctime,lay.n,length(lay.E),nB,nnz(K),ex,ez,digest)
            flush(io)
        end
    end
end
println("TREE_ABLATION_REPLAY_PASS output=",out)
