using Test,Random,LinearAlgebra
include(joinpath(@__DIR__,"unit_bootstrap.jl"))
include(joinpath(@__DIR__,"cache_core.jl"))
Random.seed!(29527);n=7
function tiny()
 PintElem(Matrix(8.0I,n,n),randn(n,n),Matrix(9.0I,n,n),randn(n),randn(n))
end
for val in (false,true)
 x=tiny();y=val ? PintVal(Matrix(10.0I,n,n),randn(n)) : tiny()
 fresh=_pp_combine!(deepcopy(x),y,PairWork(n));out,c=populate_merge!(deepcopy(x),y,PairWork(n))
 @test c!==nothing
 for f in (hasproperty(fresh,:V) ? (:V,:v) : (:Maa,:Mab,:Mbb,:ma,:mb))
  @test exactbits(getproperty(out,f),getproperty(fresh,f))
 end
 z=deepcopy(x);z.ma .+= 0.1;z.mb .-=0.2
 yy=deepcopy(y);(val ? yy.v : yy.ma).+=0.3
 a=hit_merge!(deepcopy(z),yy,PairWork(n),c);b=_pp_combine!(z,yy,PairWork(n))
 for f in (hasproperty(b,:V) ? (:V,:v) : (:ma,:mb));@test exactbits(getproperty(a,f),getproperty(b,f));end
end
x=tiny();y=tiny();x.Mbb[1,1]=-100.0
out,c=populate_merge!(deepcopy(x),y,PairWork(n));ref=_pp_combine!(deepcopy(x),y,PairWork(n))
@test c===nothing
for f in (:Maa,:Mab,:Mbb,:ma,:mb);@test exactbits(getproperty(out,f),getproperty(ref,f));end
p=(k=4,j=0,mu=1.0,context="pinned",proof=:pinned_unshifted_tp_pre,stage_keys=["a"],terminal_key="t")
q=merge(p,(j=1,mu=0.2));@test eligible(p,q)
for delta in ((k=5,),(context="new_model",),(proof=:unknown,),(stage_keys=["changed_factor"],),(terminal_key="other",))
 @test !eligible(p,merge(q,delta))
end
@test !eligible(p,q;complete=false)
# Full68-node closure test: hit intermediates need no matrix materialization.
leaves=[tiny() for t in 1:24];term=PintVal(zeros(n,n),zeros(n))
xs=Any[deepcopy(x) for x in leaves];push!(xs,deepcopy(term));cache=Dict{Any,Any}()
for (depth,ops) in enumerate(rounds(24)),(t,mid) in ops
 xs[t],entry=populate_merge!(xs[t],xs[mid],PairWork(n));cache[(depth,t,mid)]=entry
end
@test length(cache)==68&&count(c->c.V!==nothing,values(cache))==24
newleaves=deepcopy(leaves)
for e in newleaves;e.ma .+=0.11;e.mb .-=0.07;end
actual=Any[deepcopy(x) for x in newleaves];push!(actual,deepcopy(term))
ref=deepcopy(actual)
for (depth,ops) in enumerate(rounds(24)),(t,mid) in ops
 actual[t]=hit_merge!(actual[t],actual[mid],PairWork(n),cache[(depth,t,mid)])
 ref[t]=_pp_combine!(ref[t],ref[mid],PairWork(n))
end
for t in 1:25
 @test exactbits(actual[t].V,ref[t].V)&&exactbits(actual[t].v,ref[t].v)
end
@test sum(length,rounds(24))==68
storage=zeros(3n*n+2n);slot=shared_slot_view(storage,n)
slot.ma .=1;slot.mb .=2
@test storage[3n*n+1:3n*n+n]==ones(n)
@test storage[3n*n+n+1:end]==fill(2.0,n)
@test all(iszero,storage[1:3n*n])
println("CACHE_UNIT_PASS exact_kernel elem_value current_vectors fallback invalidation global68")
