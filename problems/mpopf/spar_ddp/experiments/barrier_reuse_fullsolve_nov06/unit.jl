using Test,Random,LinearAlgebra
include(joinpath(@__DIR__,"unit_bootstrap.jl"))
include(joinpath(@__DIR__,"right_kernel.jl"))
include(joinpath(@__DIR__,"schedule.jl"))
include(joinpath(@__DIR__,"cache_kernel.jl"))
Random.seed!(31196);n=7
same(a,b)=reinterpret(UInt64,vec(a))==reinterpret(UInt64,vec(b))
tiny()=PintElem(Matrix(8.0I,n,n),randn(n,n),Matrix(9.0I,n,n),randn(n),randn(n))
for N in (7,24,96,144,288)
 leaves=[tiny() for _ in 1:N];term=PintVal(zeros(n,n),zeros(n))
 xs=Any[deepcopy(leaves)...,deepcopy(term)];ref=deepcopy(xs);cache=Dict{Any,Any}()
 for (d,ops) in enumerate(replay_rounds(N,2)),(t,m) in ops
  xs[t],c=br_populate!(xs[t],xs[m],PairWork(n));cache[(d,t,m)]=c
  ref[t]=right_combine!(ref[t],ref[m],PairWork(n))
  for f in (xs[t] isa PintVal ? (:V,:v) : (:Maa,:Mab,:Mbb,:ma,:mb))
   @test same(getproperty(xs[t],f),getproperty(ref[t],f))
  end
 end
 @test length(cache)==sum(length,replay_rounds(N,2))
 @test count(c->c.V!==nothing,values(cache))==N
 for x in leaves;x.ma .+=0.17;x.mb .-=0.09;end
 xs=Any[deepcopy(leaves)...,deepcopy(term)];ref=deepcopy(xs)
 for (d,ops) in enumerate(replay_rounds(N,2)),(t,m) in ops
  xs[t]=br_hit!(xs[t],xs[m],PairWork(n),cache[(d,t,m)])
  ref[t]=right_combine!(ref[t],ref[m],PairWork(n))
 end
 for t in 1:N+1
  @test same(xs[t].V,ref[t].V);@test same(xs[t].v,ref[t].v)
 end
 # Re-population writes existing storage and still matches right kernel.
 x=tiny();y=tiny();_,c=br_populate!(deepcopy(x),y,PairWork(n))
 _,cc=br_populate!(deepcopy(x),y,PairWork(n),c);@test cc===c
 @test br_bytes(c)==3n*n*8
end
x=tiny();y=tiny();x.Mbb[1,1]=-100
z,c=br_populate!(deepcopy(x),y,PairWork(n));ref=right_combine!(deepcopy(x),y,PairWork(n))
@test c===nothing
for f in (:Maa,:Mab,:Mbb,:ma,:mb);@test same(getproperty(z,f),getproperty(ref,f));end
for reason in (:accepted,:model,:terminal,:regularization,:path,:unknown,:state,:factor)
 g=BRGate();br_event!(g,:solve_start);@test !br_request!(g,4,0,1.0,(:pinned,),true)
 br_event!(g,reason);@test !br_request!(g,4,1,0.2,(:pinned,),true)
end
g=BRGate();br_event!(g,:solve_start);@test !br_request!(g,4,0,1.0,(:pinned,),true)
br_event!(g,:barrier);@test br_request!(g,4,1,0.2,(:pinned,),true)
@test !br_request!(g,4,2,0.04,(:pinned,),true)
for (k,j,mu,ctx,supported) in ((5,1,.2,(:pinned,),true),(4,2,.2,(:pinned,),true),(4,1,1.,(:pinned,),true),(4,1,.2,(:changed,),true),(4,1,.2,(:pinned,),false))
 local g=BRGate();br_request!(g,4,0,1.0,(:pinned,),true);br_event!(g,:barrier)
 @test !br_request!(g,k,j,mu,ctx,supported)
end
println("BR_UNIT_PASS grouped_right_population_hits_full_T7_T24_T96_T144_T288_buffer_reuse_fallback_epoch_negatives")
