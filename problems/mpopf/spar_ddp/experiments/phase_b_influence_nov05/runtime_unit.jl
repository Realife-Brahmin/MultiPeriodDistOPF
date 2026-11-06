# Same pinned production methods, tiny synthetic inputs only.
using Test,Random,LinearAlgebra
include(joinpath(@__DIR__,"lib","replay.jl"))
include(joinpath(@__DIR__,"influence_core.jl"))
Random.seed!(9272);n=7
leaves=[PintElem(Matrix(8.0I,n,n),0.1randn(n,n),Matrix(9.0I,n,n),randn(n),randn(n)) for t in 1:24]
xs=Any[deepcopy(e) for e in leaves];push!(xs,PintVal(zeros(n,n),zeros(n)));w=PairWork(n)
for op in scan_ops(24)
 _,t,mid=op
 if xs[mid] isa PintElem
  z=merge_inputs(xs[t],xs[mid]);actual=input_combine!(deepcopy(xs[t]),xs[mid],PairWork(n),z)
  original=_pp_combine!(deepcopy(xs[t]),xs[mid],PairWork(n))
  @test all(f->getproperty(actual,f)==getproperty(original,f),fields_of(actual))
 end
 xs[t]=_pp_combine!(xs[t],xs[mid],w)
end
ref=production_values(leaves,PintVal(zeros(n,n),zeros(n)))
@test all(t->xs[t].V==ref[t].V&&xs[t].v==ref[t].v,1:25)
println("PRODUCTION_RUNTIME_UNIT_PASS exact_no_change_kernel full_global_association")
