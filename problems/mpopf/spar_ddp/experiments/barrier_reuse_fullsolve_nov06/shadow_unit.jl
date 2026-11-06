# Actual diagnostic reference must be independent and leave shared leaves untouched.
using Test,Random,LinearAlgebra,SHA,Serialization
include(joinpath(@__DIR__,"unit_bootstrap.jl"))
include(joinpath(@__DIR__,"right_kernel.jl"))
include(joinpath(@__DIR__,"schedule.jl"))
include(joinpath(@__DIR__,"cache_kernel.jl"))
const BR_DEBUG=true;const BR_DEBUG_ROWS=Ref(0)
const _PP_SLOTS=PintElem[]
study_scan_rounds(N,grouped)=replay_rounds(N,grouped ? 2 : 1)
study_thread_check()=(@assert BLAS.get_num_threads()==2;2)
BLAS.set_num_threads(2)
const SOURCE_AST=Meta.parseall(read(joinpath(@__DIR__,"integration.jl"),String))
wanted=Set((:br_shadow_worker,:br_check,:br_digest));loaded=Set{Symbol}()
for ex in SOURCE_AST.args
 ex isa Expr||continue
 if ex.head==:function && ex.args[1] isa Expr && ex.args[1].head==:call
  name=ex.args[1].args[1]
  if name in wanted;Core.eval(Main,ex);push!(loaded,name);end
 end
end
@test wanted==loaded
Random.seed!(31724);n=7;N=24
for _ in 1:N;push!(_PP_SLOTS,PintElem(Matrix(8.0I,n,n),randn(n,n),Matrix(9.0I,n,n),randn(n),randn(n)));end
before=deepcopy(_PP_SLOTS)
ref=Any[deepcopy(_PP_SLOTS)...,PintVal(zeros(n,n),zeros(n))]
for ops in replay_rounds(N,2),(t,mid) in ops;ref[t]=right_combine!(ref[t],ref[mid],PairWork(n));end
out=br_shadow_worker(N,n)
@test length(out.nodes)==52 && length(out.values)==25
for t in 1:N,field in (:Maa,:Mab,:Mbb,:ma,:mb)
 @test getproperty(_PP_SLOTS[t],field)==getproperty(before[t],field)
end
for t in 1:N+1,field in (:V,:v)
 @test reinterpret(UInt64,vec(getproperty(out.values[t],field)))==reinterpret(UInt64,vec(getproperty(ref[t],field)))
end
println("BR_SHADOW_UNIT_PASS actual_worker_reference_full_grouped_graph_all_outputs_bitexact_no_leaf_mutation")
