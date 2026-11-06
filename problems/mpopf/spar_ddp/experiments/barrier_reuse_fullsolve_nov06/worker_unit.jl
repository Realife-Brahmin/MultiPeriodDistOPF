# Exercise the actual integrated worker admission and invalidation methods.
using Test,LinearAlgebra,Random
const TEST_STUDY=get(ENV,"BR_STUDY_ROOT",@__DIR__)
include(joinpath(TEST_STUDY,"unit_bootstrap.jl"))
include(joinpath(TEST_STUDY,"right_kernel.jl"))
include(joinpath(TEST_STUDY,"cache_kernel.jl"))
const BR_ON=true;const BR_DEBUG=false;const BR_GATE=BRGate()
const BR_ENTRIES=Dict{Tuple{Int,Int,Int},BRMatrixEntry}()
const BR_MODE=Ref(:fresh);const BR_CONTEXT=Ref{Any}(nothing)
const BR_EXPECTED=Ref(0);const BR_VALID=Ref(false);const BR_FAILED=Ref(false)
const BR_LOCAL_INVALID=Ref(false);const BR_STAGE_KEYS=Dict{Int,String}()
const BR_REF_RULES=Dict{Int,Any}()
const BR_DEBUG_ROWS=Ref(0)
const INTEGRATION_AST=Meta.parseall(read(joinpath(TEST_STUDY,"integration.jl"),String))
wanted=Set((:br_release!,:br_worker_start,:br_force_population,:br_worker_finish,:br_leaf!,:barrier_cache_invalidate!))
loaded=Set{Symbol}()
for ex in INTEGRATION_AST.args
 ex isa Expr||continue
 if ex.head in (:function,:(=)) && ex.args[1] isa Expr && ex.args[1].head==:call
  name=ex.args[1].args[1]
  if name in wanted;Core.eval(Main,ex);push!(loaded,name);end
 end
end
@test loaded==wanted
Random.seed!(9471);n=7
x=PintElem(Matrix(8.0I,n,n),randn(n,n),Matrix(9.0I,n,n),randn(n),randn(n))
y=deepcopy(x);_,entry=br_populate!(deepcopy(x),y,PairWork(n))
function seed!()
 br_release!();BR_ENTRIES[(1,1,2)]=deepcopy(entry)
 BR_VALID[]=true;BR_EXPECTED[]=1;BR_CONTEXT[]=:same;BR_LOCAL_INVALID[]=false
end
seed!();@test br_worker_start(true,:same,1,false)&&BR_MODE[]===:hit
for reason in (:accepted,:state,:matrix,:model,:factor,:terminal,:regularization,:path,:unknown)
 seed!();barrier_cache_invalidate!(reason)
 @test !br_worker_start(true,:same,1,false)
 @test BR_MODE[]===:populate && !BR_LOCAL_INVALID[]
end
seed!();@test !br_worker_start(true,:changed,1,false);@test BR_MODE[]===:populate
seed!();empty!(BR_ENTRIES);@test !br_worker_start(true,:same,1,false)
seed!();@test !br_worker_start(true,:same,1,true);@test isempty(BR_ENTRIES)
seed!();br_worker_start(true,:same,1,false);barrier_cache_invalidate!(:matrix)
@test_throws AssertionError br_leaf!(nothing,1,n)
# On admission failure, fresh population recomputes changed matrix work;
# it may reuse storage, but cannot retain the old factor/block contents.
for f in (:Maa,:Mab,:Mbb)
 changed=deepcopy(x);getproperty(changed,f)[1,1]+=0.01
 actual,c=br_populate!(deepcopy(changed),y,PairWork(n),deepcopy(entry))
 ref=right_combine!(deepcopy(changed),y,PairWork(n))
 for out in (:Maa,:Mab,:Mbb,:ma,:mb)
  @test reinterpret(UInt64,vec(getproperty(actual,out)))==reinterpret(UInt64,vec(getproperty(ref,out)))
 end
end
br_release!();@test isempty(BR_ENTRIES)&&BR_EXPECTED[]==0&&!BR_VALID[]
println("BR_WORKER_UNIT_PASS actual_admission_state_matrix_model_factor_terminal_reg_path_unknown_context_missing_reset_late_invalidation_changed_matrix_fresh")
