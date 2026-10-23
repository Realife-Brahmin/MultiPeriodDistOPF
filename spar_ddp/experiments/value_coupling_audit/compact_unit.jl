using LinearAlgebra,Random,Serialization
BLAS.set_num_threads(1)
root=ENV["VALUE_REPLAY_NUMERICAL_ROOT"]
include(joinpath(root,"spar_ddp/examples/power_system/pint_scan.jl"))
include(joinpath(@__DIR__,"hierarchy.jl"))
include(joinpath(@__DIR__,"followup_blocks.jl"))
include(joinpath(@__DIR__,"compact_kernel.jl"))
Random.seed!(20261023);n=96
# Interleaving deliberately stresses ORIGINAL-coordinate upper completion.
ids=vcat(collect(1:2:n),collect(2:2:n))
function manual_tree(v)
    length(v)<=6 && return coupling_leaf(v)
    k=length(v)÷2;l=manual_tree(v[1:k]);r=manual_tree(v[k+1:end])
    CouplingTree(vcat(l.ids,r.ids),l,r,Int[])
end
tree=manual_tree(ids)
function unit_elem()
    U=randn(n,3);H=20Matrix{Float64}(I,n,n)+U*U'
    H[2,1]+=1e-8
    B=randn(n,2)*randn(2,n)+.1Matrix{Float64}(I,n,n)
    PintElem(H,B,2H,randn(n),randn(n))
end
x=unit_elem();y=unit_elem();p=construct_blocks(x,tree,1e-12);q=construct_blocks(y,tree,1e-12)
s=Dict{String,Int}();c=ck_merge(p,q;tol=1e-12,stats=s)
X=randn(n,3)
@assert norm(ck_apply(p.reps[2],X)-x.Mab*X)/norm(x.Mab*X)<1e-10
@assert norm(ck_apply(p.reps[2],X;trans=true)-x.Mab'*X)/norm(x.Mab'*X)<1e-10
KU=Matrix(Symmetric(x.Mbb+y.Maa,:U))
@assert norm(ck_apply(c.K,X)-KU*X)/norm(KU*X)<1e-10
@assert norm(ck_solve(c.factor,X)-KU\X)/norm(KU\X)<1e-9
ref=pint_combine(x,y)
for f in (:Maa,:Mab,:Mbb)
    # Allocating pint_combine may handle raw diagonal skew differently; compare
    # explicit upper-completed production formulas for the diagonal blocks.
    B=x.Mab;C=y.Mab
    A=f==:Maa ? Matrix(Symmetric(x.Maa,:U))-B*(KU\B') : f==:Mab ? -B*(KU\C) : Matrix(Symmetric(y.Mbb,:U))-C'*(KU\C)
    @assert norm(ck_action(c,f,X)-A*X)/max(1,norm(A*X))<1e-9
end
g=x.mb+y.ma
@assert norm(c.ma-(x.ma-x.Mab*(KU\g)))<1e-9
@assert norm(c.mb-(y.mb-y.Mab'*(KU\g)))<1e-9
@assert s["max_local_dense_dim"]<=CK_LOCAL_DENSE_CAP
@assert s["interleaved_upper_nodes"]>0
println("COMPACT_UNIT_PASS interleaved_original_upper=true raw_asymmetry=true matvec_transpose=true Woodbury_solve=true Schur_blocks_vectors=true max_local_dense=$(s["max_local_dense_dim"])")
