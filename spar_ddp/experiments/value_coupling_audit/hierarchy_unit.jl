using LinearAlgebra,SparseArrays,Random,Serialization
unit_root=normpath(get(ENV,"VALUE_REPLAY_NUMERICAL_ROOT",joinpath(@__DIR__,"../../..")))
include(joinpath(unit_root,"spar_ddp/examples/power_system/pint_scan.jl"))
include(joinpath(@__DIR__,"hierarchy.jl"))
Random.seed!(271828);n=32
data=Dict(:Nset=>collect(1:n+1),:Lset=>[(i,i+1) for i in 1:n],:Bset=>collect(2:n+1))
tree=topology_hierarchy(data;maxleaf=4)
@assert sort(vcat(coupling_groups(tree)...))==collect(1:n)
@assert maximum(length,coupling_groups(tree))<=4
H=Matrix{Float64}(I,2n,2n)
for g in coupling_groups(tree)
    ids=vcat(g,n.+g);A=randn(length(ids),length(ids));H[ids,ids]+=A'A
end
U=randn(2n,3);H+=U*U'
e=PintElem(H[1:n,1:n],H[1:n,n+1:end],H[n+1:end,n+1:end],randn(n),randn(n))
x=retain_coupling(e,tree;tol=1e-8)
@assert x.relative_frobenius_error<=1e-8
@assert x.elem.ma==e.ma && x.elem.mb==e.mb
@assert any(!r.dense_fallback && r.rank>0 for r in x.reports)
@assert norm(x.elem.Mab)>0
exact=retain_coupling(e,tree;tol=0.0)
@assert norm(joint_value(exact.elem)-joint_value(e))/norm(joint_value(e))<1e-12
asymmetric=deepcopy(e);asymmetric.Maa[1,2]+=1e-3
for tol in (0.0,1e-12)
    fallback=retain_coupling(asymmetric,tree;tol=tol)
    @assert fallback.dense_raw_fallback && fallback.relative_frobenius_error==0.0
    @assert fallback.numeric_payload_bytes==8*(3n*n+2n)
    @assert all(getproperty(fallback.elem,f)==getproperty(asymmetric,f) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
end
println("COUPLING_HIERARCHY_UNIT_PASS states_retained=true within_groups_full=true nonzero_cross_coupling_retained=true frobenius_budget=true")
if length(ARGS)==2
    payload=deserialize(ARGS[1]);model=deserialize(ARGS[2]);hier=topology_hierarchy(model;maxleaf=8)
    for t in (1,length(payload.elems)),tol in (0.0,1e-8,1e-5)
        raw=payload.elems[t];z=retain_coupling(raw,hier;tol=tol)
        @assert z.relative_frobenius_error<=tol+1e-12
        @assert z.elem.ma==raw.ma && z.elem.mb==raw.mb
        z.dense_raw_fallback && @assert all(getproperty(z.elem,f)==getproperty(raw,f) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
        println("REAL_LEAF_BUDGET_PASS stage=$t tol=$tol error=$(z.relative_frobenius_error) raw_fallback=$(z.dense_raw_fallback) symmetry_defect=$(z.symmetry_defect_relative)")
    end
end
