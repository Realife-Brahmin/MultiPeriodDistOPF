using LinearAlgebra,SparseArrays,Random
include(joinpath(@__DIR__,"../../examples/power_system/pint_scan.jl"))
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
println("COUPLING_HIERARCHY_UNIT_PASS states_retained=true within_groups_full=true nonzero_cross_coupling_retained=true frobenius_budget=true")
