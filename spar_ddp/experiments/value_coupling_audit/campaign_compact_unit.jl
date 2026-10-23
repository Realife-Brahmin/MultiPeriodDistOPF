include(joinpath(@__DIR__,"direct_initial_unit.jl"))
include(joinpath(@__DIR__,"campaign_compact.jl"))
# SPD four-period model, raw asymmetric diagonal lower entries, interleaved
# spatial tree, distinct directed cross matrices. No data/model simplification.
function campaign_unit_leaf()
    A=25Matrix{Float64}(I,n,n)+.01X*X';A[2,1]+=1e-9
    B=.03Matrix{Float64}(I,n,n)+.001X*X'
    PintElem(A,B,2A,randn(n),randn(n))
end
es=[campaign_unit_leaf() for _ in 1:4]
ps=[construct_blocks(e,tree,1e-12) for e in es]
T=3Matrix{Float64}(I,n,n);terminal=(rep=block_payload(T,tree,1e-12,Any[],:T),v=randn(n))
stats=Dict{String,Int}();values=campaign_build(ps,terminal;stats=stats)
rv=PintVal(T,terminal.v);refs=Vector{PintVal}(undef,5);refs[5]=rv
for t in 4:-1:1
    e=es[t];K=Matrix(Symmetric(e.Mbb+refs[t+1].V,:U));local B=e.Mab
    refs[t]=PintVal(Matrix(Symmetric(e.Maa,:U))-B*(K\B'),e.ma-B*(K\(e.mb+refs[t+1].v)))
end
for t in 1:4
    @assert norm(campaign_action(values[t],X)-refs[t].V*X)/norm(refs[t].V*X)<1e-9
    @assert norm(values[t].v-refs[t].v)/max(1,norm(refs[t].v))<1e-9
end
@assert [length(c.g) for c in values[2:4]]==[3n,2n,n]
@assert all(isfinite,campaign_conditions(values))
println("CAMPAIGN_COMPACT_UNIT_PASS all_four_values=true retained_growth=true raw_upper=true independent_dense_formula=true")
