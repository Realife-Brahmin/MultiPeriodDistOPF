using LinearAlgebra,Random,Serialization
include(joinpath(ENV["VALUE_REPLAY_NUMERICAL_ROOT"],"spar_ddp/examples/power_system/pint_scan.jl"))
include(joinpath(@__DIR__,"hierarchy.jl"))
include(joinpath(@__DIR__,"followup_blocks.jl"))
Random.seed!(1223);n=24
tree=topology_hierarchy(Dict(:Nset=>collect(1:n),:Lset=>[(i,i+1) for i in 1:n-1],:Bset=>collect(1:n));maxleaf=4)
X=randn(n,2);A=10Matrix{Float64}(I,n,n)+X*X';B=randn(n,2)*randn(2,n)
A[1,2]+=1e-8 # Preserve real raw skew, not an artificial symmetric projection.
e=PintElem(copy(A),B,3A,randn(n),randn(n))
for tol in (0.0,1e-6,1e-10)
    p=construct_blocks(e,tree,tol);m=materialize_blocks(p,e,tol)
    @assert m.error<=tol+32eps(Float64)
    @assert m.elem.ma==e.ma && m.elem.mb==e.mb
    for f in (:Maa,:Mab,:Mbb),g in coupling_groups(tree)
        @assert getproperty(m.elem,f)[g,g]==getproperty(e,f)[g,g]
    end
    if tol==0
        @assert all(getproperty(m.elem,f)==getproperty(e,f) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
        @assert m.bytes==8*(3n*n+2n)
    end
end
@assert controlled_tol(1e-6,1e8)≈1e-10
@assert controlled_tol(1e-6,Inf)==0.0
@assert elimination_condition([1.0 0;1e-8 2])>=1
println("BLOCK_CONTROL_UNIT_PASS raw_skew=true local_blocks_exact=true vectors_exact=true zero_tol_bit_exact=true")
