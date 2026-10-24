using LinearAlgebra,Random,Serialization
include(joinpath(ENV["VALUE_REPLAY_NUMERICAL_ROOT"],"spar_ddp/examples/power_system/pint_scan.jl"))
include(joinpath(get(ENV,"VALUE_FOLLOWUP_BASE",@__DIR__),"hierarchy.jl"))
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

# Exercise the actual measurement function with tiny mocks, without loading an
# optimizer/model. Regression: nested norm accumulator must not overwrite path
# or policy variables captured from the enclosing Julia function.
const DIR_SCALED=1e-7
const DIR_RELATIVE=1e-5
emit(io,x...)=println(io,join(x,','))
battery_path(elems::Vector{Int},values,f)=[(dx=[0.0,1.0],dp=[2.0,3.0]),(dx=[4.0,5.0],dp=[6.0,7.0])]
policy(packet,future,mu)=packet.rule
directions(rule,dx)=rule
source_text=read(joinpath(@__DIR__,"followup.jl"),String)
lo=first(findfirst("function measure_dirs(",source_text))
hi=first(findfirst("\nfunction followup_scan(",source_text))-1
include_string(Main,source_text[lo:hi],"measure_dirs_regression")
rule=(du=[1.0,2.0],dphi=[3.0,4.0],dzl=[5.0,6.0],dzu=[7.0,8.0])
packet=(f=[1.0,1.0],probes=[(dx=[0.0,0.0],),(dx=[0.01,0.02],)],rule=rule)
measured=measure_dirs(IOBuffer(),"unit",[1,2],[0,0,0],[0,0,0],Dict(1=>packet),1.0,Dict(1=>rule))
@assert measured.pass && measured.scaled==0 && measured.relative==0
println("DIRECTION_MEASUREMENT_REGRESSION_PASS path_and_full_fields=true")
