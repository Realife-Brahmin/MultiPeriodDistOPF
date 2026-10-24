using LinearAlgebra,SparseArrays,Random
const PS=joinpath(@__DIR__,"../../examples/power_system")
include(joinpath(PS,"pint_scan.jl"))
s=read(joinpath(PS,"paired_pass.jl"),String)
include_string(Main,s[first(findfirst("struct PairWork",s)):first(findfirst("# Collector time",s))-1])
include(joinpath(PS,"value_group_diagnostic.jl"))
Random.seed!(2718)
function elem(n)
 A=randn(2n,2n);H=A'A+10I
 PintElem(H[1:n,1:n],H[1:n,n+1:end],H[n+1:end,n+1:end],randn(n),randn(n))
end
function project(e)
 x=deepcopy(e);for f in (:Maa,:Mab,:Mbb);vg_zero_cross!(getproperty(x,f));end;x
end
rel(x,y)=maximum(norm(getproperty(x,f)-getproperty(y,f))/max(1,norm(getproperty(y,f))) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
d=Dict(:Bset=>[2,3,4,5],:substationBus=>1,:parent=>Dict(2=>1,3=>1,4=>2,5=>3))
ENV["FILTERDDP_VALUE_GROUPING"]="chunks16";vg_install(d)
x=elem(4);y=elem(4)
@assert rel(vg_combine!(deepcopy(x),y),pint_combine(x,y))<1e-12
ENV["FILTERDDP_VALUE_GROUPING"]="root";vg_install(d)
@assert length(_VG_GROUPS[])==2
@assert rel(vg_combine!(deepcopy(x),y),pint_combine(project(x),project(y)))<1e-12
v=PintVal(Matrix{Float64}(I,4,4),randn(4))
a=vg_combine!(deepcopy(x),v);b=pint_combine(project(x),v)
@assert norm(a.V-b.V)<1e-12 && norm(a.v-b.v)<1e-12
for file in ("procs_pass.jl","value_group_diagnostic.jl")
 Meta.parseall(read(joinpath(PS,file),String))
end
println("VALUE_GROUP_UNIT_PASS one_group_dense=true two_group_projected_reference=true closing=true source_parse=true")
