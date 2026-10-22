using Serialization,LinearAlgebra,SparseArrays,Printf,SharedArrays
include(joinpath(@__DIR__,"../../examples/power_system/pint_scan.jl"))
include(joinpath(@__DIR__,"../../examples/power_system/frozen_value_capture.jl"))
ENV["FILTERDDP_FROZEN_MIDDLE_K"]="7"
state=Dict{Symbol,Any}(:tags=>Set{String}())
options=(optimality_tolerance=1e-5,κ_ϵ=10.0)
ocp=(N=24,control_limits=(nl=2,nu=2))
point(k,j,mu,err,cs)= (k=k,j=j,μ=mu,dual_inf=err,cs_inf_0=cs,primal_inf=err,cs_inf_μ=err)
seen=String[]
for d in (point(0,0,1.0,50.0,1.0),point(2,0,1.0,1.0,1.0),point(2,1,.2,50.0,1.0),point(7,4,.001,.1,.1),point(14,8,1e-6,1e-6,1e-6))
    tags=fc_select!(state,d,options,ocp);append!(seen,tags);union!(state[:tags],tags)
end
@assert Set(seen)==Set(("early","middle","late","barrier_before","barrier_after"))
@assert count(==("middle"),seen)==count(==("late"),seen)==1
Meta.parseall(read(joinpath(@__DIR__,"../../examples/power_system/procs_pass.jl"),String))
println("FROZEN_CAPTURE_UNIT_PASS middle_j4=true late_j8=true transition_pair=true")
