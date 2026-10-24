using Serialization, LinearAlgebra
root=normpath(joinpath(@__DIR__,"../../.."))
for path in ("spar_ddp/experiments/battery_count/prepare.jl","spar_ddp/experiments/battery_count/validate_physical.jl","spar_ddp/examples/power_system/centralized_ipopt_matched.jl","spar_ddp/examples/power_system/procs_pass.jl","spar_ddp/examples/power_system/ieee123c_filterddp.jl")
    Meta.parseall(read(joinpath(root,path),String));println("PARSE_PASS ",path)
end
include(joinpath(root,"spar_ddp/examples/power_system/control_layout.jl"))
include(joinpath(root,"spar_ddp/examples/power_system/loadflow_start.jl"))
d=deserialize(ARGS[1]); idx,nu=control_layout(d);u=zeros(nu)
loadflow_start!(u,d,idx,1)
v,p,q,ell,_=loadflow_initial_values(d,1)
@assert u[idx.v]==v && u[idx.P]==p && u[idx.Q]==q && u[idx.ell]==ell
@assert u[idx.pb]==zeros(length(d[:Bset])) && u[idx.qnorm]==zeros(length(d[:Dset]))
println("MATCHED_PHYSICAL_START_PASS exact array equality")
