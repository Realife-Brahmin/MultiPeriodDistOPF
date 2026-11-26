using Serialization, LinearAlgebra, Printf
C=ARGS[1]
include(joinpath(C,"source/spar_ddp/examples/power_system/control_layout.jl"))
rel(a,b)=norm(a-b)/max(norm(a),norm(b),1.0)
println("system,repeat,mask,relative_x,relative_u,relative_battery_dispatch,relative_import,finite,status")
for sys in ("ieee123C_1ph","ieee2522C_1ph","large10kC_1ph")
    base=deserialize(joinpath(C,"results","$(sys)_T96_r1_m7","solution.jls"))
    data=deserialize(joinpath(C,"source/spar_ddp/results/network_filterddp/network_data_$(sys)_T96_periodic.jls"))
    idx,_=control_layout(data)
    for rep in 1:2, mask in 0:7
        s=deserialize(joinpath(C,"results","$(sys)_T96_r$(rep)_m$(mask)","solution.jls"))
        @assert s[:system]==base[:system] && s[:T]==base[:T]
        finite=all(v->all(isfinite,v),vcat(s[:x],s[:u]))
        @assert finite && s[:status]==0
        dx=maximum(rel.(base[:x],s[:x])); du=maximum(rel.(base[:u],s[:u]))
        dp=maximum(rel(a[idx.pb],b[idx.pb]) for (a,b) in zip(base[:u],s[:u]))
        ds=maximum(abs(a[idx.ps]-b[idx.ps])/max(abs(a[idx.ps]),abs(b[idx.ps]),1.0) for (a,b) in zip(base[:u],s[:u]))
        @printf("%s,%d,%d,%.12e,%.12e,%.12e,%.12e,%s,%d\n",sys,rep,mask,dx,du,dp,ds,string(finite),s[:status])
    end
end
