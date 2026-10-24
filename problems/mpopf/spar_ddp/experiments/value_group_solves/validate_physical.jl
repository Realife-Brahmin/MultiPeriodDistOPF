using Serialization, LinearAlgebra, Printf
root=abspath(ARGS[1]); T=parse(Int,ARGS[2]); sol=deserialize(ARGS[3])
system=ARGS[4]
cb,gamma,tol=system=="ieee123C_1ph" ? (1e-3,1.75e4,1e-6) : system=="ieee2522C_1ph" ? (1.9e-3,3.29e4,1e-5) : system=="large10kC_1ph" ? (2.4e-5,420.0,1e-5) : error("Unknown system")
include(joinpath(root,"spar_ddp/examples/power_system/control_layout.jl"))
d=deserialize(joinpath(root,"spar_ddp/results/network_filterddp/network_data_$(system)_T$(T)_periodic.jls"))
idx,nu=control_layout(d);@assert length(sol[:u])==T
bp=Dict(j=>k for (k,j) in enumerate(d[:Nset]));lp=Dict(e=>k for (k,e) in enumerate(d[:Lset]))
bpos=Dict(j=>k for (k,j) in enumerate(d[:Bset]));dp=Dict(j=>k for (k,j) in enumerate(d[:Dset]))
energy=[d[:B0_pu][j] for j in d[:Bset]]; initial=copy(energy)
vmax=0.0;objective=0.0
for t in 1:T
 u=sol[:u][t];@assert length(u)==nu && all(isfinite,u)
 global vmax=max(vmax,abs(u[idx.ps]-sum(u[idx.P[lp[e]]] for e in d[:L1set])),abs(u[idx.qs]-sum(u[idx.Q[lp[e]]] for e in d[:L1set])),abs(u[idx.v[bp[d[:substationBus]]]]-1.05^2),-u[idx.ps])
 if haskey(sol,:x); global vmax=max(vmax,norm(sol[:x][t]-energy,Inf)); end
 for j in d[:Nm1set]
  e=(d[:parent][j],j);k=lp[e];pL=j in d[:NLset] ? d[:p_L_pu][j,t] : 0.0;qL=j in d[:NLset] ? d[:q_L_pu][j,t] : 0.0
  pD=j in d[:Dset] ? d[:p_D_pu][j,t] : 0.0;pb=haskey(bpos,j) ? u[idx.pb[bpos[j]]] : 0.0
  qD=haskey(dp,j) ? sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))*u[idx.qnorm[dp[j]]] : 0.0
  p=sum((u[idx.P[lp[(j,c)]]] for c in d[:children][j]);init=0.0)-u[idx.P[k]]+d[:rdict_pu][e]*u[idx.ell[k]]-pb-pD+pL
  q=sum((u[idx.Q[lp[(j,c)]]] for c in d[:children][j]);init=0.0)-u[idx.Q[k]]+d[:xdict_pu][e]*u[idx.ell[k]]-qD+qL
  global vmax=max(vmax,abs(p),abs(q))
 end
 for (k,e) in enumerate(d[:Lset])
  vi=u[idx.v[bp[e[1]]]];vj=u[idx.v[bp[e[2]]]];p=u[idx.P[k]];q=u[idx.Q[k]];ell=u[idx.ell[k]];r=d[:rdict_pu][e];x=d[:xdict_pu][e]
  global vmax=max(vmax,abs(vj-vi+2*(r*p+x*q)-(r^2+x^2)*ell),p^2+q^2-vi*ell,-ell)
 end
 for (k,j) in enumerate(d[:Nset]);global vmax=max(vmax,d[:Vminpu][j]^2-u[idx.v[k]],u[idx.v[k]]-d[:Vmaxpu][j]^2);end
 for k in idx.qnorm;global vmax=max(vmax,abs(u[k])-1);end
 for (k,j) in enumerate(d[:Bset])
  pb=u[idx.pb[k]];energy[k]-=d[:delta_t_h]*pb
  global vmax=max(vmax,abs(pb)-d[:P_B_R_pu][j],d[:soc_min][j]*d[:B_R_pu][j]-energy[k],energy[k]-d[:soc_max][j]*d[:B_R_pu][j])
 end
 global objective+=d[:LoadShapeCost][t]*d[:kVA_B]*d[:delta_t_h]*u[idx.ps]+cb*d[:kVA_B]^2*d[:delta_t_h]*sum(abs2,u[idx.pb])
end
objective+=gamma*sum(abs2,energy-initial)
@printf("VALUE_GROUP_PHYSICAL original_max=%.12e objective=%.12f\n",vmax,objective)
@assert isfinite(objective) && vmax<=tol
