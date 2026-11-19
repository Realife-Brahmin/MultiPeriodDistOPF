# Pack the same original-model checks into concrete numeric arrays once.
# No constraint, tolerance, or diagnostic frequency is changed.
function typed_physical_layout(d)
    key=(d[:C_B],d[:T])
    cache=get!(d,:_typed_diagnostic_cache,Dict{Any,Any}())
    return get!(cache,key) do
        old=physical_layout(d); (;ix,bp,ep,bat,der,root,loads)=old
        ns=d[:Nm1set]; es=d[:Lset]; bs=d[:Bset]; T=d[:T]
        incoming=Int[ep[(d[:parent][j],j)] for j in ns]
        children=[Int[ep[(j,k)] for k in d[:children][j]] for j in ns]
        pl=Float64[j in loads ? d[:p_L_pu][j,t] : 0.0 for j in ns,t in 1:T]
        ql=Float64[j in loads ? d[:q_L_pu][j,t] : 0.0 for j in ns,t in 1:T]
        pd=Float64[haskey(der,j) ? d[:p_D_pu][j,t] : 0.0 for j in ns,t in 1:T]
        dq=Float64[haskey(der,j) ? sqrt(max(0.0,d[:S_D_R][j]^2-pd[n,t]^2)) : 0.0 for (n,j) in enumerate(ns),t in 1:T]
        (;ix,T,incoming,children,pl,ql,pd,dq,root=bp[root],
            rootedges=Int[ep[e] for e in d[:L1set]],
            battery=Int[get(bat,j,0) for j in ns],der=Int[get(der,j,0) for j in ns],
            resistance=Float64[d[:rdict_pu][e] for e in es],reactance=Float64[d[:xdict_pu][e] for e in es],
            parent=Int[bp[e[1]] for e in es],child=Int[bp[e[2]] for e in es],
            vmin=Float64[d[:Vminpu][j]^2 for j in d[:Nset]],vmax=Float64[d[:Vmaxpu][j]^2 for j in d[:Nset]],
            b0=Float64[d[:B0_pu][j] for j in bs],pmax=Float64[d[:P_B_R_pu][j] for j in bs],
            bmin=Float64[d[:soc_min][j]*d[:B_R_pu][j] for j in bs],
            bmax=Float64[d[:soc_max][j]*d[:B_R_pu][j] for j in bs],
            dt=Float64(d[:delta_t_h]),base=Float64(d[:kVA_B]),cb=Float64(d[:C_B]),price=Float64[d[:LoadShapeCost][t] for t in 1:T])
    end
end

function typed_physical_metrics(us,bb,c)
    ix=c.ix
    maxeq=0.0; maxbound=0.0; maxcone=0.0; gap=0.0; loss_gap=0.0; obj=0.0
    recovered_energy=copy(c.b0)
    for t in 1:c.T
        u=us[t]; p=u[ix.P]; q=u[ix.Q]; ell=u[ix.ell]; v=u[ix.v]; pb=u[ix.pb]
        maxeq=max(maxeq,abs(u[ix.ps]-sum(p[n] for n in c.rootedges)),
            abs(u[ix.qs]-sum(q[n] for n in c.rootedges)),abs(v[c.root]-1.05^2))
        maxbound=max(maxbound,-u[ix.ps])
        for j in eachindex(c.incoming)
            n=c.incoming[j]; bn=c.battery[j]; dn=c.der[j]
            pbat=bn==0 ? 0.0 : pb[bn]
            qd=dn==0 ? 0.0 : u[ix.qnorm[dn]]*c.dq[j,t]
            maxeq=max(maxeq,abs(sum(p[k] for k in c.children[j];init=0.0)-p[n]+c.resistance[n]*ell[n]-pbat-c.pd[j,t]+c.pl[j,t]),
                abs(sum(q[k] for k in c.children[j];init=0.0)-q[n]+c.reactance[n]*ell[n]-qd+c.ql[j,t]))
        end
        for n in eachindex(c.resistance)
            i=c.parent[n]; j=c.child[n]; r=c.resistance[n]; x=c.reactance[n]
            maxeq=max(maxeq,abs(v[j]-v[i]+2*(r*p[n]+x*q[n])-(r^2+x^2)*ell[n]))
            cg=v[i]*ell[n]-p[n]^2-q[n]^2
            maxcone=max(maxcone,-cg); gap=max(gap,abs(cg))
            loss_gap+=r*cg/v[i]*c.base*c.dt
            maxbound=max(maxbound,-ell[n])
        end
        for n in eachindex(c.vmin)
            maxbound=max(maxbound,c.vmin[n]-v[n],v[n]-c.vmax[n])
        end
        for n in eachindex(c.b0)
            recovered_energy[n]-=c.dt*pb[n]
            prev=t==1 ? c.b0[n] : bb[n,t-1]
            maxeq=max(maxeq,abs(bb[n,t]-prev+c.dt*pb[n]))
            maxbound=max(maxbound,abs(pb[n])-c.pmax[n],c.bmin[n]-bb[n,t],bb[n,t]-c.bmax[n],
                c.bmin[n]-recovered_energy[n],recovered_energy[n]-c.bmax[n])
        end
        maxbound=max(maxbound,maximum(abs.(u[ix.qnorm]);init=0.0)-1)
        obj+=c.price[t]*c.base*c.dt*u[ix.ps]+c.cb*c.base^2*c.dt*sum(abs2,pb)
    end
    return (;objective=obj,original_primal=max(maxeq,maxbound,maxcone),maxeq,maxbound,maxcone,soc_equality_gap=gap,fictitious_loss_kwh=loss_gap)
end
