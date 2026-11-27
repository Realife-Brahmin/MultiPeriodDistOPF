# How constrained are the test cases at the optimum, and how tight is the
# relaxation? (reviewer items R1-3, AE-1, R5-2). Reads an instance and a saved
# dispatch and counts active limits with explicit thresholds.
#
#   julia binding_stats.jl out.csv folder [folder ...]
# each folder holds instance.jls and central_dispatch.jls (and optionally
# converged_dispatch.jls); one row is written per dispatch.
include("model.jl")
const TOL=1e-6   # p.u. (p.u. squared for voltage and the cone)

function stats(d,sol)
    ix,_=control_layout(d); T=d[:T]; ns=d[:Nset]; es=d[:Lset]; bs=d[:Bset]; ds=d[:Dset]; root=d[:substationBus]
    bp=Dict(j=>n for (n,j) in enumerate(ns))
    v_low=0; v_high=0; periods_v=0; pb_lim=0; soc_low=0; soc_high=0; q_lim=0; ps_zero=0
    cone_max=0.0; cone_loose=0; cone_rel_max=0.0; loss=0.0; fict=0.0
    for t in 1:T
        u=sol[:u][t]; v=u[ix.v]; any_v=false
        for (n,j) in enumerate(ns)
            j==root && continue
            lo=v[n]-d[:Vminpu][j]^2<=TOL; hi=d[:Vmaxpu][j]^2-v[n]<=TOL
            v_low+=lo; v_high+=hi; any_v|=lo||hi
        end
        periods_v+=any_v
        for (n,j) in enumerate(bs)
            pb_lim+=abs(u[ix.pb[n]])>=d[:P_B_R_pu][j]-TOL
            soc_low+=sol[:B][n,t]-d[:soc_min][j]*d[:B_R_pu][j]<=TOL
            soc_high+=d[:soc_max][j]*d[:B_R_pu][j]-sol[:B][n,t]<=TOL
        end
        q_lim+=count(x->abs(x)>=1-TOL,u[ix.qnorm]); ps_zero+=u[ix.ps]<=TOL
        for (n,e) in enumerate(es)
            p=u[ix.P[n]]; q=u[ix.Q[n]]; ell=u[ix.ell[n]]; vi=v[bp[e[1]]]
            slack=vi*ell-p^2-q^2; cone_max=max(cone_max,abs(slack)); cone_loose+=slack>TOL
            cone_rel_max=max(cone_rel_max,abs(slack)/max(vi*ell,1e-8))
            loss+=d[:rdict_pu][e]*ell; fict+=d[:rdict_pu][e]*max(slack,0.0)/vi
        end
    end
    nbus=(length(ns)-1)*T; nbat=length(bs)*T
    return (;T,buses=length(ns),batteries=length(bs),bus_periods=nbus,v_low,v_high,periods_with_binding_voltage=periods_v,
        battery_periods=nbat,battery_power_at_limit=pb_lim,soc_at_min=soc_low,soc_at_max=soc_high,
        inverter_periods=length(ds)*T,inverter_q_at_limit=q_lim,periods_with_zero_import=ps_zero,
        branch_periods=length(es)*T,cone_max_gap=cone_max,cone_gaps_above_tol=cone_loose,cone_max_relative_gap=cone_rel_max,
        fictitious_loss_share=fict/max(loss,1e-12))
end

function main()
    out=ARGS[1]; first=true
    open(out,"w") do io
        for folder in ARGS[2:end], name in ("central","converged")
            file=joinpath(folder,"$(name)_dispatch.jls"); isfile(file) || continue
            d=deserialize(joinpath(folder,"instance.jls")); sol=deserialize(file); s=stats(d,sol)
            first && println(io,"case,dispatch,",join(keys(s),',')); first=false
            println(io,basename(folder),",",name,",",join(values(s),',')); flush(io)
            println("STATS ",basename(folder)," ",name," v_low=",s.v_low," v_high=",s.v_high," of ",s.bus_periods," | pb_at_limit=",s.battery_power_at_limit," of ",s.battery_periods,
                " | cone_max_gap=",s.cone_max_gap," loose=",s.cone_gaps_above_tol," of ",s.branch_periods," | fictitious_loss_share=",s.fictitious_loss_share)
        end
    end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
