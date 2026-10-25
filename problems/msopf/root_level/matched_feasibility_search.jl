# Fixed-angle continuation for difficult matched cases. Physical feasibility is
# checked independently; a local positive-slack optimum is not a global proof.
const SEARCH_ARGS=copy(ARGS)
ENV["MATCHED_LIBRARY_ONLY"]="1"
include(joinpath(@__DIR__,"matched_secondary_study.jl"))
const SEARCH_LABEL=SEARCH_ARGS[1]
const SEARCH_CASE=Symbol(SEARCH_ARGS[2])
const SEARCH_ID=parse(Int,SEARCH_ARGS[3])
s,names,deck=matched_network(SEARCH_LABEL)
packet=deserialize(joinpath(MATCH_ROOT,SEARCH_LABEL,"state.jls"))
u=matchdirs(length(names))[SEARCH_ID]
at(r)=Dict(names[j] => (j==1 ? 0.0 : r*u[j-1]) for j in eachindex(names))
const TRIALS=joinpath(ROOT8500,"oracle_trials.csv")
const HEADER="configuration,case,direction,trial,radius_deg,status,iterations,seconds,replay_ok,physical_feasible,local_outside,dP_kW,dQ_kvar,dV_pu,dB_kW,vmin,vmax,min_import_kW,terminal_kWh,soc_violation_kWh,energy_balance_kWh,battery_power_violation_kW,pv_q_violation_kvar,violation"
function oracle(r,start,trial)
    x=solve_multiperiod(s.net,at(r);problem(s,SEARCH_CASE)...,objective=:violation,start,no_backflow=true,max_iter=300)
    if isempty(x.V_pu)
        appendrow(TRIALS,HEADER,[SEARCH_LABEL,SEARCH_CASE,SEARCH_ID,trial,r,x.status,x.iterations,x.time,false,false,false,fill(NaN,13)...])
        return (;x,feasible=false,outside=false)
    end
    rp=replay_opendss(s.net,x,at(r);load_factor=s.lam,pvs=s.P,pv_factor=s.F)
    dg=diagnostics8500(s,x)
    good=rp.converged && rp.pq_band && rp.dP<1e-3 && rp.dQ<1e-3 && rp.dV<1e-7 && rp.dB<1e-3 && max(dg.terminal,dg.socviol,dg.dyn,dg.bpow,dg.qviol)<.01
    feasible=good && dg.vmin>=.95-1e-6 && dg.vmax<=1.05+1e-6 && dg.minp>=-.01
    # A failed status never supplies the upper side of a local sign bracket.
    outside=good && !feasible && x.status in RAY_OK && (dg.minp<-.1 || dg.vmin<.95-1e-5 || dg.vmax>1.05+1e-5)
    appendrow(TRIALS,HEADER,[SEARCH_LABEL,SEARCH_CASE,SEARCH_ID,trial,r,x.status,x.iterations,x.time,good,feasible,outside,rp.dP,rp.dQ,rp.dV,rp.dB,dg.vmin,dg.vmax,dg.minp,dg.terminal,dg.socviol,dg.dyn,dg.bpow,dg.qviol,x.violation])
    say("ORACLE r=$r status=$(x.status) replay=$good feasible=$feasible local_outside=$outside")
    return (;x,feasible,outside)
end
# The idle trajectory is physically admissible for cyclic dispatch (zero battery
# power and constant SOC); this avoids importing a cost-induced cycling pattern.
xlo=packet.x0s[:nobatt];lo=0.0;hi=NaN;xhi=xlo;method="feasible_lower_bound";trial=0
pred,_=predicted_reach(packet.sens,Dict(zip(names[2:end],u));day=SEARCH_CASE===:cyclic)
r=clamp(isfinite(pred) && pred>0 ? .85pred : 1.0,.05,RMAX8500)
for _ in 1:6
    trial+=1;q=oracle(r,xlo,trial)
    if q.feasible
        lo=r;xlo=q.x
        r>=RMAX8500 && (method="at_limit";break)
        r=min(RMAX8500,1.3r)
    elseif q.outside
        hi=r;xhi=q.x;break
    else
        # One continuation retry can complete a nearly converged outside solve.
        trial+=1;qq=oracle(r,isempty(q.x.V_pu) ? xlo : q.x,trial)
        if qq.feasible
            lo=r;xlo=qq.x;r=min(RMAX8500,1.3r)
        elseif qq.outside
            hi=r;xhi=qq.x;break
        else
            break
        end
    end
end
if isfinite(hi)
    for _ in 1:10
        hi-lo<=.05 && break
        r=(lo+hi)/2;trial+=1;q=oracle(r,xlo,trial)
        if q.feasible
            lo=r;xlo=q.x
        elseif q.outside
            hi=r;xhi=q.x
        else
            trial+=1;qq=oracle(r,isempty(q.x.V_pu) ? xlo : q.x,trial)
            if qq.feasible;lo=r;xlo=qq.x
            elseif qq.outside;hi=r;xhi=qq.x
            else;break
            end
        end
    end
    method=hi-lo<=.05 ? "oracle_bracket" : "oracle_coarse_bracket"
end
dg=diagnostics8500(s,xlo);imports=Dict(k=>minimum(xlo.P_subs_kW[k]) for k in names)
bind=dg.minp<=1 ? argmin(imports) : (dg.vmin<.9501 || dg.vmax>1.0499 ? "voltage" : "unknown")
appendrow(EDGES8500,"configuration,case,T,direction,u1,u2,u3,predicted_deg,measured_deg,upper_deg,method,binding,min_import_kW,vmin,vmax,imports",[SEARCH_LABEL,SEARCH_CASE,24,SEARCH_ID,[j<=length(u) ? u[j] : 0.0 for j in 1:3]...,pred,lo,hi,method,bind,dg.minp,dg.vmin,dg.vmax,join(["$k=$(imports[k])" for k in names],';')])
say("COMPLETED oracle $SEARCH_LABEL $SEARCH_CASE $SEARCH_ID: $lo to $hi ($method)")
