# July split implementation: primal-first phase, then gated residual balancing.
# Residuals use the same owned/unique RMS definitions as the historical stop.
function legacy_rho_update(rho,r,s,k,seen;eps_pri=1e-4,tau_decrease=5.0)
    isfinite(tau_decrease) && tau_decrease>1 || error("Invalid rho decrease factor")
    k%5==0 || return rho,seen
    seen |= r<=eps_pri
    next=rho
    if !seen
        r>eps_pri && (next=min(1e6,2rho))
    elseif r>5s
        next=min(1e6,2rho)
    elseif s>5r && r<=eps_pri
        next=max(1.0,rho/tau_decrease)
    end
    return next,seen
end

function legacy_rho_watchdog(rho,r,count;eps_pri=1e-4)
    if r>2eps_pri
        count+=1
    elseif r<0.5eps_pri
        count=0
    end
    if count>=15
        return min(1e6,2.5rho),0
    end
    return rho,count
end
