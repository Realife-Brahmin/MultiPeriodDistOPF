# July split implementation: primal-first phase, then gated residual balancing.
# Residuals use the same owned/unique RMS definitions as the historical stop.
function legacy_rho_update(rho,r,s,k,seen;eps_pri=1e-4)
    k%5==0 || return rho,seen
    seen |= r<=eps_pri
    next=rho
    if !seen
        r>eps_pri && (next=min(1e6,2rho))
    elseif r>5s
        next=min(1e6,2rho)
    elseif s>5r && r<=eps_pri
        next=max(1.0,rho/5)
    end
    return next,seen
end
