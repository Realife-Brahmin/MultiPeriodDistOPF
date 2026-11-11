# Historical tadmm_socp.jl metric: one owned energy copy per battery/time.
# User-selected primal threshold 1e-4, not the later relaxed July defaults.
function historical_rms(copies,z,zold,windows,rho)
    nb,T=size(z)
    owned2=0.0
    for t in 1:T
        a=findfirst(==(t),windows[t])
        a===nothing && error("Missing owned-period energy copy")
        for n in 1:nb
            owned2+=(copies[t][n,a]-z[n,t])^2
        end
    end
    return sqrt(owned2/(nb*T)),rho*sqrt(sum(abs2,z.-zold)/(nb*T))
end

function historical_stop(r,s,objectives)
    stagnant=false
    if length(objectives)>=5
        recent=@view objectives[end-4:end]
        hi=maximum(recent)
        stagnant=hi>0 && (hi-minimum(recent))/abs(hi)<1e-3
    end
    return (;met=r<=1e-4 && (s<=1e-2 || stagnant),stagnant)
end
