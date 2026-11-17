# Feeders directly below the substation and batteries per feeder, per deck
# (what FILTERDDP_VALUE_BLOCKS can split the battery value function into).
const REPO = pwd()
include(joinpath(REPO, "test_networks", "parse_opendss.jl"))
for (sys, kv) in (("large10kC_1ph", 12.47), ("ieee2522C_1ph", 7.2), ("ieee123C_1ph", 2.4018),
                  ("ieee9500xf_1ph", 12.47), ("smartds_p1u_xf_1ph", 12.47))
    d = try
        redirect_stdout(devnull) do
            parse_system_from_dss(sys, 2; kV_B=kv)
        end
    catch e
        println(rpad(sys, 20), " parse failed: ", sprint(showerror, e)[1:min(end, 200)]); continue
    end
    sub = d[:substationBus]; par = d[:parent]
    bat = collect(d[:Bset])
    tops = Int[]
    for b in bat
        b == sub && continue
        k = b
        while par[k] != sub; k = par[k]; end
        push!(tops, k)
    end
    nchild = count(j -> get(par, j, 0) == sub, collect(keys(par)))
    cnt = sort([count(==(k), tops) for k in unique(tops)]; rev=true)
    nB = length(tops)
    work = isempty(cnt) ? 0.0 : sum(c^3 for c in cnt) / max(nB, 1)^3
    println(rpad(sys, 20), " buses=", length(d[:Nset]), " batteries=", nB, " substation children=", nchild,
            " feeders with batteries=", length(cnt), " largest blocks=", cnt[1:min(6, end)],
            " | battery-block work vs dense: ", round(work; sigdigits=3))
end
