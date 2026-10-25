# Voltage-banded reference matching the paper; boundary measurements are independent.
const REFERENCE_ARGS=copy(ARGS)
ENV["SECONDARY_LIBRARY_ONLY"]="1"
include(joinpath(@__DIR__,"secondary_8500_study.jl"))
function matched_network(label)
    n=parse(Int,string(last(label)));names=["subs$i" for i in 1:n]
    if startswith(label,"sf")
        deck="smartds_p1u_xf_1ph";extra=variant_switches(deck,"Master_4sub_typeIII.dss")
    else
        deck=startswith(label,"secondary") ? "ieee9500xf_1ph" : "ieee9500_1ph"
        extra=["Line.$t.enabled=true" for t in TIES_9500]
        if n==4
            # Preserve the historical 2705 fourth-source site and its positive delivery impedance.
            append!(extra,["New Line.poi4 phases=1 bus1=4s bus2=l2804256 r1=1e-6 x1=1e-6 c1=0 c0=0",
             "New Vsource.subs4 phases=1 bus1=4s basekv=7.1996 pu=1.05 angle=0 R1=0.104186 X1=0.617339 R0=0.104186 X0=0.617339"])
        end
    end
    net,pvs=load_sources(names;system=deck,extra)
    # Same delivery-point voltage contract in both representations, including direct-MV loads.
    bs=union(customer_buses(net),Set(l.bus for l in net.loads))
    s=setup(net,pvs;T=24,feas=1e-7,band_buses=bs)
    @assert Set(x.name for x in net.sources)==Set(names)
    say("$label: $(length(net.buses)) buses, $(length(bs)) delivery buses, $(length(net.batteries)) batteries; T=24")
    return s,names,deck
end
s,names,deck=matched_network(REFERENCE_ARGS[1])
sens=sensitivity8500(s,REFERENCE_ARGS[1],names;V_band=(.95,1.05),h=.1)
sens===nothing && error("Voltage-banded reference failed")
println("COMPLETED voltage-banded reference at h=0.1 degrees, ",REFERENCE_ARGS[1])
