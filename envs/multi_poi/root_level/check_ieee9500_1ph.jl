# Check rawData/ieee9500_1ph against the original IEEE 9500-node model in OpenDSS, both at
# nominal load with no DERs and the ties open: each substation's import (the original's
# three-phase import / 3) and the primary voltages, bus by bus (the original's mean phase
# magnitude). Run from the repo root after convert_ieee9500_1ph.jl:
#   julia --project=envs/msopf envs/multi_poi/root_level/check_ieee9500_1ph.jl [variant ...]
# Variants of the original, to separate what the per-phase equivalent changes from what the
# dropped regulators and capacitors change: "asis" (taps solved, capacitors in), "flat"
# (every regulator at tap 1, capacitors out, substations at 1.05 pu behind the transformer).
const VARIANTS = copy(ARGS)
include(joinpath(@__DIR__, "convert_ieee9500_1ph.jl"))   # defines solve_original, bus, ...

"Original: per-phase import of each substation and mean phase voltage (pu) per primary bus."
function original_state(variant)
    solve_original()
    if variant == "flat"
        for nm in ODD.Transformers.AllNames()
            ODD.Transformers.Name(nm); occursin(r"reg", lowercase(nm)) && (ODD.Transformers.Tap(1.0))
        end
        ODD.Text.Command("BatchEdit RegControl..* enabled=false")
        ODD.Text.Command("BatchEdit CapControl..* enabled=false")
        ODD.Text.Command("BatchEdit Capacitor..* enabled=false")
        ODD.Text.Command("Set controlmode=off")
        ODD.Text.Command("Solve")
    end
    imp = Float64[]
    for i in 1:3
        ODD.Circuit.SetActiveElement("Transformer.HVMV69_11Sub$i")
        n = ODD.CktElement.NumConductors(); pw = ODD.CktElement.Powers()
        push!(imp, -sum(real(pw[n+1:2n])) / 3)          # out of the low-side terminal
    end
    v = Dict{String,Float64}()
    for b in ODD.Circuit.AllBusNames()
        primary(b) || continue
        m = ODD.Bus.puVmagAngle()[1:2:end]
        v[lowercase(b)] = sum(m) / length(m)
    end
    return imp, v
end

"The 1-phase deck, snapshot at nominal load, no DERs, ties open."
function deck_state()
    ODD.Text.Command("Clear")
    ODD.Text.Command("Redirect \"$(joinpath(OUT, "Master.dss"))\"")
    ODD.Text.Command("BatchEdit PVSystem..* enabled=false"); ODD.Text.Command("BatchEdit Storage..* enabled=false")
    ODD.Text.Command("Set mode=Snapshot"); ODD.Text.Command("Solve")
    ODD.Solution.Converged() || error("deck did not converge")
    imp = Float64[]
    for i in 1:3
        ODD.Circuit.SetActiveElement("Vsource.subs$i"); push!(imp, -real(ODD.CktElement.Powers()[1]))
    end
    v = Dict{String,Float64}()
    for b in ODD.Circuit.AllBusNames()
        ODD.Circuit.SetActiveBus(String(b)); v[lowercase(b)] = ODD.Bus.VMagAngle()[1] / (KV_LN * 1e3)
    end
    return imp, v
end


dimp, dv = deck_state()
for variant in (isempty(VARIANTS) ? ["asis", "flat"] : VARIANTS)
    oimp, ov = original_state(variant)
    diffs = [dv[b] - ov[b] for b in keys(ov) if haskey(dv, b)]
    @printf("%-5s imports original %s kW, deck %s kW; voltage original %.4f-%.4f, deck %.4f-%.4f pu; deck - original: mean %+.4f, max |.| %.4f pu over %d buses\n",
            variant, join([@sprintf("%.0f", x) for x in oimp], "/"), join([@sprintf("%.0f", x) for x in dimp], "/"),
            minimum(values(ov)), maximum(values(ov)), minimum(values(dv)), maximum(values(dv)),
            sum(diffs) / length(diffs), maximum(abs, diffs), length(diffs))
end

# The network conversion alone, without the sources and regulators: the voltage drop from each
# feeder's head (its substation bus) to every bus of that feeder, deck against the original
# with its regulators at tap 1 and capacitors out ("flat").
c = convert(); fd = feeders(c)
oimp, ov = original_state("flat")
for i in 1:3
    h = c.subs[i]
    d = [(dv[h] - dv[b]) - (ov[h] - ov[b]) for b in keys(ov) if haskey(dv, b) && get(fd, b, 0) == i]
    od = maximum(ov[h] - ov[b] for b in keys(ov) if haskey(dv, b) && get(fd, b, 0) == i)
    dd = maximum(dv[h] - dv[b] for b in keys(dv) if get(fd, b, 0) == i)
    @printf("feeder %d: largest drop from the head, original %.4f pu, deck %.4f pu; deck - original per bus: mean %+.4f, max |.| %.4f pu\n",
            i, od, dd, sum(d) / length(d), maximum(abs, d))
end
