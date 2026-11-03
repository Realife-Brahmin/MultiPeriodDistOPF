# Check each converted deck (test_networks/<system>_1ph, or CONVERT_OUT/<system>_1ph) against its
# original in OpenDSS, at nominal load with no DERs and ties open:
#   * each substation's import: the original's three-phase import / 3, out of the substation
#     transformer's MV terminal, against the deck's Vsource;
#   * the voltage drop from each substation bus to every bus it feeds, per unit (the original's
#     mean phase magnitude), deck against the original.
# The original is solved as the deck represents it: regulator taps 1 and their controls off,
# capacitors out (in the deck too), the loads the converter drops (other feeders) out, SMART-DS loads at
# the converter's peak scale, no DERs. The upstream
# source impedance is in the original only, so absolute voltages differ; drops from the head
# do not depend on it.
#   julia --project=problems/msopf problems/msopf/root_level/check_1ph.jl [system ...]
include(joinpath(@__DIR__, "convert_to_1ph.jl"))

function original_flat(c)
    cfg = c.cfg
    compile_original(cfg)
    if cfg.drop_loads !== nothing
        for nm in real_names(ODD.Loads.AllNames())
            occursin(cfg.drop_loads, lowercase(nm)) && ODD.Text.Command("Load.$nm.enabled=false")
        end
    end
    scale = peak_scale(cfg)
    if scale !== nothing
        for nm in real_names(ODD.Loads.AllNames())
            ODD.Loads.Name(nm); kw, kvar = ODD.Loads.kW(), ODD.Loads.kvar(); k = scale[lowercase(nm)]
            ODD.Loads.kW(k * kw); ODD.Loads.kvar(k * kvar)
        end
    end
    for nm in real_names(ODD.Transformers.AllNames())
        ODD.Circuit.SetActiveElement("Transformer.$nm"); bs = bus.(ODD.CktElement.BusNames())
        length(bs) == 2 || continue
        k1, k2 = kvll(bs[1]), kvll(bs[2])
        ismv(k1) && ismv(k2) && abs(k1 - k2) / k1 < 0.02 && (ODD.Transformers.Name(nm); ODD.Transformers.Tap(1.0))
    end
    for cls in ("RegControl", "CapControl", "Capacitor")
        any(startswith(lowercase(e), lowercase(cls) * ".") for e in ODD.Circuit.AllElementNames()) &&
            ODD.Text.Command("BatchEdit $cls..* enabled=false")
    end
    ODD.Text.Command("Set controlmode=off"); ODD.Solution.Solve()
    ODD.Solution.Converged() || @warn "$(c.name): original (flat) did not converge"
    imp = Dict{String,Float64}()
    for nm in real_names(ODD.Transformers.AllNames())
        ODD.Circuit.SetActiveElement("Transformer.$nm"); bs = bus.(ODD.CktElement.BusNames())
        ks = kvll.(bs); m = findall(ismv, ks)
        (length(m) == 1 && any(>(MV[2]), ks)) || continue
        n = ODD.CktElement.NumConductors(); pw = ODD.CktElement.Powers()
        w = m[1]; P = -sum(real(pw[(w-1)*n+1:w*n])) / 3      # out of the MV terminal, per phase
        imp[bs[w]] = get(imp, bs[w], 0.0) + P
    end
    v = Dict{String,Float64}()
    for b in ODD.Circuit.AllBusNames()
        ODD.Circuit.SetActiveBus(String(b)); ismv(ODD.Bus.kVBase() * sqrt(3)) || continue
        m = ODD.Bus.puVmagAngle()[1:2:end][[n in 1:3 for n in ODD.Bus.Nodes()]]   # phases only, not a neutral node
        isempty(m) || (v[lowercase(b)] = sum(m) / length(m))
    end
    return imp, v
end

function deck_flat(c, dir)
    ODD.Text.Command("Clear"); ODD.Text.Command("Redirect \"$(joinpath(dir, "Master.dss"))\"")
    ODD.Text.Command("BatchEdit PVSystem..* enabled=false"); ODD.Text.Command("BatchEdit Storage..* enabled=false")
    any(startswith(lowercase(e), "capacitor.") for e in ODD.Circuit.AllElementNames()) && ODD.Text.Command("BatchEdit Capacitor..* enabled=false")
    ODD.Text.Command("Set mode=Snapshot"); ODD.Text.Command("Set maxiterations=100"); ODD.Solution.Solve()
    conv = ODD.Solution.Converged()
    imp = Float64[]
    for i in eachindex(c.subs)
        ODD.Circuit.SetActiveElement("Vsource.subs$i"); push!(imp, -real(ODD.CktElement.Powers()[1]))
    end
    base = c.cfg.base_kV / sqrt(3) * 1e3
    v = Dict{String,Float64}()
    for b in ODD.Circuit.AllBusNames()
        ODD.Circuit.SetActiveBus(String(b)); v[lowercase(b)] = ODD.Bus.VMagAngle()[1] / base
    end
    return conv, imp, v
end

function check(name)
    c = convert(name)
    dir = get(ENV, "CONVERT_OUT", "") == "" ? joinpath(REPO, "test_networks", deckname(c)) : joinpath(ENV["CONVERT_OUT"], deckname(c))
    conv, dimp, dv = deck_flat(c, dir)
    oimp, ov = original_flat(c)
    rows = String[]
    say(s) = (println(s); push!(rows, s))
    vb = [dv[b] for b in c.buses if haskey(dv, b)]
    say(@sprintf("%s: deck %s; voltage %.4f-%.4f pu at nominal load, no DERs (original, flat: %.4f-%.4f)", name, conv ? "converged" : "DID NOT CONVERGE", minimum(vb), maximum(vb), minimum(values(ov)), maximum(values(ov))))
    for (i, s) in enumerate(c.subs)
        o = sum((p for (b, p) in oimp if b in s.mv); init = 0.0)
        h = s.mv[1]                                  # the distribution-side bus of its transformer
        bs = [b for (b, f) in c.feeder if f == i && haskey(ov, b) && haskey(dv, b)]
        d = [(dv[h] - dv[b]) - (ov[h] - ov[b]) for b in bs]
        say(@sprintf("  subs%d: import original %.1f kW, deck %.1f kW (%+.1f%%); largest drop from the head original %.4f, deck %.4f pu; per bus deck - original mean %+.4f, max |.| %.4f pu over %d buses",
                     i, o, dimp[i], 100 * (dimp[i] - o) / o, maximum(ov[h] - ov[b] for b in bs), maximum(dv[h] - dv[b] for b in bs),
                     sum(d) / length(d), maximum(abs, d), length(bs)))
    end
    open(io -> foreach(r -> println(io, r), rows), joinpath(dir, "check.txt"), "w")
end

for name in (isempty(ARGS) ? sort(collect(keys(SYSTEMS))) : ARGS)
    check(name)
end
