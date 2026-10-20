# Convert a published OpenDSS distribution system into this repo's single-phase deck format
# (rawData/<name>_1ph): one voltage level, 1-phase Lines with series R + jX only, constant-PQ
# loads, 1-phase Vsources behind each substation transformer, our own PV and batteries. The
# full-angle MSOPF model reads only those, so everything else is turned into them:
#   julia --project=envs/msopf envs/multi_poi/root_level/convert_to_1ph.jl <system ...>
# Systems are listed in SYSTEMS below; the originals stay unmodified (rawData/ieee9500_source,
# test_systems/). Steps:
#   * Compile the original from its master, skipping commands that run studies or open windows;
#     built-in generators, PV and storage are left out (we place our own).
#   * Kept network = every bus between 1 and 35 kV line-to-line (the primary, "MV"). Above it
#     (sub-transmission, transmission) is not modeled; below it (customer secondaries) each load is
#     moved to the MV bus of its service transformer, whose losses are dropped.
#   * One voltage base for the whole deck, the MV level the substations feed (`base_kV`). An
#     impedance on another MV level is referred to it by (V_base / V_level)^2, so per-unit
#     voltages and powers are unchanged.
#   * Per-phase equivalent: loads and capacitors / 3; a k-phase section gets (3/k) x its
#     positive-sequence impedance (k = 3: mean self minus mean mutual; k = 1: self), which keeps
#     the voltage drop of a section carrying S on its k phases equal to the original's.
#   * Closed switches and voltage regulators (transformers with both windings at one voltage) are
#     not in the model: each becomes a short line of 1e-3 + j1e-3 ohm between its two buses, kept
#     apart (user, 2026-10-20) and grouped under a comment in BranchData.dss; regulator taps are 1.
#     A transformer between two MV levels becomes a line with its leakage impedance (%R of both
#     windings + j%X on its rating), referred to the base. A series reactor becomes a line.
#   * Open switches become lines with enabled=false (Ties.dss), closed by a study.
#   * Each transformer from above MV into the MV network is a substation transformer: a line from
#     its high-voltage bus with its leakage impedance in ohm at the distribution-side voltage
#     (parallel units combined). Each distinct high-voltage bus carries one stiff 1-phase Vsource
#     at 1.05 pu (user, 2026-10-20); nothing upstream of it is modeled.
#   * Loads that stand for other feeders lumped at the substation bus (`drop_loads`) are left
#     out: they are not part of the modeled network and would set the substation's import.
#   * Buses no source reaches after all this are dropped, with their loads, and reported.
#   * DERs (user's rule, 2026-10-20): about 10% of the buses get a battery and about 10% a PV unit,
#     one third of each set shared, on buses with load, sized in proportion to the bus load so
#     each fleet is 40% of the total load. PV kVA = 1.2 Pmpp; batteries 4 h, kVA = 1.2 kWrated,
#     62.5% stored, 30% reserve (as ieee2522C_1ph). Fixed seed per system.
using OpenDSSDirect, Printf, Random, LinearAlgebra
const ODD = OpenDSSDirect

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TS = joinpath(REPO, "test_systems")
const MV = (1.0, 35.0)                      # kV line-to-line kept as the primary network
const DER_SHARE, DER_FRACTION, OVERLAP = 0.40, 0.10, 1 / 3
const ZSHORT = complex(1e-3, 1e-3)           # ohm: a closed switch or a regulator, which the model does not represent
const ZFLOOR = 1e-4                           # ohm: a line below this gets ZSHORT

# Commands that run a study or open a window, skipped in every master.
const SKIP = r"(?i)^\s*(solve|plot|show|export|visualize|buscoords|latlongcoords|calcv\w*\s*$|set\s+mode|set\s+mark|summary|interpolate)\b"

const SYSTEMS = Dict(
    "ieee9500" => (; master = joinpath(REPO, "rawData", "ieee9500_source", "base", "Master-bal-initial-config.dss"), base_kV = 12.47, seed = 9500,
                    skip = r"(?i)^\s*(redirect\s+(generators|energystorage|pv_)|setkvbase|generator\.|new energymeter)", drop_loads = nothing),
    "ieee8500" => (; master = joinpath(TS, "ieee8500", "Master.dss"), base_kV = 12.47, seed = 8500, skip = nothing, drop_loads = nothing),
    "epri_ckt5" => (; master = joinpath(TS, "epri_test_circuits", "ckt5", "Master_ckt5.dss"), base_kV = 12.47, seed = 5, skip = nothing, drop_loads = nothing),
    "epri_ckt7" => (; master = joinpath(TS, "epri_test_circuits", "ckt7", "Master_ckt7.dss"), base_kV = 12.47, seed = 7, skip = nothing,
                     drop_loads = r"^256\d\d_[abc]$"),                       # 13 feeders lumped on bus ckt7 (Substation_ckt7.dss)
    "epri_ckt24" => (; master = joinpath(TS, "epri_test_circuits", "ckt24", "master_ckt24.dss"), base_kV = 34.5, seed = 24, skip = nothing,
                      drop_loads = r"^other_feeders$"),
    "epri_j1" => (; master = joinpath(TS, "epri_test_circuits", "epri_dpv", "J1", "Master_noPV.dss"), base_kV = 12.47, seed = 101, skip = nothing,
                   drop_loads = r"^aggregate_load$"),
    "epri_k1" => (; master = joinpath(TS, "epri_test_circuits", "epri_dpv", "K1", "Master_NoPV.dss"), base_kV = 13.2, seed = 111, skip = nothing,
                   drop_loads = r"^k_22_aggregate$"),
    "epri_m1" => (; master = joinpath(TS, "epri_test_circuits", "epri_dpv", "M1", "Master_NoPV.dss"), base_kV = 12.47, seed = 131, skip = nothing,
                   drop_loads = r"^aggregated$"),
    "smartds_p1u" => (; master = joinpath(TS, "smartds_sfo_p1u", "Master.dss"), base_kV = 12.47, seed = 2018, skip = nothing, drop_loads = nothing,
                       flat_shapes = joinpath(TS, "smartds_sfo_p1u", "LoadShapes.dss")))

"Compile the original. SMART-DS load shapes read annual profiles not downloaded: flat 1.0 in memory."
function compile_original(cfg)
    cmds = String[]
    for l in readlines(cfg.master)
        l = strip(split(split(l, '!')[1], "//")[1]); isempty(l) && continue
        startswith(l, "~") ? (cmds[end] *= " " * l[2:end]) : push!(cmds, l)
    end
    flat = haskey(cfg, :flat_shapes) ?
        ["New Loadshape.$(m[1]) npts=1 interval=1 mult=(1) qmult=(1)" for m in eachmatch(r"(?i)New\s+Loadshape\.(\S+)", read(cfg.flat_shapes, String))] : nothing
    ODD.Text.Command("Clear"); ODD.Text.Command("cd \"$(dirname(cfg.master))\"")
    for c in cmds
        (occursin(SKIP, c) || (cfg.skip !== nothing && occursin(cfg.skip, c))) && continue
        if flat !== nothing && occursin(r"(?i)^redirect\s+loadshapes\.dss", c)
            foreach(ODD.Text.Command, flat)
        else
            ODD.Text.Command(c)
        end
    end
    for cls in ("Generator", "PVSystem", "Storage")
        any(startswith(lowercase(e), lowercase(cls) * ".") for e in ODD.Circuit.AllElementNames()) &&
            ODD.Text.Command("BatchEdit $cls..* enabled=false")
    end
    ODD.Text.Command("Calcvoltagebases"); ODD.Text.Command("Set mode=snapshot")
    ODD.Text.Command("Set maxiterations=100"); ODD.Text.Command("Set maxcontroliter=200")
    ODD.Solution.Solve()
end

bus(ref) = lowercase(String(split(ref, '.')[1]))
kvll(b) = (ODD.Circuit.SetActiveBus(String(b)); ODD.Bus.kVBase() * sqrt(3))
ismv(k) = MV[1] < k <= MV[2]
clean(s) = replace(lowercase(s), r"[^a-z0-9_\-]" => "_")
real_names(v) = filter(n -> lowercase(n) != "none", v)

struct DSU; p::Dict{String,String}; end
DSU() = DSU(Dict{String,String}())
root(d::DSU, x) = (r = get!(d.p, x, x); r == x ? x : (d.p[x] = root(d, r)))
function join!(d::DSU, a, b)
    ra, rb = root(d, a), root(d, b)
    ra != rb && (d.p[rb] = ra)
    return nothing
end

"Per-phase-equivalent series impedance (ohm, at the element's own voltage) of the active 2-terminal element, from YPrim."
function series_z()
    n = ODD.CktElement.NumConductors(); k = ODD.CktElement.NumPhases()
    Y = ODD.CktElement.YPrim(); Y = Y isa AbstractMatrix ? Y : reshape(Y, 2n, 2n)
    Z = inv(-Y[1:n, n+1:2n])
    n > k && (Z = Z[1:k, 1:k] - Z[1:k, k+1:n] * (Z[k+1:n, k+1:n] \ Z[k+1:n, 1:k]))
    z1 = k == 1 ? Z[1, 1] : sum(diag(Z)) / k - (sum(Z) - sum(diag(Z))) / (k * (k - 1))
    return 3 / k * z1
end

"Leakage impedance of the active 2-winding transformer, per-phase equivalent, in ohm at line-to-line voltage `kv_ref`."
function transformer_z(nm, kv_ref)
    ODD.Transformers.Name(nm)
    ODD.Transformers.Wdg(1.0); r1 = ODD.Transformers.R(); kva = ODD.Transformers.kVA()
    ODD.Transformers.Wdg(2.0); r2 = ODD.Transformers.R()
    zpu = complex(r1 + r2, ODD.Transformers.Xhl()) / 100
    k = ODD.CktElement.NumPhases()
    return 3 / k * zpu * (kv_ref / sqrt(3))^2 * 1000 / (kva / k)
end

function convert(name)
    cfg = SYSTEMS[name]
    compile_original(cfg)
    vref = cfg.base_kV
    kv = Dict(lowercase(b) => kvll(b) for b in ODD.Circuit.AllBusNames())
    mv = Set(b for (b, k) in kv if ismv(k))
    refer(z, b) = z * (vref / kv[b])^2                     # ohm at bus b's level -> ohm at the base
    dsu = DSU(); foreach(b -> root(dsu, b), mv)
    branches, ties, notes = [], [], String[]
    nsw = 0; nreg = 0; nstep = 0; nreac = 0
    for nm in real_names(ODD.Lines.AllNames())
        ODD.Lines.Name(nm); el = "Line.$nm"; ODD.Circuit.SetActiveElement(el)
        ODD.CktElement.Enabled() || continue
        a, b = bus(ODD.Lines.Bus1()), bus(ODD.Lines.Bus2())
        (a in mv && b in mv) || continue
        if ODD.CktElement.IsOpen(1, 0) || ODD.CktElement.IsOpen(2, 0)
            push!(ties, (; name = clean(nm), a, b))
        elseif ODD.Lines.IsSwitch()
            push!(branches, (; name = clean(nm), a, b, z = ZSHORT, kind = :switch)); nsw += 1
        else
            push!(branches, (; name = clean(nm), a, b, z = refer(series_z(), a), kind = :line))
        end
    end
    for nm in real_names(ODD.Reactors.AllNames())
        ODD.Circuit.SetActiveElement("Reactor.$nm"); refs = ODD.CktElement.BusNames()
        length(refs) == 2 && all(r -> bus(r) in mv, refs) && bus(refs[1]) != bus(refs[2]) || continue
        push!(branches, (; name = clean("reactor_" * nm), a = bus(refs[1]), b = bus(refs[2]), z = refer(series_z(), bus(refs[1])), kind = :reactor)); nreac += 1
    end
    lv = DSU()
    for nm in real_names(ODD.Lines.AllNames())
        ODD.Lines.Name(nm); a, b = bus(ODD.Lines.Bus1()), bus(ODD.Lines.Bus2())
        (a in mv || b in mv || !(haskey(kv, a) && haskey(kv, b))) || join!(lv, a, b)
    end
    lvtop = Dict{String,String}(); subx = []
    for nm in real_names(ODD.Transformers.AllNames())
        ODD.Circuit.SetActiveElement("Transformer.$nm"); ODD.CktElement.Enabled() || continue
        bs = bus.(ODD.CktElement.BusNames()); ks = [kv[b] for b in bs]
        m = [b for b in bs if b in mv]
        if length(m) == length(bs) && length(bs) == 2
            if abs(ks[1] - ks[2]) / ks[1] < 0.02          # regulator: merge, tap 1
                push!(branches, (; name = clean("reg_" * nm), a = bs[1], b = bs[2], z = ZSHORT, kind = :regulator)); nreg += 1
            else                                            # step between MV levels: a line
                push!(branches, (; name = clean("xfmr_" * nm), a = bs[1], b = bs[2], z = transformer_z(nm, vref), kind = :transformer)); nstep += 1
            end
        elseif length(m) == 1 && any(>(MV[2]), ks)          # substation: HV into MV
            length(bs) == 2 || push!(notes, "substation transformer $nm has $(length(bs)) windings; windings 1-2 used")
            hv = bs[findfirst(b -> kv[b] > MV[2], bs)]
            push!(subx, (; name = clean("xfmr_" * nm), a = hv, b = m[1], z = transformer_z(nm, vref), kind = :substation))
        elseif length(m) == 1                               # service: MV into secondaries
            foreach(x -> x in mv || (lvtop[root(lv, x)] = m[1]), bs)
        end
    end
    append!(branches, subx)                                 # substation transformers: lines from the source bus
    load = Dict{String,ComplexF64}(); dropped = 0.0; nload = 0
    for nm in real_names(ODD.Loads.AllNames())
        ODD.Circuit.SetActiveElement("Load.$nm"); ODD.CktElement.Enabled() || continue
        ODD.Loads.Name(nm); s = complex(ODD.Loads.kW(), ODD.Loads.kvar())
        if cfg.drop_loads !== nothing && occursin(cfg.drop_loads, lowercase(nm))
            dropped += real(s); continue
        end
        b = bus(ODD.CktElement.BusNames()[1])
        pb = b in mv ? b : get(lvtop, root(lv, b), nothing)
        pb === nothing && error("load $nm at $b ($(round(get(kv, b, NaN); digits = 3)) kV): no service transformer found")
        r = root(dsu, pb); load[r] = get(load, r, 0im) + s / 3; nload += 1
    end
    caps = Dict{String,Float64}()
    for nm in real_names(ODD.Capacitors.AllNames())
        ODD.Circuit.SetActiveElement("Capacitor.$nm"); b = bus(ODD.CktElement.BusNames()[1])
        b in mv || continue
        ODD.Capacitors.Name(nm); r = root(dsu, b); caps[r] = get(caps, r, 0.0) + ODD.Capacitors.kvar() / 3
    end
    # merged branches: drop self-loops, combine parallels
    par = Dict{Tuple{String,String},Vector{Any}}()
    for l in branches
        a, b = root(dsu, l.a), root(dsu, l.b); a == b && continue
        push!(get!(par, a < b ? (a, b) : (b, a), Any[]), l)
    end
    br = sort([(; name = ls[1].name, a = k[1], b = k[2], z = 1 / sum(1 / l.z for l in ls), n = length(ls), kind = ls[1].kind) for (k, ls) in par]; by = x -> x.name)
    # Lines of almost no impedance (bus bars, jumpers; 1e-8 ohm on some decks) get the short-line value.
    br = [abs(x.z) < ZFLOOR && x.kind in (:line, :reactor) ? merge(x, (; z = ZSHORT, kind = :short)) : x for x in br]
    # reach from the sources; what is not reached is dropped
    adj = Dict{String,Vector{String}}()
    for x in br
        push!(get!(adj, x.a, String[]), x.b); push!(get!(adj, x.b, String[]), x.a)
    end
    srcs = sort(unique([x.a for x in subx]))                # one source per high-voltage bus
    feeder = Dict{String,Int}()
    for (i, s) in enumerate(srcs)
        q = [s]; feeder[s] = i
        while !isempty(q)
            x = popfirst!(q)
            for y in get(adj, x, String[])
                haskey(feeder, y) || (feeder[y] = i; push!(q, y))
            end
        end
    end
    allroots = Set(root(dsu, b) for b in mv)
    lost = setdiff(allroots, keys(feeder))
    lostload = sum((real(get(load, b, 0im)) for b in lost); init = 0.0)
    br = filter(x -> haskey(feeder, x.a) && haskey(feeder, x.b), br)
    load = Dict(b => s for (b, s) in load if haskey(feeder, b))
    caps = Dict(b => q for (b, q) in caps if haskey(feeder, b))
    ties = [(; t..., a = root(dsu, t.a), b = root(dsu, t.b)) for t in ties if haskey(feeder, root(dsu, t.a)) && haskey(feeder, root(dsu, t.b))]
    stats = (; nsw, nreg, nstep, nreac, nload, dropped, lost = length(lost), lostload, levels = sort(unique(round.([kv[b] for b in mv]; digits = 2))))
    return (; name, cfg, buses = Set(keys(feeder)), branches = br, ties, load, caps,
            subs = [(; bus = s, mv = sort(unique([x.b for x in subx if x.a == s])), z = 1 / sum(1 / x.z for x in subx if x.a == s)) for s in srcs], feeder, stats, notes)
end

function place_ders(c)
    rng = MersenneTwister(c.cfg.seed)
    subset = Set([[s.bus for s in c.subs]; reduce(vcat, [s.mv for s in c.subs])])
    cands = sort([b for (b, s) in c.load if real(s) > 0 && !(b in subset)])
    n = round(Int, DER_FRACTION * (length(c.buses) - length(c.subs))); nboth = round(Int, OVERLAP * n)   # buses below the sources
    n = min(n, (length(cands) + nboth) ÷ 2)
    pick = shuffle(rng, cands)[1:(2n - nboth)]
    pvb, btb = sort([pick[1:nboth]; pick[nboth+1:n]]), sort([pick[1:nboth]; pick[n+1:end]])
    total = sum(real, values(c.load))
    size(bs) = (w = sum(real(c.load[b]) for b in bs); Dict(b => DER_SHARE * total * real(c.load[b]) / w for b in bs))
    return size(pvb), size(btb), nboth
end

f3(x) = @sprintf("%.6g", x)
deckname(c) = c.name * "_1ph"

function write_deck(c, pv, bt; out = joinpath(REPO, "rawData", deckname(c)))
    mkpath(out)
    kvln = c.cfg.base_kV / sqrt(3); kvs = f3(round(kvln; digits = 4))
    src = relpath(c.cfg.master, REPO)
    hdr = "! $(deckname(c)): single-phase per-phase equivalent of $(src) (unmodified original).\n" *
          "! Written by envs/multi_poi/root_level/convert_to_1ph.jl -- do not edit by hand.\n"
    w(body, f) = open(io -> (print(io, hdr); body(io)), joinpath(out, f), "w")
    w("BranchData.dss") do io
        println(io, "! Per-phase equivalent, referred to $(c.cfg.base_kV) kV: a k-phase section has (3/k) x its positive-sequence")
        println(io, "! impedance. Series R + jX only (c1=0 c0=0). Grouped by what each line is in the original.")
        heads = [:substation => "! Substation transformers: a line with the leakage impedance (%R of both windings + j%X on its rating) in ohm at the distribution-side voltage, from the source bus (parallel units combined).",
                 :line => "! Lines of the original.",
                 :transformer => "! Transformers between MV levels in the original: a line with the leakage impedance (%R of both windings + j%X).",
                 :reactor => "! Series reactors of the original.",
                 :switch => "! Closed switches in the original. A switch is not in the model: a short line of 1e-3 + j1e-3 ohm.",
                 :regulator => "! Voltage regulators in the original. A regulator is not in the model (tap 1): a short line of 1e-3 + j1e-3 ohm.",
                 :short => "! Lines of the original below 1e-4 ohm (bus bars, jumpers): 1e-3 + j1e-3 ohm, to keep the admittance matrix well scaled."]
        for (kind, head) in heads
            sel = filter(x -> x.kind == kind, c.branches); isempty(sel) && continue
            println(io, head)
            for x in sel
                println(io, "New Line.$(x.name) Phases=1 Bus1=$(x.a).1 Bus2=$(x.b).1 r1=$(f3(real(x.z))) x1=$(f3(imag(x.z))) c1=0 c0=0 length=1 units=none")
            end
        end
    end
    w("Ties.dss") do io
        println(io, "! The original's normally-open switches, disabled. Close one with Line.<name>.enabled=true; each is given")
        println(io, "! 1e-3 + j1e-3 ohm, as the other decks' shortest lines, to keep the admittance matrix well scaled.")
        for t in c.ties
            println(io, "New Line.$(t.name) Phases=1 Bus1=$(t.a).1 Bus2=$(t.b).1 r1=0.001 x1=0.001 c1=0 c0=0 length=1 units=none enabled=false")
        end
    end
    w("Loads.dss") do io
        println(io, "! Each customer load summed onto the MV bus of its service transformer, divided by 3.")
        for b in sort(collect(keys(c.load)))
            s = c.load[b]; real(s) > 0 || continue
            println(io, "New load.L_$b Bus=$b.1 Phases=1 Model=1 kV=$kvs kW=$(f3(real(s))) kVAr=$(f3(imag(s))) Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault")
        end
    end
    w("Capacitor.dss") do io
        println(io, "! The original's capacitor banks (kvar / 3), not redirected by Master.dss.")
        for b in sort(collect(keys(c.caps)))
            println(io, "New generator.$(b)C Phases=1 Bus1=$b.1 kV=$kvs Model=1 kW=1e-06 kVAr=$(f3(c.caps[b]))")
        end
    end
    w("PVSystem.dss") do io
        @printf(io, "! %d PV units on load buses, sized in proportion to the bus load; fleet Pmpp = %.0f %% of total load.\n", length(pv), 100DER_SHARE)
        for b in sort(collect(keys(pv)))
            println(io, "New PVsystem.PV_$b irrad=1.0 Phases=1 Bus1=$b.1 kV=$kvs kVA=$(f3(1.2pv[b])) Pmpp=$(f3(pv[b])) %cutin=0.001 %cutout=0.001")
        end
    end
    w("Storage.dss") do io
        @printf(io, "! %d batteries on load buses, sized in proportion to the bus load; fleet rating = %.0f %% of total load, 4 h.\n", length(bt), 100DER_SHARE)
        for b in sort(collect(keys(bt)))
            println(io, "New Storage.BESS_$b Phases=1 Bus1=$b.1 kV=$kvs kVA=$(f3(1.2bt[b])) kWrated=$(f3(bt[b])) kWhrated=$(f3(4bt[b])) %stored=62.5 %reserve=30 %EffCharge=95 %EffDischarge=95 %idlingkW=0 Vminpu=0.95 Vmaxpu=1.05 DispMode=External")
        end
    end
    open(joinpath(out, "Master.dss"), "w") do io
        print(io, "Clear\n", hdr)
        println(io, "New Circuit.$(deckname(c)) basekv=$kvs Bus1=$(c.subs[1].bus) pu=1.05 phases=1 R1=0 X1=0.0001 R0=0 X0=0.0001")
        println(io, "VSource.source.enabled=false  //Disable the default VSource created when declaring the circuit.")
        println(io, "! One stiff source per substation, on the high-voltage bus of its transformer(s) (original bus name, carried at the")
        println(io, "! deck's $(c.cfg.base_kV) kV base). The transformers are lines in BranchData.dss; nothing upstream is modeled.")
        for (i, s) in enumerate(c.subs)
            println(io, "New Vsource.subs$i Phases=1 Bus1=$(s.bus) basekv=$kvs pu=1.05 angle=0.0 R1=0 X1=0.00000001 R0=0 X0=0.00000001")
        end
        println(io, "Redirect BranchData.dss\nRedirect Ties.dss")
        print(io, """
Redirect ../LoadShapeDefault.dss
Redirect ../LoadShapePVDefault.dss
Redirect Loads.dss     ! Balanced Loads
// Redirect  Capacitor.dss
Redirect PVSystem.dss
Redirect Storage.dss

Set VoltageBases = [$(c.cfg.base_kV)]
Set mode = Daily
Set stepsize = 1h
Set number = 1

CalcVoltageBases
""")
    end
    return out
end

function report(c, pv, bt, nboth)
    s = c.stats
    @printf("%s: MV levels %s kV -> %d buses, %d lines (%d transformer steps, %d reactors, %d closed switches and %d regulators as short lines); %d ties\n",
            c.name, join(s.levels, "/"), length(c.buses), length(c.branches), s.nstep, s.nreac, s.nsw, s.nreg, length(c.ties))
    tot = sum(values(c.load))
    @printf("  load (per phase) %.1f kW %.1f kvar on %d buses; by substation %s kW; dropped as other feeders %.0f kW (3-phase)\n",
            real(tot), imag(tot), length(c.load),
            join([@sprintf("%.0f", sum((real(v) for (b, v) in c.load if c.feeder[b] == i); init = 0.0)) for i in eachindex(c.subs)], " / "), s.dropped)
    for (i, x) in enumerate(c.subs)
        @printf("  subs%d at %s, feeding %s: transformer %.4f + j%.4f ohm\n", i, x.bus, join(x.mv, ", "), real(x.z), imag(x.z))
    end
    s.lost > 0 && @printf("  ** %d buses no source reaches were dropped, with %.1f kW of load (per phase)\n", s.lost, s.lostload)
    zs = sort(abs.(getfield.(c.branches, :z)))
    @printf("  |z| per line: min %.2e, median %.2e, max %.2e ohm\n", zs[1], zs[end÷2], zs[end])
    @printf("  DERs: %d PV, %d batteries, %d buses with both; each fleet %.0f kW\n", length(pv), length(bt), nboth, sum(values(pv)))
    foreach(n -> println("  note: ", n), c.notes)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    for name in (isempty(ARGS) ? sort(collect(keys(SYSTEMS))) : ARGS)
        c = convert(name); pv, bt, nboth = place_ders(c)
        report(c, pv, bt, nboth)
        println("  wrote ", write_deck(c, pv, bt; out = get(ENV, "CONVERT_OUT", "") == "" ? joinpath(REPO, "rawData", deckname(c)) : joinpath(ENV["CONVERT_OUT"], deckname(c))))
    end
end
