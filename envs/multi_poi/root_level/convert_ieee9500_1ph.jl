# Convert the IEEE 9500-node test system (PNNL / GridAPPS-D CIMHub, branch feature/SETO,
# ieee9500/base; kept unmodified in rawData/ieee9500_source) into rawData/ieee9500_1ph, a
# single-phase per-phase equivalent of its 12.47 kV primary network in the format of the other
# decks, with this repo's own PV and batteries. Run from the repo root:
#   julia --project=envs/msopf envs/multi_poi/root_level/convert_ieee9500_1ph.jl
#
# The original has three feeders, each fed by its own 69/12.47 kV substation, joined by
# normally-open switches. Conversion, in order:
#   * Solve the original (balanced loads, its normal switch configuration, built-in DERs left
#     out: we place our own) for the operating point the checks compare against.
#   * Primary network = every bus at 12.47 kV (7.2 kV line-to-neutral). Customer loads sit
#     behind service transformers and triplex secondaries; each is moved to the primary bus
#     of its service transformer. Secondary and service-transformer losses are dropped.
#   * Per-phase equivalent: the three-phase system carries three times what the one-phase
#     model carries, so every load and capacitor is divided by 3. A k-phase line section gets
#     impedance (3/k) * z, with z its positive-sequence impedance (k = 3: mean self minus
#     mean mutual; k = 2 the same over two phases; k = 1 its self impedance). The factor 3/k
#     keeps the voltage drop of a section carrying S, spread over its k phases, equal to the
#     original's, and its loss at one third of the original's.
#   * Closed switches, and the 1-phase regulators (feeder heads and in-line), are merged into
#     the buses they join. Regulator taps are therefore 1, and capacitors are written to
#     Capacitor.dss but not redirected (as in the other decks): the full-angle model reads
#     lines, loads, sources and our DERs only.
#   * The normally-open switches are written as lines with enabled=false (Ties.dss), so a
#     study closes a tie with `Line.<name>.enabled=true`.
#   * Each substation becomes a 1-phase Vsource behind its 69/12.47 kV transformer impedance
#     (20 MVA, X = 7.94 %, R = 2 x 0.67 %: Zs = 0.104 + j0.617 ohm at 12.47 kV), at 1.05 pu
#     like the other decks; the 69 kV and 115 kV networks upstream are not modeled.
#   * DERs (the user's rule, 2026-10-20): about 10 % of the primary buses get a battery and
#     about 10 % a PV unit, partly overlapping (one third of each set shared); units go on
#     buses with load, sized in proportion to that bus's load so that each fleet's rating is
#     40 % of the total load (inside the 30-50 % the user uses). PV kVA = 1.2 Pmpp; batteries
#     4 h, kVA = 1.2 kWrated, 62.5 % stored, 30 % reserve, as in ieee2522C_1ph. Fixed seed.
using OpenDSSDirect, Printf, Random, LinearAlgebra
const ODD = OpenDSSDirect

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SRC = joinpath(REPO, "rawData", "ieee9500_source", "base")
const OUT = joinpath(REPO, "rawData", "ieee9500_1ph")
const KV_LN = 7.2
const ZS = complex(0.0134, 0.0794) * 12.47^2 / 20     # substation transformer, ohm at 12.47 kV
const DER_SHARE, DER_FRACTION, OVERLAP, SEED = 0.40, 0.10, 1 / 3, 9500

primary(b) = (ODD.Circuit.SetActiveBus(String(b)); abs(ODD.Bus.kVBase() - KV_LN) < 0.2)
bus(ref) = lowercase(String(split(ref, '.')[1]))
nodes(ref) = (p = split(ref, '.')[2:end]; isempty(p) ? [1, 2, 3] : [parse(Int, n) for n in p])

"Solve the original without its built-in DERs, run line by line from its folder."
function solve_original()
    skip = r"(?i)^\s*(redirect\s+(generators|energystorage|pv_)|setkvbase|generator\.|//|!|$|latlongcoords|buscoords|set\s+mark|new energymeter)"
    cmds = String[]
    for l in readlines(joinpath(SRC, "Master-bal-initial-config.dss"))
        occursin(skip, l) && continue
        l = strip(split(l, '!')[1]); isempty(l) && continue
        startswith(l, "~") ? (cmds[end] *= " " * l[2:end]) : push!(cmds, l)
    end
    ODD.Text.Command("Clear"); ODD.Text.Command("cd \"$SRC\"")
    foreach(ODD.Text.Command, cmds)
    ODD.Text.Command("Solve")
    ODD.Solution.Converged() || error("the original did not converge")
end

"Union-find over bus names."
struct DSU; p::Dict{String,String}; end
DSU() = DSU(Dict{String,String}())
root(d::DSU, x) = (r = get!(d.p, x, x); r == x ? x : (d.p[x] = root(d, r)))
function join!(d::DSU, a, b)
    ra, rb = root(d, a), root(d, b)
    ra != rb && (d.p[rb] = ra)
    return nothing
end

"Series impedance of a line in ohms, per phase conductor, from YPrim (Kron-reducing extra conductors)."
function line_z()
    n = ODD.CktElement.NumConductors(); k = ODD.CktElement.NumPhases()
    Y = ODD.CktElement.YPrim()
    Y = Y isa AbstractMatrix ? Y : reshape(Y, 2n, 2n)
    Z = inv(-Y[1:n, n+1:2n])
    if n > k                                  # neutral grounded at both ends
        Z = Z[1:k, 1:k] - Z[1:k, k+1:n] * (Z[k+1:n, k+1:n] \ Z[k+1:n, 1:k])
    end
    z1 = k == 1 ? Z[1, 1] : sum(diag(Z)) / k - (sum(Z) - sum(diag(Z))) / (k * (k - 1))
    return 3 / k * z1, k
end

function convert()
    solve_original()
    allbus = lowercase.(ODD.Circuit.AllBusNames())
    prim = Set(b for b in allbus if primary(b))
    # 1. merges: closed switches and regulators between primary buses
    dsu = DSU(); foreach(b -> root(dsu, b), prim)
    lines = []; ties = []
    for nm in ODD.Lines.AllNames()
        ODD.Lines.Name(nm); el = "Line.$nm"; ODD.Circuit.SetActiveElement(el)
        ODD.CktElement.Enabled() || continue
        b1, b2 = ODD.Lines.Bus1(), ODD.Lines.Bus2()
        (bus(b1) in prim && bus(b2) in prim) || continue
        isopen = ODD.CktElement.IsOpen(1, 0) || ODD.CktElement.IsOpen(2, 0)
        if isopen
            push!(ties, (; name = lowercase(nm), a = bus(b1), b = bus(b2), k = ODD.Lines.Phases()))
        elseif ODD.Lines.IsSwitch()
            join!(dsu, bus(b1), bus(b2))
        else
            z, k = line_z(); push!(lines, (; name = lowercase(nm), a = bus(b1), b = bus(b2), z, k))
        end
    end
    nreg = 0
    for nm in ODD.Transformers.AllNames()
        ODD.Circuit.SetActiveElement("Transformer.$nm"); refs = ODD.CktElement.BusNames()
        all(r -> bus(r) in prim, refs) || continue
        join!(dsu, bus(refs[1]), bus(refs[2])); nreg += 1
    end
    # 2. loads: each to the primary bus of its service transformer
    lv = DSU()
    for nm in ODD.Lines.AllNames()
        ODD.Lines.Name(nm); a, b = bus(ODD.Lines.Bus1()), bus(ODD.Lines.Bus2())
        (a in prim || b in prim) || join!(lv, a, b)
    end
    lvtop = Dict{String,String}()
    for nm in ODD.Transformers.AllNames()
        ODD.Circuit.SetActiveElement("Transformer.$nm"); bs = bus.(ODD.CktElement.BusNames())
        p = filter(in(prim), bs); s = filter(!in(prim), bs)
        length(p) == 1 && !isempty(s) && foreach(x -> (lvtop[root(lv, x)] = p[1]), s)
    end
    load = Dict{String,ComplexF64}()
    for nm in ODD.Loads.AllNames()
        ODD.Loads.Name(nm); ODD.Circuit.SetActiveElement("Load.$nm"); b = bus(ODD.CktElement.BusNames()[1])
        pb = b in prim ? b : get(lvtop, root(lv, b), nothing)
        pb === nothing && error("load $nm at $b: no service transformer found")
        load[root(dsu, pb)] = get(load, root(dsu, pb), 0im) + complex(ODD.Loads.kW(), ODD.Loads.kvar()) / 3
    end
    caps = Dict{String,Float64}()
    for nm in ODD.Capacitors.AllNames()
        ODD.Capacitors.Name(nm); ODD.Circuit.SetActiveElement("Capacitor.$nm")
        b = root(dsu, bus(ODD.CktElement.BusNames()[1]))
        caps[b] = get(caps, b, 0.0) + ODD.Capacitors.kvar() / 3
    end
    # 3. substations: the low-side bus of each 69/12.47 kV transformer
    subs = String[]
    for i in 1:3
        ODD.Circuit.SetActiveElement("Transformer.HVMV69_11Sub$i")
        push!(subs, root(dsu, bus(ODD.CktElement.BusNames()[2])))
    end
    # merged lines: drop self-loops, combine parallels
    par = Dict{Tuple{String,String},Vector{Any}}()
    for l in lines
        a, b = root(dsu, l.a), root(dsu, l.b); a == b && continue
        push!(get!(par, a < b ? (a, b) : (b, a), Any[]), l)
    end
    branches = [(; name = ls[1].name, a = key[1], b = key[2], z = 1 / sum(1 / l.z for l in ls), n = length(ls))
                for (key, ls) in par]
    sort!(branches; by = x -> x.name)
    buses = Set(root(dsu, b) for b in prim)
    return (; buses, branches, ties = [merge(t, (; a = root(dsu, t.a), b = root(dsu, t.b))) for t in ties],
            load, caps, subs, nreg, nlines = length(lines))
end

"Feeder of every bus: breadth-first from each substation over the branches (ties open)."
function feeders(c)
    adj = Dict{String,Vector{String}}()
    for br in c.branches
        push!(get!(adj, br.a, String[]), br.b); push!(get!(adj, br.b, String[]), br.a)
    end
    f = Dict{String,Int}()
    for (i, s) in enumerate(c.subs)
        q = [s]; f[s] = i
        while !isempty(q)
            x = popfirst!(q)
            for y in get(adj, x, String[])
                haskey(f, y) || (f[y] = i; push!(q, y))
            end
        end
    end
    return f
end

function place_ders(c)
    rng = MersenneTwister(SEED)
    cands = sort([b for (b, s) in c.load if real(s) > 0 && !(b in c.subs)])
    n = round(Int, DER_FRACTION * length(c.buses)); nboth = round(Int, OVERLAP * n)
    pick = shuffle(rng, cands)[1:(2n - nboth)]
    both, ponly, bonly = pick[1:nboth], pick[nboth+1:n], pick[n+1:end]
    pvb, btb = sort([both; ponly]), sort([both; bonly])
    total = sum(real, values(c.load))
    size(bs) = (w = sum(real(c.load[b]) for b in bs); Dict(b => DER_SHARE * total * real(c.load[b]) / w for b in bs))
    return size(pvb), size(btb), length(both)
end

f3(x) = @sprintf("%.6g", x)

function write_deck(c, pv, bt)
    mkpath(OUT)
    hdr = "! ieee9500_1ph: single-phase per-phase equivalent of the IEEE 9500-node test system\n" *
          "! (PNNL / GridAPPS-D CIMHub; original in ../ieee9500_source, BSD-style license there).\n" *
          "! Written by envs/multi_poi/root_level/convert_ieee9500_1ph.jl -- do not edit by hand.\n"
    open(joinpath(OUT, "BranchData.dss"), "w") do io
        print(io, hdr, "! Per-phase equivalent: a k-phase section has (3/k) x its positive-sequence impedance.\n",
              "! Lines are series R + jX only (c1=0 c0=0), as in the other decks.\n")
        for br in c.branches
            println(io, "New Line.$(br.name) Phases=1 Bus1=$(br.a).1 Bus2=$(br.b).1 r1=$(f3(real(br.z))) x1=$(f3(imag(br.z))) c1=0 c0=0 length=1 units=none")
        end
    end
    open(joinpath(OUT, "Ties.dss"), "w") do io
        print(io, hdr, "! The original's normally-open switches, disabled. Close one with Line.<name>.enabled=true.\n",
              "! A switch has (almost) no impedance; each is given the 1e-3 + j1e-3 ohm of the other decks'\n",
              "! shortest lines rather than its 1e-6 ohm, to keep the admittance matrix well scaled.\n")
        for t in c.ties
            println(io, "New Line.$(t.name) Phases=1 Bus1=$(t.a).1 Bus2=$(t.b).1 r1=0.001 x1=0.001 c1=0 c0=0 length=1 units=none enabled=false")
        end
    end
    open(joinpath(OUT, "Loads.dss"), "w") do io
        print(io, hdr, "! Each customer load summed onto the primary bus of its service transformer, divided by 3.\n")
        for b in sort(collect(keys(c.load)))
            s = c.load[b]; real(s) > 0 || continue
            println(io, "New load.L_$b Bus=$b.1 Phases=1 Model=1 kV=$KV_LN kW=$(f3(real(s))) kVAr=$(f3(imag(s))) Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault")
        end
    end
    open(joinpath(OUT, "Capacitor.dss"), "w") do io
        print(io, hdr, "! The original's capacitor banks (kvar / 3), not redirected by Master.dss.\n")
        for b in sort(collect(keys(c.caps)))
            println(io, "New generator.$(b)C Phases=1 Bus1=$b.1 kV=$KV_LN Model=1 kW=1e-06 kVAr=$(f3(c.caps[b]))")
        end
    end
    open(joinpath(OUT, "PVSystem.dss"), "w") do io
        print(io, hdr, @sprintf("! %d PV units on load buses, sized in proportion to the bus load; fleet Pmpp = %.0f %% of total load.\n", length(pv), 100DER_SHARE))
        for b in sort(collect(keys(pv)))
            println(io, "New PVsystem.PV_$b irrad=1.0 Phases=1 Bus1=$b.1 kV=$KV_LN kVA=$(f3(1.2pv[b])) Pmpp=$(f3(pv[b])) %cutin=0.001 %cutout=0.001")
        end
    end
    open(joinpath(OUT, "Storage.dss"), "w") do io
        print(io, hdr, @sprintf("! %d batteries on load buses, sized in proportion to the bus load; fleet rating = %.0f %% of total load, 4 h.\n", length(bt), 100DER_SHARE))
        for b in sort(collect(keys(bt)))
            println(io, "New Storage.BESS_$b Phases=1 Bus1=$b.1 kV=$KV_LN kVA=$(f3(1.2bt[b])) kWrated=$(f3(bt[b])) kWhrated=$(f3(4bt[b])) %stored=62.5 %reserve=30 %EffCharge=95 %EffDischarge=95 %idlingkW=0 Vminpu=0.95 Vmaxpu=1.05 DispMode=External")
        end
    end
    open(joinpath(OUT, "Master.dss"), "w") do io
        print(io, "Clear\n", hdr, """
New Circuit.ieee9500_1ph basekv=$KV_LN Bus1=s1 pu=1.05 phases=1 R1=0 X1=0.0001 R0=0 X0=0.0001
VSource.source.enabled=false  //Disable the default VSource created when declaring the circuit.
! Each substation is an EMF behind its 69/12.47 kV transformer (20 MVA, X = 7.94 %, R = 2 x 0.67 %):
! Zs = $(f3(real(ZS))) + j$(f3(imag(ZS))) ohm. A 1-phase Vsource's impedance is (2 Z1 + Z0)/3, so Z1 = Z0 = Zs.
""")
        for (i, s) in enumerate(c.subs)
            println(io, "New Vsource.subs$i Phases=1 Bus1=s$i basekv=$KV_LN pu=1.05 angle=0.0 R1=$(f3(real(ZS))) X1=$(f3(imag(ZS))) R0=$(f3(real(ZS))) X0=$(f3(imag(ZS)))")
        end
        println(io, "Redirect BranchData.dss\nRedirect Ties.dss\n! Substation connection lines")
        for (i, s) in enumerate(c.subs)
            println(io, "New Line.poi$i Phases=1 Bus1=s$i Bus2=$s.1 r1=1e-06 x1=1e-06")
        end
        print(io, """
Redirect ../LoadShapeDefault.dss
Redirect ../LoadShapePVDefault.dss
Redirect Loads.dss     ! Balanced Loads
// Redirect  Capacitor.dss
Redirect PVSystem.dss
Redirect Storage.dss

Set VoltageBases = [12.47]
Set mode = Daily
Set stepsize = 1h
Set number = 1

CalcVoltageBases
""")
    end
end

function main()
    c = convert(); f = feeders(c)
    @printf("original: %d primary lines (%d kept after merging), %d regulators merged; %d buses after merging\n",
            c.nlines, length(c.branches), c.nreg, length(c.buses))
    @printf("parallel lines combined: %d; unreached buses: %d\n", count(b -> b.n > 1, c.branches), count(b -> !haskey(f, b), c.buses))
    zs = sort(abs.(getfield.(c.branches, :z)))
    @printf("|z| per branch: min %.2e, 1%% %.2e, median %.2e, max %.2e ohm; R/X median %.2f\n",
            zs[1], zs[max(1, length(zs) ÷ 100)], zs[end÷2], zs[end], sort([real(b.z) / imag(b.z) for b in c.branches])[end÷2])
    tot = sum(values(c.load))
    @printf("load (per phase): %.1f kW %.1f kvar on %d buses; by feeder %s kW\n", real(tot), imag(tot), count(s -> real(s) > 0, values(c.load)),
            join([@sprintf("%.0f", sum(real(s) for (b, s) in c.load if get(f, b, 0) == i)) for i in 1:3], " / "))
    for t in c.ties
        @printf("  tie %-16s %d-phase  feeder %s -- feeder %s\n", t.name, t.k, get(f, t.a, 0), get(f, t.b, 0))
    end
    pv, bt, nboth = place_ders(c)
    @printf("DERs: %d PV, %d batteries, %d buses with both; PV %.0f kW, batteries %.0f kW / %.0f kWh\n",
            length(pv), length(bt), nboth, sum(values(pv)), sum(values(bt)), 4sum(values(bt)))
    write_deck(c, pv, bt)
    println("wrote ", OUT)
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
