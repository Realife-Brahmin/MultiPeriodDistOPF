# Compile and solve each downloaded test system once (snapshot, its own master file) and write
# test_networks/sources/<folder>/summary.csv plus a combined test_networks/sources/summary.csv. Nothing on disk is
# changed; what is changed in memory is listed per system below and in each folder's README.
#   julia --project=problems/msopf tools/summarize_test_systems.jl
# Uses OpenDSSDirect.jl (DSS C-API, the same engine OpenDSSDirect.py wraps).
using OpenDSSDirect, Printf
const ODD = OpenDSSDirect
const ROOT = joinpath(normpath(joinpath(@__DIR__, "..")), "test_networks", "sources")

# Commands of a master that run a study or open a window rather than build the circuit.
const SKIP = r"(?i)^\s*(solve|plot|show|export|visualize|buscoords|latlongcoords|calcv\w*\s*$|set\s+mode)\b"

"Master commands, run one by one from its folder; `swap` replaces a Redirect target by commands."
function compile(master; swap = Dict{String,Vector{String}}())
    cmds = String[]
    for l in readlines(master)
        l = strip(split(split(l, '!')[1], "//")[1]); isempty(l) && continue
        if startswith(l, "~")
            cmds[end] *= " " * l[2:end]
        else
            push!(cmds, l)
        end
    end
    ODD.Text.Command("Clear")
    ODD.Text.Command("cd \"$(dirname(master))\"")
    for c in cmds
        occursin(SKIP, c) && continue
        m = match(r"(?i)^redirect\s+(\S+)", c)
        if m !== nothing && haskey(swap, lowercase(m[1]))
            foreach(ODD.Text.Command, swap[lowercase(m[1])])
        else
            ODD.Text.Command(c)
        end
    end
    ODD.Text.Command("Calcvoltagebases")
    ODD.Text.Command("Set mode=snapshot")
    ODD.Text.Command("Set maxcontroliter=200")
    ODD.Text.Command("Set maxiterations=100")      # IEEE8500 needs more than the default 15 (its Run file sets 20)
    ODD.Solution.Solve()
end

"SMART-DS: its LoadShapes.dss reads annual profiles not downloaded; each shape as a flat 1.0 instead."
function flat_shapes(file)
    names = [m[1] for m in eachmatch(r"(?i)New\s+Loadshape\.(\S+)", read(file, String))]
    return ["New Loadshape.$n npts=1 interval=1 mult=(1) qmult=(1)" for n in names]
end

function summary(system, circuit, master; kw...)
    t = time()
    err = ""
    try
        compile(master; kw...)
    catch e
        err = replace(sprint(showerror, e), r"\s+" => " ")[1:min(end, 200)]
    end
    conv = err == "" && ODD.Solution.Converged()
    nlines = ODD.Lines.Count()
    sw = 0; openl = 0
    for nm in ODD.Lines.AllNames()
        nm == "NONE" && continue
        ODD.Lines.Name(nm); ODD.Lines.IsSwitch() && (sw += 1)
        ODD.Circuit.SetActiveElement("Line.$nm")
        (ODD.CktElement.IsOpen(1, 0) || ODD.CktElement.IsOpen(2, 0)) && (openl += 1)
    end
    loadkw = 0.0
    for nm in ODD.Loads.AllNames()
        nm == "NONE" && continue
        ODD.Loads.Name(nm); loadkw += ODD.Loads.kW()
    end
    losses = real(ODD.Circuit.Losses()) / 1e3
    src = -real(ODD.Circuit.TotalPower())
    kvb = sort(unique(round.(ODD.Settings.VoltageBases(); digits = 3)))
    row = (; system, circuit, master = relpath(master, ROOT), buses = ODD.Circuit.NumBuses(), nodes = ODD.Circuit.NumNodes(),
           lines = nlines, loads = ODD.Loads.Count(), transformers = ODD.Transformers.Count(), regulators = ODD.RegControls.Count(),
           capacitors = ODD.Capacitors.Count(), switches = sw, open_lines = openl, vsources = ODD.Vsources.Count(),
           converged = conv ? "yes" : "no", iterations = ODD.Solution.Iterations(), total_load_kW = round(loadkw; digits = 1),
           source_kW = round(src; digits = 1), losses_kW = round(losses; digits = 1), kV_bases = join(kvb, " "),
           seconds = round(time() - t; digits = 1), error = err)
    @printf("%-20s %-28s buses %6d  loads %5d  regs %3d  caps %3d  switches %4d (open %2d)  vsources %d  %s  load %9.1f kW  losses %7.1f kW%s\n",
            system, circuit, row.buses, row.loads, row.regulators, row.capacitors, row.switches, row.open_lines, row.vsources,
            conv ? "converged" : "NOT converged", loadkw, losses, err == "" ? "" : "  ERROR " * err)
    return row
end

const EPRI = joinpath(ROOT, "epri_test_circuits")
const SDS = joinpath(ROOT, "smartds_sfo_p1u")
cases = Any[("epri_test_circuits", "ckt5", joinpath(EPRI, "ckt5", "Master_ckt5.dss"), (;)),
         ("epri_test_circuits", "ckt7", joinpath(EPRI, "ckt7", "Master_ckt7.dss"), (;)),
         ("epri_test_circuits", "ckt24", joinpath(EPRI, "ckt24", "master_ckt24.dss"), (;)),
         ("epri_test_circuits", "J1 (no PV)", joinpath(EPRI, "epri_dpv", "J1", "Master_noPV.dss"), (;)),
         ("epri_test_circuits", "K1 (no PV)", joinpath(EPRI, "epri_dpv", "K1", "Master_NoPV.dss"), (;)),
         ("epri_test_circuits", "M1 (no PV)", joinpath(EPRI, "epri_dpv", "M1", "Master_NoPV.dss"), (;)),
         ("ieee8500", "IEEE8500 (balanced)", joinpath(ROOT, "ieee8500", "Master.dss"), (;)),
         ("ieee8500", "IEEE8500 (unbalanced)", joinpath(ROOT, "ieee8500", "Master-unbal.dss"), (;))]
shapes = Dict("loadshapes.dss" => flat_shapes(joinpath(SDS, "LoadShapes.dss")))
push!(cases, ("smartds_sfo_p1u", "P1U full region", joinpath(SDS, "Master.dss"), (; swap = shapes)))
for d in sort(filter(x -> isdir(joinpath(SDS, x)) && startswith(x, "p1uhs"), readdir(SDS)))
    own = joinpath(SDS, d, "LoadShapes.dss")
    sw = Dict("loadshapes.dss" => isfile(own) ? flat_shapes(own) : shapes["loadshapes.dss"])
    push!(cases, ("smartds_sfo_p1u", "substation $d", joinpath(SDS, d, "Master.dss"), (; swap = sw)))
end

rows = [summary(s, c, m; kw...) for (s, c, m, kw) in cases]
csv(r) = join([occursin(',', string(v)) ? "\"" * string(v) * "\"" : string(v) for v in values(r)], ',')
hdr = join(string.(keys(rows[1])), ',')
for s in unique(r.system for r in rows)
    open(io -> (println(io, hdr); foreach(r -> println(io, csv(r)), filter(r -> r.system == s, rows))), joinpath(ROOT, s, "summary.csv"), "w")
end
open(io -> (println(io, hdr); foreach(r -> println(io, csv(r)), rows)), joinpath(ROOT, "summary.csv"), "w")
println("wrote summary.csv files under ", ROOT)
