#!/usr/bin/env julia
#
# check_voltage_bases.jl -- diagnose per-unit voltage bases in the OpenDSS decks.
#
# WHY THIS EXISTS
# ---------------
# `Set VoltageBases = [...]` is specified in kV LINE-TO-LINE. `CalcVoltageBases`
# then assigns each bus a kVBase by matching against that list, dividing by
# sqrt(3) for buses with fewer than three phases. For a 1-phase Vsource,
# `basekv` is already LINE-TO-NEUTRAL. So a correct 1-phase deck must declare
#
#     Set VoltageBases = [ basekv * sqrt(3) ]
#
# If a deck instead sets VoltageBases equal to the L-N `basekv`, every bus base
# comes out sqrt(3) too small and every per-unit voltage reads sqrt(3) too high.
# The solved voltages in VOLTS are unaffected -- `Vsource.basekv` sets the real
# source voltage correctly. Only the per-unit REPORTING base is wrong.
#
# THE TEST
# --------
# Each Vsource is commanded to a known per-unit value (`pu=` in the deck). By
# construction its own bus therefore sits at exactly that per-unit voltage. If
# OpenDSS reports a different per-unit value at that bus, the bus kVBase -- not
# the power flow -- is wrong. The ratio is printed; sqrt(3) = 1.7321.
#
# Run from the repo root:
#     julia --project=envs/multi_poi envs/multi_poi/check_voltage_bases.jl

using OpenDSSDirect
using Printf

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SQRT3 = sqrt(3)
const TOL = 0.01            # relative tolerance on commanded vs reported pu

const SYSTEMS = ["ads10A_1ph", "small2poi_1ph", "ieee123C_1ph",
                 "ieee123_5poi_1ph", "ieee2522C_1ph", "large10kC_1ph"]

"Strip any node qualifier: \"12.1.2\" -> \"12\"."
strip_nodes(bus) = String(split(bus, '.')[1])   # ODD's SetActiveBus accepts String only, not SubString

"Compile and solve one system; return its Vsource rows plus circuit-wide stats."
function probe(system::AbstractString)
    master = joinpath(REPO_ROOT, "rawData", system, "Master.dss")
    isfile(master) || error("Master.dss not found: $master")

    OpenDSSDirect.Text.Command("Clear")
    OpenDSSDirect.Text.Command("Redirect \"$master\"")
    OpenDSSDirect.Solution.Solve()

    converged = OpenDSSDirect.Solution.Converged()

    rows = NamedTuple[]
    for name in OpenDSSDirect.Vsources.AllNames()
        OpenDSSDirect.Vsources.Name(name)
        OpenDSSDirect.Circuit.SetActiveElement("Vsource.$name")

        # A disabled Vsource carries no bus in the solved circuit; skip it.
        buses = OpenDSSDirect.CktElement.BusNames()
        isempty(buses) && continue
        OpenDSSDirect.CktElement.Enabled() || continue

        basekv       = OpenDSSDirect.Vsources.BasekV()     # L-N for phases=1
        commanded_pu = OpenDSSDirect.Vsources.PU()         # what the deck asked for

        bus = strip_nodes(first(buses))
        OpenDSSDirect.Circuit.SetActiveBus(bus)
        kvbase      = OpenDSSDirect.Bus.kVBase()           # what CalcVoltageBases assigned
        reported_pu = first(OpenDSSDirect.Bus.puVmagAngle())

        push!(rows, (; name, bus, basekv, kvbase, commanded_pu, reported_pu,
                       ratio = commanded_pu == 0 ? NaN : reported_pu / commanded_pu))
    end

    pu = filter(>(0.01), OpenDSSDirect.Circuit.AllBusMagPu())
    return (; converged, rows, vmin = minimum(pu), vmax = maximum(pu))
end

function main()
    @printf "%-18s %-9s %9s %9s %7s %8s %8s  %s\n" "system" "vsource" "basekv" "kVBase" "cmd_pu" "got_pu" "ratio" "verdict"
    println("-"^104)

    suspect = String[]

    for system in SYSTEMS
        local r
        try
            r = probe(system)
        catch err
            @printf "%-18s  ERROR: %s\n" system sprint(showerror, err)
            continue
        end

        r.converged || @printf "%-18s  WARNING: power flow did NOT converge\n" system

        bad_here = false
        for (i, row) in enumerate(r.rows)
            ok = isfinite(row.ratio) && abs(row.ratio - 1) <= TOL
            bad_here |= !ok
            verdict = ok ? "ok" :
                      isapprox(row.ratio, SQRT3; rtol = 0.01) ? "OFF by sqrt(3)" :
                      @sprintf("OFF by %.4f", row.ratio)
            @printf "%-18s %-9s %9.4f %9.4f %7.4f %8.4f %8.4f  %s\n" (i == 1 ? system : "") row.name row.basekv row.kvbase row.commanded_pu row.reported_pu row.ratio verdict
        end

        @printf "%-18s %-9s %9s %9s %7s %8.4f %8.4f  (circuit vmin/vmax pu)\n" "" "" "" "" "" r.vmin r.vmax
        bad_here && push!(suspect, system)
        println()
    end

    println("="^104)
    if isempty(suspect)
        println("All decks report per-unit voltages consistent with their commanded Vsource setpoints.")
        println("No change to `Set VoltageBases` is warranted.")
    else
        println("Decks whose reported per-unit voltages disagree with their commanded setpoints:")
        println()
        for system in suspect
            master = joinpath(REPO_ROOT, "rawData", system, "Master.dss")
            current = strip(join(filter(l -> occursin(r"^\s*Set\s+VoltageBases"i, l),
                                        readlines(master)), " | "))
            OpenDSSDirect.Text.Command("Clear")
            OpenDSSDirect.Text.Command("Redirect \"$master\"")
            OpenDSSDirect.Solution.Solve()
            OpenDSSDirect.Vsources.First()
            suggested = OpenDSSDirect.Vsources.BasekV() * SQRT3
            @printf "  %-18s  now: %-32s  should be: Set VoltageBases = [%.4g]\n" system current suggested
        end
        println()
        println("Solved voltages in volts are NOT affected -- `Vsource.basekv` is correct in")
        println("every deck. Only the per-unit reporting base is wrong, which matters wherever")
        println("per-unit voltages are compared or checked against limits.")
    end
end

main()
