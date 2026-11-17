# Can one substation carry a whole multi-source deck as a radial network?
#
#   julia --project=problems/mpopf/spar_ddp problems/mpopf/spar_ddp/experiments/scaling_effort_nov16/single_source_scan.jl <deck> [load scale]
#
# The deck's areas are the parts the substations feed with every tie open (the
# ties are the original normally-open switches, Ties.dss). For each substation
# s: every other source off, and ties closed one at a time to reach the other
# areas, each time the tie (into an area not yet reached) whose power flow
# keeps the lowest bus voltage highest; a tie inside an area or into a reached
# one would close a loop and is never used. Power flow at nominal load (times
# the scale) with the PV and the batteries off: the worst case the OPF faces at
# full load, before any dispatch helps. Prints, per substation, the lowest and
# highest voltage and the ties used; the DDP band is 0.95-1.05 pu.
using OpenDSSDirect, Printf
const ODD = OpenDSSDirect
const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
deck = ARGS[1]
scale = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.0
master = joinpath(REPO, "test_networks", deck, "Master.dss")
bus(s) = lowercase(split(s, '.')[1])

function compile!(extra::Vector{String})
    ODD.Text.Command("clear")
    ODD.Text.Command("redirect [$master]")
    for c in extra; ODD.Text.Command(c); end
    ODD.Text.Command("BatchEdit PVSystem..* enabled=false")
    ODD.Text.Command("BatchEdit Storage..* enabled=false")
    scale == 1.0 || ODD.Text.Command("Set LoadMult=$scale")
    ODD.Text.Command("Set mode=snapshot")
    ODD.Text.Command("Set maxiterations=100")
    ODD.Text.Command("Solve")
    ok = ODD.Solution.Converged()
    v = ODD.Circuit.AllBusMagPu()
    return ok, (ok ? minimum(v) : NaN), (ok ? maximum(v) : NaN), -real(first(ODD.Circuit.TotalPower()))
end

# topology with every tie open
compile!(String[])
tienames = String[]
for l in eachline(joinpath(dirname(master), "Ties.dss"))
    m = match(r"New\s+Line\.(\S+)"i, l); isnothing(m) || push!(tienames, lowercase(m.captures[1]))
end
adj = Dict{String, Vector{String}}()
for name in ODD.Lines.AllNames()
    lowercase(name) in tienames && continue
    ODD.Lines.Name(name); ODD.CktElement.Enabled() || continue
    bs = ODD.CktElement.BusNames(); a, b = bus(bs[1]), bus(bs[2])
    push!(get!(adj, a, String[]), b); push!(get!(adj, b, String[]), a)
end
srcs = Pair{String, String}[]
for n in ODD.Vsources.AllNames()
    ODD.Circuit.SetActiveElement("Vsource." * n)
    ODD.CktElement.Enabled() && push!(srcs, n => bus(ODD.CktElement.BusNames()[1]))
end
area = Dict{String, Int}()
for (k, (_, b)) in enumerate(srcs)
    haskey(area, b) && continue
    stack = [b]; area[b] = k
    while !isempty(stack)
        x = pop!(stack)
        for y in get(adj, x, String[]); haskey(area, y) || (area[y] = k; push!(stack, y)); end
    end
end
nbus = length(ODD.Circuit.AllBusNames())
@printf("%s: %d buses, %d in source areas, %d sources, %d ties, load scale %.2f\n", deck, nbus, length(area), length(srcs), length(tienames), scale)
ties = Tuple{String, Int, Int}[]
for t in tienames
    ODD.Circuit.SetActiveElement("Line." * t)
    bs = ODD.CktElement.BusNames(); a, b = get(area, bus(bs[1]), 0), get(area, bus(bs[2]), 0)
    push!(ties, (t, a, b))
end
inter = [x for x in ties if x[2] != x[3] && x[2] > 0 && x[3] > 0]
println("  ties between areas: ", length(inter), " (", join(sort(unique(["$(min(a,b))-$(max(a,b))" for (_, a, b) in inter])), " "), ")")
for (k, (s, sb)) in enumerate(srcs)
    off = ["Vsource.$n.enabled=false" for (n, _) in srcs if n != s]
    reached = Set([area[sb]]); closed = String[]
    allareas = Set(area[b] for (_, b) in srcs)
    best = (false, NaN, NaN, NaN)
    while reached != allareas
        cand = [x for x in inter if (x[2] in reached) != (x[3] in reached)]
        isempty(cand) && break
        best = (false, -Inf, NaN, NaN); bt = ""
        for (t, a, b) in cand
            r = compile!(vcat(off, ["Line.$c.enabled=true" for c in vcat(closed, [t])]))
            r[1] && r[2] > best[2] && (best = r; bt = t)
        end
        bt == "" && break
        push!(closed, bt); x = only(y for y in inter if y[1] == bt); push!(reached, x[2], x[3])
    end
    reached == allareas || (println("  ", s, ": cannot reach every area without loops"); continue)
    isempty(closed) && (best = compile!(off))
    @printf("  %-6s at %-28s lowest %.4f highest %.4f pu, import %.0f kW, %s; ties closed: %s\n", s, sb, best[2], best[3], best[4],
            best[1] && best[2] >= 0.95 && best[3] <= 1.05 ? "WITHIN 0.95-1.05" : "outside the band", join(closed, " "))
end
