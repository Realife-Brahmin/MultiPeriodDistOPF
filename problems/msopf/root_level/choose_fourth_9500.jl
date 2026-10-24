# Where to add a fourth substation to ieee9500_1ph (the five ties between feeders closed): a bus about
# equally far, electrically, from the deck's three substations, and far from each, so that the four-
# substation region is not a thin wedge (as on large10kD, whose added substations sat on one feeder).
# Distance = the Thevenin impedance between two buses through the network's series impedances. Candidates
# must sit on a three-phase section of the original model (not a single-phase lateral).
#   julia --project=problems/msopf problems/msopf/root_level/choose_fourth_9500.jl
empty!(ARGS); push!(ARGS, "none")
include(joinpath(@__DIR__, "multiperiod_study.jl"))
using LinearAlgebra, Statistics
tc = THREE["three_ieee9500"]
net, _ = load_sources(["subs1", "subs2", "subs3"]; system = tc.system, extra = ["Line.$t.enabled=true" for t in tc.ties])
idx = Dict(b => i for (i, b) in enumerate(net.buses)); n = length(net.buses)
L = zeros(ComplexF64, n, n)
for l in net.lines
    y, i, j = -l.Y[1, 2], idx[l.from], idx[l.to]
    L[i, i] += y; L[j, j] += y; L[i, j] -= y; L[j, i] -= y
end
Z = zeros(ComplexF64, n, n); Z[1:n-1, 1:n-1] = inv(L[1:n-1, 1:n-1])
zeff(a, b) = abs(Z[a, a] + Z[b, b] - 2Z[a, b])
subs = [idx[s.bus] for s in net.sources]
say("between the substations: ", join([@sprintf("%.2f", zeff(subs[i], subs[j])) for i in 1:3 for j in i+1:3], ", "), " ohm")
# three-phase buses of the original: both ends of a three-phase line in the original deck
three = Set{String}()
for l in eachline(joinpath(REPO_ROOT, "test_networks", "ieee9500_source", "base", "LinesSwitchesGeometry.dss"))
    m = match(r"(?i)^New Line\.\S+\s+bus1=(\S+)\s+bus2=(\S+)(.*)$", l); m === nothing && continue
    occursin(r"(?i)phases=[12]", m[3]) && continue
    occursin(r"\.\d", m[1]) && !occursin(r"\.1\.2\.3", m[1]) && continue
    push!(three, lowercase(split(m[1], '.')[1])); push!(three, lowercase(split(m[2], '.')[1]))
end
rows = []
for (b, i) in idx
    (b in three) || continue
    d = [zeff(i, s) for s in subs]
    push!(rows, (; b, d, spread = (maximum(d) - minimum(d)) / mean(d), mind = minimum(d)))
end
dmax = maximum(r.mind for r in rows)
# score: equal distances first, then far from the nearest substation
sort!(rows; by = r -> r.spread - 0.5 * r.mind / dmax)
say("candidates (three-phase buses): distance to substations 1, 2, 3 (ohm), spread")
for r in rows[1:12]
    say(@sprintf("  %-22s %6.2f %6.2f %6.2f   spread %4.1f%%", r.b, r.d..., 100r.spread))
end
# central candidates: about as far from each substation as the substations are from one another
dsub = [zeff(subs[i], subs[j]) for i in 1:3 for j in i+1:3]
central = sort([r for r in rows if minimum(r.d) >= 0.75minimum(dsub) && maximum(r.d) <= 1.25maximum(dsub)]; by = r -> r.spread)
say("central candidates (each distance within 0.75-1.25x the substations' own separations):")
for r in central[1:min(12, end)]
    say(@sprintf("  %-22s %6.2f %6.2f %6.2f   spread %4.1f%%", r.b, r.d..., 100r.spread))
end
