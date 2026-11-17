# Feeder structure of a multi-source deck: the network without its ties, split
# by substation; per substation, the branches directly below it and the
# batteries in each. Usage: julia feeder_blocks_dss.jl <deck dir>
using OpenDSSDirect
const ODD = OpenDSSDirect
dir = ARGS[1]
ODD.Text.Command("clear")
ODD.Text.Command("compile [$(joinpath(dir, "Master.dss"))]")
bus(s) = lowercase(split(s, '.')[1])
tie = Set{String}()
if isfile(joinpath(dir, "Ties.dss"))
    for l in eachline(joinpath(dir, "Ties.dss"))
        m = match(r"New\s+(\w+\.\S+)"i, l); isnothing(m) || push!(tie, lowercase(m.captures[1]))
    end
end
adj = Dict{String, Vector{String}}()
addedge(a, b) = (push!(get!(adj, a, String[]), b); push!(get!(adj, b, String[]), a))
ODD.Circuit.SetActiveClass("Line")  # touch
for name in ODD.PDElements.AllNames()
    lowercase(name) in tie && continue
    ODD.Circuit.SetActiveElement(name)
    bs = ODD.CktElement.BusNames()
    length(bs) >= 2 || continue
    a, b = bus(bs[1]), bus(bs[2]); a == b || addedge(a, b)
end
srcs = String[]
for n in ODD.Vsources.AllNames()
    ODD.Vsources.Name(n); ODD.Circuit.SetActiveElement("Vsource." * n)
    push!(srcs, bus(ODD.CktElement.BusNames()[1]))
end
bats = String[]
for n in ODD.Storages.AllNames()
    ODD.Circuit.SetActiveElement("Storage." * n); push!(bats, bus(ODD.CktElement.BusNames()[1]))
end
println("deck ", basename(dir), ": buses=", length(adj), " batteries=", length(bats), " sources=", length(srcs), " ties excluded=", length(tie))
for s in unique(srcs)
    haskey(adj, s) || (println("  source $s: not connected"); continue)
    # the path from the source down to the first bus with more than one way on
    root = s; prev = ""
    while true
        nb = [x for x in adj[root] if x != prev]
        length(nb) == 1 || break
        prev, root = root, nb[1]
    end
    owner = Dict{String, String}(root => "")
    for c in adj[root]
        c == prev && continue
        stack = [c]; owner[c] = c
        while !isempty(stack)
            x = pop!(stack)
            for y in adj[x]
                haskey(owner, y) || y == prev || (owner[y] = c; push!(stack, y))
            end
        end
    end
    cnt = Dict{String, Int}()
    nb = 0
    for b in bats
        o = get(owner, b, nothing); isnothing(o) && continue
        nb += 1; o == "" || (cnt[o] = get(cnt, o, 0) + 1)
    end
    v = sort(collect(values(cnt)); rev=true)
    tot = sum(v; init=0)
    work = tot == 0 ? 0.0 : sum(c^3 for c in v) / tot^3
    println("  source ", s, ": first branching bus ", root, " (", s == root ? "at the source" : "below a single line", "), branches=",
            count(!=(prev), adj[root]), ", batteries reached=", nb, ", feeders with batteries=", length(v),
            ", blocks=", v[1:min(8, end)], ", battery-block work vs dense=", round(work; sigdigits=3))
end
