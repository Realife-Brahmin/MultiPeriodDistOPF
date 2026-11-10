# Current ratings for a deck that publishes none (user, 2026-11-03: ieee123). Two modes:
#   julia problems/msopf/root_level/uniform_ratings.jl <deck> blanket <factor> <envelope.csv> [<envelope.csv> ...]
#     every line at FACTOR x the largest current in any validated dispatch (the envelopes of a Farkas campaign), so the
#     rating cannot cut an OPF edge found without it;
#   julia problems/msopf/root_level/uniform_ratings.jl <deck> tiers
#     ieee123 only: each line at its IEEE 123-node conductor's ampacity x k/3 (the per-phase equivalent carries 1/3 of
#     each load, so a k-phase section carries k/3 of the real phase current). The conductor class is read from the
#     line's X/R, which separates the benchmark's configurations exactly (55 / 3 / 55 / 5 lines, as in its line data):
#       X/R 2.05-2.06: configs 1-6, 3-phase 336,400 26/7 ACSR, 530 A (Kersting)            -> 530 A
#       X/R 2.24-2.25: configs 7-8, 2-phase 336,400 26/7 ACSR, 530 A                        -> 353.3 A
#       X/R 1.00-1.01: configs 9-11, 1-phase 1/0 ACSR, 230 A (Kersting)                     -> 76.7 A
#       X/R 0.48:      config 12, 3-phase 1/0 AA 1/3 concentric-neutral cable, 195 A
#                      (direct buried, 90 C; catalog listings give 195-225 A)                -> 195 A
#     Zero-impedance elements (switches, the regulator stand-ins, the substations' poi lines) take the conductor in
#     series with them, as branch_ratings.jl's switches: at an end with exactly one other branch, that branch's rating;
#     otherwise the smaller end's largest rated branch.
# The old per-line envelope limits went down to the 1e-4 pu floor on dead or unloaded spurs (poi2/3/5, l36, l102, ...),
# which no reader would take for a thermal rating. Writes test_networks/<deck>/BranchRatings.csv in branch_ratings.jl's
# format (limit_pu is what the campaigns read).
using Printf
deck, mode = ARGS[1], ARGS[2]
root = normpath(joinpath(@__DIR__, "..", "..", ".."))
master = read(joinpath(root, "test_networks", deck, "Master.dss"), String)
kvb = parse(Float64, match(r"(?i)New\s+Vsource\.\S+.*?basekv=([0-9.]+)", master)[1])
ibase = 1000.0 / kvb                          # A per pu: 1000 kVA base over the sources' line-to-neutral kV (the solver's base)
lines = NamedTuple[]
for f in ("BranchData.dss", "Master.dss"), l in eachline(joinpath(root, "test_networks", deck, f))
    m = match(r"(?i)^\s*New\s+Line\.(\S+)\s.*?Bus1=(\S+)\s+Bus2=(\S+)\s+r1=(\S+)\s+x1=(\S+)", l)
    m === nothing || push!(lines, (; name = lowercase(m[1]), a = m[2], b = m[3], r = parse(Float64, m[4]), x = parse(Float64, m[5])))
end
A = Dict{String,Float64}(); src = Dict{String,String}()
if mode == "blanket"
    factor = parse(Float64, ARGS[3])
    peak, where = let p = 0.0, w = ""
        for f in ARGS[4:end], l in readlines(f)[2:end]
            c = split(l, ','); c[2] == "line" || continue
            v = parse(Float64, c[3]); v > p && (p = v; w = "$(c[1]) in $(basename(dirname(f)))/$(basename(f))")
        end
        p, w
    end
    for l in lines; A[l.name] = factor * peak * ibase; src[l.name] = @sprintf("%.2f x the largest validated current (%.4f pu at %s)", factor, peak, where); end
elseif mode == "tiers"
    for l in lines
        l.r < 1e-5 && continue
        q = l.x / l.r
        t = 2.0 <= q < 2.15 ? (530.0, "3-phase 336,400 26/7 ACSR, 530 A") : 2.15 <= q < 2.4 ? (530.0 * 2 / 3, "2-phase 336,400 26/7 ACSR, 530 A x 2/3") :
            0.9 <= q < 1.1 ? (230.0 / 3, "1-phase 1/0 ACSR, 230 A x 1/3") : 0.4 <= q < 0.6 ? (195.0, "3-phase 1/0 AA concentric-neutral cable, 195 A") :
            error(@sprintf("line %s: X/R %.3f matches no IEEE 123 configuration", l.name, q))
        A[l.name] = t[1]; src[l.name] = t[2]
    end
    nbr = Dict{String,Vector{Int}}()
    for (i, l) in enumerate(lines); push!(get!(nbr, l.a, Int[]), i); push!(get!(nbr, l.b, Int[]), i); end
    # A chain of zero-impedance lines (poi1-L1 at the head, L13-L14 at 11-126-16) is one element: its rated neighbours at
    # each bus of the chain; a bus with exactly one gives the conductor in series, else the smaller bus's largest. A bus
    # with none (a substation's own bus) is skipped, so a poi line takes its feeder end's conductor.
    for (i, l) in enumerate(lines)
        haskey(A, l.name) && continue
        chain = Set([i]); buses = Set([l.a, l.b]); grow = true
        while grow
            grow = false
            for bus in collect(buses), j in get(nbr, bus, Int[])
                j in chain || haskey(A, lines[j].name) || (push!(chain, j); push!(buses, lines[j].a, lines[j].b); grow = true)
            end
        end
        rated = [[A[lines[j].name] for j in get(nbr, bus, Int[]) if !(j in chain)] for bus in buses]
        rated = filter(!isempty, rated)
        isempty(rated) && error("zero-impedance chain at $(l.name) touches no rated line")
        pass = [r[1] for r in rated if length(r) == 1]
        v = isempty(pass) ? minimum(maximum.(rated)) : maximum(pass)
        for j in chain; A[lines[j].name] = v; src[lines[j].name] = "zero impedance: the conductor in series"; end
    end
    missing_ = [l.name for l in lines if !haskey(A, l.name)]
    isempty(missing_) || error("no rating for ", join(missing_, ", "))
else
    error("mode must be blanket or tiers")
end
open(joinpath(root, "test_networks", deck, "BranchRatings.csv"), "w") do io
    println(io, "name,kind,members,phases,rated_by,original,limit_A,limit_pu,source")
    for l in lines
        @printf(io, "%s,line,%s,1,%s,NaN,%.4f,%.8f,%s\n", l.name, l.name, mode, A[l.name], A[l.name] / ibase, src[l.name])
    end
end
cnt = Dict{String,Int}(); for l in lines; cnt[src[l.name]] = get(cnt, src[l.name], 0) + 1; end
@printf("%s (%s, base %.4f kV, %.3f A/pu): %d lines\n", deck, mode, kvb, ibase, length(lines))
for (k, v) in sort(collect(cnt)); @printf("  %3d  %s\n", v, k); end
