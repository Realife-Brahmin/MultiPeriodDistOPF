# One blanket current rating for every line of a deck that publishes none (user, 2026-11-03: ieee123). The rating is
# FACTOR x the largest current in any validated dispatch of the study (the envelopes of its Farkas campaigns), so it
# cannot cut an OPF edge found without it; the old per-line envelope limits went down to the 1e-4 pu floor on dead or
# unloaded spurs (poi2/3/5, l36, l102, ...), which no reader would take for a thermal rating.
#   julia problems/msopf/root_level/uniform_ratings.jl <deck> <factor> <envelope.csv> [<envelope.csv> ...]
# Writes test_networks/<deck>/BranchRatings.csv in branch_ratings.jl's format (limit_pu is what the campaigns read).
using Printf
deck, factor = ARGS[1], parse(Float64, ARGS[2])
root = normpath(joinpath(@__DIR__, "..", "..", ".."))
peak, at = 0.0, ""
for f in ARGS[3:end], l in readlines(f)[2:end]
    c = split(l, ','); c[2] == "line" || continue
    v = parse(Float64, c[3]); v > peak && (global peak = v; global at = "$(c[1]) in $(basename(dirname(f)))/$(basename(f))")
end
names = String[]
for f in ("BranchData.dss", "Master.dss"), l in eachline(joinpath(root, "test_networks", deck, f))
    m = match(r"(?i)^\s*New\s+Line\.(\S+)", l); m === nothing || push!(names, lowercase(m[1]))
end
kvb = parse(Float64, match(r"(?i)New\s+Vsource\.\S+.*?basekv=([0-9.]+)", read(joinpath(root, "test_networks", deck, "Master.dss"), String))[1])
ibase = 1000.0 / kvb                          # A per pu: 1000 kVA base over the sources' line-to-neutral kV (the solver's base)
lim = factor * peak
open(joinpath(root, "test_networks", deck, "BranchRatings.csv"), "w") do io
    println(io, "name,kind,members,phases,rated_by,original,limit_A,limit_pu,source")
    for n in names
        @printf(io, "%s,line,%s,1,blanket,NaN,%.4f,%.8f,%.2f x the largest validated current (%.4f pu at %s)\n", n, n, lim * ibase, lim, factor, peak, at)
    end
end
@printf("%s: %d lines at %.6f pu = %.2f A (%.2f x %.6f pu, %s)\n", deck, length(names), lim, lim * ibase, factor, peak, at)
