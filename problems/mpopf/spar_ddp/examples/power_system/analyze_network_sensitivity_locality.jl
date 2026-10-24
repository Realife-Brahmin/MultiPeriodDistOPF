using LinearAlgebra
using Printf
using Serialization
using SparseArrays
using Statistics

length(ARGS) >= 3 || error("usage: capture.jls network_data.jls output_prefix")
capture_path, data_path, output_prefix = ARGS[1:3]

capture = deserialize(capture_path)
beta = capture.beta
omega = capture.omega
capture = nothing
GC.gc()

data = deserialize(data_path)
buses = data[:Nset]
lines = data[:Lset]
batteries = data[:Bset]
ders = data[:Dset]
root = data[:substationBus]
nx = length(batteries)
size(beta, 2) == nx || error("capture and network data have different state dimensions")

adjacency = Dict(bus => Int[] for bus in buses)
for (i, j) in lines
    push!(adjacency[i], j)
    push!(adjacency[j], i)
end

function distances_from(source)
    distance = Dict(source => 0)
    queue = Vector{Int}(undef, length(buses))
    queue[1] = source
    first = 1
    last = 1
    while first <= last
        bus = queue[first]
        first += 1
        for neighbor in adjacency[bus]
            if !haskey(distance, neighbor)
                distance[neighbor] = distance[bus] + 1
                last += 1
                queue[last] = neighbor
            end
        end
    end
    distance
end

# Match control_layout in ieee123c_filterddp.jl.
N, L, B, D = length(buses), length(lines), length(batteries), length(ders)
k = Ref(0)
take(n) = (r = (k[] + 1):(k[] + n); k[] += n; r)
idx = (ps=first(take(1)), qs=first(take(1)), P=take(L), Q=take(L),
       v=take(N), ell=take(L), pb=take(B), qnorm=take(D),
       soc_slack=take(L), energy_slack=take(B))

line_buses = [j for (_, j) in lines]
control_families = [
    ("P", idx.P, line_buses), ("Q", idx.Q, line_buses),
    ("v", idx.v, buses), ("ell", idx.ell, line_buses),
    ("pb", idx.pb, batteries), ("qnorm", idx.qnorm, ders),
    ("soc_slack", idx.soc_slack, line_buses),
    ("energy_slack", idx.energy_slack, batteries),
]

# Constraint order in build_model: P balance, Q balance, voltage drop, SOC,
# root voltage, then one energy equation per battery.
balance_buses = vcat([root], data[:Nm1set])
row = Ref(0)
rows(n) = (r = (row[] + 1):(row[] + n); row[] += n; r)
constraint_families = [
    ("P_balance", rows(N), balance_buses),
    ("Q_balance", rows(N), balance_buses),
    ("voltage_drop", rows(L), line_buses),
    ("soc", rows(L), line_buses),
    ("root_voltage", rows(1), [root]),
    ("energy", rows(B), batteries),
]
row[] == size(omega, 1) || error("constraint row mapping does not match omega")

radii = [0, 1, 2, 4, 8, 16, 32, 64, typemax(Int)]
fraction_rows = NamedTuple[]
radius_rows = NamedTuple[]

function analyze_family!(map_name, matrix, family_name, indices, associated_buses)
    fractions = [Float64[] for _ in radii]
    radius90 = Int[]
    radius99 = Int[]
    for state_index in 1:nx
        distance = distances_from(batteries[state_index])
        values = @view matrix[indices, state_index]
        weights = abs2.(values)
        total = sum(weights)
        total <= eps(Float64) && continue
        entry_distances = [distance[bus] for bus in associated_buses]
        for (rindex, radius) in enumerate(radii)
            push!(fractions[rindex], sum((weights[i] for i in eachindex(weights)
                if entry_distances[i] <= radius); init=0.0) / total)
        end
        order = sortperm(entry_distances)
        cumulative = 0.0
        found90 = false
        for i in order
            cumulative += weights[i] / total
            if !found90 && cumulative >= 0.90
                push!(radius90, entry_distances[i])
                found90 = true
            end
            if cumulative >= 0.99
                push!(radius99, entry_distances[i])
                break
            end
        end
    end
    for (rindex, radius) in enumerate(radii)
        values = fractions[rindex]
        isempty(values) && continue
        push!(fraction_rows, (map=map_name, family=family_name,
            radius=radius == typemax(Int) ? "all" : string(radius),
            mean=mean(values), median=median(values),
            p10=quantile(values, 0.10), p90=quantile(values, 0.90),
            columns=length(values)))
    end
    if !isempty(radius90)
        push!(radius_rows, (map=map_name, family=family_name,
            target="90pct", median=median(radius90), p10=quantile(radius90, 0.10),
            p90=quantile(radius90, 0.90), maximum=maximum(radius90)))
        push!(radius_rows, (map=map_name, family=family_name,
            target="99pct", median=median(radius99), p10=quantile(radius99, 0.10),
            p90=quantile(radius99, 0.90), maximum=maximum(radius99)))
    end
end

for (name, indices, associated) in control_families
    analyze_family!("beta", beta, name, indices, associated)
end
for (name, indices, associated) in constraint_families
    analyze_family!("omega", omega, name, indices, associated)
end

open(output_prefix * "_distance_energy.csv", "w") do io
    println(io, "map,family,radius,mean,median,p10,p90,columns")
    for r in fraction_rows
        @printf(io, "%s,%s,%s,%.12g,%.12g,%.12g,%.12g,%d\n",
            r.map, r.family, r.radius, r.mean, r.median, r.p10, r.p90, r.columns)
    end
end
open(output_prefix * "_locality_radius.csv", "w") do io
    println(io, "map,family,target,median_distance,p10_distance,p90_distance,max_distance")
    for r in radius_rows
        @printf(io, "%s,%s,%s,%.6g,%.6g,%.6g,%d\n",
            r.map, r.family, r.target, r.median, r.p10, r.p90, r.maximum)
    end
end

for r in radius_rows
    r.target == "90pct" && @printf(
        "LOCALITY map=%s family=%s radius90_median=%.1f radius90_p90=%.1f radius90_max=%d\n",
        r.map, r.family, r.median, r.p90, r.maximum)
end
