# Merge the multi_windows_<dirset>.csv rows of a folder of Kamiak job outputs into
# results/multiperiod/multi_windows_<dirset>.csv: a later job (higher Slurm id) replaces an earlier one's row
# for the same (system, sources, T, case, direction). A row from a job that ran before capped solves were
# retried (no "retrying" support in its code: its log never prints "capped solve") and that parked a solve
# along that row's direction is flagged "+parked_oldcode" in the method column: the old bracket read a
# capped solve as infeasible, so such an edge may be low.
#   julia --project=envs/msopf envs/multi_poi/root_level/merge_multi.jl <jobs folder> [dirset] [stage]
using Printf, LinearAlgebra
const RES = joinpath(@__DIR__, "..", "results", "multiperiod")
dir = ARGS[1]; dset = length(ARGS) >= 2 ? ARGS[2] : "fib78"; stage = length(ARGS) >= 3 ? ARGS[3] : ""
fib_dirs(n) = [(z = 1 - 2(i + 0.5) / n; r = sqrt(1 - z^2); ph = i * pi * (3 - sqrt(5)); [r * cos(ph), r * sin(ph), z]) for i in 0:n-1]
DIRS = dset == "fib78" ? fib_dirs(78) : error("only fib78 so far")
jobid(d) = parse(Int, match(r"_(\d+)$", d)[1])
jobs = sort(filter(d -> isdir(joinpath(dir, d)) && occursin(r"_\d+$", d), readdir(dir)); by = jobid)
header = ""; rows = Dict{String,Vector{String}}(); src = Dict{String,String}()
for j in jobs
    f = joinpath(dir, j, "multi_windows_$dset.csv"); isfile(f) || continue
    ls = readlines(f); global header = ls[1]
    for l in ls[2:end]
        r = String.(split(l, ',')); k = join(r[1:5], ',')
        rows[k] = r; src[k] = j
    end
end
# old-code jobs that parked solves: flag the directions they parked along
nflag = 0
for j in jobs
    p = joinpath(dir, j, "parked_solves.csv"); isfile(p) || continue
    log = joinpath(dir, j, "run.log")
    (isfile(log) && occursin("capped solve", read(log, String))) && continue      # new code: retried, resolved or flagged itself
    for l in readlines(p)[2:end]
        r = split(l, ','); case = r[2]
        ang = Dict(split(kv, '=')[1] => parse(Float64, split(kv, '=')[2]) for kv in split(r[5], ';'))
        v = [ang[k] for k in sort(collect(keys(ang))) if ang[k] != 0 || k != first(sort(collect(keys(ang))))]
        length(v) == 3 || continue
        u = v / norm(v); i = argmax([dot(u, d) for d in DIRS])
        for (k, row) in rows
            if row[4] == case && parse(Int, row[5]) == i && src[k] == j && !occursin("parked_oldcode", row[17])
                row[17] *= "+parked_oldcode"; global nflag += 1
            end
        end
    end
end
out = joinpath(RES, "multi_windows_$dset.csv")
old = isfile(out) ? [String.(split(l, ',')) for l in readlines(out)[2:end]] : Vector{String}[]
keep = [r for r in old if !haskey(rows, join(r[1:5], ','))]                         # rows of other systems stay
allrows = sort([keep; collect(values(rows))]; by = r -> (r[1], r[2], r[4], parse(Int, r[5])))
open(io -> (println(io, header); foreach(r -> println(io, join(r, ',')), allrows)), out, "w")
@printf("merged %d rows from %d jobs into %s (%d rows in all); %d flagged +parked_oldcode\n", length(rows), length(jobs), basename(out), length(allrows), nflag)
