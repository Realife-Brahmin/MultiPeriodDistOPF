# Drift of the stage KKT matrix over several iterations -- what a factorization
# frozen for p iterations would be solving against. Reads the kept per-iteration
# captures (run_kkt_evolution.sh keeps ieee123's) and reports, per stage and
# block, ||K_{k+p} - K_k||_F / ||K_k||_F for p = 1, 2, 5, 10, with the barrier
# parameter at both ends.
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_evolution_drift.jl \
#         <capture_dir> <out.csv>

using Serialization, SparseArrays, LinearAlgebra, Printf

capdir, outcsv = ARGS[1], ARGS[2]
files = filter(f -> occursin(r"^iter\d+_stage\d+\.jls$", f), readdir(capdir))
key(f) = (m = match(r"iter(\d+)_stage(\d+)", f); (parse(Int, m[1]), parse(Int, m[2])))
caps = Dict{Tuple{Int,Int}, Any}()
for f in files
    c = deserialize(joinpath(capdir, f))
    nu = c.nu
    H = c.K[1:nu, 1:nu]
    bar = c.Sigma_L .+ c.Sigma_U
    caps[key(f)] = (mu = c.barrier_mu, K = c.K, Hd = Vector(diag(H)), bar = bar,
                    J = c.K[nu+1:end, 1:nu])
end
iters = sort(unique(first.(keys(caps)))); stages = sort(unique(last.(keys(caps))))
rel(a, b) = (n = norm(b); n == 0 ? 0.0 : norm(a - b) / n)
open(outcsv, "w") do io
    println(io, "stage,iter,gap,mu_start,mu_end,K,H_barrier,H_nonbarrier,J")
    for s in stages, k in iters, p in (1, 2, 5, 10)
        (haskey(caps, (k, s)) && haskey(caps, (k + p, s))) || continue
        a, b = caps[(k + p, s)], caps[(k, s)]
        @printf(io, "%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s, k, p, b.mu, a.mu,
                rel(a.K, b.K), rel(a.bar, b.bar), rel(a.Hd .- a.bar, b.Hd .- b.bar), rel(a.J, b.J))
    end
end
println("wrote $outcsv")
