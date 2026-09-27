# When is the diagonal stage Hessian a safe approximation? Measures, on
# exact-Hessian captures, how diagonal the battery-power block of the stage
# Hessian is and what its diagonal is made of.
#
# With linear battery dynamics (f_u = -dt on the battery powers) the exact
# block is
#     H_pb = diag(cost) + diag(Sigma_pb) + dt^2 * Vxx
#   cost   : 2 C_B S^2 dt (+ 2 gamma dt^2 at the last stage), constant
#   Sigma  : interior-point barrier terms of the battery-power bounds
#   Vxx    : the incoming value-function Hessian -- the only off-diagonal part,
#            and exactly what the diagonal approximation drops.
# Per capture (stage, iteration) it reports
#   diag_energy        ||diag(H_pb)||^2 / ||H_pb||_F^2
#   dom_min/median     row dominance |H_jj| / sum_{k != j} |H_jk|
#   frac_nondominant   share of rows with dominance < 1
#   rho                ||D^{-1/2} E D^{-1/2}||_2, E = H_pb - D, D = diag(H_pb):
#                      rho < 1 means H_pb is positive definite and within a
#                      factor (1 +- rho) of its diagonal in every direction
#   lam_min_H / min_diag   smallest eigenvalue of H_pb against smallest pivot of D
#   share_cost/barrier/vxx the median share of each part in the diagonal
#
# Follows a capture directory like kkt_evolution_analysis.jl (done-file,
# KKT_EVOL_DELETE=1 deletes processed captures).
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/hessian_dominance_analysis.jl \
#         <capture_dir> <network_data.jls> <system> <cb_label> <out.csv> [done_file]

using Serialization, SparseArrays, LinearAlgebra, Printf, Statistics

capdir, datafile, system, cblabel, outcsv = ARGS[1:5]
donefile = length(ARGS) >= 6 ? ARGS[6] : ""
delete_after = get(ENV, "KKT_EVOL_DELETE", "0") == "1"
BLAS.set_num_threads(1)

d = deserialize(datafile)
N, L, B, D = length(d[:Nset]), length(d[:Lset]), length(d[:Bset]), length(d[:Dset])
dt = d[:delta_t_h]
pb0 = 2 + L + L + N + L                # ps, qs, P, Q, v, ell precede pb (control_layout)
pb = (pb0 + 1):(pb0 + B)

io = open(outcsv, isfile(outcsv) ? "a" : "w")
filesize(outcsv) == 0 && println(io, "system,cb,stage,iter,mu,B,diag_energy,dom_min,dom_median," *
    "frac_nondominant,rho,lam_min_H,min_diag,share_cost_median,share_barrier_median,share_vxx_median")

function analyse(c)
    Hb = Matrix(c.K[pb, pb])
    Hb = (Hb + Hb') / 2
    dg = diag(Hb)
    E = Hb - Diagonal(dg)
    off = vec(sum(abs, E; dims=2))
    dom = [off[j] == 0 ? Inf : abs(dg[j]) / off[j] for j in eachindex(dg)]
    rho = if all(>(0), dg)
        s = 1 ./ sqrt.(dg)
        maximum(abs, eigvals(Symmetric(s .* E .* s')))
    else
        Inf
    end
    bar = (c.Sigma_L .+ c.Sigma_U)[pb]
    vxx = dt^2 .* diag(c.Vxx_incoming)
    cost = dg .- bar .- vxx
    share(x) = median(x ./ dg)
    @printf(io, "%s,%s,%d,%d,%.6e,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
            system, cblabel, c.stage, c.iteration, c.barrier_mu, B,
            sum(abs2, dg) / sum(abs2, Hb), minimum(dom), median(dom), count(<(1), dom) / B,
            rho, eigmin(Symmetric(Hb)), minimum(dg), share(cost), share(bar), share(vxx))
    flush(io)
end

parse_name(f) = (m = match(r"iter(\d+)_stage(\d+)\.jls$", f); m === nothing ? nothing :
                 (parse(Int, m[1]), parse(Int, m[2])))
done = Set{String}(); sizes = Dict{String, Int}(); n = 0
while true
    finished = !isempty(donefile) && isfile(donefile)
    files = sort([f for f in readdir(capdir) if parse_name(f) !== nothing && !(f in done)], by = parse_name)
    ready = String[]
    for f in files
        sz = filesize(joinpath(capdir, f))
        if finished || isempty(donefile) ||
           (get(sizes, f, -1) == sz && sz > 0 && time() - mtime(joinpath(capdir, f)) > 5)
            push!(ready, f)
        end
        sizes[f] = sz
    end
    for f in ready
        path = joinpath(capdir, f)
        analyse(deserialize(path))
        push!(done, f); global n += 1
        delete_after && rm(path)
    end
    (isempty(donefile) || finished) && isempty(setdiff(Set(files), Set(ready))) && break
    isempty(ready) && sleep(3)
end
close(io)
println("HESSIAN_DOMINANCE processed=$n captures -> $outcsv")
