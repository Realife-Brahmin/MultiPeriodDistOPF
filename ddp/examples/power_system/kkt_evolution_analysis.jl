# How much does FilterDDP's stage KKT matrix change between iterations, and in
# which blocks? (Question from the meeting with A. of 2026-09-26: reusing a stale
# factorization diverged except on the easiest problems.)
#
# Reads the per-iteration, per-stage captures written by DDP4OPF's periodic
# capture hook (FILTERDDP_PERIODIC_CAPTURE_DIR, iterNNNN_stageSSS.jls) and, for
# each stage, compares every capture with the previous one of the same stage.
#
#   K = [H cu'; cu 0]
#   H block, per variable group (ps qs P Q v ell pb qnorm soc_slack energy_slack),
#     and its diagonal split into
#       barrier    Sigma_L + Sigma_U           (interior-point bound terms)
#       Vxx        dt^2 * diag(V_xx incoming)  (value-function curvature, pb only)
#       rest       diag(H) - barrier - Vxx     (constraint curvature from the
#                  SOCP multipliers, the constant cost terms, the diagonal floor)
#     plus the off-diagonal part of H (nonzero only with the exact Hessian).
#   J = cu, per constraint-row group (P balance, Q balance, voltage drop, SOCP,
#     root voltage, energy). Every group but SOCP is linear, so it cannot move.
#
# Output is long format, one row per (stage, iteration pair, block):
#   rel_change   ||new - old||_F / ||old||_F
#   share_dK2    that block's share of ||K_new - K_old||_F^2
#   max_log10    largest |log10(new/old)| over entries nonzero in both
#   frac_10x     share of those entries that changed by more than 10x
#
# With a done-file argument it follows a capture directory while FilterDDP is
# still writing it: a file is processed once its size has been stable for a few
# seconds, or once the done-file exists, and (KKT_EVOL_DELETE=1) deleted after,
# so large systems never hold more than a few snapshots on disk.
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_evolution_analysis.jl \
#         <capture_dir> <network_data.jls> <system> <out.csv> [done_file]

using Serialization, SparseArrays, LinearAlgebra, Printf, Statistics

capdir, datafile, system, outcsv = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
donefile = length(ARGS) >= 5 ? ARGS[5] : ""
delete_after = get(ENV, "KKT_EVOL_DELETE", "0") == "1"

d = deserialize(datafile)
N, L, B, D = length(d[:Nset]), length(d[:Lset]), length(d[:Bset]), length(d[:Dset])
dt = d[:delta_t_h]
# Control layout of ieee123c_filterddp.jl (control_layout) and its constraint rows.
function ranges(sizes)
    k = 0; out = Pair{String, UnitRange{Int}}[]
    for (name, n) in sizes
        push!(out, name => (k+1):(k+n)); k += n
    end
    return out, k
end
ugroups, nu_expect = ranges(["ps" => 1, "qs" => 1, "P" => L, "Q" => L, "v" => N, "ell" => L,
                             "pb" => B, "qnorm" => D, "soc_slack" => L, "energy_slack" => B])
cgroups, nc_expect = ranges(["Pbal" => N, "Qbal" => N, "vdrop" => L, "SOCP" => L, "rootV" => 1,
                             "energy" => B])
pbrange = last(ugroups[findfirst(p -> first(p) == "pb", ugroups)])

relchg(new, old) = (n = norm(old); n == 0 ? (norm(new) == 0 ? 0.0 : Inf) : norm(new - old) / n)
function ratio_stats(new, old)
    both = findall(i -> old[i] != 0 && new[i] != 0, eachindex(old))
    isempty(both) && return (NaN, NaN)
    r = abs.(log10.(abs.(new[both]) ./ abs.(old[both])))
    return (maximum(r), count(>(1.0), r) / length(r))
end

function slim(c)
    nu, nc = c.nu, c.nc
    (nu == nu_expect && nc == nc_expect) || error("layout mismatch: nu=$nu/$nu_expect nc=$nc/$nc_expect")
    K = c.K
    H = K[1:nu, 1:nu]
    Hd = Vector(diag(H))
    Hoff = H - spdiagm(0 => Hd)
    bar = c.Sigma_L .+ c.Sigma_U
    vxx = zeros(nu); vxx[pbrange] .= dt^2 .* diag(c.Vxx_incoming)
    return (iter = c.iteration, stage = c.stage, mu = c.barrier_mu, K = K, Hd = Hd, Hoff = Hoff,
            J = K[nu+1:end, 1:nu], bar = bar, vxx = vxx, rest = Hd .- bar .- vxx)
end

io = open(outcsv, isfile(outcsv) ? "a" : "w")
filesize(outcsv) == 0 && println(io,
    "system,stage,iter_prev,iter,mu_prev,mu,block,rel_change,share_dK2,max_log10,frac_10x")
function row(s, o, block, new, old, dK2)
    rc = relchg(new, old)
    part = norm(new - old)^2
    ml, f10 = new isa AbstractVector ? ratio_stats(new, old) :
              ratio_stats(nonzeros(new), nonzeros(old))
    @printf(io, "%s,%d,%d,%d,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e\n", system, s.stage, o.iter, s.iter,
            o.mu, s.mu, block, rc, dK2 == 0 ? 0.0 : part / dK2, ml, f10)
end

function compare(s, o)
    dK2 = norm(s.K - o.K)^2
    row(s, o, "K", s.K, o.K, dK2)
    row(s, o, "H.diag", s.Hd, o.Hd, dK2)
    nnz(s.Hoff) + nnz(o.Hoff) > 0 && row(s, o, "H.offdiag", s.Hoff, o.Hoff, dK2)
    # The Jacobian appears twice in K (cu and cu'), hence 2x its squared change.
    Jd2 = 2 * norm(s.J - o.J)^2
    @printf(io, "%s,%d,%d,%d,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e\n", system, s.stage, o.iter, s.iter,
            o.mu, s.mu, "J", relchg(s.J, o.J), dK2 == 0 ? 0.0 : Jd2 / dK2,
            ratio_stats(nonzeros(s.J), nonzeros(o.J))...)
    for (name, part) in (("H.barrier", :bar), ("H.Vxx", :vxx), ("H.rest", :rest))
        row(s, o, name, getfield(s, part), getfield(o, part), dK2)
    end
    for (g, r) in ugroups
        row(s, o, "H.diag:" * g, s.Hd[r], o.Hd[r], dK2)
        any(!iszero, o.bar[r]) && row(s, o, "H.barrier:" * g, s.bar[r], o.bar[r], dK2)
    end
    for (g, r) in cgroups
        Jn, Jo = s.J[r, :], o.J[r, :]
        rc = relchg(Jn, Jo)
        ml, f10 = ratio_stats(nonzeros(Jn), nonzeros(Jo))
        @printf(io, "%s,%d,%d,%d,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e\n", system, s.stage, o.iter, s.iter,
                o.mu, s.mu, "J:" * g, rc, dK2 == 0 ? 0.0 : 2 * norm(Jn - Jo)^2 / dK2, ml, f10)
    end
    flush(io)
end

prev = Dict{Int, Any}()
done = Set{String}()
sizes = Dict{String, Int}()
parse_name(f) = (m = match(r"iter(\d+)_stage(\d+)\.jls$", f); m === nothing ? nothing :
                 (parse(Int, m[1]), parse(Int, m[2])))
nprocessed = 0
while true
    finished = !isempty(donefile) && isfile(donefile)
    files = sort([f for f in readdir(capdir) if parse_name(f) !== nothing && !(f in done)],
                 by = f -> parse_name(f))
    ready = String[]
    for f in files
        sz = filesize(joinpath(capdir, f))
        if finished || isempty(donefile) || (get(sizes, f, -1) == sz && sz > 0 &&
                                             time() - mtime(joinpath(capdir, f)) > 5)
            push!(ready, f)
        end
        sizes[f] = sz
    end
    for f in ready
        path = joinpath(capdir, f)
        s = slim(deserialize(path))
        haskey(prev, s.stage) && prev[s.stage].iter < s.iter && compare(s, prev[s.stage])
        prev[s.stage] = s
        push!(done, f); global nprocessed += 1
        delete_after && rm(path)
    end
    (isempty(donefile) || finished) && isempty(setdiff(Set(files), Set(ready))) && break
    isempty(ready) && sleep(3)
end
close(io)
println("KKT_EVOLUTION processed=$nprocessed captures -> $outcsv")
