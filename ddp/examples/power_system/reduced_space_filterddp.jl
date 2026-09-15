# Reduced-space MPOPF solved by the REAL FilterDDP solver.
#
# The network is eliminated: FilterDDP optimises only over battery quantities,
# and every algebraic network variable is recovered by an inner single-period
# OPF that runs inside the stage-cost callback.
#
#   state    x = B^t                    (nx = nB)
#   control  u = [P_B^t ; s_E^t]        (nu = 2 nB)
#   dynamics B^t = B^{t-1} - dt * P_B^t              (linear)
#   stage    l_t(x,u) = Phi_t(P_B^t) + C_B * S^2 * dt * sum (P_B^t)^2
#   stage eq c_t(x,u) = x - dt*P_B - Bmin - s_E = 0  (SOC box via slack)
#
# where Phi_t(P_B) is the optimal value of the inner network OPF at fixed P_B.
# Its gradient is FREE: the real-power balance duals are exactly dPhi/dP_B
# (verified to ~1e-9 in probe_reduced_value_function.jl, and re-verified here).
#
# NOTHING IN THE PRODUCTION SOLVER IS MODIFIED. FilterDDP's Objective/Dynamics/
# EqualityConstraints are plain structs over six callables; the Symbolics-based
# constructors are only a convenience. We build the structs directly from
# hand-written closures, which is also the only option available -- an Ipopt
# solve is not automatically differentiable, so the derivatives have to be
# supplied analytically (gradient) and by finite differences (curvature).
#
# CURVATURE is the one genuinely approximate ingredient, so it is a switch:
#
#   hessian=exact         d2Phi/dP_B2 by forward differences on the DUALS.
#                         Costs nB extra inner solves per stage per backward
#                         pass. Affordable only for small nB.
#   hessian=frozen        compute d2Phi ONCE per stage at the initial iterate
#                         and reuse it for every later iteration. Justified by
#                         measurement, not convenience: on ieee123 the Hessian
#                         moves 1.1% between the initial point and the optimum,
#                         6.3% at a 70%-charge outlier and 2.6% between stages.
#                         Costs nB solves per stage TOTAL rather than per
#                         backward pass.
#   hessian=battery_only  drop d2Phi entirely, keep only the battery term's
#                         2*C_B*S^2*dt*I. Zero extra solves. Since Phi is convex
#                         this UNDER-estimates curvature, which is the safe
#                         direction for a trust-region/filter method but may
#                         cost iterations. Whether it converges decides whether
#                         the approach scales to ieee2522 and large10k.
#
# Infeasible trial dispatches are NOT silently relaxed: the hard inner solve is
# tried first, and only if it fails does the adaptive L1 penalty supply a finite
# value and gradient so the outer method can retreat. Every such event is
# counted and reported -- a run with a nonzero count did NOT solve the original
# problem in the ordinary sense, and says so.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#     ddp/examples/power_system/reduced_space_filterddp.jl [system] [T] [hessian]

using FilterDDP
using LinearAlgebra
using Random
using Printf
using Serialization
using SparseArrays

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))
include(joinpath(@__DIR__, "inner_opf_persistent.jl"))
include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))   # control_layout, for the reference

# ---------------------------------------------------------------- Phi cache --
# FilterDDP evaluates l, lu and luu at the same iterate repeatedly within one
# iteration. Without memoisation every one of those is a fresh Ipopt solve.
mutable struct PhiStats
    inner_solves::Int
    hessian_solves::Int
    cache_hits::Int
    infeasible_events::Int
    inner_time::Float64
end
PhiStats() = PhiStats(0, 0, 0, 0, 0.0)

struct PhiCache
    data::Dict
    t::Int
    pin::Union{Nothing,PersistentInner}
    stats::PhiStats
    slots::Vector{Tuple{Vector{Float64},NamedTuple}}   # small exact-match cache
    hslots::Vector{Tuple{Vector{Float64},Matrix{Float64}}}
    frozen::Base.RefValue{Union{Nothing,Matrix{Float64}}}
    maxslots::Int
end
PhiCache(data, t, stats; maxslots=6, persistent::Bool=true) =
    PhiCache(data, t, persistent ? persistent_inner(data, t) : nothing, stats, Tuple{Vector{Float64},NamedTuple}[],
             Tuple{Vector{Float64},Matrix{Float64}}[],
             Ref{Union{Nothing,Matrix{Float64}}}(nothing), maxslots)

"""
Value and gradient of `Phi_t` at `pb`. Hard solve first; adaptive penalty only
as a fallback, and the fallback is counted.
"""
function phi(cache::PhiCache, pb::Vector{Float64})
    for (k, v) in cache.slots
        k == pb && (cache.stats.cache_hits += 1; return v)
    end
    t0 = time()
    r = cache.pin === nothing ? inner_opf(cache.data, cache.t, pb) :
                                solve_at!(cache.pin, pb)
    cache.stats.inner_solves += 1
    out = if r.feasible
        (; value = r.substation_cost, grad = copy(r.lambda_bal),
           feasible = true, violation = 0.0)
    else
        cache.stats.infeasible_events += 1
        a = inner_opf_adaptive(cache.data, cache.t, pb)
        cache.stats.inner_solves += a.rounds
        if a.res.feasible
            (; value = a.res.substation_cost + a.rho * a.res.total_violation,
               grad = copy(a.res.lambda_bal), feasible = false,
               violation = a.res.total_violation)
        else
            (; value = 1e12, grad = zeros(length(pb)), feasible = false,
               violation = Inf)
        end
    end
    cache.stats.inner_time += time() - t0
    push!(cache.slots, (copy(pb), out))
    length(cache.slots) > cache.maxslots && popfirst!(cache.slots)
    return out
end

"""
Randomised low-rank Hessian (Nystrom). `H*d` for ANY direction costs exactly one
inner solve -- perturb along `d`, difference the gradients -- so a randomised
range finder applies directly and needs `rank+oversample` solves instead of `nB`.

Nystrom rather than a generic randomised SVD because it needs only the one
sketch `Y = H*Omega` (no second pass), and it is valid precisely because `Phi` is
convex: every eigenvalue of the measured Hessian is positive (1.19 .. 2482 on
ieee123), so `H` is PSD and `H ~ Y (Omega'Y)^-1 Y'` is the right construction.

The step `h` is larger here than for coordinate differences on purpose. A
unit-norm random direction spreads the perturbation over all `nB` coordinates,
so each one moves by only `h/sqrt(nB)`; at `h = 1e-6` on large10k that is `3e-8`
per coordinate, below the level at which the duals are trustworthy. Since only
~6% accuracy is needed (measured: `frozen` tolerates it), trading truncation
error for noise is the right direction.
"""
function phi_hessian_lowrank(cache::PhiCache, pb::Vector{Float64};
                             rank::Int=25, oversample::Int=10,
                             h::Float64=1e-4, seed::Int=20260914,
                             floor_eigs::Bool=true)
    n = length(pb)
    l = min(rank + oversample, n)
    rng = MersenneTwister(seed + cache.t)
    g0 = phi(cache, pb).grad
    Om = zeros(n, l); Y = zeros(n, l)
    t0 = time()
    for i in 1:l
        d = randn(rng, n); d ./= norm(d)
        Om[:, i] = d
        pert = pb .+ h .* d
        r = cache.pin === nothing ? inner_opf(cache.data, cache.t, pert) :
                                    solve_at!(cache.pin, pert)
        cache.stats.hessian_solves += 1
        Y[:, i] = r.feasible ? (r.lambda_bal .- g0) ./ h : zeros(n)
    end
    cache.stats.inner_time += time() - t0
    C = Symmetric(0.5 .* (Om' * Y .+ (Om' * Y)'))
    H = Y * (pinv(Matrix(C), 1e-8) * Y')
    H = 0.5 .* (H .+ H')

    # Eigenvalue FLOOR. Without it the truncated directions are handed curvature
    # ~0, and since the battery term contributes only 2.24 against d2Phi
    # eigenvalues reaching 2482 (C_B = 1.4e-07 in the network cases, the tADMM
    # regime -- d2Phi IS the curvature model here), there is nothing to bound the
    # Newton step in exactly the directions the sketch never measured. That is
    # what makes plain Nystrom stall rather than merely slow down. Flooring at
    # the smallest CAPTURED eigenvalue is deliberately conservative: it
    # over-states curvature in the unexplored subspace, which shortens steps
    # instead of letting them run away.
    if floor_eigs
        E = eigen(Symmetric(H))
        pos = filter(>(0.0), E.values)
        sigma = isempty(pos) ? 1.0 : minimum(pos)
        return E.vectors * Diagonal(max.(E.values, sigma)) * E.vectors'
    end
    return H
end

"""
`d2Phi/dP_B2` by forward differences on the analytic gradient: nB extra solves,
then symmetrised. Forward (not central) differences because the gradient is
already exact to solver tolerance, so the error is dominated by the step, and
nB solves is already the affordability ceiling.
"""
function phi_hessian(cache::PhiCache, pb::Vector{Float64}; h::Float64=1e-6)
    for (k, v) in cache.hslots
        k == pb && (cache.stats.cache_hits += 1; return v)
    end
    n = length(pb)
    g0 = phi(cache, pb).grad
    H = zeros(n, n)
    t0 = time()
    pert = copy(pb)
    for b in 1:n
        pert[b] = pb[b] + h
        r = cache.pin === nothing ? inner_opf(cache.data, cache.t, pert) :
                                    solve_at!(cache.pin, pert)
        cache.stats.hessian_solves += 1
        H[:, b] = r.feasible ? (r.lambda_bal .- g0) ./ h : zeros(n)
        pert[b] = pb[b]
    end
    cache.stats.inner_time += time() - t0
    H = 0.5 .* (H .+ H')       # Phi is C^2 where smooth; symmetrise FD noise
    push!(cache.hslots, (copy(pb), H))
    length(cache.hslots) > 2 && popfirst!(cache.hslots)
    return H
end

# ------------------------------------------------------- reduced-space OCP --
function reduced_stage_objective(data::Dict, t::Int, nB::Int, stats::PhiStats,
                                 hessian_mode::Symbol)
    nx = nB; nu = 2nB
    dt = data[:delta_t_h]; pbase = data[:kVA_B]
    cb = data[:C_B] * pbase^2 * dt
    cache = PhiCache(data, t, stats;
                     persistent = get(ENV, "REDUCED_PERSISTENT", "1") == "1")

    l = function (x, u)
        pb = u[1:nB]
        [phi(cache, pb).value + cb * sum(abs2, pb)]
    end
    lx = (x, u) -> zeros(nx)
    lu = function (x, u)
        pb = u[1:nB]
        g = zeros(nu)
        g[1:nB] .= phi(cache, pb).grad .+ 2cb .* pb
        g
    end
    lxx = (x, u) -> zeros(nx, nx)
    lux = (x, u) -> spzeros(nu, nx)
    luu = function (x, u)
        H = zeros(nu, nu)
        pb = u[1:nB]
        if hessian_mode === :exact
            H[1:nB, 1:nB] .= phi_hessian(cache, pb)
        elseif hessian_mode === :frozen
            if cache.frozen[] === nothing
                cache.frozen[] = phi_hessian(cache, pb)
            end
            H[1:nB, 1:nB] .= cache.frozen[]
        elseif hessian_mode === :lowrank
            if cache.frozen[] === nothing
                cache.frozen[] = phi_hessian_lowrank(cache, pb;
                    rank = parse(Int, get(ENV, "REDUCED_HESS_RANK", "25")),
                    oversample = parse(Int, get(ENV, "REDUCED_HESS_OVERSAMPLE", "10")),
                    h = parse(Float64, get(ENV, "REDUCED_HESS_STEP", "1e-4")),
                    floor_eigs = get(ENV, "REDUCED_HESS_FLOOR", "1") == "1")
            end
            H[1:nB, 1:nB] .= cache.frozen[]
        end
        for b in 1:nB
            H[b, b] += 2cb
        end
        H
    end
    return FilterDDP.Objective{nx,nu,typeof(l),typeof(lx),typeof(lu),
                               typeof(lxx),typeof(lux),typeof(luu)}(
        l, lx, lu, lxx, lux, luu), cache
end

function build_reduced_ocp(data::Dict; hessian_mode::Symbol=:exact)
    Bset = data[:Bset]; nB = length(Bset)
    nx = nB; nu = 2nB; nc = nB
    T = size(data[:p_L_pu], 2)
    dt = data[:delta_t_h]

    stats = PhiStats()
    stage_objs = Any[]; caches = PhiCache[]
    for t in 1:T
        o, c = reduced_stage_objective(data, t, nB, stats, hessian_mode)
        push!(stage_objs, o); push!(caches, c)
    end

    # dynamics: B^t = B^{t-1} - dt * P_B^t   (linear -> zero second derivatives)
    f  = (x, u) -> x .- dt .* u[1:nB]
    fx = (x, u) -> Matrix{Float64}(I, nx, nx)
    fu = function (x, u)
        J = zeros(nx, nu)
        for b in 1:nB; J[b, b] = -dt; end
        J
    end
    fxx = (x, u, w) -> zeros(nx, nx)
    fux = (x, u, w) -> spzeros(nu, nx)
    fuu = (x, u, w) -> spzeros(nu, nu)
    dyn = FilterDDP.Dynamics{nx,nu,typeof(f),typeof(fx),typeof(fu),
                             typeof(fxx),typeof(fux),typeof(fuu)}(f, fx, fu, fxx, fux, fuu)

    # stage equality: x - dt*P_B - Bmin - s_E = 0  (same SOC box as full space)
    emin = Float64[data[:soc_min][j] * data[:B_R_pu][j] for j in Bset]
    c  = (x, u) -> [x[b] - dt * u[b] - emin[b] - u[nB + b] for b in 1:nB]
    cx = function (x, u)
        J = spzeros(nc, nx); for b in 1:nB; J[b, b] = 1.0; end; J
    end
    cu = function (x, u)
        J = spzeros(nc, nu)
        for b in 1:nB; J[b, b] = -dt; J[b, nB + b] = -1.0; end
        J
    end
    cxx = (x, u, p) -> zeros(nx, nx)
    cux = (x, u, p) -> spzeros(nu, nx)
    cuu = (x, u, p) -> spzeros(nu, nu)
    con = FilterDDP.EqualityConstraints{nx,nu,nc,typeof(c),typeof(cx),typeof(cu),
                                        typeof(cxx),typeof(cux),typeof(cuu)}(
        c, cx, cu, cxx, cux, cuu)

    lower = fill(-Inf, nu); upper = fill(Inf, nu)
    for (b, j) in enumerate(Bset)
        lower[b] = -data[:P_B_R_pu][j]; upper[b] = data[:P_B_R_pu][j]
        width = (data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j]
        lower[nB + b] = 0.0; upper[nB + b] = width
    end
    limits = ControlLimits(lower, upper)

    ocp = build_ocp(T, stage_objs[1], stage_objs[end], dyn, con, limits;
                    stage_objectives = stage_objs, stage_constraints = nothing)
    return (; ocp, nx, nu, nc, nB, T, stats, caches, emin)
end

# ----------------------------------------------------------------- driver ---
function main(args = ARGS)
    system = length(args) >= 1 ? args[1] : "ieee123C_1ph"
    T      = length(args) >= 2 ? parse(Int, args[2]) : 3
    hmode  = Symbol(length(args) >= 3 ? args[3] : "exact")

    datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                        "network_data_$(system)_T$(T).jls")
    data = deserialize(datafile)
    # Opt-in C_B override. C_B is the battery cycling cost and it sets how much
    # perfectly-conditioned damping (2*C_B*S^2*dt*I) sits under d2Phi in the
    # curvature model, so it governs how forgiving the problem is of curvature
    # error. Changing it changes the PROBLEM, so any comparison must re-run the
    # full-space reference at the same value.
    if haskey(ENV, "REDUCED_CB")
        data[:C_B] = parse(Float64, ENV["REDUCED_CB"])
        @printf("C_B OVERRIDE: %.6g (battery diagonal 2*C_B*S^2*dt = %.6g)
",
                data[:C_B], 2 * data[:C_B] * data[:kVA_B]^2 * data[:delta_t_h])
    end
    dt = data[:delta_t_h]; pbase = data[:kVA_B]

    t_build = time()
    R = build_reduced_ocp(data; hessian_mode = hmode)
    @printf("reduced OCP: nx=%d nu=%d nc=%d  (full-space nu would be %d)\n",
            R.nx, R.nu, R.nc,
            2 + 4length(data[:Lset]) + length(data[:Nset]) + length(data[:Dset]) + 2R.nB)
    @printf("build: %.2f s, hessian=%s\n", time() - t_build, hmode)

    max_iter = parse(Int, get(ENV, "FILTERDDP_MAX_ITERATIONS", "200"))
    opt_tol  = parse(Float64, get(ENV, "FILTERDDP_OPTIMALITY_TOLERANCE", "1e-7"))
    solver = Solver(R.ocp; options = Options{Float64}(
        verbose = get(ENV, "FILTERDDP_QUIET", "0") != "1",
        optimality_tolerance = opt_tol, max_iterations = max_iter))

    x0 = Float64[data[:B0_pu][j] for j in data[:Bset]]
    u0 = zeros(R.nu)
    for b in 1:R.nB; u0[R.nB + b] = x0[b] - R.emin[b]; end
    ubar = [copy(u0) for _ in 1:T]

    t1 = time()
    status = solve!(solver, x0, ubar)
    wall = time() - t1
    @printf("\nsolve: %.2f s, iterations=%d, status=%s\n", wall, solver.data.k, string(status))
    @printf("residuals: primal=%.6e dual=%.6e complementarity=%.6e\n",
            solver.data.primal_inf, solver.data.dual_inf, solver.data.cs_inf_0)
    s = R.stats
    @printf("inner solves=%d  hessian solves=%d  cache hits=%d  inner time=%.1f s (%.0f%% of wall)\n",
            s.inner_solves, s.hessian_solves, s.cache_hits, s.inner_time,
            100 * s.inner_time / max(wall, 1e-9))
    @printf("INFEASIBLE trial dispatches (penalty fallback used): %d\n", s.infeasible_events)

    xr, ur = get_trajectory(solver)
    cb = data[:C_B] * pbase^2 * dt
    obj_red = 0.0
    for t in 1:T
        pb = ur[t][1:R.nB]
        r = inner_opf(data, t, pb)
        obj_red += (r.feasible ? r.substation_cost : NaN) + cb * sum(abs2, pb)
    end
    @printf("reduced-space objective: %.10f USD\n", obj_red)

    # ---- compare against the full-space FilterDDP solution, when present ----
    reffile = joinpath(REPO, "ddp", "results", "network_filterddp",
                       "filterddp_solution_$(system)_T$(T).jls")
    if isfile(reffile)
        ref = deserialize(reffile)
        idx, _ = control_layout(data)
        obj_full = 0.0
        dpb = 0.0; dB = 0.0
        for t in 1:T
            pbf = ref[:u][t][idx.pb]
            obj_full += data[:LoadShapeCost][t] * pbase * dt * ref[:u][t][idx.ps] +
                        cb * sum(abs2, pbf)
            dpb = max(dpb, maximum(abs, ur[t][1:R.nB] .- pbf))
            dB  = max(dB,  maximum(abs, xr[t] .- ref[:x][t]))
        end
        @printf("full-space  objective: %.10f USD  (%d iterations)\n",
                obj_full, ref[:iterations])
        @printf("objective gap: %.6e USD (rel %.3e)\n",
                abs(obj_red - obj_full), abs(obj_red - obj_full) / abs(obj_full))
        @printf("max |dP_B| vs full space: %.6e pu (%.4f kW)\n", dpb, dpb * pbase)
        @printf("max |dB|   vs full space: %.6e pu\n", dB)
    else
        println("no full-space reference at $(basename(reffile)); objective reported unchecked")
    end

    outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
    serialize(joinpath(outdir, "reduced_filterddp_$(system)_T$(T)_$(hmode).jls"),
              Dict(:system => system, :T => T, :hessian => string(hmode),
                   :status => string(status), :iterations => solver.data.k,
                   :wall => wall, :objective => obj_red,
                   :inner_solves => s.inner_solves, :hessian_solves => s.hessian_solves,
                   :infeasible_events => s.infeasible_events,
                   :x => xr, :u => ur))
    return status
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
