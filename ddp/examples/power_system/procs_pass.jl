# One process per period: the pairwise backward pass of paired_pass.jl and the
# forward pass of threaded_pass.jl with every period in a Julia process of its
# own, on one machine.
#
# The one-worker-per-period figures so far are a count taken from a sequential
# run (FILTERDDP_PAIRED_SWEEP=count), and the threaded pass shares one heap, so
# a collection stops every period. Here a period has its own memory and its own
# collector, and the pass is timed by a clock.
#
#   coordinator (process 1): the solver's loop, the filter and the barrier
#       update, as in any other run. It holds no stage data and does no stage
#       work; it sends each period the scalars of the pass and adds up what
#       comes back.
#   worker of period t: builds that period's model only, keeps its trajectory
#       element, its factors and its stored policy, and runs steps 1 and 3 of
#       the pass, its combinations in the rounds, and its part of the forward
#       pass. The routines are those of threaded_pass.jl and paired_pass.jl.
#
# What passes between processes:
#   - the stages' functions G_t (three n_B x n_B blocks and two vectors each)
#     live in shared memory, one segment per period, mapped by every process.
#     A combination reads its right operand there and writes over its left
#     operand, as with FILTERDDP_PAIRED_INPLACE=1; in a round the left operands
#     are distinct and none of them is a right operand. One machine only.
#   - step 1 returns the stage's scalars and two n_B vectors (for the multiplier
#     recursion of the dual residual, done by the coordinator); step 3 receives
#     one n_B vector and returns scalars.
#   - forward pass: the battery states of all periods come from a chain of
#     T-1 messages of n_B numbers (the recursion of threaded_affine_directions,
#     one hop per period), then every period computes its directions and ratio
#     test at once, and every trial point at once. Only scalars come back.
#   - the trajectory itself comes back once, when the solve stops.
#
# Differences from the counted run:
#   - no factor reuse after a barrier update (FILTERDDP_REUSE_BARRIER): that
#     pass is sequential by construction, and here a full pass at the new
#     barrier parameter is one parallel pass. Same linear systems.
#   - the states of the forward pass come from the recursion, as in the
#     threaded pass, so the iterates agree with the sequential run's to
#     rounding, not bit for bit.
#
# A failed stage (a configuration outside what _tp_pre supports, a failed
# factorization, a non-finite number) stops the solve with status 1: the
# coordinator has no trajectory to fall back on.
#
# Fewer workers than periods (FILTERDDP_PROCS_WORKERS): a worker then owns
# several periods and does them one after another.

using Distributed
using SharedArrays

isdefined(@__MODULE__, :paired_backward_pass!) || include(joinpath(@__DIR__, "paired_pass.jl"))

# ------------------------------------------------------------------ worker --

mutable struct ProcsWorker
    solver::Any
    ocp::Any
    data::Any
    idx::Any
    nx::Int
    nu::Int
    nc::Int
    S::Dict{Int, Any}            # phase-1 record of each owned period
    dir::Dict{Int, Any}          # directions of each owned period
    shared::Vector{Any}          # keeps the shared segments mapped
end
const _PW = Ref{Any}(nothing)

_pw_gc() = Base.gc_num().total_time / 1e9

function procs_worker_setup(system::String, T::Int, own::Vector{Int}, shared::Vector)
    BLAS.set_num_threads(1); _BLAS_DEFAULT[] = 1
    ENV["FILTERDDP_PAIRED_INPLACE"] = "1"
    t0 = time()
    ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
    data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                                "network_data_$(system)_T$(T)$(ptag).jls"))
    haskey(ENV, "REDUCED_CB") && (data[:C_B] = battery_cb(system, ENV["REDUCED_CB"]))
    gammaT = terminal_soc_soft() ? gamma_terminal(system) : 0.0
    built = redirect_stdout(devnull) do
        build_model(data; gamma=gammaT, stages=own)
    end
    ocp, idx = built[1], built[2]
    solver = Solver(ocp; options=Options{Float64}())
    nx, nu, nc = typeof(solver).parameters[2], typeof(solver).parameters[3], typeof(solver).parameters[4]
    # storage of the periods of other workers: not needed here
    none = DDP4OPF.TrajectoryElement{Float64, nx, nu, nc}(Float64[], Float64[], Float64[], Float64[], Float64[])
    for t in 1:T
        t in own && continue
        solver.nominal[t] = none; solver.current[t] = none
        r = solver.update[t]
        r.α = Float64[]; r.ψ = Float64[]; r.χl = Float64[]; r.χu = Float64[]; r.Σ_L = Float64[]; r.Σ_U = Float64[]
    end
    # every period's function, in the shared segments
    n = nx
    empty!(_PP_SLOTS); empty!(_PP_WORK)
    for t in 1:T
        p = pointer(sdata(shared[t]))
        blk = i -> unsafe_wrap(Array, p + (i - 1) * n * n * sizeof(Float64), (n, n))
        vec = i -> unsafe_wrap(Array, p + (3n * n + (i - 1) * n) * sizeof(Float64), (n,))
        push!(_PP_SLOTS, PintElem(blk(1), blk(2), blk(3), vec(1), vec(2)))
    end
    push!(_PP_WORK, PairWork(n))
    _PP_END[] = PintVal(zeros(n, n), zeros(n))
    _TP_LAYOUT[] = nothing
    _TP_LAST_DIAG[] = nothing
    _PW[] = ProcsWorker(solver, ocp, data, idx, nx, nu, nc, Dict{Int, Any}(), Dict{Int, Any}(), collect(Any, shared))
    return (setup_s=time() - t0, rss_mib=Sys.maxrss() / 2^20)
end

# Step 1 of period t. elem: the period's trajectory element at the start of a
# solve; advance: the last trial point was accepted.
function pw_step1(t::Int, k::Int, μ::Float64, elem, advance::Bool)
    W = _PW[]; s = W.solver
    t0 = time_ns(); g0 = _pw_gc()
    isnothing(elem) || (s.nominal[t] = elem)
    advance && (s.nominal[t] = s.current[t])
    S = _tp_pre(s, W.ocp, s.nominal, t, k, μ, false, W.nx, W.nu, W.nc, W.data, W.idx, true)
    isnothing(S) && return nothing
    isnothing(_pp_elem(S, t, W.nx)) && return nothing
    W.S[t] = S
    return (lval=S.lval, c1=S.c1, cinf=S.cinf, logsum=S.logsum, cs0=S.cs0, csμ=S.csμ, cϕ=S.cϕ,
            lx=Vector{Float64}(S.lx), cxϕ=Vector{Float64}(S.cxϕ),
            el=(time_ns() - t0) / 1e9, gc=_pw_gc() - g0)
end

# One combination: period t's function with the one that starts at period mid
# (yval: that one is already a value function; last: it is the zero value
# function after the last period).
function pw_combine(t::Int, mid::Int, yval::Bool, last::Bool)
    t0 = time_ns(); g0 = _pw_gc()
    x = _PP_SLOTS[t]
    if last
        _pp_combine!(x, _PP_END[], _PP_WORK[1])
    elseif yval
        _pp_combine!(x, PintVal(_PP_SLOTS[mid].Maa, _PP_SLOTS[mid].ma), _PP_WORK[1])
    else
        _pp_combine!(x, _PP_SLOTS[mid], _PP_WORK[1])
    end
    return (el=(time_ns() - t0) / 1e9, gc=_pw_gc() - g0)
end

# Step 3 of period t, with the value function of period t+1 (in its segment).
function pw_step3(t::Int, μ::Float64, λ::Vector{Float64}, last::Bool)
    W = _PW[]; s = W.solver
    t0 = time_ns(); g0 = _pw_gc()
    S = W.S[t]
    Vn = last ? _PP_END[] : PintVal(_PP_SLOTS[t+1].Maa, _PP_SLOTS[t+1].ma)
    R = _pp_solve(S, Vn.v, Vn.V)
    isnothing(R) && return nothing
    P = _tp_post(s, S, R, t, μ, λ, W.nx, W.nu, W.nc, false, true, W.ocp.control_limits)
    isnothing(P) && return nothing
    return (dual=P.dual, zl_sum=P.zl_sum, zu_sum=P.zu_sum, ϕ1=P.ϕ1, e1=P.e1, e2=P.e2,
            el=(time_ns() - t0) / 1e9, gc=_pw_gc() - g0)
end

# The state after period t from the state before it, by the battery rows of the
# feedback (the recursion of threaded_affine_directions).
function pw_state_step(t::Int, x::Vector{Float64})
    W = _PW[]; s = W.solver
    nom = s.nominal[t]; rule = s.update[t]; pol = rule.factor_policy
    (isnothing(pol) || isempty(pol.βB)) && error("period $t has no stored battery rows")
    rows = pol.active_rows
    ur = copy(nom.u)
    @views ur[rows] .+= rule.α[rows] .+ pol.βB * (x .- nom.x)
    return Vector{Float64}(W.ocp.dynamics.f(x, ur))
end

# Period t's directions and its ratio test (the loop body of threaded_affine_directions).
function pw_directions(t::Int, x::Vector{Float64}, τ::Float64, split::Bool)
    W = _PW[]; s = W.solver; cl = W.ocp.control_limits
    t0 = time_ns(); g0 = _pw_gc()
    nom = s.nominal[t]; rule = s.update[t]
    δx = x - nom.x
    βδx, ωδx = DDP4OPF.policy_actions!(rule, δx)
    du = rule.α .+ βδx
    dϕ = rule.ψ .+ ωδx
    dzl = rule.χl .- rule.Σ_L .* βδx
    dzu = rule.χu .+ rule.Σ_U .* βδx
    γ = 1.0; γd = 1.0
    @inbounds for i in eachindex(du)
        if cl.maskl[i]
            du[i] < 0 && (γ = min(γ, τ * (nom.u[i] - cl.l[i]) / -du[i]))
            if dzl[i] < 0
                g = τ * nom.zl[i] / -dzl[i]
                split ? (γd = min(γd, g)) : (γ = min(γ, g))
            end
        end
        if cl.masku[i]
            du[i] > 0 && (γ = min(γ, τ * (cl.u[i] - nom.u[i]) / du[i]))
            if dzu[i] < 0
                g = τ * nom.zu[i] / -dzu[i]
                split ? (γd = min(γd, g)) : (γ = min(γ, g))
            end
        end
    end
    W.dir[t] = (dx=δx, du=du, dϕ=dϕ, dzl=dzl, dzu=dzu)
    return (γ=γ, γd=γd, el=(time_ns() - t0) / 1e9, gc=_pw_gc() - g0)
end

# Period t's part of a trial point (the loop body of threaded_affine_rollout!).
function pw_rollout(t::Int, γ::Float64, γz::Float64, μ::Float64)
    W = _PW[]; s = W.solver; ocp = W.ocp; cl = ocp.control_limits
    t0 = time_ns(); g0 = _pw_gc()
    nom = s.nominal[t]; dir = W.dir[t]
    x = nom.x .+ γ .* dir.dx
    u = nom.u .+ γ .* dir.du
    ϕ = nom.ϕ .+ γ .* dir.dϕ
    zl = nom.zl .+ γz .* dir.dzl
    zu = nom.zu .+ γz .* dir.dzu
    s.current[t] = DDP4OPF.TrajectoryElement{Float64, W.nx, W.nu, W.nc}(x, u, ϕ, zl, zu)
    c = DDP4OPF.stage_con(ocp, t).c(x, u)
    c1 = norm(c, 1); cϕ = dot(c, ϕ)
    ul = u - cl.l
    uu = cl.u - u
    bad = any((ul .<= 0) .* cl.maskl) || any((uu .<= 0) .* cl.masku) ||
          any((zl .<= 0) .* cl.maskl) || any((zu .<= 0) .* cl.masku)
    logl = 0.0; logu = 0.0; lv = 0.0
    if !bad
        logl = dot(log.(ul), cl.maskl)
        logu = dot(log.(uu), cl.masku)
        lv = t == ocp.N ? ocp.term_objective.l(x, u)[1] : DDP4OPF.stage_obj(ocp, t).l(x, u)[1]
    end
    return (c1=c1, cϕ=cϕ, logl=logl, logu=logu, lv=lv, bad=bad, el=(time_ns() - t0) / 1e9, gc=_pw_gc() - g0)
end

function pw_nominal(t::Int, advance::Bool)
    s = _PW[].solver
    advance && (s.nominal[t] = s.current[t])
    return s.nominal[t]
end

pw_memory() = Sys.maxrss() / 2^20

# ------------------------------------------------------------- coordinator --

mutable struct ProcsMain
    pids::Vector{Int}
    owner::Vector{Int}           # the process of each period
    shared::Vector{Any}
    last_k::Int                  # iteration of the last backward pass
    # The state before the first period. solve! copies the coordinator's unused
    # `current` over its `nominal` after every accepted step, so the
    # coordinator's trajectory is not to be read between the start and the end.
    x1::Vector{Float64}
end
const _PM = Ref{Any}(nothing)

# f(t) on the process of every period t, all at once.
function _pm_all(f, M::ProcsMain, N::Int)
    out = Vector{Any}(undef, N)
    @sync for t in 1:N
        @async out[t] = f(t, M.owner[t])
    end
    return out
end

function procs_backward_pass!(solver::DDP4OPF.Solver{T, nx, nu, nc}, ocp, traj, d, options) where {T, nx, nu, nc}
    M = _PM[]::ProcsMain
    N = ocp.N
    μ = Float64(d.μ)
    cl = ocp.control_limits
    t_start = time_ns()
    fresh = d.k == 0 && d.j == 0                   # a new solve: the workers take the starting point
    advance = !fresh && d.k != M.last_k            # the last trial point was accepted
    M.last_k = d.k
    fresh && (M.x1 = Vector{Float64}(traj[1].x))

    # 1. every period: its stage data and its function of the two end states
    r1 = _pm_all((t, w) -> remotecall_fetch(pw_step1, w, t, d.k, μ, fresh ? traj[t] : nothing, advance), M, N)
    if any(isnothing, r1)
        d.status = 1
        @printf("FILTERDDP_PROCS_FAIL iteration=%d step=1 periods=%s\n", d.k, string(findall(isnothing, r1)))
        return true
    end
    t_1 = time_ns()

    # 2. the rounds
    isval = falses(N + 1); isval[N+1] = true
    rounds = 0; merges = 0; pair_sum = 0.0; pair_max = 0.0; pair_gc = 0.0
    n = N + 1; h = 1
    while h < n
        pairs = Tuple{Int, Int}[]
        for i in 1:2h:n
            mid = i + h
            mid > n && continue
            for t in i:mid-1; push!(pairs, (t, mid)); end
        end
        res = Vector{Any}(undef, length(pairs))
        @sync for (q, (t, mid)) in enumerate(pairs)
            @async res[q] = remotecall_fetch(pw_combine, M.owner[t], t, mid, isval[mid], mid == N + 1)
        end
        for (t, mid) in pairs; isval[t] = isval[mid]; end
        rounds += 1; merges += length(pairs)
        pair_sum += sum(r.el for r in res); pair_max += maximum(r.el for r in res); pair_gc += sum(r.gc for r in res)
        h *= 2
    end
    t_2 = time_ns()

    # the multiplier recursion of the dual residual (vectors only; f_x is the identity)
    λs = Vector{Vector{Float64}}(undef, N)
    λ = zeros(nx)
    for t in N:-1:1
        λs[t] = λ
        λ = r1[t].lx + λ
        λ = λ + r1[t].cxϕ
    end

    # 3. every period with the value function after it
    r3 = _pm_all((t, w) -> remotecall_fetch(pw_step3, w, t, μ, λs[t], t == N), M, N)
    if any(isnothing, r3)
        d.status = 1
        @printf("FILTERDDP_PROCS_FAIL iteration=%d step=3 periods=%s\n", d.k, string(findall(isnothing, r3)))
        return true
    end
    t_3 = time_ns()

    # the pass's scalars, in the order of the sequential pass
    d.status = 0
    d.barrier_lagrangian_curr = 0.0
    d.primal_1_curr = 0.0
    d.primal_inf = 0.0
    d.cs_inf_μ = 0.0
    d.cs_inf_0 = 0.0
    d.objective = 0.0
    d.dual_inf = 0.0
    d.expected_change_L = 0.0
    ϕ_norm = 0.0
    z_norm = 0.0
    for t in N:-1:1
        s, p = r1[t], r3[t]
        d.objective += s.lval
        d.primal_1_curr += s.c1
        d.primal_inf = max(d.primal_inf, s.cinf)
        d.barrier_lagrangian_curr -= μ * s.logsum
        d.cs_inf_0 = max(d.cs_inf_0, s.cs0)
        d.cs_inf_μ = max(d.cs_inf_μ, s.csμ)
        d.barrier_lagrangian_curr += s.cϕ
        d.dual_inf = max(d.dual_inf, p.dual)
        z_norm += p.zl_sum
        z_norm += p.zu_sum
        ϕ_norm += p.ϕ1
        d.expected_change_L += p.e1
        d.expected_change_L += p.e2
    end
    ni = (cl.nl + cl.nu) * N
    scaling_dual = max(options.s_max, (ϕ_norm + z_norm) / max(ni + nc * N, 1.0)) / options.s_max
    scaling_cs = max(options.s_max, z_norm / max(ni, 1.0)) / options.s_max
    d.dual_inf /= scaling_dual
    d.cs_inf_0 /= scaling_cs
    d.cs_inf_μ /= scaling_cs
    d.barrier_lagrangian_curr += d.objective
    d.reg_last = 0.0

    # Wall times of the three steps as the coordinator sees them, and beside
    # them what the periods' own clocks say: the slowest period (max) and all
    # of them (sum). The difference between a step's wall time and its slowest
    # period is the cost of running it as separate processes.
    @printf("FILTERDDP_PROCS_PASS iteration=%d barrier_iteration=%d stages=%d workers=%d rounds=%d combinations=%d step1_s=%.6f rounds_s=%.6f step3_s=%.6f total_s=%.6f pre_max_s=%.6f pre_sum_s=%.6f pair_max_s=%.6f pair_sum_s=%.6f post_max_s=%.6f post_sum_s=%.6f gc_sum_s=%.6f objective=%.12e primal_inf=%.6e dual_inf=%.6e cs_inf=%.6e mu=%.6e\n",
            d.k, d.j, N, length(M.pids), rounds, merges, (t_1 - t_start) / 1e9, (t_2 - t_1) / 1e9, (t_3 - t_2) / 1e9,
            (time_ns() - t_start) / 1e9, maximum(r.el for r in r1), sum(r.el for r in r1), pair_max, pair_sum,
            maximum(r.el for r in r3), sum(r.el for r in r3),
            sum(r.gc for r in r1) + pair_gc + sum(r.gc for r in r3),
            d.objective, d.primal_inf, d.dual_inf, d.cs_inf_0, μ)
    return true
end

function procs_affine_directions(solver, ocp, d, τ::T) where T
    M = _PM[]::ProcsMain
    N = ocp.N
    t0 = time_ns()
    xs = Vector{Vector{Float64}}(undef, N)
    xs[1] = M.x1                                          # the initial state: the same in every iterate
    for t in 1:N-1
        xs[t+1] = remotecall_fetch(pw_state_step, M.owner[t], t, xs[t])
    end
    t1 = time_ns()
    split = DDP4OPF._split_step()
    r = _pm_all((t, w) -> remotecall_fetch(pw_directions, w, t, xs[t], Float64(τ), split), M, N)
    γmax = T(minimum(x.γ for x in r))
    γdual = split ? T(minimum(x.γd for x in r)) : T(NaN)
    @printf("FILTERDDP_PROCS_DIRECTIONS iteration=%d chain_s=%.6f directions_s=%.6f dir_max_s=%.6f dir_sum_s=%.6f gc_sum_s=%.6f gamma_max=%.6e\n",
            d.k, (t1 - t0) / 1e9, (time_ns() - t1) / 1e9, maximum(x.el for x in r), sum(x.el for x in r),
            sum(x.gc for x in r), γmax)
    E = Vector{Vector{T}}()
    return DDP4OPF.AffineDirections{T}(E, E, E, E, E, γmax, γdual)
end

function procs_affine_rollout!(solver, ocp, d, dir, γ::T) where T
    M = _PM[]::ProcsMain
    N = ocp.N
    μ = Float64(d.μ)
    γz = isnan(dir.γdual) ? γ : dir.γdual
    t0 = time_ns()
    r = _pm_all((t, w) -> remotecall_fetch(pw_rollout, w, t, Float64(γ), Float64(γz), μ), M, N)
    d.status = 0
    d.primal_1_next = T(0.0)
    d.barrier_lagrangian_next = T(0.0)
    bad = false
    for t in 1:N                                   # sums in the order of the sequential rollout
        d.primal_1_next += r[t].c1
        d.barrier_lagrangian_next += r[t].cϕ
        if r[t].bad
            d.status = 2
            bad = true
            break
        end
        d.barrier_lagrangian_next -= μ * r[t].logl
        d.barrier_lagrangian_next -= μ * r[t].logu
        d.barrier_lagrangian_next += r[t].lv
    end
    @printf("FILTERDDP_PROCS_TRIAL iteration=%d gamma=%.6e wall_s=%.6f trial_max_s=%.6f trial_sum_s=%.6f gc_sum_s=%.6f rejected=%d\n",
            d.k, γ, (time_ns() - t0) / 1e9, maximum(x.el for x in r), sum(x.el for x in r), sum(x.gc for x in r), bad)
    return true
end

# When the solve stops: the trajectory back to the coordinator's solver.
function procs_solve_exit(solver)
    M = _PM[]::ProcsMain
    N = solver.ocp.N
    advance = solver.data.k != M.last_k            # the last accepted point was not followed by a pass
    els = _pm_all((t, w) -> remotecall_fetch(pw_nominal, w, t, advance), M, N)
    for t in 1:N; solver.nominal[t] = els[t]; end
    M.last_k = solver.data.k
    @printf("FILTERDDP_PROCS_MEMORY coordinator_maxrss_mib=%.1f worker_maxrss_mib_max=%.1f worker_maxrss_mib_sum=%.1f\n",
            Sys.maxrss() / 2^20, maximum(remotecall_fetch(pw_memory, p) for p in M.pids),
            sum(remotecall_fetch(pw_memory, p) for p in M.pids))
    return nothing
end

# Call after the workers exist and have this file and the driver loaded.
function procs_install(system::String, T::Int, nx::Int)
    pids = workers()
    (isempty(pids) || pids == [1]) && error("procs_pass.jl needs worker processes (addprocs)")
    owner = [pids[(t - 1) % length(pids) + 1] for t in 1:T]
    len = 3nx * nx + 2nx
    shared = Any[SharedArray{Float64}((len,); pids=vcat(myid(), pids)) for _ in 1:T]
    t0 = time()
    info = Vector{Any}(undef, length(pids))
    @sync for (i, p) in enumerate(pids)
        @async info[i] = remotecall_fetch(procs_worker_setup, p, system, T, findall(==(p), owner), shared)
    end
    _PM[] = ProcsMain(pids, owner, shared, -1, Float64[])
    DDP4OPF.BACKWARD_PASS_HOOK[] = procs_backward_pass!
    DDP4OPF.AFFINE_DIRECTIONS_HOOK[] = procs_affine_directions
    DDP4OPF.AFFINE_ROLLOUT_HOOK[] = procs_affine_rollout!
    DDP4OPF.SOLVE_EXIT_HOOK[] = procs_solve_exit
    @printf("PROCS_SWEEP installed: %d workers for %d periods (at most %d periods a worker), setup %.1f s (slowest worker %.1f s), worker memory %.0f MiB at most, shared %.1f MiB\n",
            length(pids), T, maximum(count(==(p), owner) for p in pids), time() - t0,
            maximum(x.setup_s for x in info), maximum(x.rss_mib for x in info), T * len * 8 / 2^20)
    return nothing
end
