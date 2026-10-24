# FILTERDDP_REUSE_BARRIER=1: no second factorization after a barrier update.
#
# When the barrier problem is solved to its tolerance, solve! lowers mu and
# calls backward_pass! again at the SAME trajectory. Nothing in a stage matrix
# depends on mu: the barrier terms of the Hessian are z ./ slack, and the
# value Hessian is built from the stage matrices alone. Only the right-hand
# side of the feedforward column changes (through mu * (1/uu - 1/ul)), and
# with it the step, the value gradient and the bound-multiplier updates.
#
# So the repeated pass needs no derivative evaluation, no assembly, no
# factorization and no feedback rows: every stage keeps its factor (it is
# already held by the stage's factor-backed policy, with B and c_x), and the
# pass is one solve per stage with that factor plus vector work. The same
# linear systems are solved; the step is the one the full pass would compute.
#
# Applies to the battery-block path with identity state dynamics and no
# regularisation (what the MPOPF runs use); otherwise the full pass runs.
# The repeated pass prints no FILTERDDP_FEASIBILITY / FILTERDDP_ORIGINAL_PRIMAL
# line: the iterate is the one of the pass before it.

_reuse_barrier() = get(ENV, "FILTERDDP_REUSE_BARRIER", "0") != "0"
reset_phase_times!(F) = nothing

struct ReuseStage{T}
    c::Vector{T}
    lu::Vector{T}
    cuϕ::Vector{T}
    lx::Vector{T}
    cxϕ::Vector{T}
    fu::SparseMatrixCSC{T, Int}
    cϕ::T                         # dot(c, ϕ)
    logsum::T                     # dot(log.(ul), maskl) + dot(log.(uu), masku)
end

struct ReusePass{T}
    solver_id::UInt
    k::Int
    diag::Bool                    # the Hessian model the factors were built with
    stages::Vector{ReuseStage{T}}
    primal_1::T
    primal_inf::T
    objective::T
    dual_inf::T                   # already scaled
    cs_inf_0::T                   # already scaled
    scaling_cs::T
end

const _REUSE_CACHE = Ref{Any}(nothing)
reuse_reset!() = (_REUSE_CACHE[] = nothing)

function _reuse_stage(c, lu_, cu, ϕ, lx, cx, fu, u, cl)
    ul = u - cl.l; uu = cl.u - u
    return ReuseStage{Float64}(copy(c), Vector(lu_), cu' * ϕ, Vector(lx), cx' * ϕ, sparse(fu), dot(c, ϕ),
                               dot(log.(ul), cl.maskl) + dot(log.(uu), cl.masku))
end

# The cached pass, if this call repeats it at the same trajectory.
function _reuse_candidate(solver, data)
    _reuse_barrier() || return nothing
    rc = _REUSE_CACHE[]
    (rc isa ReusePass && rc.solver_id == objectid(solver) && rc.k == data.k &&
        rc.diag == _diag_hessian_enabled()) || return nothing
    all(t -> !isnothing(solver.update[t].factor_policy), eachindex(rc.stages)) || return nothing
    return rc
end

# One-worker-per-period accounting of a stage of the repeated pass: its barrier
# terms and first network solve need nothing from the next period, its
# multiplier updates and second network solve nothing but its own battery
# solution; the rest (right-hand side, battery block, value gradient) is the
# sequential sweep.
function _backward_pass_reuse!(solver::Solver{T, nx, nu, nc, nux, ncx}, ocp::OCP{T, nx, nu, nc},
        traj::Vector{TrajectoryElement{T, nx, nu, nc}}, data::SolverData{T}, rc::ReusePass) where {T, nx, nu, nc, nux, ncx}
    μ = data.μ
    cl = ocp.control_limits
    parsim = _parsim()
    timing_diagnostic = get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"
    pass_start_ns = time_ns()
    data.status = 0
    data.barrier_lagrangian_curr = T(0.0)
    data.cs_inf_μ = T(0.0)
    data.expected_change_L = T(0.0)
    data.primal_1_curr = rc.primal_1
    data.primal_inf = rc.primal_inf
    data.objective = rc.objective
    data.dual_inf = rc.dual_inf
    data.cs_inf_0 = rc.cs_inf_0
    V̂x = zeros(T, nx)
    rhs = Vector{T}(undef, nu + nc)
    ps = zeros(5)                       # pre_sum, pre_max, seq_sum, post_sum, post_max
    for t = ocp.N:-1:1
        t0 = time_ns()
        st = rc.stages[t]
        u, zl, zu = traj[t].u, traj[t].zl, traj[t].zu
        rule = solver.update[t]
        fp = rule.factor_policy
        F = fp.factor

        ul = u - cl.l
        uu = cl.u - u
        data.barrier_lagrangian_curr -= μ * st.logsum
        cs_l = ul .* zl
        cs_u = uu .* zu
        data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_l .- μ) .* cl.maskl, Inf))
        data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_u .- μ) .* cl.masku, Inf))
        inv_ul = inv.(ul) .* cl.maskl
        inv_uu = inv.(uu) .* cl.masku
        t1 = time_ns()

        Qû = st.lu + st.fu' * V̂x + μ .* (inv_uu - inv_ul)
        data.barrier_lagrangian_curr += st.cϕ
        Qû = Qû + st.cuϕ
        @views @. rhs[1:nu] = -Qû
        @views @. rhs[nu+1:end] = -st.c
        reset_phase_times!(F)
        ldiv!(F, rhs)
        all(isfinite, rhs) || (data.status = 1)
        copyto!(rule.α, @view rhs[1:nu])
        copyto!(rule.ψ, @view rhs[nu+1:end])
        α = rule.α
        ψ = rule.ψ
        Vx_inc = fp.B_active' * α[fp.active_rows] + fp.cx' * ψ
        V̂x = st.lx + Vx_inc + V̂x
        V̂x = V̂x + st.cxϕ
        t2 = time_ns()

        @. rule.χl = inv_ul * μ - zl - rule.Σ_L * α
        @. rule.χu = inv_uu * μ - zu + rule.Σ_U * α
        data.expected_change_L += dot(Qû, α)
        data.expected_change_L += dot(st.c, ψ)
        t3 = time_ns()

        pt = stage_phase_times(F)
        pre_s = (t1 - t0) / 1e9 + (isnothing(pt) ? 0.0 : pt[3])
        post_s = (t3 - t2) / 1e9 + (isnothing(pt) ? 0.0 : pt[5])
        seq_s = (t3 - t0) / 1e9 - pre_s - post_s
        ps[1] += pre_s; ps[2] = max(ps[2], pre_s); ps[3] += seq_s; ps[4] += post_s; ps[5] = max(ps[5], post_s)
        parsim && @printf("FILTERDDP_PARSIM_STAGE iteration=%d stage=%d pre_s=%.6f seq_s=%.6f post_s=%.6f reuse=1\n",
                          data.k, t, pre_s, seq_s, post_s)
        data.status == 1 && break
    end
    data.cs_inf_μ /= rc.scaling_cs
    data.barrier_lagrangian_curr += data.objective
    data.reg_last = 0.0
    parsim && @printf(
        "FILTERDDP_PARSIM_BACKWARD iteration=%d stages=%d pre_sum_s=%.9f pre_max_s=%.9f seq_sum_s=%.9f post_sum_s=%.9f post_max_s=%.9f\n",
        data.k, ocp.N, ps[1], ps[2], ps[3], ps[4], ps[5])
    timing_diagnostic && @printf("FILTERDDP_REUSE_PASS iteration=%d barrier_iteration=%d stages=%d total_s=%.9f mu=%.9e\n",
                                 data.k, data.j, ocp.N, (time_ns() - pass_start_ns) / 1e9, μ)
    return nothing
end
