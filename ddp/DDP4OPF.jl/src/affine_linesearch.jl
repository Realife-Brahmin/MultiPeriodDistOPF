# Line search for AFFINE dynamics (FILTERDDP_AFFINE_LINESEARCH=1).
#
# The trial point of the forward pass is
#     u_t = ū_t + γ α_t + β_t (x_t - x̄_t),      x_{t+1} = f(x_t, u_t),
# and likewise for the multipliers. When f is affine (the battery dynamics
# B^{t+1} = B^t - Δt P_B are), the whole trial trajectory is affine in the
# step size γ:
#     (x, u, ϕ, zl, zu)_t(γ) = nominal_t + γ (dx, du, dϕ, dzl, dzu)_t,
# with the directions given by ONE rollout at γ = 1. Two consequences:
#
#   1. the largest step that keeps every control and bound multiplier inside
#      the fraction-to-boundary rule is a ratio test, as in Ipopt, instead of
#      "halve the step until the trial passes";
#   2. every trial of the line search is an axpy and a constraint evaluation:
#      the stage systems are solved once per iteration, not once per trial.
#
# The acceptance tests (filter, switching, Armijo) are not touched. This is a
# change of algorithm (the accepted step is no longer a power of two), so the
# iterates differ from the halving line search.

_affine_linesearch() = get(ENV, "FILTERDDP_AFFINE_LINESEARCH", "0") != "0"
# FILTERDDP_SPLIT_STEP=1 (with the affine line search): separate step sizes for
# the controls and for the bound multipliers, as in Ipopt. The controls,
# states and equality multipliers take the line-search step, limited by the
# control bounds only; the bound multipliers take their own largest step.
# Without it one step size serves both and the smaller limit wins.
_split_step() = get(ENV, "FILTERDDP_SPLIT_STEP", "0") != "0"

struct AffineDirections{T}
    dx::Vector{Vector{T}}
    du::Vector{Vector{T}}
    dϕ::Vector{Vector{T}}
    dzl::Vector{Vector{T}}
    dzu::Vector{Vector{T}}
    γmax::T
    γdual::T                    # NaN: the multipliers take the line-search step
end

function affine_directions(solver, ocp, data, τ::T) where T
    N = ocp.N
    cl = ocp.control_limits
    dx = Vector{Vector{T}}(undef, N); du = Vector{Vector{T}}(undef, N)
    dϕ = Vector{Vector{T}}(undef, N); dzl = Vector{Vector{T}}(undef, N); dzu = Vector{Vector{T}}(undef, N)
    γmax = one(T); γdual = one(T)
    split = _split_step()
    lim_t = 0; lim_i = 0; lim_kind = "none"
    x = solver.nominal[1].x
    for t in 1:N
        nom = solver.nominal[t]; rule = solver.update[t]
        δx = x - nom.x
        βδx, ωδx = policy_actions!(rule, δx)        # views into the policy's buffer: copy below
        dx[t] = δx
        du[t] = rule.α .+ βδx
        dϕ[t] = rule.ψ .+ ωδx
        dzl[t] = rule.χl .- rule.Σ_L .* βδx
        dzu[t] = rule.χu .+ rule.Σ_U .* βδx
        d, el, eu = du[t], dzl[t], dzu[t]
        @inbounds for i in eachindex(d)
            if cl.maskl[i]
                if d[i] < 0
                    g = τ * (nom.u[i] - cl.l[i]) / -d[i]
                    g < γmax && (γmax = g; lim_t = t; lim_i = i; lim_kind = "control_lower")
                end
                if el[i] < 0
                    g = τ * nom.zl[i] / -el[i]
                    if split
                        γdual = min(γdual, g)
                    elseif g < γmax
                        γmax = g; lim_t = t; lim_i = i; lim_kind = "multiplier_lower"
                    end
                end
            end
            if cl.masku[i]
                if d[i] > 0
                    g = τ * (cl.u[i] - nom.u[i]) / d[i]
                    g < γmax && (γmax = g; lim_t = t; lim_i = i; lim_kind = "control_upper")
                end
                if eu[i] < 0
                    g = τ * nom.zu[i] / -eu[i]
                    if split
                        γdual = min(γdual, g)
                    elseif g < γmax
                        γmax = g; lim_t = t; lim_i = i; lim_kind = "multiplier_upper"
                    end
                end
            end
        end
        if t < N
            u1 = nom.u .+ d
            xnext = solver.ocp.dynamics.f(x, u1)
            if t == 1                                # the premise: f affine along the step
                mid = solver.ocp.dynamics.f(nom.x .+ T(0.5) .* δx, nom.u .+ T(0.5) .* d)
                ref = T(0.5) .* (solver.ocp.dynamics.f(nom.x, nom.u) .+ xnext)
                norm(mid - ref, Inf) <= 1e-10 * (1 + norm(ref, Inf)) ||
                    error("FILTERDDP_AFFINE_LINESEARCH needs affine dynamics")
            end
            x = xnext
        end
    end
    _ftb_diagnostic() && @printf("FILTERDDP_STEP_LIMIT iteration=%d gamma_max=%.6e stage=%d kind=%s index=%d gamma_dual=%.6e\n",
                                  data.k, γmax, lim_t, lim_kind, lim_i, split ? γdual : γmax)
    return AffineDirections{T}(dx, du, dϕ, dzl, dzu, γmax, split ? γdual : T(NaN))
end

function affine_rollout!(solver::Solver{T, nx, nu, nc, nux, ncx}, ocp::OCP{T, nx, nu, nc},
            data::SolverData{T}, dir::AffineDirections{T}, γ::T) where {T, nx, nu, nc, nux, ncx}
    μ = data.μ
    cl = ocp.control_limits
    data.status = 0
    data.primal_1_next = T(0.0)
    data.barrier_lagrangian_next = T(0.0)
    γz = isnan(dir.γdual) ? γ : dir.γdual
    for t = 1:ocp.N
        nom = solver.nominal[t]
        x = nom.x .+ γ .* dir.dx[t]
        u = nom.u .+ γ .* dir.du[t]
        ϕ = nom.ϕ .+ γ .* dir.dϕ[t]
        zl = nom.zl .+ γz .* dir.dzl[t]
        zu = nom.zu .+ γz .* dir.dzu[t]
        solver.current[t] = TrajectoryElement{T, nx, nu, nc}(x, u, ϕ, zl, zu)

        if nc > 0
            c = stage_con(ocp, t).c(x, u)
            data.primal_1_next += norm(c, 1)
            data.barrier_lagrangian_next += dot(c, ϕ)
        end

        # The ratio test already enforces the fraction-to-boundary rule; only
        # rounding can leave a slack or multiplier non-positive.
        ul = u - cl.l
        uu = cl.u - u
        if any((ul .<= 0) .* cl.maskl) || any((uu .<= 0) .* cl.masku) ||
           any((zl .<= 0) .* cl.maskl) || any((zu .<= 0) .* cl.masku)
            data.status = 2
            return
        end

        data.barrier_lagrangian_next -= μ * dot(log.(ul), cl.maskl)
        data.barrier_lagrangian_next -= μ * dot(log.(uu), cl.masku)
        if t == ocp.N
            data.barrier_lagrangian_next += ocp.term_objective.l(x, u)[1]
        else
            data.barrier_lagrangian_next += stage_obj(ocp, t).l(x, u)[1]
        end
    end
end
