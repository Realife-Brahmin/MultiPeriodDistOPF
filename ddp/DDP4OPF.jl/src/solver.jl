struct TrajectoryElement{T, nx, nu, nc}
    x::Vector{T}
    u::Vector{T}
    ϕ::Vector{T}
    zl::Vector{T}
    zu::Vector{T}
    function TrajectoryElement{T, nx, nu, nc}(x::AbstractVector, u::AbstractVector,
            ϕ::AbstractVector, zl::AbstractVector, zu::AbstractVector) where {T, nx, nu, nc}
        new{T, nx, nu, nc}(Vector{T}(x), Vector{T}(u), Vector{T}(ϕ),
            Vector{T}(zl), Vector{T}(zu))
    end
end

function TrajectoryElement(T::T_, nx::Int, nu::Int, nc::Int) where T_
    TrajectoryElement{T, nx, nu, nc}(
        zeros(T, nx), zeros(T, nu), zeros(T, nc), zeros(T, nu), zeros(T, nu)
        )
end

mutable struct FactorBackedPolicy{T}
    factor::Any
    active_rows::Vector{Int}
    B_active::Matrix{T}
    cx::SparseMatrixCSC{T, Int}
    rhs::Vector{T}
    Bδx::Vector{T}
end

mutable struct UpdateRule{T, nx, nu, nc, nux, ncx}
    α::Vector{T}
    ψ::Vector{T}
    β::Matrix{T}
    ω::Matrix{T}
    χl::Vector{T}
    χu::Vector{T}
    Σ_L::Vector{T}
    Σ_U::Vector{T}
    factor_policy::Union{Nothing, FactorBackedPolicy{T}}
    function UpdateRule{T, nx, nu, nc, nux, ncx}(α::Vector{T},
            ψ::Vector{T}, β::Matrix{T}, ω::Matrix{T},
            χl::Vector{T}, χu::Vector{T}, Σ_L::Vector{T},
            Σ_U::Vector{T}) where {T, nx, nu, nc, nux, ncx}
        new{T, nx, nu, nc, nux, ncx}(α, ψ, β, ω, χl, χu, Σ_L, Σ_U, nothing)
    end
    function UpdateRule{T, nx, nu, nc, nux, ncx}(α::AbstractVector,
            ψ::AbstractVector, β::AbstractMatrix, ω::AbstractMatrix,
            χl::AbstractVector, χu::AbstractVector, Σ_L::AbstractVector,
            Σ_U::AbstractVector) where {T, nx, nu, nc, nux, ncx}
        UpdateRule{T, nx, nu, nc, nux, ncx}(Vector{T}(α), Vector{T}(ψ),
            Matrix{T}(β), Matrix{T}(ω), Vector{T}(χl), Vector{T}(χu),
            Vector{T}(Σ_L), Vector{T}(Σ_U))
    end
end

function policy_actions!(rule::UpdateRule{T}, δx::AbstractVector{T}) where T
    policy = rule.factor_policy
    if isnothing(policy)
        return rule.β * δx, rule.ω * δx
    end

    nu = length(rule.α)
    fill!(policy.rhs, zero(T))
    mul!(policy.Bδx, policy.B_active, δx)
    @views @. policy.rhs[policy.active_rows] = -policy.Bδx
    mul!(@view(policy.rhs[nu+1:end]), policy.cx, δx)
    @views policy.rhs[nu+1:end] .*= -one(T)
    ldiv!(policy.factor, policy.rhs)
    return @view(policy.rhs[1:nu]), @view(policy.rhs[nu+1:end])
end

function UpdateRule(T::T_, nx::Int, nu::Int, nc::Int) where T_
    nux = nu * nx
    ncx = nc * nx
    factor_backed = get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1"
    β = factor_backed ? zeros(T, 0, 0) : zeros(T, nu, nx)
    ω = factor_backed ? zeros(T, 0, 0) : zeros(T, nc, nx)
    UpdateRule{T, nx, nu, nc, nux, ncx}(
        zeros(T, nu), zeros(T, nc), β, ω,
        zeros(T, nu), zeros(T, nu), zeros(T, nu), zeros(T, nu)
    )
end

struct Solver{T, nx, nu, nc, nux, ncx}
    # Deliberately erase the many generated-function type parameters here.
    # Keeping them in Solver's type makes inference recurse through very large
    # symbolic expressions on network-scale models.
    ocp::OCP
    nominal::Vector{TrajectoryElement{T, nx, nu, nc}}
    current::Vector{TrajectoryElement{T, nx, nu, nc}}
    update::Vector{UpdateRule{T, nx, nu, nc, nux, ncx}}
	data::SolverData{T}
    options::Options{T}
end

function Solver(ocp::OCP{T, nx, nu, nc, F1, F2, F3, F4, F5, F6, C1, C2, C3, C4, C5, C6, OS1, OS2, OS3, OS4, OS5, OS6, OT1, OT2, OT3, OT4, OT5, OT6};
        options::Union{Options{T}, Nothing}=nothing) where {T, nx, nu, nc, F1, F2, F3, F4, F5, F6, C1, C2, C3, C4, C5, C6, OS1, OS2, OS3, OS4, OS5, OS6, OT1, OT2, OT3, OT4, OT5, OT6}
    nominal = [TrajectoryElement(T, nx, nu, nc) for _ = 1:ocp.N]
    current = [TrajectoryElement(T, nx, nu, nc) for _ = 1:ocp.N]
    update = [UpdateRule(T, nx, nu, nc) for _ = 1:ocp.N]
    data = solver_data(T)
    options = isnothing(options) ? Options{T}() : options
	return Solver{T, nx, nu, nc, nx * nu, nx * nc}(
        ocp, nominal, current, update, data, options)
end

function get_trajectory(solver::Solver{T, nx, nu, nc, nux, ncx}
        ) where {T, nx, nu, nc, nux, ncx}
    x = [nom.x for nom in solver.nominal]
    u = [nom.u for nom in solver.nominal]
	return x, u
end

function initialize_trajectory!(solver::Solver{T, nx, nu, nc, nux, ncx},
        u::AbstractVector{<:AbstractVector}, x1::AbstractVector) where {T, nx, nu, nc, nux, ncx}
    options = solver.options
    κ_1 = options.κ_1
    κ_2 = options.κ_2

    cl = solver.ocp.control_limits
    x = Vector{T}(x1)
    for t = 1:solver.ocp.N
        ūt = zeros(T, nu)
        mask_lo = cl.maskl .* .!cl.masku
        ūt = ūt + max.(u[t],  options.κ_1 .* max.(cl.l, 1.0) + cl.l) .* mask_lo

        mask_up = cl.masku .* .!cl.maskl
        ūt = ūt + min.(u[t], -options.κ_1 .* max.(cl.u, 1.0) + cl.u) .* mask_up

        mask_both = cl.masku .* cl.maskl
        u1 = cl.l + min.(κ_1 * max.(1.0, abs.(cl.l)), κ_2 .* (cl.u - cl.l))
        u2 = cl.u - min.(κ_1 * max.(1.0, abs.(cl.u)), κ_2 .* (cl.u - cl.l))
        ūt = ūt + min.(max.(u[t], u1), u2) .* mask_both

        mask_none = .!cl.masku .* .!cl.maskl
        ūt = ūt + u[t] .* mask_none

        ϕ = zeros(T, nc)
        tmp = ones(T, nu)
        zl = Vector{T}(tmp .* cl.maskl)
        zu = Vector{T}(tmp .* cl.masku)

        solver.nominal[t] = TrajectoryElement{T, nx, nu, nc}(x, ūt, ϕ, zl, zu)
        
        if t < solver.ocp.N
            x = solver.ocp.dynamics.f(x, ūt)
        end
    end
end

function update_nominal_trajectory!(solver::Solver{T, nx, nu, nc, nux, ncx}) where
        {T, nx, nu, nc, nux, ncx}
    for t = 1:solver.ocp.N
        solver.nominal[t] = solver.current[t]
    end
    return nothing
end

function get_feedback(solver::Solver{T, nx, nu, nc, nux, ncx},
        t::Int) where {T, nx, nu, nc, nux, ncx}
    α = solver.update[t].α
    rule = solver.update[t]
    ū = solver.nominal[t].u
    x̄ = solver.nominal[t].x
    f = x -> begin
        βδx, _ = policy_actions!(rule, x - x̄)
        ū + α + βδx
    end
    return f
end
