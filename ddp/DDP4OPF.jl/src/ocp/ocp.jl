struct OCP{T, nx, nu, nc, F1, F2, F3, F4, F5, F6, C1, C2, C3, C4, C5, C6, OS1, OS2, OS3, OS4, OS5, OS6, OT1, OT2, OT3, OT4, OT5, OT6}
    N::Int64
    stage_objective::Objective{nx, nu, OS1, OS2, OS3, OS4, OS5, OS6}
    term_objective::Objective{nx, nu, OT1, OT2, OT3, OT4, OT5, OT6}
    dynamics::Dynamics{nx, nu, F1, F2, F3, F4, F5, F6}
    constraints::EqualityConstraints{nx, nu, nc, C1, C2, C3, C4, C5, C6}
    control_limits::ControlLimits{T, nu}
    # --- PATCH: optional per-stage data (see ddp/patches/per_stage_data.patch) ---
    # Remark 1 of the global-convergence paper states that the objective,
    # dynamics and constraint functions "can in general be time-varying"; the
    # package only ever exposed one of each for the whole horizon. These two
    # fields carry a distinct Objective / EqualityConstraints per stage when a
    # problem needs one. Both default to `nothing`, in which case every call
    # site falls back to the single shared function and the original code path
    # is bit-for-bit unchanged.
    stage_objectives::Union{Nothing, Vector}
    stage_constraints::Union{Nothing, Vector}
end

"""Objective in force at stage `t` (the shared one unless per-stage data was given)."""
@inline stage_obj(ocp::OCP, t::Int) =
    isnothing(ocp.stage_objectives) ? ocp.stage_objective : ocp.stage_objectives[t]

"""Equality constraints in force at stage `t`."""
@inline stage_con(ocp::OCP, t::Int) =
    isnothing(ocp.stage_constraints) ? ocp.constraints : ocp.stage_constraints[t]

function build_ocp(N::Int64, stage_objective::Objective{nx, nu, OS1, OS2, OS3, OS4, OS5, OS6},
        term_objective::Objective{nx, nu, OT1, OT2, OT3, OT4, OT5, OT6}, dynamics::Dynamics{nx, nu, F1, F2, F3, F4, F5, F6},
        constraints::EqualityConstraints{nx, nu, nc, C1, C2, C3, C4, C5, C6}, control_limits::ControlLimits{T, nu};
        stage_objectives::Union{Nothing, Vector} = nothing,
        stage_constraints::Union{Nothing, Vector} = nothing) where
                        {T<:Real, nx, nu, nc, F1, F2, F3, F4, F5, F6,
                        C1, C2, C3, C4, C5, C6, OS1, OS2, OS3, OS4, OS5, OS6,
                        OT1, OT2, OT3, OT4, OT5, OT6}
    l = ones(T, nx)
    u = floatmax(T) .* l
    l = -floatmax(T) .* l
    control_limits_ = isnothing(control_limits) ? ControlLimits(l, u) : control_limits
    constraints_ = isnothing(constraints) ? EqualityConstraints(nx, nu) : constraints
    return OCP{T, nx, nu, nc, F1, F2, F3, F4, F5, F6, C1, C2, C3, C4, C5, C6, OS1, OS2, OS3, OS4, OS5, OS6, OT1, OT2, OT3, OT4, OT5, OT6}(
        N, stage_objective, term_objective, dynamics, constraints_, control_limits_,
        stage_objectives, stage_constraints)
end
