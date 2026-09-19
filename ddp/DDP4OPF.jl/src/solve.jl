function solve!(solver::Solver{T, nx, nu, nc, nux, ncx},
        x1::AbstractVector, u::AbstractVector{<:AbstractVector}; kwargs...) where {T, nx, nu, nc, nux, ncx}
    initialize_trajectory!(solver, u, x1)
    status = solve!(solver; kwargs...)
    return status
end

# ------------------------------------------------------ near-optimal stop --
# Standing rule (user, 2026-09-19): a run is scored by the time until its iterate
# is FEASIBLE (primal infeasibility <= FILTERDDP_NEAR_OPT_PRIMAL, default 1e-6)
# AND its objective is within FILTERDDP_NEAR_OPT_GAP (default 0.5%) of the
# centralized reference FILTERDDP_NEAR_OPT_REFERENCE, which is always solved
# first. The strict-tolerance tail (dual residual, mu -> 1e-8) does not count --
# as with tADMM, the dual residual can keep a run going long after the answer is
# in hand. Off unless the reference is set, so default behaviour is unchanged.
# Stops with status 9 and prints one FILTERDDP_NEAR_OPT line.
function _near_opt_settings()
    haskey(ENV, "FILTERDDP_NEAR_OPT_REFERENCE") || return nothing
    ref = parse(Float64, ENV["FILTERDDP_NEAR_OPT_REFERENCE"])
    gap = parse(Float64, get(ENV, "FILTERDDP_NEAR_OPT_GAP", "0.005"))
    primal = parse(Float64, get(ENV, "FILTERDDP_NEAR_OPT_PRIMAL", "1e-6"))
    isfinite(ref) || error("FILTERDDP_NEAR_OPT_REFERENCE must be finite")
    0.0 <= gap < 1.0 || error("FILTERDDP_NEAR_OPT_GAP must lie in [0, 1)")
    primal >= 0.0 || error("FILTERDDP_NEAR_OPT_PRIMAL must be nonnegative")
    return (ref = ref, gap = gap, primal = primal)
end

function solve!(solver::Solver{T, nx, nu, nc, nux, ncx,}) where {T, nx, nu, nc, nux, ncx,}
    (solver.options.verbose && solver.data.k==0) && solver_info()
    solve_start_ns = time_ns()
    near_opt = _near_opt_settings()

	ocp = solver.ocp
    options = solver.options
	data = solver.data
    
    reset!(data)
    reset_filter!(data)

    data.μ = options.μ_init
    ni = (ocp.control_limits.nl + ocp.control_limits.nu) * ocp.N
    timing_diagnostic = get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"

    while data.k < options.max_iterations
        backward_start_ns = time_ns()
        backward_pass!(solver, ocp, solver.nominal, data, options, verbose=options.verbose)
        backward_s = (time_ns() - backward_start_ns) / 1e9
        timing_diagnostic && flush(stdout)
        data.status != 0 && break

        # near-optimal stop: feasible and within the gap of the centralized reference
        if near_opt !== nothing && data.primal_inf <= near_opt.primal &&
                abs(data.objective - near_opt.ref) <=
                    near_opt.gap * max(abs(near_opt.ref), eps(Float64))
            timing_diagnostic && @printf(
                "FILTERDDP_ITER_TIMING iteration=%d barrier_iteration=%d backward_s=%.9f forward_s=0.000000000 total_s=%.9f outcome=near_optimal step_size=%.9e backtracks=%d mu=%.9e\n",
                data.k, data.j, backward_s, backward_s, data.step_size, data.l, data.μ)
            options.verbose && iteration_status(data, options)
            @printf("FILTERDDP_NEAR_OPT iteration=%d elapsed_s=%.3f objective=%.12f reference=%.12f rel_gap=%.3e primal_inf=%.3e dual_inf=%.3e mu=%.3e\n",
                data.k, (time_ns() - solve_start_ns) / 1e9, data.objective, near_opt.ref,
                abs(data.objective - near_opt.ref) / max(abs(near_opt.ref), eps(Float64)),
                data.primal_inf, data.dual_inf, data.μ)
            flush(stdout)
            data.status = 9
            break
        end

        # check (outer) overall problem convergence
        # check (inner) barrier problem convergence and update barrier parameter if so
        opt_err_μ = max(data.dual_inf, data.cs_inf_μ, data.primal_inf)
        opt_err_0 = max(data.dual_inf, data.cs_inf_0, data.primal_inf)

        if opt_err_0 < options.optimality_tolerance
            timing_diagnostic && @printf(
                "FILTERDDP_ITER_TIMING iteration=%d barrier_iteration=%d backward_s=%.9f forward_s=0.000000000 total_s=%.9f outcome=converged step_size=%.9e backtracks=%d mu=%.9e\n",
                data.k, data.j, backward_s, backward_s, data.step_size, data.l, data.μ)
            timing_diagnostic && flush(stdout)
            break
        end
        if opt_err_μ <= options.κ_ϵ * data.μ && ni > 0 && data.μ > options.optimality_tolerance / 10.0
            timing_diagnostic && @printf(
                "FILTERDDP_ITER_TIMING iteration=%d barrier_iteration=%d backward_s=%.9f forward_s=0.000000000 total_s=%.9f outcome=barrier_update step_size=%.9e backtracks=%d mu=%.9e\n",
                data.k, data.j, backward_s, backward_s, data.step_size, data.l, data.μ)
            timing_diagnostic && flush(stdout)
            data.μ = max(options.optimality_tolerance / 10.0, min(options.κ_μ * data.μ, data.μ ^ options.θ_μ))
            reset_filter!(data)
            data.j += 1
            continue
        end
        
        options.verbose && iteration_status(data, options)

        forward_start_ns = time_ns()
        forward_pass!(solver, ocp, data, options, verbose=options.verbose)
        forward_s = (time_ns() - forward_start_ns) / 1e9
        timing_diagnostic && @printf(
            "FILTERDDP_ITER_TIMING iteration=%d barrier_iteration=%d backward_s=%.9f forward_s=%.9f total_s=%.9f outcome=%s step_size=%.9e backtracks=%d mu=%.9e\n",
            data.k, data.j, backward_s, forward_s, backward_s + forward_s,
            data.status == 0 ? "accepted" : "forward_failed", data.step_size, data.l, data.μ)
        timing_diagnostic && flush(stdout)
        data.status != 0 && break
        
        update_nominal_trajectory!(solver)
        (!data.armijo_passed && !data.switching) && update_filter!(data, options)
        data.barrier_lagrangian_curr = data.barrier_lagrangian_next
        data.primal_1_curr = data.primal_1_next
        
        data.k += 1
    end
    
    data.k == options.max_iterations && (data.status = 8)
    options.verbose && iteration_status(data, options)
    options.verbose && on_exit(data)
    return data.status
end

function update_filter!(data::SolverData{T}, options::Options{T}) where T
    new_filter_pt = [(1. - options.γ_θ) * data.primal_1_curr,
                        data.barrier_lagrangian_curr - options.γ_L * data.primal_1_curr]
    push!(data.filter, new_filter_pt)
end

function reset_filter!(data::SolverData{T}) where T
    empty!(data.filter)
    push!(data.filter, [data.max_primal_1, T(-Inf)])
    data.status = 0
end
