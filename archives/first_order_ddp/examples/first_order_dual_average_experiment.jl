# Compare the user's original first-order DDP dual update with a rolling mean
# of the current raw dual vector and the two preceding raw dual vectors.
#
# Run from the repository root:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/first_order_dual_average_experiment.jl

include("user_ddp.jl")

function detected_period(history; max_period=20, atol=1e-10)
    n = length(history)
    for period in 1:min(max_period, n ÷ 3)
        matched = true
        for offset in 0:(3 * period - 1)
            a = history[n - offset]
            b = history[n - offset - period]
            if !(isapprox(a.B, b.B; atol=atol, rtol=0) &&
                 isapprox(a.mu, b.mu; atol=atol, rtol=0))
                matched = false
                break
            end
        end
        matched && return period
    end
    return 0
end

T = 6
max_iter = 300
reference = centralized_reference(T)
reference === nothing && error("centralized T=6 reference is missing")
output = joinpath(@__DIR__, "..", "..", "results", "copper_plate",
                  "first_order_dual_average.csv")

open(output, "w") do io
    println(io, "variant,dual_window,sweeps,converged,cycle_period,J_min,J_max," *
                "objective_gap_min,objective_gap_max,max_PB_gap_min,max_PB_gap_max," *
                "state_step_max,dual_step_max")
    for (variant, window) in (("original", 1), ("three_dual_mean", 3))
        result = ddp_solve(T; max_iter=max_iter, tol=1e-8, verbose=false,
                           dual_average_window=window)
        period = detected_period(result.history)
        cycle_length = period == 0 ? min(20, length(result.history)) : period
        cycle = result.history[(end - cycle_length + 1):end]
        objective_gaps = abs.([h.J for h in cycle] .- reference.J)
        pb_gaps = [maximum(abs.(h.P_B .- reference.pb)) for h in cycle]
        converged = result.history[end].err < 1e-8
        row = (
            variant, window, length(result.history), converged, period,
            minimum(h.J for h in cycle), maximum(h.J for h in cycle),
            minimum(objective_gaps), maximum(objective_gaps),
            minimum(pb_gaps), maximum(pb_gaps),
            maximum(h.err_state for h in cycle),
            maximum(h.err_dual for h in cycle),
        )
        println(io, join(row, ','))
        println(join(row, ','))
    end
end
