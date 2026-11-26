# Reproducible T=3 diagnostic for the user's first-order DDP sweep.

include("user_ddp.jl")
include("mp_optimality.jl")

cost, demand = tadmm_tutorial3()
result = ddp_solve(cost, demand; max_iter=2000, tol=1e-10, verbose=false)
reference = solve_explicit(3; bounded=true, cost=cost, pL=demand)

println("T=3 adjusted price = ", cost)
println("T=3 adjusted demand = ", demand)
println("sweeps = ", length(result.history))
println("last 10 (iteration, residual, objective):")
for row in result.history[max(1, end-9):end]
    println((row.k, row.err, row.J))
end
println("centralized objective = ", reference.objective)
println("last objective gap = ", abs(result.history[end].J - reference.objective))
println("maximum P_B error = ", maximum(abs.(result.P_B .- reference.P_B)))

# Paper data: generated here so figures and tables cannot drift from the run.
figure_dir = normpath(joinpath(@__DIR__, "..", "..", "paper", "figures"))
mkpath(figure_dir)

open(joinpath(figure_dir, "user_ddp_t3_inputs.csv"), "w") do io
    println(io, "t,cost,load_factor")
    for t in eachindex(cost)
        println(io, "$(t),$(cost[t]),$(demand[t] / P_RATED)")
    end
end

open(joinpath(figure_dir, "user_ddp_t3_objective.csv"), "w") do io
    println(io, "iteration,objective,centralized")
    for row in result.history[1:min(20, end)]
        println(io, "$(row.k),$(row.J),$(reference.objective)")
    end
end

open(joinpath(figure_dir, "user_ddp_t3_mu_table.tex"), "w") do io
    println(io, raw"\begin{tabular}{cS[table-format=1.5]S[table-format=1.5]S[table-format=1.5]}")
    println(io, raw"\toprule")
    println(io, raw"{$k$} & {$\lambda_{\mathrm{soc}}^{k,1}$} & {$\lambda_{\mathrm{soc}}^{k,2}$} & {$\lambda_{\mathrm{soc}}^{k,3}$} " * "\\\\")
    println(io, raw"\midrule")
    for row in result.history[15:20]
        @printf(io, "%d & %.5f & %.5f & %.5f \\\\\n", row.k, row.mu...)
    end
    println(io, raw"\bottomrule")
    println(io, raw"\end{tabular}")
end
