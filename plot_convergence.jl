# Extract convergence data and plot using UnicodePlots fallback
ENV["GKSwstype"] = "nul"

import Pkg
env_path = joinpath(@__DIR__, "envs", "tadmm")
Pkg.activate(env_path)

# Load necessary packages for deserialization
using JuMP
using Ipopt
using Serialization
using Printf
using Plots

# Try GR with png workstation type
ENV["GKSwstype"] = "png"
gr()

data_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData", "large10kC_1ph_T4")
println("Loading solutions from: $data_dir")

sol_tadmm = deserialize(joinpath(data_dir, "sol_socp_tadmm.jls"))
sol_bf = deserialize(joinpath(data_dir, "sol_socp_bf.jls"))
println("Solutions loaded successfully")

# Extract convergence history
hist = sol_tadmm[:convergence_history]
obj_history = hist[:obj_history]
r_norm_history = hist[:r_norm_history]
s_norm_history = hist[:s_norm_history]
ρ_history = get(hist, :ρ_history, Float64[])
bf_obj = sol_bf[:objective]

n_iter = length(obj_history)

# Save convergence data to CSV
csv_path = joinpath(data_dir, "convergence", "convergence_data.csv")
mkpath(dirname(csv_path))
open(csv_path, "w") do f
    println(f, "iteration,objective,r_norm,s_norm,rho")
    for k in 1:n_iter
        rho_val = k <= length(ρ_history) ? ρ_history[k] : NaN
        @printf(f, "%d,%.6f,%.10e,%.10e,%.2f\n", k, obj_history[k], r_norm_history[k], s_norm_history[k], rho_val)
    end
end
println("Convergence data saved to: $csv_path")

# Print data for viewing
println("\n=== CONVERGENCE DATA ===")
println("BF Optimal: \$$(round(bf_obj, digits=2))")
println()
@printf("%-5s  %-15s  %-12s  %-12s  %-8s\n", "k", "Objective", "‖r‖", "‖s‖", "ρ")
@printf("%-5s  %-15s  %-12s  %-12s  %-8s\n", "---", "----------", "------", "------", "---")
for k in 1:n_iter
    rho_val = k <= length(ρ_history) ? @sprintf("%.1f", ρ_history[k]) : "N/A"
    @printf("%-5d  \$%-14.2f  %-12.4e  %-12.4e  %s\n", k, obj_history[k], r_norm_history[k], s_norm_history[k], rho_val)
end

# Try plotting
println("\nAttempting plot save...")
iterations = 1:n_iter
eps_pri = 8e-5
eps_dual = 3e-4

p1 = plot(iterations, obj_history,
    xlabel="Iteration (k)", ylabel="Objective [\$]",
    title="Objective Function", lw=3, color=:dodgerblue,
    markershape=:circle, markersize=5, label="tADMM",
    dpi=150, size=(800, 400))
hline!(p1, [bf_obj], color=:red, lw=2, linestyle=:dash,
    label="BF = \$$(round(bf_obj, digits=2))")

p2 = plot(iterations, r_norm_history,
    xlabel="Iteration (k)", ylabel="‖r‖",
    title="Primal Residual", lw=3, color=:darkgreen,
    markershape=:square, markersize=4, yscale=:log10,
    label="‖r‖", dpi=150, size=(800, 400))
hline!(p2, [eps_pri], color=:red, lw=2, linestyle=:dash,
    label="ε_pri=$(eps_pri)")

# Filter out zero values for log scale
s_nonzero = [(k, s) for (k, s) in zip(iterations, s_norm_history) if s > 0]
s_iters = [x[1] for x in s_nonzero]
s_vals = [x[2] for x in s_nonzero]

p3 = plot(s_iters, s_vals,
    xlabel="Iteration (k)", ylabel="‖s‖",
    title="Dual Residual", lw=3, color=:darkorange,
    markershape=:diamond, markersize=4, yscale=:log10,
    label="‖s‖", dpi=150, size=(800, 400))
hline!(p3, [eps_dual], color=:red, lw=2, linestyle=:dash,
    label="ε_dual=$(eps_dual)")

p_combined = plot(p1, p2, p3, layout=(3,1), size=(900, 1200), dpi=150,
    bottom_margin=5Plots.mm, left_margin=5Plots.mm)

outfile = joinpath(data_dir, "convergence", "tadmm_convergence_socp.png")
try
    savefig(p_combined, outfile)
    fsize = filesize(outfile)
    println("Plot saved: $outfile ($fsize bytes)")
catch e
    println("savefig failed: $e")
    println("Convergence data is available in CSV: $csv_path")
end
