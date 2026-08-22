# Emits figures/input_curves.csv for input_curves_fig.tex.
#
# The solver receives one load and one price value per stage.  For display only,
# this script evaluates the analytic profile between those stage samples and
# marks the actual values supplied to the solver on top.
#
# Run:  julia --startup-file=no ddp/paper/figures/make_input_curves.jl

# The per-period data are the analytic benchmark profiles at T = 6.
const PL_RATED_ = 2.0
const PL = PL_RATED_ .* (0.8 .+ 0.2 .* (sin.(range(0, 2pi, length = 6) .- 0.8) .+ 1) ./ 2)
const C  = 0.08 .+ 0.12 .* (sin.(range(0, 2pi, length = 6)) .+ 1) ./ 2

# Display conversions, so the figure carries the same quantities as the
# paper's input-profile plot: a dimensionless load factor and a price in
# cents/kWh.
const PL_RATED = 2.0     # rated load [p.u.]; λ = p_L / PL_RATED lands in (0,1]
const CENTS = 100.0      # with P_BASE = 1 kW and Δt = 1 h, 1 p.u.h = 1 kWh,
                         # so c [$/p.u.h] * 100 = cents/kWh

profile_load(τ) = 0.8 + 0.2 * (sin(2pi * (τ - 1) / 5 - 0.8) + 1) / 2
profile_cost(τ) = 0.08 + 0.12 * (sin(2pi * (τ - 1) / 5) + 1) / 2

out = joinpath(@__DIR__, "input_curves.csv")
open(out, "w") do io
    println(io, "tau,lambda,cents")
    for k = 0:400
        τ = 1.0 + 5.0 * k / 400
        println(io, τ, ",", profile_load(τ), ",", profile_cost(τ) * CENTS)
    end
end
println("wrote ", out)

# node file, plotted as markers on top of the curves
outn = joinpath(@__DIR__, "input_nodes.csv")
open(outn, "w") do io
    println(io, "tau,lambda,cents")
    for t = 1:length(PL)
        println(io, t, ",", PL[t] / PL_RATED, ",", C[t] * CENTS)
    end
end
println("wrote ", outn)
