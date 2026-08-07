# Emits figures/input_curves.csv for input_curves_fig.tex.
#
# The benchmark's per-period data are not sampled from an analytic shape; they
# are given at the nodes and then represented inside the OCP by the Lagrange
# interpolants through those nodes (see the model section, eq. (4)), because
# FilterDDP admits one stage cost and one constraint function for the whole
# horizon. The input figure therefore plots the interpolants themselves, which
# is what the solver actually evaluates, with the nodes marked on top.
#
# Sampling here rather than writing the polynomial into the .tex avoids relying
# on pgfplots' limited-precision arithmetic for a degree-5 interpolant.
#
# Run:  julia --startup-file=no ddp/paper/figures/make_input_curves.jl

# The per-period data are the tADMM profiles at T = 6, identical to
# envs/tadmm/root_level/config.jl and to ddp/examples/power_system/tadmm_profiles.jl.
const PL_RATED_ = 2.0
const PL = PL_RATED_ .* (0.8 .+ 0.2 .* (sin.(range(0, 2pi, length = 6) .- 0.8) .+ 1) ./ 2)
const C  = 0.08 .+ 0.12 .* (sin.(range(0, 2pi, length = 6)) .+ 1) ./ 2

# Display conversions, so the figure carries the same quantities as the tADMM
# paper's input-profile plot: a dimensionless load factor and a price in
# cents/kWh.
const PL_RATED = 2.0     # rated load [p.u.]; λ = p_L / PL_RATED lands in (0,1]
const CENTS = 100.0      # with P_BASE = 1 kW and Δt = 1 h, 1 p.u.h = 1 kWh,
                         # so c [$/p.u.h] * 100 = cents/kWh

"Lagrange interpolant through (1,v[1]), ..., (n,v[n]) evaluated at τ."
function lagrange(v, τ)
    n = length(v)
    acc = 0.0
    for i = 1:n
        term = v[i]
        for j = 1:n
            j == i && continue
            term *= (τ - j) / (i - j)
        end
        acc += term
    end
    return acc
end

# sanity: the interpolant must reproduce the data exactly at the nodes
for t = 1:length(PL)
    @assert abs(lagrange(PL, float(t)) - PL[t]) < 1e-12
    @assert abs(lagrange(C, float(t)) - C[t]) < 1e-12
end
println("node reproduction OK (max error < 1e-12)")

out = joinpath(@__DIR__, "input_curves.csv")
open(out, "w") do io
    println(io, "tau,lambda,cents")
    for k = 0:400
        τ = 1.0 + 5.0 * k / 400          # τ ∈ [1, 6], the domain the model uses
        println(io, τ, ",", lagrange(PL, τ) / PL_RATED, ",", lagrange(C, τ) * CENTS)
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

# how far the interpolant strays between nodes -- worth knowing, since the
# degree grows with T and this is the device that does not scale
mx_l = maximum(abs(lagrange(PL, 1.0 + 5.0k/400)) for k = 0:400)
mx_c = maximum(abs(lagrange(C, 1.0 + 5.0k/400)) for k = 0:400)
println("max |p_L(τ)| on [1,6] = ", round(mx_l, digits = 4),
        "   (node max ", maximum(PL), ")")
println("max |c(τ)|   on [1,6] = ", round(mx_c, digits = 4),
        "   (node max ", maximum(C), ")")
