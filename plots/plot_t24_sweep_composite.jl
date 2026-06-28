#!/usr/bin/env julia
# Build a 2x3 composite of the 5 ieee2552C_1ph T=24 rho-sweep convergence PNGs.
# Output: pngs_review/ieee2552C_1ph/T024_rho_sweep_composite.png

using FileIO, Images, Plots
Plots.default(fontfamily = "Computer Modern")

const REVIEW_DIR = joinpath(@__DIR__, "pngs_review", "ieee2552C_1ph")
const RHOS = [4000, 8000, 16000, 24000, 32000]

panels = []
for r in RHOS
    fname = joinpath(REVIEW_DIR, "T024_rho$(lpad(r,6,'0'))_convergence.png")
    img = load(fname)
    p = plot(img;
        axis = false, grid = false, showaxis = false, ticks = nothing,
        title = "rho_init = $r", titlefontsize = 11,
    )
    push!(panels, p)
end
# Sixth slot: blank placeholder so the layout stays 2x3.
push!(panels, plot(legend = false, axis = false, grid = false,
                   showaxis = false, ticks = nothing, framestyle = :none))

fig = plot(panels...; layout = (2, 3), size = (1650, 1500), dpi = 150,
           plot_title = "ieee2552C  T=24  rho sweep (eps_pri=1e-4, eps_dual=5e-4)",
           plot_titlefontsize = 14,
           left_margin = 2Plots.mm, right_margin = 2Plots.mm,
           top_margin  = 2Plots.mm, bottom_margin = 2Plots.mm)

out = joinpath(REVIEW_DIR, "T024_rho_sweep_composite.png")
savefig(fig, out)
println("Saved: ", out)
