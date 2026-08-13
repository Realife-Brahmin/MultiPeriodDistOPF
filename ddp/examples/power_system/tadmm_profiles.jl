# The tADMM load and price profiles, shared by every copper-plate script.
#
# These mirror envs/tadmm/root_level/config.jl EXACTLY, so the copper-plate
# benchmark and the tADMM experiments are driven by the same time series:
#
#     LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2pi, length=T) .- 0.8) .+ 1) ./ 2
#     LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2pi, length=T)) .+ 1) ./ 2
#
# Note these RESAMPLE with T -- they are not a fixed 24-point series that gets
# sliced. tadmm_cost(3) is not tadmm_cost(6)[1:3]. Any code that used to write
# `PL6[1:T]` must call the functions instead.
#
# Two properties of the tADMM formulas that matter here:
#
#   1. `range(0, 2pi, length=T)` includes BOTH endpoints, and sin(0) = sin(2pi),
#      so c^1 == c^T for every T. The price series is periodic with a duplicated
#      endpoint rather than a clean cycle of T distinct points.
#
#   2. CHECK THE PRICE BEFORE ADOPTING A NEW T. At T = 3 the sample points are
#      0, pi and 2pi, where sin vanishes at all three, so the price comes out
#      CONSTANT: c = [0.14, 0.14, 0.14]. Such an instance has no arbitrage signal
#      at all -- the optimal P_B is identical in every period and set purely by
#      the terminal penalty -- which makes it a poor storage benchmark however
#      well the solver behaves on it. T = 3 was used and then dropped for exactly
#      this reason; the committed results are T = 6 only. Any small T is worth
#      inspecting the same way before it is used.
#
# The load profile is phase-shifted by -0.8 rad relative to the price. That
# offset is what keeps the two series from being collinear: Pearson r is 0.644 at
# T = 6 and 0.682 at T = 24, against 0.997 for the hand-picked series this
# replaced. Real day-ahead price/load correlation sits in much the same range.

"""Load factor lambda^t, in [0.8, 1.0]. Identical to tADMM's `LoadShapeLoad`."""
tadmm_load(T::Int) = T == 1 ? [0.97174] :
    0.8 .+ 0.2 .* (sin.(range(0, 2pi, length = T) .- 0.8) .+ 1) ./ 2

"""Energy price c^t in \$/kWh, in [0.08, 0.20]. Identical to tADMM's `LoadShapeCost`."""
tadmm_cost(T::Int) = T == 1 ? [0.14] :
    0.08 .+ 0.12 .* (sin.(range(0, 2pi, length = T)) .+ 1) ./ 2

"""Rated load, p.u. The copper-plate demand is the load factor against this."""
const P_RATED = 2.0

"""Demand p_L^t in p.u. = rated load x load factor."""
tadmm_pL(T::Int) = P_RATED .* tadmm_load(T)
