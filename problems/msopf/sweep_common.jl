# sweep_common.jl -- helpers shared by the angle-sweep drivers in root_level/: labels, the
# window in which both substations import, and the plot style (the reference palette on the
# light surface: categorical slots 1 and 2, validated, for the two substations).

using Plots
using Printf

"\"subs3\" -> \"Subs 3\", \"grid1\" -> \"Subs 1\"."
label(name) = (m = match(r"(\d+)$", name)) === nothing ? name : "Subs " * m.captures[1]
"\"subs3\" -> \"δ3\"."
dlabel(name) = (m = match(r"(\d+)$", name)) === nothing ? "δ_" * name : "δ" * m.captures[1]

"Voltage at every load, per unit of that load's own kV."
load_V(net, x) = [abs(x.V_pu[d.bus]) * net.kV_base / d.kV for d in net.loads]

"Does every load sit inside the (Vminpu, Vmaxpu) band given for it (by load name)?"
function in_band(net, band, x)
    for (d, v) in zip(net.loads, load_V(net, x))
        lo, hi = band[d.name]
        lo <= v <= hi || return false
    end
    return true
end

"Network losses in kW: everything the substations deliver beyond the load (batteries add
their output: charging is extra load, discharging covers some)."
losses(x, P_load) = sum(values(x.P_subs_kW)) + sum(values(x.P_B_kW); init = 0.0) - P_load

"The same result with every value NaN: for a point with no power-flow solution."
blank(x) = merge(x, (; P_subs_kW = Dict(k => NaN for k in keys(x.P_subs_kW)),
                     Q_subs_kvar = Dict(k => NaN for k in keys(x.Q_subs_kvar)),
                     P_B_kW = Dict(k => NaN for k in keys(x.P_B_kW)),
                     V_pu = Dict(k => complex(NaN, NaN) for k in keys(x.V_pu))))

ipopt_status(x) = string(x.status) * (x.iterations === missing ? "" : " ($(x.iterations) it)")

"Indices (l, r) of the unbroken run of `true` in `ok` that contains index i."
function stretch(ok, i)
    l, r = i, i
    while l > 1 && ok[l-1]
        l -= 1
    end
    while r < length(ok) && ok[r+1]
        r += 1
    end
    return l, r
end

"""
    import_window(dd, Pa, Pb)

The interval around dd = 0 in which both substations import (Pa, Pb > 0; NaN where there
is no solution), each edge placed where the P that ends it crosses zero (linear
interpolation). `nothing` if they do not both import at 0; an edge is `nothing` if the
sweep, or the solutions, end first.
"""
function import_window(dd, Pa, Pb)
    both = [isfinite(a) && isfinite(b) && a > 0 && b > 0 for (a, b) in zip(Pa, Pb)]
    i0 = findfirst(==(0), dd)
    (i0 === nothing || !both[i0]) && return nothing
    l, r = stretch(both, i0)
    function edge(i, j)                          # i inside the window, j just outside
        1 <= j <= length(dd) && isfinite(Pa[j]) && isfinite(Pb[j]) || return nothing
        P = Pa[j] <= 0 ? Pa : Pb                 # whichever stopped importing
        return dd[i] + (dd[j] - dd[i]) * P[i] / (P[i] - P[j])
    end
    return (edge(l, l - 1), edge(r, r + 1))
end

"Last finite entry of a vector (end labels sit beside it), NaN if none."
lastfinite(v) = (i = findlast(isfinite, v)) === nothing ? NaN : v[i]

"Integer tick label with thousands separators: -12500 -> \"-12,500\"."
function with_commas(v)
    s = string(round(Int, abs(v)))
    s = reverse(join((join(c) for c in Iterators.partition(reverse(s), 3)), ","))
    return (round(Int, v) < 0 ? "-" : "") * s
end

const SURFACE = colorant"#fcfcfb"
const INK, INK2, MUTED = colorant"#0b0b0b", colorant"#52514e", colorant"#898781"
const GRIDC, AXISC, WASH = colorant"#e1e0d9", colorant"#c3c2b7", colorant"#f0efec"
const SUB_COLORS = (colorant"#2a78d6", colorant"#eb6834")   # substation a, substation b

"Apply the shared plot defaults."
plot_style!() = default(; fontfamily = "sans-serif", background_color = SURFACE,
                        foreground_color_axis = AXISC, foreground_color_border = AXISC,
                        foreground_color_text = INK2, foreground_color_guide = INK2,
                        foreground_color_title = INK, gridcolor = GRIDC, gridalpha = 1.0,
                        gridlinewidth = 1, gridstyle = :solid, titlefontsize = 11,
                        guidefontsize = 10, tickfontsize = 9, legendfontsize = 9,
                        legend_foreground_color = GRIDC)
