# Where the T = 24 relaxation of the dispatched day cheats (2026-10-31). DIAGNOSTIC, NOT A CERTIFICATE.
# Solves the CD relaxation (envelope current limits, cuts as farkas_bracketing.jl) at one unresolved point beyond the
# dispatched reach and measures, per line and hour, how far its moments are from a real power flow:
#   fictitious loss  R_l (E|I|^2 - |E I|^2) in kW: losses booked beyond what the relaxation's own currents imply;
#   product angle    arg W_ij, W_ij = V_i conj(V_j) from the product variables, against the first-moment angle;
#   loop mismatch    sum of product angles around each fundamental cycle (tie loops), which is 0 for a real flow;
#   source paths     product angles summed along the tree path between source terminals, against the first moments.
# Lines on a cycle (closed ties) are told apart from radial lines, to see whether the gap sits on the loops.
# FB_CUTS picks the cuts (farkas_bracketing.jl; kcl0,jbound by default, lsec the line power secants).
# FB_LDL picks Clarabel's linear solver (qdldl by default; mkl = MKL Pardiso, needs Pardiso in the environment, threads
# from MKL_NUM_THREADS; mosek = the MOSEK solver instead of Clarabel, needs MosekTools in the environment and a licence
# in ~/mosek, threads from FB_MOSEK_THREADS); FB_GAPDIR names the output folder (relaxation_gap by default).
#   FB_SYSTEM=sfp1u_four FB_DIR=24 FB_RATIO=1.1 FB_IDLE=... MULTIPERIOD_OUT_DIR=<campaign folder> julia ... relaxation_gap_diagnose.jl
ENV["FB_PHASE"] = "none"                             # farkas_dispatched.jl: load its setup, run none of its phases
include(joinpath(@__DIR__, "farkas_dispatched.jl"))
const LDL = Symbol(get(ENV, "FB_LDL", "qdldl"))
LDL === :mkl && @eval using Pardiso                  # loads Clarabel's Pardiso extension
LDL === :mosek && @eval using MosekTools
const SOLVER = LDL === :mosek ?
    (; optimizer = Mosek.Optimizer, attrs = [("MSK_IPAR_NUM_THREADS", parse(Int, get(ENV, "FB_MOSEK_THREADS", "1")))]) :
    (; optimizer = Clarabel.Optimizer, attrs = [CLARABEL_DEFAULTS; ("direct_solve_method", LDL)])

const GD = joinpath(OUT, get(ENV, "FB_GAPDIR", "relaxation_gap")); mkpath(GD)
i = parse(Int, ENV["FB_DIR"]); d = DIRS[i]; ratio = parse(Float64, get(ENV, "FB_RATIO", "1.1"))
rs = parse(Float64, split(readlines(joinpath(OUT, "fb_AD_$(d.id).csv"))[2], ',')[4]); r = ratio * rs
amp = Dict{String,Float64}(); senv = Dict{String,Float64}()
for l in readlines(joinpath(OUT, "fb_envelope_disp.csv"))[2:end]
    c = split(l, ','); c[2] == "line" ? (amp[String(c[1])] = parse(Float64, c[4])) : (senv[String(c[1])] = parse(Float64, c[4]))
end
kwL = merge(kwb, (; ampacity_pu = amp, source_ampacity_pu = senv))
say(@sprintf("%s: %.4f x reach %.5f = %.5f deg; cuts %s; linear solver %s, %s threads", d.label, ratio, rs, r, join(CUTS, ","), LDL,
             LDL === :mosek ? get(ENV, "FB_MOSEK_THREADS", "1") : get(ENV, "MKL_NUM_THREADS", "1")))
# FB_LIFTCHECK=1: a validated dispatched day at 98% of the reach (cold solve), lifted into the relaxation with these
# cuts: every constraint must hold for it (a cut that a real point violates is not valid).
if get(ENV, "FB_LIFTCHECK", "0") == "1"
    tl = time(); m98, x98 = dmargin(d, 0.98rs, nothing)
    x98 === nothing && error("the 98% point does not validate")
    R98 = relax_margin(net, angles(d, 0.98rs), kwL; k = KV, cuts = CUTS)
    lv = lifted_violation(R98, net, x98, kwL; margin_kW = m98)
    say(@sprintf("  lifted check at 98%% (margin %.1f kW): worst violation %.2e at %s, box %.2e, lsec cuts %.2e over %d; %.0f s",
                 m98, lv.worst, lv.at, lv.box_worst, lv.lsec_worst, length(R98.lseccons), time() - tl))
    R98 = nothing; x98 = nothing; GC.gc()
end
t0 = time(); c = certify(net, angles(d, r), kwL; k = KV, cuts = CUTS, solver_attributes = SOLVER.attrs, optimizer = SOLVER.optimizer)
its = try barrier_iterations(c.R.model) catch; -1 end; tcl = try solve_time(c.R.model) catch; NaN end
say(@sprintf("  relaxation: %s, U %+.1f kW, fictitious loss %.1f kW (mean over hours), rank ratio %.2e, %d lsec cuts, %.0f s (solver %.0f s, %d iterations), max RSS %.1f GB",
             c.status, c.U, c.fictitious_loss_kW, c.rank_ratio, length(c.R.lseccons), time() - t0, tcl, its, Sys.maxrss() / 2^30))
if isfinite(c.farkas.bound); say(@sprintf("  Farkas bound %.4e (> 0 proves infeasibility)", c.farkas.bound)); end
has_values(c.R.model) || (say("  no primal point to diagnose (the relaxation is infeasible or failed)"); exit(0))
R = c.R; S = net.S_base_kVA; T = kwL.T; v(x) = value(x)

# spanning tree (BFS from the first source terminal); every non-tree line closes one fundamental cycle
adj = Dict(b => Tuple{String,String,Int}[] for b in net.buses)          # (neighbour, line, +1 if b is the from-end)
for l in net.lines; push!(adj[l.from], (l.to, l.name, 1)); push!(adj[l.to], (l.from, l.name, -1)); end
root = net.sources[1].bus; par = Dict(root => ("", "", 0)); depth = Dict(root => 0); q = [root]; intree = Set{String}()
while !isempty(q)
    b = popfirst!(q)
    for (nb, ln, sg) in adj[b]
        haskey(par, nb) && continue
        par[nb] = (b, ln, -sg); depth[nb] = depth[b] + 1; push!(intree, ln); push!(q, nb)   # sg seen from nb
    end
end
"Oriented lines (name, +1 when traversed from -> to) on the tree path from a to b."
function treepath(a, b)
    up, down = Tuple{String,Int}[], Tuple{String,Int}[]
    while a != b
        if depth[a] >= depth[b]; (p, ln, sg) = par[a]; push!(up, (ln, sg)); a = p          # a -> parent
        else; (p, ln, sg) = par[b]; push!(down, (ln, -sg)); b = p end                       # parent -> b, reversed later
    end
    return vcat(up, reverse(down))
end
lines = Dict(l.name => l for l in net.lines)
cycles = [(l.name, vcat([(l.name, 1)], treepath(l.to, l.from))) for l in net.lines if !(l.name in intree)]
onloop = Set{String}(); for (_, cy) in cycles, (ln, _) in cy; push!(onloop, ln); end
say(@sprintf("  %d lines, %d fundamental cycles, %d lines on a cycle", length(net.lines), length(cycles), length(onloop)))

g(s, n, t) = v(R.Lv[s][n, t])
thW(n, t) = atan(g(:Lfe, n, t) - g(:Lef, n, t), g(:Lee, n, t) + g(:Lff, n, t))           # arg V_i conj(V_j), products
thV(n, t) = (l = lines[n]; angle(complex(v(R.e[l.from, t]), v(R.f[l.from, t]))) - angle(complex(v(R.e[l.to, t]), v(R.f[l.to, t]))))
wrap(x) = mod(x + pi, 2pi) - pi
fict(n, t) = real(-net.Y_base / lines[n].Y[1, 2]) * (g(:lrlr, n, t) + g(:lili, n, t) - g(:lr, n, t)^2 - g(:li, n, t)^2) * S

F = Dict(l.name => [fict(l.name, t) for t in 1:T] for l in net.lines)
open(joinpath(GD, "gap_hours_$(d.id).csv"), "w") do io
    println(io, "hour,fict_loop_kW,fict_radial_kW,battery_net_kW,max_cycle_mismatch_deg")
    for t in 1:T
        fl = sum(F[n][t] for n in onloop; init = 0.0); fr = sum(F[l.name][t] for l in net.lines) - fl
        pb = sum(v(R.P_B[kk, t]) for kk in axes(R.P_B, 1)) * S
        mm = maximum((abs(rad2deg(wrap(sum(sg * thW(ln, t) for (ln, sg) in cy)))) for (_, cy) in cycles); init = 0.0)
        @printf(io, "%d,%.3f,%.3f,%.1f,%.5f\n", t, fl, fr, pb, mm)
    end
end
open(joinpath(GD, "gap_cycles_$(d.id).csv"), "w") do io
    println(io, "closing_line,length,max_mismatch_deg,hour_of_max,fict_day_kWh")
    for (cl, cy) in cycles
        mm = [abs(rad2deg(wrap(sum(sg * thW(ln, t) for (ln, sg) in cy)))) for t in 1:T]
        @printf(io, "%s,%d,%.5f,%d,%.3f\n", cl, length(cy), maximum(mm), argmax(mm), sum(sum(F[ln]) for (ln, _) in cy) * 24 / T)
    end
end
srcs = net.sources
open(joinpath(GD, "gap_sources_$(d.id).csv"), "w") do io
    println(io, "from,to,path_lines,max_product_minus_moment_deg,hour_of_max")
    for a in 1:length(srcs), b in a+1:length(srcs)
        p = treepath(srcs[a].bus, srcs[b].bus)
        dd = [abs(rad2deg(wrap(sum(sg * (thW(ln, t) - thV(ln, t)) for (ln, sg) in p)))) for t in 1:T]
        @printf(io, "%s,%s,%d,%.5f,%d\n", srcs[a].name, srcs[b].name, length(p), maximum(dd), argmax(dd))
    end
end
top = sort([(sum(F[l.name]) * 24 / T, l.name) for l in net.lines]; rev = true)[1:min(200, end)]
open(joinpath(GD, "gap_lines_$(d.id).csv"), "w") do io
    println(io, "line,on_loop,fict_day_kWh,fict_peak_kW,hour_of_peak,max_angle_gap_deg")
    for (fd, n) in top
        ag = maximum(abs(rad2deg(wrap(thW(n, t) - thV(n, t)))) for t in 1:T)
        @printf(io, "%s,%s,%.3f,%.3f,%d,%.5f\n", n, n in onloop, fd, maximum(F[n]), argmax(F[n]), ag)
    end
end
tot = sum(sum(x) for x in values(F)) * 24 / T; lp = sum(sum(F[n]) for n in onloop; init = 0.0) * 24 / T
say(@sprintf("  fictitious loss over the day %.1f kWh: %.1f%% on loop lines (%d of %d lines)", tot, 100lp / max(tot, eps()), length(onloop), length(net.lines)))
LAT = lateral_boxes(net, kwL, R.bounds, Dict(s.bus => s for s in net.sources); ideal = get(kwL, :ideal_sources, false)).lat
la = sum((sum(F[n]) for n in keys(LAT)); init = 0.0) * 24 / T
say(@sprintf("  %.1f%% on lines with a source-free side (%d lines, where :lsec applies); %.1f kWh elsewhere", 100la / max(tot, eps()), length(LAT), tot - la))
say("  top lines: ", join((@sprintf("%s%s %.1f", n, n in onloop ? "*" : "", fd) for (fd, n) in top[1:min(10, end)]), ", "))
say(@sprintf("  done in %.0f s; files in %s", time() - t0, GD))
