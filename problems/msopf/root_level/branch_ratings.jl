# Current ratings of every line of a converted deck (user, 2026-11-02: ampacity limits in the MSOPF campaign), as a
# limit in the deck's own per-phase-equivalent units, and a check of those units against the original 3-phase deck.
#   julia --project=problems/msopf problems/msopf/root_level/branch_ratings.jl ieee9500xf smartds_p1u_xf
# Writes test_networks/<deck>/BranchRatings.csv and results/ratings_nov02/<deck>_check.csv.
#
# Per-phase equivalent (convert_to_1ph.jl): every load and capacitor / 3, a k-phase section (3/k) x its positive-sequence
# impedance, everything referred to the base voltage. A section carrying S on its k phases carries S/3 in the model, so
#   line:        model current = real phase current x k/3 x (kV of the section / base kV); limit = normamps x k/3 x kV/base
#   transformer: rated S (kVA, winding 1) on its k phases is S/3 in the model; limit = kVA / (sqrt(3) x base kV LL) A,
#                whatever k and the voltage level
#   parallel units combined into one line: sum of their limits; a tie made of switches and a cable: the smallest
#   no rating in the original (switches and jumpers given by impedance alone, reactors, a cable whose conductor data
#   carry no normamps): the conductor it continues, i.e. at each end the largest rated line, then the smaller end.
# Service transformers are rated here, but the campaign leaves them unlimited: they carry only their own customers'
# load, which no dispatch or angle changes, and at the datasets' peak some exceed nameplate (up to 1.5x on sf6937).
# The check solves the original (3-phase, its own controls, our load scale) and the converted deck (snapshot, no PV or
# batteries) and compares each rated line's loading, current / rating: a factor-of-k slip in the units shows as a ratio
# of 3 or 1/3 on single-phase sections.
include(joinpath(@__DIR__, "convert_to_1ph.jl"))
using Statistics, Logging
global_logger(ConsoleLogger(stderr, Logging.Error))      # ODD warns on every "?" query that returns text

"Text index of the original's .dss files: linecodes and conductors with normamps, geometries' conductors, lines with normamps."
function rating_index(dir)
    defs = String[]
    for (r, _, fs) in walkdir(dir), f in fs
        endswith(lowercase(f), ".dss") || continue
        for l in eachline(joinpath(r, f))
            l = strip(split(split(l, '!')[1], "//")[1]); isempty(l) && continue
            startswith(l, "~") && !isempty(defs) ? (defs[end] *= " " * l[2:end]) : push!(defs, l)
        end
    end
    lc, cond, geo, ln = Set{String}(), Set{String}(), Dict{String,Vector{String}}(), Set{String}()
    for d in defs
        m = match(r"(?i)^new\s+(linecode|wiredata|cndata|tsdata|linegeometry|line)\.(\S+)", d); m === nothing && continue
        cls, nm = lowercase(m[1]), lowercase(m[2]); has = occursin(r"(?i)\bnormamps\s*=", d)
        cls == "linecode" && has && push!(lc, nm)
        cls in ("wiredata", "cndata", "tsdata") && has && push!(cond, nm)
        cls == "line" && has && push!(ln, nm)
        if cls == "linegeometry"
            ws = [lowercase(x[1]) for x in eachmatch(r"(?i)\b(?:wire|cncable|tscable)\s*=\s*([^\s\]]+)", d)]
            for x in eachmatch(r"(?i)\b(?:wires|cncables|tscables)\s*=\s*[\[\(\"']([^\]\)\"']*)", d); append!(ws, lowercase.(split(x[1]))); end
            geo[nm] = ws
        end
    end
    return (; lc, cond, geo, ln)
end

q(cmd) = strip(ODD.Text.Command("? " * cmd))

function ratings(name)
    c = convert(name); cfg = c.cfg; vref = cfg.base_kV
    Ib = 1000 / (vref / sqrt(3))                                    # A per pu: S_base 1000 kVA at the base phase voltage
    idx = rating_index(dirname(cfg.master))
    lines = Dict(clean(n) => n for n in real_names(ODD.Lines.AllNames()))
    xf = Dict{String,String}()
    for n in real_names(ODD.Transformers.AllNames()), p in ("svc_", "xfmr_", "reg_"); xf[p * clean(n)] = n; end
    # one original element: (limit in model amperes or nothing, how, the original's figure, phases)
    function member(m)
        if haskey(lines, m)
            n = lines[m]; ODD.Circuit.SetActiveElement("Line.$n"); k = ODD.CktElement.NumPhases()
            lcode, g = lowercase(q("Line.$n.linecode")), lowercase(q("Line.$n.geometry"))
            rated = lowercase(n) in idx.ln || (!isempty(lcode) && lcode in idx.lc) ||
                    (!isempty(g) && haskey(idx.geo, g) && !isempty(idx.geo[g]) && all(w -> w in idx.cond, idx.geo[g]))
            rated || return (nothing, "line without normamps", NaN, k)
            # the section's own voltage level, from an energized end: a normally-open tie's buses are dead and OpenDSS gives
            # them the lowest base on its list (0.12 kV on SMART-DS, which made tie limits 1/104 of their rating, 2026-11-02)
            ODD.Lines.Name(n); a = ODD.Lines.NormAmps()
            ends = filter(b -> (ODD.Circuit.SetActiveBus(b); maximum(abs, ODD.Bus.Voltages(); init = 0.0) > 1e-9), bus.([ODD.Lines.Bus1(), ODD.Lines.Bus2()]))
            kvb = isempty(ends) ? vref : kvll(ends[1])
            return (a * k / 3 * kvb / vref, "normamps", a, k)
        elseif haskey(xf, m)
            n = xf[m]; ODD.Transformers.Name(n); ODD.Transformers.Wdg(1.0); s = ODD.Transformers.kVA()
            ODD.Circuit.SetActiveElement("Transformer.$n")
            return (s / (sqrt(3) * vref), "kva", s, ODD.CktElement.NumPhases())
        end
        return (nothing, startswith(m, "reactor_") ? "reactor" : "not found", NaN, 0)
    end
    rows = []
    for (set, istie) in ((c.branches, false), (c.ties, true))
        for x in set
            ms = [member(m) for m in x.members]
            lim = any(isnothing(r[1]) for r in ms) ? nothing : istie ? minimum(r[1] for r in ms) : sum(r[1] for r in ms)
            push!(rows, Dict{Symbol,Any}(:name => x.name, :kind => istie ? "tie" : String(x.kind), :a => x.a, :b => x.b, :members => join(x.members, "|"),
                             :how => join(unique(r[2] for r in ms), "|"), :orig => join([@sprintf("%.6g", r[3]) for r in ms], "|"),
                             :phases => join(unique(r[4] for r in ms), "|"), :A => lim, :source => lim === nothing ? "" : "original"))
        end
    end
    # unrated (a switch or jumper in series with a conductor): at each end the largest rated line there (the conductor
    # it continues; laterals and service transformers at the same bus are smaller), then the smaller of its two ends,
    # an end with no rated line left out; repeated for chains of unrated elements. (The smallest neighbour of all,
    # used first, gave trunk switches a lateral's or a service transformer's 1-3 A, 2026-11-02.)
    at = Dict{String,Vector{Int}}()
    for (i, r) in enumerate(rows); push!(get!(at, r[:a], Int[]), i); push!(get!(at, r[:b], Int[]), i); end
    for _ in 1:100
        changed = false
        for (i, r) in enumerate(rows)
            r[:A] === nothing || continue
            side(b) = (v = [rows[j][:A] for j in at[b] if j != i && rows[j][:A] !== nothing && rows[j][:kind] ∉ ("service", "tie")]; isempty(v) ? nothing : maximum(v))
            nb = filter(!isnothing, [side(r[:a]), side(r[:b])])
            isempty(nb) && continue
            r[:A] = minimum(nb); r[:source] = "inherited"; changed = true
        end
        changed || break
    end
    out = joinpath(REPO, "test_networks", deckname(c), "BranchRatings.csv")
    open(out, "w") do io
        println(io, "name,kind,members,phases,rated_by,original,limit_A,limit_pu,source")
        for r in sort(rows, by = r -> r[:name])
            @printf(io, "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", r[:name], r[:kind], r[:members], r[:phases], r[:how], r[:orig],
                    r[:A] === nothing ? "" : @sprintf("%.6g", r[:A]), r[:A] === nothing ? "" : @sprintf("%.6g", r[:A] / Ib), r[:source])
        end
    end
    cnt(f) = count(f, rows)
    @printf("%s: %d lines and ties; rated in the original %d, inherited %d, none %d; base %.2f A per pu\n", deckname(c), length(rows),
            cnt(r -> r[:source] == "original"), cnt(r -> r[:source] == "inherited"), cnt(r -> r[:A] === nothing), Ib)
    for h in sort(unique(r[:how] for r in rows if r[:source] != "original"))
        @printf("  not rated in the original, %s: %d\n", h, cnt(r -> r[:how] == h && r[:source] != "original"))
    end
    println("  wrote ", out)
    return c, rows, Ib
end

"Loading (current / rating) of each rated single-element line, in the original 3-phase deck and in the converted deck."
function check(c, rows, Ib)
    cfg = c.cfg; scale = peak_scale(cfg)
    compile_original(cfg)
    if scale !== nothing                                            # the converted deck carries the dataset's peak scale
        for n in real_names(ODD.Loads.AllNames())
            ODD.Loads.Name(n); k = get(scale, lowercase(n), 1.0); ODD.Loads.kW(ODD.Loads.kW() * k); ODD.Loads.kvar(ODD.Loads.kvar() * k)
        end
        ODD.Solution.Solve()
    end
    lines = Dict(clean(n) => n for n in real_names(ODD.Lines.AllNames()))
    orig = Dict{String,Tuple{Float64,Float64,Int,Float64}}()   # name -> (largest phase current A, normamps, phases, mean phase current A)
    for r in rows
        r[:source] == "original" && r[:kind] in ("line", "short") && !occursin('|', r[:members]) || continue
        n = lines[r[:members]]; ODD.Circuit.SetActiveElement("Line.$n"); k = ODD.CktElement.NumPhases()
        I = ODD.CktElement.CurrentsMagAng()[1:2:2k]; ODD.Lines.Name(n)
        orig[r[:name]] = (maximum(I), ODD.Lines.NormAmps(), k, sum(I) / k)
    end
    deck = joinpath(REPO, "test_networks", deckname(c))
    ODD.Text.Command("Clear"); ODD.Text.Command("Redirect \"$(joinpath(deck, "Master.dss"))\"")
    for cls in ("PVSystem", "Storage")
        any(startswith(lowercase(e), lowercase(cls) * ".") for e in ODD.Circuit.AllElementNames()) && ODD.Text.Command("BatchEdit $cls..* enabled=false")
    end
    ODD.Text.Command("Set mode=snapshot"); ODD.Solution.Solve()
    ODD.Solution.Converged() || error("converted deck did not converge")
    lim = Dict(r[:name] => r[:A] for r in rows)
    res = []
    for (nm, (Io, a, k, Iavg)) in orig
        ODD.Circuit.SetActiveElement("Line.$nm"); Im = ODD.CktElement.CurrentsMagAng()[1]
        Io > 0.05 * a || continue                                    # lightly loaded lines say little
        push!(res, (; nm, k, Io, a, Im, lim = lim[nm], lo = Io / a, lm = Im / lim[nm], ravg = Im / (Iavg * k / 3)))
    end
    mkpath(joinpath(REPO, "problems", "msopf", "results", "ratings_nov02"))
    out = joinpath(REPO, "problems", "msopf", "results", "ratings_nov02", deckname(c) * "_check.csv")
    open(out, "w") do io
        println(io, "line,phases,original_A,normamps,model_A,model_limit_A,loading_original,loading_model,ratio,model_over_mean_phase_x_k_3")
        for x in sort(res, by = x -> x.nm)
            @printf(io, "%s,%d,%.4f,%.4g,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f\n", x.nm, x.k, x.Io, x.a, x.Im, x.lim, x.lo, x.lm, x.lm / x.lo, x.ravg)
        end
    end
    for k in 1:3
        r = [x.lm / x.lo for x in res if x.k == k]; isempty(r) && continue
        q5, q95 = quantile(r, [0.05, 0.95])
        @printf("  check, %d-phase lines above 5%% loading: %d; loading model / original median %.4f (5-95%%: %.4f-%.4f)\n", k, length(r), median(r), q5, q95)
        ra = [x.ravg for x in res if x.k == k]
        # the model carries the mean of the phase currents (a balanced equivalent), so this ratio is the unit check
        @printf("    model current / (mean phase current x k/3): median %.4f (5-95%%: %.4f-%.4f)\n", median(ra), quantile(ra, [0.05, 0.95])...)
    end
    @printf("  largest loading at nominal load: original %.3f, model %.3f\n", maximum(x.lo for x in res), maximum(x.lm for x in res))
    println("  wrote ", out)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    for name in ARGS
        c, rows, Ib = ratings(name)
        check(c, rows, Ib)
    end
end
