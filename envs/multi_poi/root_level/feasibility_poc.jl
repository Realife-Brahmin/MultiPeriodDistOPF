# Proof of concept: exact AC dispatch (lower bound L) against a sparse SDP outer relaxation
# (upper bound U) for the clean maximum-margin problem at fixed angles (feasibility_certificate.jl).
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/feasibility_poc.jl debug     # T = 2, code check only
#   julia --project=envs/msopf_sdp envs/multi_poi/root_level/feasibility_poc.jl poc123    # the T = 24 experiment
# poc123: IEEE 123-node, sources 1, 3, 4, PV reactive power dispatched, batteries idle and dispatched.
# Directions in (delta3 - delta1, delta4 - delta1): 292.5 deg (the substation-4 face, measured well
# beyond the prediction) and 157.5 deg (the substation-3 face). Per battery case: equal angles,
# 98% and 102% of the measured reach in each direction (simplex_check.csv), and one point
# outside the predicted simplex but inside the measured region on the substation-4 side.
# FEAS_CUTS=rlt adds the McCormick cuts to the relaxation. Writes feasibility_poc_<label>.csv.
const POC_ARGS = copy(ARGS)
include(joinpath(@__DIR__, "..", "multiperiod_full_angle.jl"))
include(joinpath(@__DIR__, "..", "feasibility_certificate.jl"))
using Logging

# The setup helpers of multiperiod_study.jl, copied (that driver loads Plots, which envs/msopf_sdp leaves out).
const OUT_DIR = normpath(get(ENV, "MULTIPERIOD_OUT_DIR", joinpath(@__DIR__, "..", "results", "multiperiod")))
mkpath(OUT_DIR)
say(args...) = (println(args...); flush(stdout))
const RAY_OK = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
function load_sources(srcs; system = "ieee123_5poi_1ph", extra = String[])
    kw = (; sources = collect(srcs), snapshot = true, load_band = (0.5, 1.5), extra)
    pvs = with_logger(NullLogger()) do
        read_pv(system; kw...)
    end
    compile_deck(system; disable = ["PVSystem"], kw...)
    return read_network(), pvs
end
function setup(net, pvs; T, pv = "qfree", ideal = false, feas = nothing)
    lam, c, pvf = load_profile(T), price_profile(T), pv_profile(T)
    P, F = pv == "off" ? ([], zeros(T)) : (pvs, pvf)
    ipopt = feas === nothing ? (; tol = IPOPT_TOL, constr_viol_tol = IPOPT_TOL) : (; tol = feas / 10, constr_viol_tol = feas)
    return (; net, T, lam, c, P, F, q = pv == "unity" ? :unity : :optimize, CB = cb_by_rule(net), GAM = gamma_by_rule(net), ideal, feas, ipopt)
end
const AMP = haskey(ENV, "FEAS_AMP") ? parse(Float64, ENV["FEAS_AMP"]) : nothing   # FEAS_AMP=5: every line at 5 pu
function problem(s, case)
    common = (; ampacity_pu = AMP, pvs = s.P, pv_q = s.q, C_B = s.CB, ideal_sources = s.ideal, battery = case === :nobatt ? :idle : :optimize, s.ipopt...)
    return (; common..., T = s.T, load_factor = s.lam, price = s.c, pv_factor = s.F, gamma = s.GAM, terminal = case === :cyclic ? :cyclic : :soft)
end

const KV = 20000.0                         # kW of margin per pu of |V|^2, as MULTI_VMARGIN
const UTOL = 0.01                          # kW: U must be below -UTOL to call a point infeasible
const DIRS = [(th = 292.5, face = "subs4"), (th = 157.5, face = "subs3")]
const HEADER = "label,case,query,r_deg,theta_deg,d3_deg,d4_deg,L_kW,nlp_status,nlp_iter,nlp_margin_kW,orig_max_violation,orig_equations,U_kW,relax_status,relax_primal_kW,relax_dual_kW,relax_gap_kW,relax_primal_residual,psd_primal_min_eig,psd_dual_violation,dual_sign_violation,box_term_kW,rounding_kW,rank_ratio,recovered_max_violation,recovered_margin_kW,class,nlp_build_s,nlp_solve_s,relax_build_s,relax_solve_s,opendss"

function measured_reach()
    d = Dict{Tuple{String,Float64},Tuple{Float64,Float64}}()
    for l in readlines(joinpath(@__DIR__, "..", "results", "multiperiod", "simplex_check.csv"))[2:end]
        r = split(l, ','); (r[1] == "three" && r[2] == "qfree") || continue
        d[(String(r[4]), parse(Float64, r[3]))] = (parse(Float64, r[5]), parse(Float64, r[6]))
    end
    return d
end

classify(L, U) = isfinite(L) && L >= 0 ? "feasible" : (isfinite(U) && U < -UTOL ? "infeasible" : "unresolved")

function run_queries(label, s, cases, queries_for; cuts = Symbol[], replay = true)
    out = joinpath(OUT_DIR, "feasibility_poc_$(label).csv")
    open(io -> println(io, HEADER), out, "w")
    for case in cases
        kw = problem(s, case)
        prev = nothing
        for q in queries_for(case)
            ang = Dict("subs1" => 0.0, "subs3" => q.r * cosd(q.th), "subs4" => q.r * sind(q.th))
            t0 = time()
            x = solve_multiperiod(s.net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0,
                                  start = q.r == 0 ? nothing : prev)
            tn = time() - t0
            ok = x.status in RAY_OK
            v = ok ? validate_dispatch(s.net, ang, kw, x; k = KV) : nothing
            L = ok && v.equations <= 1e-6 ? v.margin_kW : NaN
            ok && (prev = x)
            c = certify(s.net, ang, kw; k = KV, cuts)
            vr = c.recovered === nothing ? nothing : validate_dispatch(s.net, ang, kw, c.recovered; k = KV)
            lv = ok && isfinite(L) ? lifted_violation(c.R, s.net, x, kw; margin_kW = max(L, -999.0)) : (; worst = NaN, at = "")
            say(@sprintf("    validity: the nonlinear solution lifted into the relaxation violates it by %.2e (%s)", lv.worst, first(lv.at, 120)))
            cls = classify(L, c.U)
            rp = ""
            if replay && ok && isfinite(L) && L >= 0 && q.name in ("equal", "98%")
                rr = replay_opendss(s.net, x, ang; load_factor = s.lam, pvs = s.P, pv_factor = s.F)
                rp = @sprintf("dP %.1e kW; dQ %.1e kvar; dV %.1e pu; converged %s", rr.dP, rr.dQ, rr.dV, rr.converged)
            end
            row = [label, string(case), q.name, @sprintf("%.6f", q.r), @sprintf("%.1f", q.th), @sprintf("%.6f", ang["subs3"]), @sprintf("%.6f", ang["subs4"]),
                   @sprintf("%.4f", L), string(x.status), string(something(x.iterations, -1)), @sprintf("%.4f", x.margin_kW),
                   ok ? @sprintf("%.2e", v.max_violation) : "NaN", ok ? @sprintf("%.2e", v.equations) : "NaN",
                   @sprintf("%.4f", c.U), string(c.status), @sprintf("%.4f", c.primal), @sprintf("%.4f", c.dualobj), @sprintf("%.2e", c.gap),
                   @sprintf("%.2e", c.primal_residual), @sprintf("%.2e", c.psd_primal_min_eig), @sprintf("%.2e", c.psd_dual_violation),
                   @sprintf("%.2e", c.dual_sign_violation), @sprintf("%.3e", c.box_term_kW), @sprintf("%.1e", c.rounding_kW), @sprintf("%.2e", c.rank_ratio),
                   vr === nothing ? "NaN" : @sprintf("%.2e", vr.max_violation), vr === nothing ? "NaN" : @sprintf("%.4f", vr.margin_kW),
                   cls, @sprintf("%.1f", x.time === nothing ? NaN : tn - x.time), @sprintf("%.1f", something(x.time, NaN)),
                   @sprintf("%.1f", c.build_s), @sprintf("%.1f", c.solve_s), rp]
            open(io -> println(io, join(row, ',')), out, "a")
            say(@sprintf("  %-6s %-9s r %.4f th %5.1f | L %+9.3f kW (%s, viol %.1e) | U %+9.3f kW (%s, gap %.1e, psd %.1e, rank %.1e) | %s | nlp %.1f s, relax %.1f+%.1f s %s",
                         case, q.name, q.r, q.th, L, x.status, ok ? v.max_violation : NaN, c.U, c.status, c.gap, c.psd_dual_violation, c.rank_ratio,
                         cls, tn, c.build_s, c.solve_s, rp))
        end
    end
    say("wrote ", out)
end

mode = isempty(POC_ARGS) ? "debug" : POC_ARGS[1]
cuts = Symbol.(split(get(ENV, "FEAS_CUTS", ""), ','; keepempty = false))
net, pvs = load_sources(["subs1", "subs3", "subs4"])
if mode == "debug"                                    # T = 2: code check, kept apart from the experiment
    s = setup(net, pvs; T = 2, pv = "qfree")
    run_queries("debug_T2", s, (:nobatt, :cyclic), case -> [(name = "equal", r = 0.0, th = 0.0), (name = "far", r = 2.0, th = 292.5)]; cuts, replay = false)
elseif mode == "probe"                                # one point, verbose: the validator's breakdown and Clarabel's log
    T = parse(Int, get(ENV, "PROBE_T", "24")); case = Symbol(get(ENV, "PROBE_CASE", "nobatt"))
    s = setup(net, pvs; T, pv = "qfree"); kw = problem(s, case)
    ang = Dict("subs1" => 0.0, "subs3" => 0.0, "subs4" => 0.0)
    x = solve_multiperiod(s.net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
    v = validate_dispatch(s.net, ang, kw, x; k = KV)
    say("nlp: ", x.status, " margin ", round(x.margin_kW, digits = 4), " kW; validated margin ", round(v.margin_kW, digits = 4), " kW; violations ", v.viol)
    opts = [(String(a), something(tryparse(Int, b), tryparse(Float64, b), b == "true" ? true : b == "false" ? false : String(b)))
            for (a, b) in (split(kv, '=') for kv in split(get(ENV, "CLARABEL_OPTS", ""), ';'; keepempty = false))]
    c = certify(s.net, ang, kw; k = KV, cuts, solver_attributes = [CLARABEL_DEFAULTS; opts])
    say(@sprintf("relaxation: %s primal %.3f dual %.3f kW, U %.3f kW, box term %.3e kW, rank %.2e, fictitious losses %.1f kW per period, max |V|^2 gap %.3e, build %.1f s, solve %.1f s",
                 c.status, c.primal, c.dualobj, c.U, c.box_term_kW, c.rank_ratio, c.fictitious_loss_kW, c.vsq_gap, c.build_s, c.solve_s))
elseif mode == "exotic"                               # other branches of the AC equations: scrambled starts
    T = parse(Int, get(ENV, "EXOTIC_T", "2")); nstart = parse(Int, get(ENV, "EXOTIC_STARTS", "12"))
    r, th = parse(Float64, get(ENV, "EXOTIC_R", "2.0")), parse(Float64, get(ENV, "EXOTIC_TH", "292.5"))
    s = setup(net, pvs; T, pv = "qfree")
    using Random
    rng = MersenneTwister(1)
    for case in (:nobatt,)
        kw = problem(s, case)
        ang = Dict("subs1" => 0.0, "subs3" => r * cosd(th), "subs4" => r * sind(th))
        base = solve_multiperiod(s.net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0)
        vb = validate_dispatch(s.net, ang, kw, base; k = KV)
        say(@sprintf("T = %d, r = %.3f deg at %.1f deg, %s: flat start: margin %+.1f kW (%s), largest line current %.2f pu", T, r, th, case,
                     vb.margin_kW, base.status, maximum(maximum(abs.(complex.(base.raw.lr[l], base.raw.li[l]))) for l in keys(base.raw.lr))))
        for k in 1:nstart
            spread = k <= nstart ÷ 2 ? pi / 4 : pi                 # half the starts mildly scrambled, half fully
            V = Dict(b => [base.V_pu[b][t] * cis(spread * (2rand(rng) - 1)) for t in 1:T] for b in s.net.buses)
            x = solve_multiperiod(s.net, ang; kw..., objective = :margin, margin_band_kw = KV, margin_reg = 0.0,
                                  start = start_from_voltages(s.net, ang, kw, V))
            ok = x.status in RAY_OK
            v = ok ? validate_dispatch(s.net, ang, kw, x; k = KV) : nothing
            Imax = ok ? maximum(maximum(abs.(complex.(x.raw.lr[l], x.raw.li[l]))) for l in keys(x.raw.lr)) : NaN
            loss = ok ? sum(sum(x.P_subs_kW[n]) for n in keys(x.P_subs_kW)) - sum(sum(base.P_subs_kW[n]) for n in keys(base.P_subs_kW)) : NaN
            say(@sprintf("  start %2d (angles +-%3.0f deg): %-22s margin %+10.1f kW, equations %.1e, largest current %8.2f pu, imports %+.0f kW vs flat",
                         k, rad2deg(spread), x.status, ok ? v.margin_kW : NaN, ok ? v.equations : NaN, Imax, loss))
        end
    end
elseif mode == "poc123"
    s = setup(net, pvs; T = 24, pv = "qfree")
    reach = measured_reach()
    function queries(case)
        qs = [(name = "equal", r = 0.0, th = 0.0)]
        for d in DIRS
            meas, pred = reach[(string(case), d.th)]
            append!(qs, [(name = "98%", r = 0.98meas, th = d.th), (name = "102%", r = 1.02meas, th = d.th)])
            d.face == "subs4" && push!(qs, (name = "beyond_pred", r = (meas + pred) / 2, th = d.th))
        end
        return qs
    end
    run_queries("ieee123_T24" * (AMP === nothing ? "" : "_amp" * string(round(Int, AMP))) * (isempty(cuts) ? "" : "_" * join(cuts, "+")), s, (:nobatt, :cyclic), queries; cuts)
else
    error("unknown mode $mode")
end
