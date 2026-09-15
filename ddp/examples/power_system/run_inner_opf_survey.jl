# Feasibility survey for the reduced-space MPOPF decomposition.
#
# Question under test:
#     does every admissible battery dispatch P_B^t admit a feasible
#     single-period network solution?
#
# For a chosen system/horizon this fixes P_B^t to a battery of structured and
# random dispatch patterns at low-, medium- and high-demand hours, solves the
# inner network OPF with IPOPT, and records everything the study needs. Any
# infeasible case is preserved as infeasible and then separately diagnosed with
# the violation-slack restoration model (never relabelled feasible).
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/run_inner_opf_survey.jl [system] [T]

using Printf
using Random
using Serialization
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))
include(joinpath(@__DIR__, "ieee123c_filterddp.jl"))   # build_model, control_layout

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 24
n_random = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 24

data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T$(T).jls"))
refpath = joinpath(REPO, "envs", "tadmm", "processedData", "$(system)_T$(T)", "sol_socp_bf.jls")
ref = isfile(refpath) ? deserialize(refpath) : nothing
Bset = data[:Bset]
nB = length(Bset)
pbr = Float64[data[:P_B_R_pu][j] for j in Bset]

# ---------------------------------------------------------------------------
# time samples: low / medium / high NET demand (load minus PV)
# ---------------------------------------------------------------------------
net = [sum(data[:p_L_pu][j, t] for j in data[:NLset]) -
       sum(data[:p_D_pu][j, t] for j in data[:Dset]) for t in 1:T]
ord = sortperm(net)
t_low  = ord[1]
t_high = ord[end]
t_med  = ord[cld(length(ord), 2)]
samples = [(t_low, "low"), (t_med, "medium"), (t_high, "high")]

# ---------------------------------------------------------------------------
# battery groups by electrical depth from the substation
# ---------------------------------------------------------------------------
root = data[:substationBus]
depth = Dict{Int,Int}(root => 0)
queue = [root]
while !isempty(queue)
    j = popfirst!(queue)
    for k in get(data[:children], j, Int[])
        depth[k] = depth[j] + 1; push!(queue, k)
    end
end
bdepth = [depth[j] for j in Bset]
med_depth = median(bdepth)
deep    = findall(>=(med_depth), bdepth)      # far from substation
shallow = findall(<(med_depth), bdepth)
q4cut   = quantile(bdepth, 0.75)
deepest = findall(>=(q4cut), bdepth)          # deepest quartile

@printf("SURVEY system=%s T=%d batteries=%d  deep=%d shallow=%d deepest_quartile=%d\n",
        system, T, nB, length(deep), length(shallow), length(deepest))
@printf("SURVEY time samples: low t=%d (net=%.4f)  medium t=%d (net=%.4f)  high t=%d (net=%.4f)\n",
        t_low, net[t_low], t_med, net[t_med], t_high, net[t_high])

# ---------------------------------------------------------------------------
# dispatch patterns.  sign convention: +P_B = discharge (injection)
# ---------------------------------------------------------------------------
function patterns_for(t::Int)
    pats = Tuple{String,Int,Vector{Float64}}[]   # (name, seed, pb)
    push!(pats, ("zero", -1, zeros(nB)))
    if ref !== nothing
        push!(pats, ("centralized_optimal", -1, Float64[ref[:P_B][j, t] for j in Bset]))
    end
    push!(pats, ("all_max_charge", -1, -copy(pbr)))
    push!(pats, ("all_max_discharge", -1, copy(pbr)))
    for (gname, gidx) in (("deep", deep), ("shallow", shallow))
        for (sname, s) in (("charge", -1.0), ("discharge", 1.0))
            v = zeros(nB); v[gidx] .= s .* pbr[gidx]
            push!(pats, ("$(gname)_max_$(sname)", -1, v))
        end
    end
    v = zeros(nB); v[deep] .= pbr[deep]; v[shallow] .= -pbr[shallow]
    push!(pats, ("opposing_deep_dis_shallow_chg", -1, copy(v)))
    v = zeros(nB); v[deep] .= -pbr[deep]; v[shallow] .= pbr[shallow]
    push!(pats, ("opposing_deep_chg_shallow_dis", -1, copy(v)))
    # structured stress directions
    v = zeros(nB); v[deepest] .= -pbr[deepest]
    push!(pats, ("struct_max_voltage_drop", -1, copy(v)))     # load concentrated at feeder end
    v = zeros(nB); v[deepest] .= pbr[deepest]
    push!(pats, ("struct_max_voltage_rise", -1, copy(v)))     # injection at feeder end
    vv = zeros(nB)
    for (rank, i) in enumerate(sortperm(bdepth))
        vv[i] = iseven(rank) ? pbr[i] : -pbr[i]
    end
    push!(pats, ("struct_alternating_circulation", -1, vv))
    # reproducible random dispatches throughout the box
    for s in 1:n_random
        rng = MersenneTwister(1000 + s)
        push!(pats, ("random", 1000 + s, [(2rand(rng) - 1) * pbr[b] for b in 1:nB]))
    end
    return pats
end

# ---------------------------------------------------------------------------
# FilterDDP cross-validation model (independent code path)
# ---------------------------------------------------------------------------
print("building FilterDDP OCP for independent validation ... "); flush(stdout)
ocp, idx, nx_, nu_, nc_ = build_model(data)
println("done")

# entering state for the (secondary) energy-row variant and for validation
x_entering(t) = t == 1 ? Float64[data[:B0_pu][j] for j in Bset] :
                (ref === nothing ? Float64[data[:B0_pu][j] for j in Bset] :
                 Float64[ref[:B][j, t-1] for j in Bset])

outdir = joinpath(REPO, "ddp", "results", "reduced_space")
mkpath(outdir)
csvpath = joinpath(outdir, "inner_opf_survey_$(system)_T$(T).csv")
io = open(csvpath, "w")
println(io, "system,horizon,time_index,demand_label,net_load_pu,price,pattern,seed," *
            "sum_pb_pu,min_pb_pu,max_pb_pu,status,feasible,iterations,objective_usd," *
            "substation_cost_usd,battery_term_usd,P_Subs_pu,Q_Subs_pu,loss_pu," *
            "reverse_export,ps_margin_pu,vmin_pu,vmax_pu,ell_max_pu,imag_max_pu," *
            "qnorm_absmax,qnorm_absmean,n_active_qnorm,n_active_vmin,n_active_vmax," *
            "min_bound_margin,filterddp_residual,validation_pass,soc_box_ok," *
            "solve_time_s,alloc_mib,rss_delta_mib")

nrows = 0; n_infeasible = 0; min_ps = Inf; worst_val = 0.0
infeasible_cases = Tuple{Int,String,Int,Vector{Float64}}[]

for (t, label) in samples
    xent = x_entering(t)
    width = [(data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j] for j in Bset]
    emin  = [data[:soc_min][j] * data[:B_R_pu][j] for j in Bset]
    for (name, seed, pb) in patterns_for(t)
        res = inner_opf(data, t, pb)
        resid = NaN; vpass = false
        if res.feasible
            resid = filterddp_residual(ocp, idx, data, t, res, pb, xent)
            vpass = resid < 1e-6
            global worst_val = max(worst_val, resid)
            global min_ps = min(min_ps, res.ps)
        else
            global n_infeasible += 1
            push!(infeasible_cases, (t, name, seed, pb))
        end
        # SOC box check is an OUTER battery question, reported alongside
        es = xent .- data[:delta_t_h] .* pb .- emin
        soc_ok = all(-1e-9 .<= es .<= width .+ 1e-9)
        @printf(io, "%s,%d,%d,%s,%.6f,%.6f,%s,%d,%.6f,%.6f,%.6f,%s,%s,%d,%.6f,%.6f,%.6e,%.6f,%.6f,%.6f,%s,%.6f,%.6f,%.6f,%.6e,%.6f,%.6f,%.6f,%d,%d,%d,%.3e,%.3e,%s,%s,%.4f,%.1f,%.1f\n",
            system, T, t, label, net[t], data[:LoadShapeCost][t], name, seed,
            sum(pb), minimum(pb), maximum(pb), res.status, res.feasible, res.iterations,
            res.objective, res.substation_cost, res.battery_term, res.ps, res.qs, res.loss,
            res.reverse_export, res.ps_margin, res.vmin, res.vmax, res.ell_max, res.imag_max,
            res.qnorm_absmax, res.qnorm_absmean,
            res.feasible ? res.n_active_qnorm : -1,
            res.feasible ? res.n_active_vmin : -1,
            res.feasible ? res.n_active_vmax : -1,
            res.min_bound_margin, resid, vpass, soc_ok,
            res.solve_time, res.alloc_mib, res.rss_delta_mib)
        global nrows += 1
    end
    @printf("  t=%2d (%-6s) done\n", t, label); flush(stdout)
end
close(io)

@printf("SURVEY rows=%d infeasible=%d min_P_Subs=%.6f pu worst_validation_residual=%.3e\n",
        nrows, n_infeasible, min_ps, worst_val)
println("SURVEY wrote=$csvpath")

# ---------------------------------------------------------------------------
# reverse-export stress: all-max-discharge at EVERY hour, closest approach to 0
# ---------------------------------------------------------------------------
expath = joinpath(outdir, "reverse_export_scan_$(system)_T$(T).csv")
open(expath, "w") do exio
    println(exio, "system,horizon,time_index,net_load_pu,total_discharge_pu,P_Subs_pu,loss_pu,feasible,status")
    worst = Inf; worst_t = -1
    for t in 1:T
        res = inner_opf(data, t, copy(pbr))
        @printf(exio, "%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%s,%s\n", system, T, t, net[t],
                sum(pbr), res.feasible ? res.ps : NaN, res.feasible ? res.loss : NaN,
                res.feasible, res.status)
        if res.feasible && res.ps < worst; worst = res.ps; worst_t = t; end
    end
    @printf("REVERSE_EXPORT closest approach to zero: P_Subs=%.6f pu at t=%d (total battery=%.4f pu)\n",
            worst, worst_t, sum(pbr))
end
println("SURVEY wrote=$expath")

# ---------------------------------------------------------------------------
# diagnostic restoration for any infeasible case (kept strictly separate)
# ---------------------------------------------------------------------------
if !isempty(infeasible_cases)
    rpath = joinpath(outdir, "infeasible_restoration_$(system)_T$(T).csv")
    open(rpath, "w") do rio
        println(rio, "system,horizon,time_index,pattern,seed,restored,status," *
                     "export_violation_pu,undervoltage_max_pu2,overvoltage_max_pu2,total_violation,ps_relaxed_pu")
        for (t, name, seed, pb) in infeasible_cases
            d = restore_feasibility(data, t, pb)
            if d.restored
                @printf(rio, "%s,%d,%d,%s,%d,%s,%s,%.6e,%.6e,%.6e,%.6e,%.6f\n",
                    system, T, t, name, seed, d.restored, d.status, d.export_violation,
                    d.undervoltage_max, d.overvoltage_max, d.total_violation, d.ps)
                @printf("RESTORE t=%d %-30s export=%.3e undervolt=%.3e overvolt=%.3e\n",
                    t, name, d.export_violation, d.undervoltage_max, d.overvoltage_max)
            else
                @printf(rio, "%s,%d,%d,%s,%d,false,%s,NaN,NaN,NaN,NaN,NaN\n",
                    system, T, t, name, seed, d.status)
            end
        end
    end
    println("SURVEY wrote=$rpath")
else
    println("SURVEY no infeasible inner solves; restoration diagnostic not required")
end
