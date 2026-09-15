# Probe the reduced stage-value function of the proposed decomposition:
#
#     Phi_t(P_B^t) = min_y { c^t * S_base * dt * P_Subs^t : network constraints }
#
# i.e. the optimal value of the inner single-period network OPF as a function of
# the fixed battery dispatch. This asks whether Phi_t is smooth, locally convex
# and approximately quadratic -- the properties a battery-only outer DDP needs --
# and whether its gradient is recoverable analytically from the inner solve
# rather than by finite differences.
#
# The analytic candidate is the dual of the real-power balance row at each
# battery bus: P_B enters that row with coefficient -1, so -lambda_bal is the
# derivative of the inner optimal value with respect to that battery's dispatch.
# Every finite-difference derivative below is checked against it.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/probe_reduced_value_function.jl [system] [T]

using Printf
using Serialization
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 24

data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
                            "network_data_$(system)_T$(T).jls"))
ref = deserialize(joinpath(REPO, "envs", "tadmm", "processedData",
                           "$(system)_T$(T)", "sol_socp_bf.jls"))
Bset = data[:Bset]; nB = length(Bset)
pbr = Float64[data[:P_B_R_pu][j] for j in Bset]

net = [sum(data[:p_L_pu][j, t] for j in data[:NLset]) -
       sum(data[:p_D_pu][j, t] for j in data[:Dset]) for t in 1:T]
ord = sortperm(net)
samples = [(ord[1], "low"), (ord[cld(length(ord), 2)], "medium"), (ord[end], "high")]

root = data[:substationBus]
depth = Dict{Int,Int}(root => 0); queue = [root]
while !isempty(queue)
    j = popfirst!(queue)
    for k in get(data[:children], j, Int[]); depth[k] = depth[j] + 1; push!(queue, k); end
end
bdepth = [depth[j] for j in Bset]

# Exclude any battery sitting on the substation bus: the root balance row
# carries no pb term, so such a battery is inert in this transcription and any
# probe direction along it is identically flat. ieee2522C_1ph has exactly one.
nonroot_set = Set(data[:Nm1set])
effective = findall(b -> Bset[b] in nonroot_set, 1:nB)
inert = setdiff(1:nB, effective)
isempty(inert) || @printf("NOTE: %d battery(ies) on the substation bus are inert in this transcription and excluded from probe directions: bus %s\n",
                          length(inert), [Bset[b] for b in inert])

med_depth = median(bdepth[effective])
deep    = [b for b in effective if bdepth[b] >= med_depth]
shallow = [b for b in effective if bdepth[b] <  med_depth]
i_shallowest = effective[argmin(bdepth[effective])]
i_deepest    = effective[argmax(bdepth[effective])]
i_median     = effective[sortperm(bdepth[effective])[cld(length(effective), 2)]]

# unit directions in battery-power space
function directions()
    ds = Tuple{String,Vector{Float64}}[]
    for (nm, i) in (("single_shallowest", i_shallowest), ("single_median_depth", i_median),
                    ("single_deepest", i_deepest))
        d = zeros(nB); d[i] = 1.0; push!(ds, (nm, d))
    end
    d = zeros(nB); d[effective] .= 1.0 / sqrt(length(effective))
    push!(ds, ("aggregate_uniform", d))
    d = zeros(nB); d[deep] .= 1.0 / sqrt(length(deep)); push!(ds, ("group_deep", d))
    d = zeros(nB); d[shallow] .= 1.0 / sqrt(length(shallow)); push!(ds, ("group_shallow", d))
    return ds
end

inbox(p) = all(-pbr .- 1e-12 .<= p .<= pbr .+ 1e-12)

outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
csvpath = joinpath(outdir, "phi_probe_$(system)_T$(T).csv")
io = open(csvpath, "w")
println(io, "system,horizon,time_index,demand_label,base,direction,h_pu," *
            "phi_base_usd,phi_plus_usd,phi_minus_usd,first_derivative,second_derivative," *
            "analytic_directional_derivative,fd_vs_analytic_rel_err," *
            "active_qnorm_base,active_qnorm_plus,active_qnorm_minus,active_changed," *
            "in_box,status_plus,status_minus")

hs = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]

for (t, label) in samples
    bases = [("zero", zeros(nB)),
             ("centralized_optimal", Float64[ref[:P_B][j, t] for j in Bset])]
    for (bname, p0) in bases
        r0 = inner_opf(data, t, p0)
        r0.feasible || (@printf("t=%d base=%s INFEASIBLE, skipped\n", t, bname); continue)
        phi0 = r0.substation_cost
        # Analytic gradient of the inner optimal value wrt each battery dispatch.
        # Verified against central differences below: with JuMP's dual sign
        # convention for this equality row, d(Phi)/d(P_B) is +lambda_bal.
        grad_analytic = r0.lambda_bal
        for (dname, d) in directions()
            dderiv_analytic = sum(grad_analytic .* d)
            for h in hs
                pp = p0 .+ h .* d; pm = p0 .- h .* d
                ok = inbox(pp) && inbox(pm)
                if !ok
                    @printf(io, "%s,%d,%d,%s,%s,%s,%.1e,%.8f,NaN,NaN,NaN,NaN,%.8f,NaN,%d,-1,-1,false,false,BOX,BOX\n",
                        system, T, t, label, bname, dname, h, phi0, dderiv_analytic,
                        r0.n_active_qnorm)
                    continue
                end
                rp = inner_opf(data, t, pp); rm = inner_opf(data, t, pm)
                if !(rp.feasible && rm.feasible)
                    @printf(io, "%s,%d,%d,%s,%s,%s,%.1e,%.8f,NaN,NaN,NaN,NaN,%.8f,NaN,%d,-1,-1,false,true,%s,%s\n",
                        system, T, t, label, bname, dname, h, phi0, dderiv_analytic,
                        r0.n_active_qnorm, rp.status, rm.status)
                    continue
                end
                d1 = (rp.substation_cost - rm.substation_cost) / (2h)
                d2 = (rp.substation_cost - 2phi0 + rm.substation_cost) / h^2
                relerr = abs(dderiv_analytic) > 1e-12 ?
                         abs(d1 - dderiv_analytic) / abs(dderiv_analytic) : NaN
                changed = (rp.n_active_qnorm != r0.n_active_qnorm) ||
                          (rm.n_active_qnorm != r0.n_active_qnorm)
                @printf(io, "%s,%d,%d,%s,%s,%s,%.1e,%.8f,%.8f,%.8f,%.8f,%.6e,%.8f,%.3e,%d,%d,%d,%s,true,%s,%s\n",
                    system, T, t, label, bname, dname, h, phi0, rp.substation_cost,
                    rm.substation_cost, d1, d2, dderiv_analytic, relerr,
                    r0.n_active_qnorm, rp.n_active_qnorm, rm.n_active_qnorm, changed,
                    rp.status, rm.status)
            end
        end
        @printf("  t=%2d (%-6s) base=%-20s phi=%.6f  done\n", t, label, bname, phi0)
        flush(stdout)
    end
end
close(io)
println("PHI_PROBE wrote=$csvpath")
