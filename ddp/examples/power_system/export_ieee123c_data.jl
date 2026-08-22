# Export the existing tADMM/OpenDSS parser output so the FilterDDP environment
# does not need to duplicate the OpenDSS dependency.
#
# Run from the repository root:
#   julia --startup-file=no --project=envs/tadmm \
#         ddp/examples/power_system/export_ieee123c_data.jl [system] [T]

using Serialization

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO, "envs", "tadmm", "parse_opendss.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
dt = 24.0 / T
load_shape = T == 1 ? [0.97174] :
    0.8 .+ 0.2 .* (sin.(range(0, 2pi, length=T) .- 0.8) .+ 1) ./ 2
cost_shape = T == 1 ? [0.14] :
    0.08 .+ 0.12 .* (sin.(range(0, 2pi, length=T)) .+ 1) ./ 2
pv_shape = let
    pv = zeros(T)
    if T >= 4
        a = max(1, round(Int, 0.25T))
        b = min(T, round(Int, 0.75T))
        pv[a:b] .= sin.(range(0, pi, length=b-a+1))
    elseif T == 3
        pv[2] = 1.0
    elseif T == 2
        pv .= 0.5
    else
        pv .= 1.0
    end
    pv
end

network_kv_base = system == "ieee2522C_1ph" ? 7.2 :
    (system == "large10kC_1ph" ? 12.47 : 2.4018)

data = parse_system_from_dss(system, T;
    LoadShapeLoad=load_shape,
    LoadShapeCost=cost_shape,
    LoadShapePV=pv_shape,
    C_B=1e-6 * minimum(cost_shape),
    delta_t_h=dt,
    kV_B=network_kv_base)

outdir = joinpath(REPO, "ddp", "results", "network_filterddp")
mkpath(outdir)
outfile = joinpath(outdir, "network_data_$(system)_T$(T).jls")
serialize(outfile, data)

println("exported=$outfile")
println("T=$T buses=$(length(data[:Nset])) branches=$(length(data[:Lset])) " *
        "batteries=$(length(data[:Bset])) pv=$(length(data[:Dset])) loads=$(length(data[:NLset]))")
