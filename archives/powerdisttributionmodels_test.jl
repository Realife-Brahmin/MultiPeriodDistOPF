using Pkg
Pkg.activate(".")
using PowerModelsDistribution
using Ipopt

pkgdir = dirname(dirname(pathof(PowerModelsDistribution)));
# run(`explorer $(pkgdir)`);
casefolder = joinpath(pkgdir, "test", "data", "opendss");
casename = "case3_unbalanced.dss";
casefile = joinpath(casefolder, casename);
# solve_mc_opf(casefile, ACPUPowerModel, Ipopt.Optimizer);

eng = parse_file(casefile);
math = transform_data_model(eng);
result = solve_mc_opf(eng, ACPUPowerModel, Ipopt.Optimizer);

using SCS
using JuMP
model = JuMP.Model(SCS.Optimizer)
result_sdp = solve_mc_opf(eng, SDPUBFPowerModel, model)