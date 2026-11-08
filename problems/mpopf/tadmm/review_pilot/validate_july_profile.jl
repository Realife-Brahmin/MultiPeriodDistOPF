# Exercise the new objective profile and the explicit provisional-cutoff path.
using Serialization, TOML, SHA
source,reference,out=ARGS
include(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))
mkpath(out)
@assert isapprox(july_battery_cb(24),8.013987384856764e-8;rtol=1e-14)
@assert july_battery_cb(1)==1.4e-7
d=deserialize(joinpath(reference,"source/problems/mpopf/spar_ddp/results/network_filterddp/network_data_ieee2522C_1ph_T24_periodic.jls"))
profiles=deepcopy((d[:LoadShapeLoad],d[:LoadShapePV],d[:LoadShapeCost]))
apply_review_cb!(d,"ieee2522C_1ph";profile="july")
@assert profiles==(d[:LoadShapeLoad],d[:LoadShapePV],d[:LoadShapeCost])
@assert d[:C_B]==july_battery_cb(24)
ENV["TADMM_LOCAL_SOLVE_POLICY"]="cap30"; ENV["TADMM_SCREEN"]="substation,ell"
s=make_model(d,[19];local_window=true,rho=3800.0)
@assert s.local_policy=="cap30"
# Force a cutoff cheaply; its primal is provisional, never a feasible-result claim.
set_optimizer_attribute(s.m,"max_iter",1)
z=[d[:B0_pu][j] for j in d[:Bset],t in 1:d[:T]]
row=stage_step!(s,d,z,zeros(length(d[:Bset]),length(s.times)),3800.0,1)
@assert row.iteration_limit_hit && row.status=="ITERATION_LIMIT" && isfinite(row.local_primal)
central=make_model(d,[19])
@assert central.local_policy=="exact"
set_optimizer_attribute(central.m,"max_iter",1)
rejected=false
try
    solve_checked!(central)
catch exception
    global rejected=occursin("ITERATION_LIMIT",sprint(showerror,exception))
end
@assert rejected
open(io->TOML.print(io,Dict("july_cb_T24"=>d[:C_B],"profiles_unchanged"=>true,
    "cutoff_preserves_status"=>true,"cutoff_local_primal"=>row.local_primal,
    "central_cutoff_rejected"=>true)),joinpath(out,"validation.toml"),"w")
println("JULY_PROFILE_VALIDATED ",out)
