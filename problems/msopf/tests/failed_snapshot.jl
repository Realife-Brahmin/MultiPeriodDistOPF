# The fixture intentionally tests the legacy elastic-band margin mode. Isolate it
# from a caller campaign's optional combined voltage/import margin setting.
withenv("MULTI_VMARGIN" => nothing) do
    include(joinpath(@__DIR__,"customer_band.jl"))
end
ENV["MULTIPERIOD_FAILED_SNAPSHOT"]="0"
x=solve_case(s,ang,:cyclic;max_iter=0)
@test x.status==MOI.ITERATION_LIMIT
@test !x.ok && isempty(x.V_pu) && isempty(x.P_subs_kW)
ENV["MULTIPERIOD_FAILED_SNAPSHOT"]="1"
z=solve_case(s,ang,:cyclic;max_iter=0)
@test z.status==MOI.ITERATION_LIMIT && !z.ok
@test !isempty(z.V_pu) && !isempty(z.P_subs_kW)
y=solve_case(s,ang,:cyclic)
@test y.ok && !isempty(y.V_pu)
println("PASS: failed-snapshot summary excludes unusable primals; successful default behavior preserved")
