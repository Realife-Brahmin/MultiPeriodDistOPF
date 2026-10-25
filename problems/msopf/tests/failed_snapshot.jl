include(joinpath(@__DIR__,"customer_band.jl"))
ENV["MULTIPERIOD_FAILED_SNAPSHOT"]="0"
x=solve_case(s,ang,:cyclic;max_iter=0)
@test x.status==MOI.ITERATION_LIMIT
@test !x.ok && isempty(x.V_pu) && isempty(x.P_subs_kW)
ENV["MULTIPERIOD_FAILED_SNAPSHOT"]="1"
y=solve_case(s,ang,:cyclic)
@test y.ok && !isempty(y.V_pu)
println("PASS: failed-snapshot summary excludes unusable primals; successful default behavior preserved")
