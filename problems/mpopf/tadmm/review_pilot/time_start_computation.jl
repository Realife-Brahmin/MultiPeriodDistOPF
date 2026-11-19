# How long does the power-flow start take? It is the per-data-set work that the
# solve clocks leave out (it sits inside the model build). One period at a time,
# after a warm-up call. usage: julia time_start_computation.jl system T
include("model.jl")
system=ARGS[1]; T=parse(Int,ARGS[2])
d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
apply_review_cb!(d,system)
loadflow_initial_values(d,1)
times=[(@elapsed loadflow_initial_values(d,t)) for t in 1:T]
@printf("START_COMPUTATION system=%s T=%d all_periods_serial_s=%.3f slowest_period_s=%.4f mean_period_s=%.4f\n",system,T,sum(times),maximum(times),sum(times)/T)
