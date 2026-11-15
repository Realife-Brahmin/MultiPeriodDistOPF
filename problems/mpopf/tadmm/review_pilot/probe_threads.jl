# Can the period OPFs be solved on several threads of ONE Julia process?
# Solves every period's model serially, then again with Threads.@threads, and
# compares the results. Ipopt with MA57 from an external HSL library is not
# documented as thread-safe, so identical answers here are necessary, not
# sufficient. Run with: julia --threads=N probe_threads.jl system T
include("storage_coordinator.jl")
system=ARGS[1]; T=parse(Int,ARGS[2])
d=deserialize(joinpath(SHARED,"..","..","results","network_filterddp","network_data_$(system)_T$(T)_periodic.jls"))
apply_review_cb!(d,system); tadmm_screen(d); physical_layout(d)
nb=length(d[:Bset]); zero_col=zeros(nb)
println("PROBE threads=",Threads.nthreads()," system=",system," T=",T)
build0=time_ns(); stages=[make_power_stage(d,t;rho=1000.0) for t in 1:T]; build_s=(time_ns()-build0)/1e9
serial=Vector{Vector{Float64}}(undef,T)
tick=time_ns()
for t in 1:T; serial[t],_,_=power_stage_step!(stages[t],d,zero_col,zero_col,1000.0,1); end
serial_s=(time_ns()-tick)/1e9
for repeat in 1:3
    threaded=Vector{Vector{Float64}}(undef,T)
    tick=time_ns()
    Threads.@threads :static for t in 1:T
        # k=1: cold start from the same initial point as the serial solve would need
        # fresh models; instead re-solve from the previous solution, as an iteration does.
        threaded[t],_,_=power_stage_step!(stages[t],d,zero_col,zero_col,1000.0,2)
    end
    threaded_s=(time_ns()-tick)/1e9
    worst=maximum(maximum(abs,threaded[t].-serial[t]) for t in 1:T)
    @printf("PROBE repeat=%d serial_s=%.2f threaded_s=%.2f speedup=%.2f max_difference=%.3e\n",repeat,serial_s,threaded_s,serial_s/threaded_s,worst)
end
@printf("PROBE_DONE build_s=%.2f maxrss_gb=%.2f\n",build_s,Sys.maxrss()/2^30)
