# Time one tADMM period solve with Gurobi under several parameter settings.
# usage: SC_STAGE_SOLVER=gurobi julia --project=conic_env gurobi_stage_bench.jl instance.jls trace.csv [periods]
# Each setting solves the same period at the penalty of a few iterations of a real run,
# with that run's battery targets replaced by zero (the time depends little on them).
using Serialization, Printf
include("storage_coordinator_actor.jl")
d=deserialize(ARGS[1])
lines=readlines(ARGS[2]); col=findfirst(==("rho"),split(lines[1],","))
trrho=[parse(Float64,split(l,",")[col]) for l in lines[2:end]]
periods=length(ARGS)>=3 ? parse.(Int,split(ARGS[3],",")) : [1,d[:T]÷2]
rhos=trrho[[1,min(10,end),min(22,end),end]]
nb=length(d[:Bset]); z=zeros(nb); u=zeros(nb)
settings=get(ENV,"BENCH_SET","1")=="3" ? [("tol1e-8",Dict())] : get(ENV,"BENCH_SET","1")=="2" ? [("tol1e-8",Dict()),
          ("tol1e-5",Dict("BarConvTol"=>1e-5,"BarQCPConvTol"=>1e-5)),
          ("tol1e-3",Dict("BarConvTol"=>1e-3,"BarQCPConvTol"=>1e-3)),
          ("scale0",Dict("ScaleFlag"=>0)),
          ("scale2",Dict("ScaleFlag"=>2)),
          ("homog1",Dict("BarHomogeneous"=>1)),
          ("tol1e-6_scale2",Dict("BarConvTol"=>1e-6,"BarQCPConvTol"=>1e-6,"ScaleFlag"=>2))] :
         [("tol1e-8",Dict()),
          ("tol1e-6",Dict("BarConvTol"=>1e-6,"BarQCPConvTol"=>1e-6)),
          ("tol1e-4",Dict("BarConvTol"=>1e-4,"BarQCPConvTol"=>1e-4)),
          ("presolve0",Dict("Presolve"=>0)),
          ("presolve1",Dict("Presolve"=>1)),
          ("order_amd",Dict("BarOrder"=>0)),
          ("order_nd",Dict("BarOrder"=>1)),
          ("tol1e-6_pre0",Dict("BarConvTol"=>1e-6,"BarQCPConvTol"=>1e-6,"Presolve"=>0)),
          ("tol1e-6_pre0_amd",Dict("BarConvTol"=>1e-6,"BarQCPConvTol"=>1e-6,"Presolve"=>0,"BarOrder"=>0))]
for t in periods
    s=make_power_stage_conic(d,t;rho=rhos[1])
    power_stage_step_conic!(s,d,z,u,rhos[1])   # compile
    base=Dict(k=>get_attribute(s.m,k) for k in ("BarConvTol","BarQCPConvTol","Presolve","BarOrder","ScaleFlag","BarHomogeneous"))
    for (name,opts) in settings
        for (k,v) in base; set_attribute(s.m,k,v); end
        for (k,v) in opts; set_attribute(s.m,k,v); end
        walls=Float64[]; solves=Float64[]; its=Int[]; objs=Float64[]
        for r in rhos, rep in 1:2
            tick=time_ns(); x,st,_=power_stage_step_conic!(s,d,z,u,r); w=(time_ns()-tick)/1e9
            rep==2 && (push!(walls,w); push!(solves,st); push!(its,MOI.get(s.m,MOI.BarrierIterations())); push!(objs,objective_value(s.m)))
        end
        @printf("BENCH period=%d setting=%-18s wall_mean=%.4f solve_mean=%.4f bar_its=%s status=%s obj=%s\n",
            t,name,sum(walls)/length(walls),sum(solves)/length(solves),join(its,"/"),termination_status(s.m),join((@sprintf("%.8f",o) for o in objs),"/"))
        flush(stdout)
    end
    # Ipopt on the same period, warm-started from its own previous solve, as tADMM runs it.
    si=make_power_stage(d,t;rho=rhos[1]); power_stage_step!(si,d,z,u,rhos[1],1)
    walls=Float64[]; its=Int[]
    for (k,r) in enumerate(rhos)
        tick=time_ns(); power_stage_step!(si,d,z,u,r,k+1); push!(walls,(time_ns()-tick)/1e9)
        push!(its,MOI.get(si.m,MOI.BarrierIterations()))
    end
    @printf("BENCH period=%d setting=%-18s wall_mean=%.4f its=%s\n",t,"ipopt_warm",sum(walls)/length(walls),join(its,"/"))
end
