# Pure saved-point evaluation: builds variable maps but never calls optimize!.
using Serialization, SHA, TOML, DelimitedFiles
source,c5,c6,out=ARGS
include(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))
cases=[joinpath(c5,"results",name) for name in
    ("large10k_T24_rho384_NO_31650191","large10k_T24_rho384_fixed_NO_31650218","large10k_T24_rho48_fixed_NO_31650238")]
append!(cases,[joinpath(c6,"results",name) for name in
    ("large10k_T24_rho48_relative_mu2_NO_31650253","large10k_T24_rho32_relative_mu2_NO_31650276")])
mkpath(out); results=Dict{String,Any}()
for case in cases
    meta=TOML.parsefile(joinpath(case,"calibration.toml"))
    ENV["TADMM_SCREEN"]=meta["screening"]
    d=deserialize(joinpath(source,"problems/mpopf/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls"))
    d[:C_B]=battery_cb("large10kC_1ph","system")
    state=deserialize(joinpath(case,"calibration_checkpoint.jls"))
    bytes2hex(sha256(read(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))))==state[:model_sha256] || error("Model mismatch")
    us=Vector{Vector{Float64}}(undef,d[:T]); ix,nu=control_layout(d)
    for t in 1:d[:T]
        s=make_model(d,[t];local_window=true,rho=state[:rho])
        vals=state[:stage_primals][t]
        length(vals)==length(s.vars) || error("Variable count mismatch")
        all(index(v).value==n for (n,v) in enumerate(s.vars)) || error("Noncontiguous variable map")
        at(v)=vals[index(v).value]
        u=zeros(nu); u[ix.ps]=at(s.ps[t]); u[ix.qs]=at(s.qs[t])
        for (n,e) in enumerate(d[:Lset]); u[ix.P[n]]=at(s.p[e,t]); u[ix.Q[n]]=at(s.q[e,t]); u[ix.ell[n]]=at(s.ell[e,t]); end
        for (n,j) in enumerate(d[:Nset]); u[ix.v[n]]=at(s.v[j,t]); end
        for (n,j) in enumerate(d[:Bset]); u[ix.pb[n]]=at(s.pb[j,t]); end
        for (n,j) in enumerate(d[:Dset])
            qm=sqrt(max(0.0,d[:S_D_R][j]^2-d[:p_D_pu][j,t]^2))
            u[ix.qnorm[n]]=qm>0 ? at(s.qd[j,t])/qm : 0.0
        end
        us[t]=u; s=nothing; GC.gc()
    end
    B=state[:z]; m=physical_metrics(d,us,B)
    trace,header=readdlm(joinpath(case,"calibration.csv"),',';header=true)
    position=Dict(String(k)=>n for (n,k) in enumerate(vec(header)))
    isapprox(m.objective,trace[end,position["objective"]];rtol=1e-13,atol=1e-9) || error("Saved objective mismatch")
    isapprox(m.original_primal,trace[end,position["original_primal"]];rtol=1e-12,atol=1e-10) || error("Saved feasibility mismatch")
    dynamic=0.0; accumulated_bound=0.0; recovered=[d[:B0_pu][j] for j in d[:Bset]]
    for t in 1:d[:T],(n,j) in enumerate(d[:Bset])
        prev=t==1 ? d[:B0_pu][j] : B[n,t-1]
        dynamic=max(dynamic,abs(B[n,t]-prev+d[:delta_t_h]*us[t][ix.pb[n]]))
        recovered[n]-=d[:delta_t_h]*us[t][ix.pb[n]]
        accumulated_bound=max(accumulated_bound,d[:soc_min][j]*d[:B_R_pu][j]-recovered[n],recovered[n]-d[:soc_max][j]*d[:B_R_pu][j])
    end
    results[basename(case)]=Dict("iteration"=>state[:iteration],"original_primal"=>m.original_primal,
        "max_equation"=>m.maxeq,"max_bound"=>m.maxbound,"max_cone"=>m.maxcone,
        "battery_dynamic_equation"=>dynamic,"integrated_power_soc_bound"=>accumulated_bound,
        "objective_trace_matches"=>true,"primal_trace_matches"=>true,"additional_optimization_solves"=>0)
end
open(io->TOML.print(io,results),joinpath(out,"components.toml"),"w")
println("SAVED_POINTS_INSPECTED ",out)
