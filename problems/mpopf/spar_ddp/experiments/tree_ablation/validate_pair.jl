using Serialization, LinearAlgebra, Printf
length(ARGS)==4 || error("usage: validate_pair.jl tree.jls umfpack.jls tree.log umfpack.log")
a,b=deserialize.(ARGS[1:2])
for key in (:system,:T);a[key]==b[key] || error("mismatched $key");end
for sol in (a,b)
    sol[:status]==0 || error("solver did not converge: $(sol[:status])")
    all(v->all(isfinite,v),vcat(sol[:x],sol[:u])) || error("nonfinite solution")
end
for (log,system) in zip(ARGS[3:4],(a[:system],b[:system]))
    text=read(log,String)
    occursin("PROCS_SWEEP installed",text) || error("not a process-parallel solve")
    occursin("FILTERDDP_PROCS_FAIL",text) && error("process pass failed")
    m=match(r"final residuals: primal=([0-9.eE+\-]+) dual=([0-9.eE+\-]+) complementarity=([0-9.eE+\-]+)",text)
    isnothing(m) && error("missing final residuals")
    tol=system=="ieee123C_1ph" ? 1e-6 : 1e-5
    maximum(parse.(Float64,m.captures))<=tol || error("final residual exceeds tolerance")
end
# A matching objective is checked independently of the trajectory comparison.
objectives=Float64[]
for log in ARGS[3:4]
    m=match(r"FilterDDP objective=([0-9.eE+\-]+)",read(log,String))
    isnothing(m) && error("missing final objective")
    push!(objectives,parse(Float64,m[1]))
end
abs(objectives[1]-objectives[2])/max(abs(objectives[1]),abs(objectives[2]),1.0)<=1e-8 || error("objective disagreement")
rel(v,w)=norm(v-w)/max(norm(v),norm(w),1.0)
eu=maximum(rel.(a[:u],b[:u]));ex=maximum(rel.(a[:x],b[:x]))
max(eu,ex)<=1e-5 || error("solutions disagree: x=$ex u=$eu")
@printf("TREE_ABLATION_PAIR_PASS system=%s T=%d max_relative_x=%.3e max_relative_u=%.3e iterations_tree=%d iterations_umfpack=%d\n",a[:system],a[:T],ex,eu,a[:iterations],b[:iterations])
