# Bounded MA57-only replay from a frozen plateau, never a fresh benchmark.
using Serialization, SHA, TOML
source,case,out=ARGS[1:3]
variants=length(ARGS)>=4 ? String.(split(ARGS[4],',')) : ["default","centered","strict_comp","adaptive"]
steps=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 8
all(v->v in ("default","centered","strict_comp","adaptive","root_bounds_removed"),variants) || error("Unknown audit variant")
include(joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl"))
BLAS.set_num_threads(1)
ispath(out) && error("Refusing existing audit directory")
mkpath(out)
cp(joinpath(case,"warmup_checkpoint.jls"),joinpath(out,"checkpoint.jls"))
state=deserialize(joinpath(out,"checkpoint.jls"))
d=deserialize(joinpath(case,"instance.jls"))
model_path=joinpath(source,"problems/mpopf/tadmm/review_pilot/model.jl")
bytes2hex(sha256(read(model_path)))==state[:model_sha256] || error("Checkpoint/model mismatch")
state[:iteration]==5000 || error("Expected the preserved 5000-iteration plateau")
ENV["TADMM_LINEAR_SOLVER"]="ma57"
lib=joinpath(homedir(),"opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so")
ENV["TADMM_IPOPT_OPTIONS"]="hsllib=$lib;linear_system_scaling=none;ma57_pivot_order=2"
open(joinpath(out,"audit.toml"),"w") do io
    TOML.print(io,Dict("checkpoint_iteration"=>state[:iteration],"rho_fixed"=>state[:rho],
        "model_sha256"=>state[:model_sha256],"driver_sha256"=>state[:driver_sha256],
        "audit_source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "checkpoint_sha256"=>bytes2hex(sha256(read(joinpath(out,"checkpoint.jls")))),
        "instance_sha256"=>bytes2hex(sha256(read(joinpath(case,"instance.jls")))),
        "hsl_sha256"=>bytes2hex(sha256(read(lib))),"linear_solver"=>"ma57",
        "note"=>"Controlled saved-state diagnosis; fixed rho, no benchmark or convergence claim"))
end

function replay_variant(d,state,variant,out; steps=8)
    dir=joinpath(out,variant); mkpath(dir)
    T=d[:T]; rho=state[:rho]; z=copy(state[:z]); dual=deepcopy(state[:dual])
    models=[make_model(d,[t];local_window=true,rho=rho) for t in 1:T]
    for (t,s) in enumerate(models)
        length(s.vars)==length(state[:stage_primals][t]) || error("Variable order/count mismatch")
        set_start_value.(s.vars,state[:stage_primals][t])
        if variant=="root_bounds_removed"
            root=d[:substationBus]
            d[:Vminpu][root]<=1.05<=d[:Vmaxpu][root] || error("Root equality outside original voltage bounds")
            for time in s.owned
                # The equality v[root,time]=1.05^2 remains; these two bounds
                # are redundant and create no change in the feasible set.
                delete_lower_bound(s.v[root,time]); delete_upper_bound(s.v[root,time])
            end
        end
        # A diagnostic cap bounds cost; a failed local solve is recorded, rejected,
        # and never advanced into the consensus iteration.
        set_optimizer_attribute(s.m,"max_iter",300)
        set_optimizer_attribute(s.m,"file_print_level",5)
        variant in ("strict_comp","adaptive") && set_optimizer_attribute(s.m,"compl_inf_tol",1e-8)
        variant=="adaptive" && set_optimizer_attribute(s.m,"mu_strategy","adaptive")
    end
    bs=d[:Bset]; nb=length(bs); counts=[count(s->t in s.times,models) for t in 1:T]
    csv=open(joinpath(dir,"replay.csv"),"w")
    println(csv,"step,objective,original_primal,r_max,s_max,rho,max_dual_sum_interior,quadratic_coefficient_error,linear_coefficient_error")
    try
        for k in 1:steps
            copies=Any[]; quad_error=0.0; linear_error=0.0
            for (t,s) in enumerate(models)
                # Query the entire previous result before changing logging or
                # objectives: optimizer-attribute edits invalidate JuMP results.
                if k>1
                    previous=value.(s.vars)
                    set_start_value.(s.vars,previous)
                end
                if variant=="centered"
                    # Same mathematical objective including the previously omitted
                    # constant; this tests setter algebra and cancellation effects.
                    @objective(s.m,Min,s.cost+rho/2*sum(
                        (s.b[j,time]-z[n,time]+dual[t][n,a])^2
                        for (n,j) in enumerate(bs),(a,time) in enumerate(s.times)))
                end
                logfile=joinpath(dir,"step_$(lpad(k,2,'0'))_period_$(lpad(t,2,'0')).log")
                set_optimizer_attribute(s.m,"output_file",logfile)
                row=try
                    # Starts were restored above; avoid querying invalidated
                    # cached results after setting the per-step log destination.
                    stage_step!(s,d,z,dual[t],rho,1)
                catch err
                    open(joinpath(dir,"failure.toml"),"w") do io
                        TOML.print(io,Dict("step"=>k,"period"=>t,
                            "status"=>string(termination_status(s.m)),"error"=>sprint(showerror,err)))
                    end
                    println("AUDIT_LOCAL_REJECTED variant=$variant step=$k period=$t status=$(termination_status(s.m))")
                    return false
                end
                objective=objective_function(s.m)
                if k==1 && t==7
                    multipliers=[Dict("magnitude"=>abs(JuMP.dual(c)),"multiplier"=>JuMP.dual(c),"constraint"=>string(c))
                        for c in all_constraints(s.m;include_variable_in_set_constraints=true)]
                    sort!(multipliers;by=r->r["magnitude"],rev=true)
                    open(joinpath(dir,"period_07_top_duals.toml"),"w") do io
                        TOML.print(io,Dict("top"=>multipliers[1:min(10,length(multipliers))]))
                    end
                end
                for (n,j) in enumerate(bs),(a,time) in enumerate(s.times)
                    quad_error=max(quad_error,abs(coefficient(objective,s.b[j,time],s.b[j,time])-rho/2))
                    linear_error=max(linear_error,abs(coefficient(objective,s.b[j,time])-rho*(dual[t][n,a]-z[n,time])))
                end
                push!(copies,row.copy_b)
            end
            quad_error<=1e-10 && linear_error<=1e-10 || error("Augmented objective coefficient mismatch")
            zold=copy(z); accumulator=zeros(size(z))
            for (t,s) in enumerate(models),(a,time) in enumerate(s.times)
                accumulator[:,time].+=copies[t][:,a].+dual[t][:,a]
            end
            for time in 1:T,(n,j) in enumerate(bs)
                z[n,time]=clamp(accumulator[n,time]/counts[time],
                    d[:soc_min][j]*d[:B_R_pu][j],d[:soc_max][j]*d[:B_R_pu][j])
            end
            rmax=0.0; sumdual=zeros(size(z))
            for (t,s) in enumerate(models),(a,time) in enumerate(s.times)
                residual=copies[t][:,a].-z[:,time]
                rmax=max(rmax,maximum(abs.(residual);init=0.0))
                dual[t][:,a].+=residual
                sumdual[:,time].+=dual[t][:,a]
            end
            # Independent consensus stationarity identity, accounting for bound
            # projection: sum of scaled duals must vanish at interior coordinates.
            interior_error=0.0
            for time in 1:T,(n,j) in enumerate(bs)
                lo=d[:soc_min][j]*d[:B_R_pu][j]; hi=d[:soc_max][j]*d[:B_R_pu][j]
                if lo+1e-10 < z[n,time] < hi-1e-10
                    interior_error=max(interior_error,abs(sumdual[n,time]))
                end
            end
            interior_error<=1e-10 || error("Consensus/dual stationarity identity failed")
            us,bb=dispatch(models,d;z=z); metrics=physical_metrics(d,us,bb)
            smax=rho*maximum(abs.(z.-zold);init=0.0)
            println(csv,join((k,metrics.objective,metrics.original_primal,rmax,smax,rho,
                interior_error,quad_error,linear_error),',')); flush(csv)
            println("AUDIT_STEP variant=$variant step=$k original=$(metrics.original_primal) rmax=$rmax smax=$smax")
            flush(stdout)
        end
    finally
        close(csv)
    end
    return true
end

for variant in variants
    replay_variant(d,state,variant,out;steps=steps)
    GC.gc()
end
println("AUDIT_FINISHED $out")
