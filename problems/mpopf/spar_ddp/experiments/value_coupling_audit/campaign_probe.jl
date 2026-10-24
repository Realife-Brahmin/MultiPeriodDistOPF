include(joinpath(@__DIR__,"direct_initial_probe.jl"))
include(joinpath(@__DIR__,"campaign_compact.jl"))
function campaign_run(cap,input,out)
    @assert !ispath(out);mkpath(out)
    meta=TOML.parsefile(joinpath(cap,"metadata.toml"));leaf=joinpath(cap,"leaves.jls")
    @assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
    saved=deserialize(leaf);raw=saved.elems;n=length(first(raw).ma)
    @assert n==1020&&length(raw)==24&&"late" in meta["tags"]&&!_paired_inplace()
    packets=Any[]
    for row in meta["stage_files"]
        f=joinpath(cap,row["file"]);@assert filesize(f)==row["bytes"]&&bytes2hex(sha256(read(f)))==row["sha256"]
        push!(packets,deserialize(f))
    end
    sort!(packets;by=p->p.stage);@assert [p.stage for p in packets]==[1,2,3,4]
    data=deserialize(input);@assert length(data[:Bset])==n
    setup=@timed topology_hierarchy(data;maxleaf=8);tree=setup.value
    full=production_values(raw,PintVal(zeros(n,n),zeros(n)))
    ord=ordinary(raw,PintVal(zeros(n,n),zeros(n)))
    oe=maximum(max(scaled_error(full[t].V,ord[t].V),scaled_error(full[t].v,ord[t].v)) for t in 1:24);@assert oe<1e-7
    tail=campaign_tail(raw);@assert scaled_error(tail.V,full[5].V)==0&&scaled_error(tail.v,full[5].v)==0
    rules=[policy(p,full[p.stage+1],saved.mu) for p in packets];probe=0.0
    for (t,p) in enumerate(packets),pr in p.probes,f in (:du,:dphi)
        probe=max(probe,scaled_error(getproperty(directions(rules[t],pr.dx),f),getproperty(pr,f)))
    end
    @assert probe<1e-7
    pathref=battery_path(raw[1:4],full,packets[1].f)
    dirs=open(joinpath(out,"directions.csv"),"w");emit(dirs,"mode,stage,probe,field,absolute,reference,scaled,relative_floor_1e_10")
    independent=(directions=[[cc_stage(p,full[p.stage+1].v,x->full[p.stage+1].V*x,saved.mu,pr.dx) for pr in p.probes] for p in packets],path=pathref)
    ig=campaign_compare(dirs,"independent_dense_operator",independent,packets,rules,pathref);@assert ig.pass
    dense_inputs=campaign_dense_inputs(packets,n);dense=production_values(dense_inputs,tail)
    dg=campaign_compare(dirs,"dense_window",campaign_dense_consume(packets,dense,dense_inputs,saved.mu),packets,rules,pathref);@assert dg.pass
    println("CAMPAIGN_REFERENCE_PASS saved=$probe independent=$(ig.scaled) dense=$(dg.scaled)");flush(stdout)
    genstats=Dict{String,Int}[];ps=Any[]
    for p in packets
        st=Dict{String,Int}();push!(ps,di_generate(p.S,tree;stats=st));push!(genstats,st)
    end
    # Match the concrete input collection type used by the timed helper.
    ps=campaign_direct_inputs(packets,tree)
    valueio=open(joinpath(out,"values.csv"),"w");emit(valueio,"phase,stage,field,absolute,reference,scaled")
    max_value=0.0;poison_pass=true
    function check(phase,t,field,A,B;count=true)
        local ae=norm(A-B);local rn=norm(B);local se=ae/max(1.0,rn)
        @assert isfinite(se);emit(valueio,phase,t,field,ae,rn,se)
        count&&(max_value=max(max_value,se));se
    end
    for t in 1:4
        poison=deepcopy(packets[t].S);fill!(poison.b0.Sh0,NaN);other=di_generate(poison,tree)
        for (j,f) in enumerate((:Maa,:Mab,:Mbb))
            A=expand_block(ps[t].reps[j],n);B=getproperty(raw[t],f)
            check("direct_leaf",t,f,A,B)
            check("raw_skew",t,f,A-A',B-B';count=false)
            poison_pass&=expand_block(other.reps[j],n)==A
        end
        for f in (:ma,:mb);check("direct_leaf",t,f,getproperty(ps[t],f),getproperty(raw[t],f));end
    end
    @assert poison_pass
    boundary=campaign_boundary(tail,tree);stats=Dict{String,Int}();values=campaign_build(ps,boundary;stats=stats)
    P=hcat(ones(n),sin.(collect(1:n)),cos.(collect(1:n)))
    # Full coordinate-basis reconstruction ONLY for validation, never execution.
    reconstructed=Matrix{Float64}[]
    for t in 1:4
        A=zeros(n,n)
        for j in 1:16:n
            ix=j:min(j+15,n);E=zeros(n,length(ix));for (k,i) in enumerate(ix);E[i,k]=1;end
            A[:,ix]=campaign_action(values[t],E)
        end
        push!(reconstructed,A)
        check("closed_value",t,:V,A,full[t].V);check("closed_value",t,:v,values[t].v,full[t].v)
        check("value_actions",t,:VtimesP,campaign_action(values[t],P),full[t].V*P)
    end
    gi=campaign_compare(dirs,"direct_all4",campaign_consume(packets,values,boundary,saved.mu),packets,rules,pathref)
    ng=campaign_compare(dirs,"negative_composed_gradient",campaign_consume(packets,values,boundary,saved.mu;corrupt=true),packets,rules,pathref)
    @assert !ng.pass
    inertio=open(joinpath(out,"inertia.csv"),"w");emit(inertio,"stage,min_eig,max_eig,threshold,negative,near_zero,reference_negative,reference_near_zero,same")
    inertia_ok=true
    for t in 1:4
        next=t==4 ? tail.V : reconstructed[t+1]
        ev=eigvals(Symmetric(raw[t].Mbb+next,:U));rv=eigvals(Symmetric(raw[t].Mbb+full[t+1].V,:U))
        threshold=n*eps(Float64)*maximum(abs,rv)
        neg=count(<(-threshold),ev);near=count(x->abs(x)<=threshold,ev)
        rneg=count(<(-threshold),rv);rnear=count(x->abs(x)<=threshold,rv)
        same=neg==rneg&&near==rnear;inertia_ok&=same
        emit(inertio,t,minimum(ev),maximum(ev),threshold,neg,near,rneg,rnear,same)
    end
    accurate=max_value<1e-7&&gi.pass&&inertia_ok
    println("CAMPAIGN_ACCURACY accurate=$accurate values=$max_value directions=$(gi.scaled) relative=$(gi.relative) inertia=$inertia_ok");flush(stdout)
    times=open(joinpath(out,"timings.csv"),"w");emit(times,"repetition,order,mode,phase,seconds,alloc_bytes,gc_s,compile_s,recompile_s")
    repeats=open(joinpath(out,"repeat_checks.csv"),"w");emit(repeats,"repetition,mode,scaled,relative,value_actions_scaled,pass")
    function record(rep,pos,mode,phase,z)
        emit(times,rep,pos,mode,phase,z.time,z.bytes,z.gctime,z.compile_time,z.recompile_time);z.value
    end
    if accurate
        # Exact helper/type paths warmed, including common incoming construction.
        wb=campaign_tail(raw);wp=campaign_direct_inputs(packets,tree);wt=campaign_boundary(wb,tree)
        wv=campaign_build(wp,wt);campaign_conditions(wv);campaign_value_actions(wv,P)
        campaign_consume(packets,wv,wt,saved.mu)
        wd=campaign_dense_inputs(packets,n);wv_d=production_values(wd,wb)
        campaign_dense_actions(wv_d,P);campaign_dense_consume(packets,wv_d,wd,saved.mu)
        for rep in 1:3,(pos,mode) in enumerate(isodd(rep) ? ("dense","direct_compact") : ("direct_compact","dense"))
            GC.gc()
            fresh_tail=record(rep,pos,mode,"incoming_dense_future_construction",@timed campaign_tail(raw))
            if mode=="dense"
                fresh=record(rep,pos,mode,"initial_dense_generation",@timed campaign_dense_inputs(packets,n))
                vv=record(rep,pos,mode,"dense_values",@timed production_values(fresh,fresh_tail))
                acts=record(rep,pos,mode,"all_value_actions",@timed campaign_dense_actions(vv,P))
                obs=record(rep,pos,mode,"all_stage_consumers",@timed campaign_dense_consume(packets,vv,fresh,saved.mu))
            else
                fresh=record(rep,pos,mode,"initial_direct_generation",@timed campaign_direct_inputs(packets,tree))
                bt=record(rep,pos,mode,"incoming_factor_discovery",@timed campaign_boundary(fresh_tail,tree))
                vv=record(rep,pos,mode,"compact_values",@timed campaign_build(fresh,bt))
                record(rep,pos,mode,"condition_probes",@timed campaign_conditions(vv))
                acts=record(rep,pos,mode,"all_value_actions",@timed campaign_value_actions(vv,P))
                obs=record(rep,pos,mode,"all_stage_consumers",@timed campaign_consume(packets,vv,bt,saved.mu))
            end
            gate=campaign_compare(dirs,"repeat$(rep)_$mode",obs,packets,rules,pathref)
            ve=maximum(scaled_error(acts[t],full[t].V*P) for t in 1:4)
            emit(repeats,rep,mode,gate.scaled,gate.relative,ve,gate.pass&&ve<1e-7);@assert gate.pass&&ve<1e-7
            record(rep,pos,mode,"cleanup_gc",@timed GC.gc())
        end
    end
    for io in (dirs,valueio,inertio,times,repeats);close(io);end
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,Dict("accurate"=>accurate,"nB"=>n,"N"=>24,"window"=>4,"generated_stages"=>[1,2,3,4],
          "generator_stats"=>genstats,"value_stats"=>stats,"value_max_scaled"=>max_value,"direction_scaled"=>gi.scaled,
          "direction_relative"=>gi.relative,"inertia_ok"=>inertia_ok,"negative_control_rejected"=>!ng.pass,
          "negative_scaled"=>ng.scaled,"poisoned_Sh0_independence_pass"=>poison_pass,
          "ordinary_error"=>oe,"saved_probe_error"=>probe,"independent_consumer_scaled"=>ig.scaled,
          "dense_window_scaled"=>dg.scaled,"initial_coefficients_bytes"=>sum(p.bytes for p in ps)+8*(boundary.rep.entries+n),
          "suffix_dimensions"=>[length(v.g) for v in values[2:4]],"suffix_coefficient_bytes"=>[8v.H.entries for v in values[2:4]],
          "object_graph_bytes_not_rss"=>Base.summarysize(values),"maxrss_all_paths_bytes"=>Sys.maxrss(),
          "topology_cold_s_not_benchmark"=>setup.time,"boundary_construction_charged_both_modes"=>true,
          "incoming_dense_tail_stages"=>20,"full_compact_horizon"=>false,"optimizer_run"=>false,
          "all_four_stage_policies_consumed"=>true,"algorithm_dense_boundary_exports"=>0,
          "true_future_given_free"=>false,"initial_factors_reused_across_iterates"=>false,
          "exact_timed_helpers_warmed"=>true,"compilation_recorded_every_phase"=>true,
          "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input)))))
    end
    println("CAMPAIGN_PROBE_COMPLETE accurate=$accurate full_horizon=false optimizer=false");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @assert length(ARGS)==3;campaign_run(ARGS...)
end
