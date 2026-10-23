include(joinpath(@__DIR__,"compact_chain_probe.jl"))
include(joinpath(@__DIR__,"direct_initial.jl"))
function di_run(cap,input,out)
    @assert !ispath(out);mkpath(out)
    meta=TOML.parsefile(joinpath(cap,"metadata.toml"));leaf=joinpath(cap,"leaves.jls")
    @assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
    for r in meta["stage_files"]
        file=joinpath(cap,r["file"]);@assert filesize(file)==r["bytes"]&&bytes2hex(sha256(read(file)))==r["sha256"]
    end
    packets=[deserialize(joinpath(cap,r["file"])) for r in meta["stage_files"]]
    available=sort([p.stage for p in packets]);@assert available==[1,12,24]
    packet=only(filter(p->p.stage==1,packets));saved=deserialize(leaf);raw=saved.elems;n=length(first(raw).ma)
    @assert n==1020&&length(raw)==24&&"late" in meta["tags"]&&!_paired_inplace()
    data=deserialize(input);tree=topology_hierarchy(data;maxleaf=8)
    full=production_values(raw,PintVal(zeros(n,n),zeros(n)));ov=ordinary(raw,PintVal(zeros(n,n),zeros(n)))
    ordinary_error=maximum(max(scaled_error(full[j].V,ov[j].V),scaled_error(full[j].v,ov[j].v)) for j in 1:24);@assert ordinary_error<1e-7
    rule=policy(packet,full[2],saved.mu);probe_error=0.0
    for pr in packet.probes,f in (:du,:dphi);probe_error=max(probe_error,scaled_error(getproperty(directions(rule,pr.dx),f),getproperty(pr,f)));end
    @assert probe_error<1e-7
    elems=raw[1:4];terminal=full[5];refpath=battery_path(elems,full,packet.f)
    dirs=open(joinpath(out,"directions.csv"),"w");emit(dirs,"mode,probe,field,absolute,reference,scaled,relative_floor_1e_10")
    reference=(directions=[cc_stage(packet,full[2].v,x->full[2].V*x,saved.mu,pr.dx) for pr in packet.probes],path=refpath)
    consumer_gate=cc_compare_dirs(dirs,"independent_dense_consumer",reference,packet,rule,refpath);@assert consumer_gate.pass
    println("DIRECT_REFERENCE_PASS stages=$(join(available,':')) consumer=$(consumer_gate.scaled)");flush(stdout)
    stats=Dict{String,Int}();direct=@timed di_generate(packet.S,tree;stats=stats);p1=direct.value
    # Poison the pre-existing dense Sh0: direct generation must NOT read it.
    poisoned=deepcopy(packet.S);fill!(poisoned.b0.Sh0,NaN)
    independent=di_generate(poisoned,tree)
    errio=open(joinpath(out,"values.csv"),"w");emit(errio,"phase,field,absolute,reference,scaled")
    maxerr=0.0
    function check(phase,f,A,B)
        local error_value=di_record_value(errio,phase,f,A,B)
        phase!="raw_skew"&&(maxerr=max(maxerr,error_value));error_value
    end
    for (j,f) in enumerate((:Maa,:Mab,:Mbb))
        A=expand_block(p1.reps[j],n);B=getproperty(raw[1],f)
        @assert check("direct_leaf1",f,A,B)<1e-7
        @assert scaled_error(expand_block(independent.reps[j],n),A)==0
        check("raw_skew",f,A-A',B-B')
    end
    for f in (:ma,:mb);@assert check("direct_leaf1",f,getproperty(p1,f),getproperty(raw[1],f))<1e-7;end
    # Leaf2:4 factors are NOT inferred from leaf1 or a different iterate. Their
    # original dense-discovery path and true future discovery remain charged.
    initial=cc_inputs(elems,terminal,tree);buildstats=Dict{String,Int}();c=cc_build(initial;stats=buildstats)
    d=cc_dense(elems,terminal)
    for (phase,a,b) in (("merge34",c.first,d.first),("merge234",c.chain,d.chain))
        e=ck_export(a) # validation only, never supplied to execution
        for f in (:Maa,:Mab,:Mbb,:ma,:mb);@assert check(phase,f,getproperty(e,f),getproperty(b,f))<1e-7;end
    end
    P=hcat(ones(n),sin.(collect(1:n)),cos.(collect(1:n)))
    @assert check("closed2",:actions,cc_value(c.value,P),full[2].V*P)<1e-7
    @assert check("closed2",:v,c.value.v,full[2].v)<1e-7
    prefixstats=Dict{String,Int}();prefix=di_prefix(p1,c,P;stats=prefixstats)
    @assert check("direct_closed1",:actions,prefix.Y,full[1].V*P)<1e-7
    @assert check("direct_closed1",:v,prefix.v,full[1].v)<1e-7
    runstats=Dict{String,Int}();z=di_consume(c,p1,packet,saved.mu;stats=runstats)
    gate=cc_compare_dirs(dirs,"direct_leaf1_pipeline",z,packet,rule,refpath);@assert gate.pass
    negative=cc_compare_dirs(dirs,"negative_composed_gradient",di_consume(c,p1,packet,saved.mu;corrupt=true),packet,rule,refpath);@assert !negative.pass
    println("DIRECT_ACCURACY pass=true values=$maxerr scaled=$(gate.scaled) relative=$(gate.relative) generator=$stats chain=$buildstats");flush(stdout)
    timings=open(joinpath(out,"timings.csv"),"w");emit(timings,"repetition,order,mode,phase,seconds,alloc_bytes,gc_s")
    checks=open(joinpath(out,"repeat_checks.csv"),"w");emit(checks,"repetition,mode,direction_scaled,direction_relative,prefix_scaled,pass")
    function record(rep,pos,mode,phase,a)
        emit(timings,rep,pos,mode,phase,a.time,a.bytes,a.gctime);a.value
    end
    # Warm every EXACT path; structural hierarchy is same fixed network.
    _pp_elem(packet.S,1,n);di_generate(packet.S,tree)
    di_dense_prefix(raw[1],d,P);di_prefix(p1,c,P)
    cc_dense_consume(d,packet,saved.mu,elems,terminal);di_consume(c,p1,packet,saved.mu)
    for rep in 1:3,(pos,mode) in enumerate(isodd(rep) ? ("dense","partial_direct") : ("partial_direct","dense"))
        GC.gc()
        if mode=="dense"
            e1=record(rep,pos,mode,"dense_leaf1_generation",@timed _pp_elem(packet.S,1,n))
            v=record(rep,pos,mode,"three_merges",@timed cc_dense(elems,terminal))
            pre=record(rep,pos,mode,"closed1_actions",@timed di_dense_prefix(e1,v,P))
            observed=record(rep,pos,mode,"directions",@timed cc_dense_consume(v,packet,saved.mu,elems,terminal))
        else
            e1=record(rep,pos,mode,"direct_leaf1_generation",@timed di_generate(packet.S,tree))
            ps=record(rep,pos,mode,"unavailable_leaf2to4_discovery",@timed [construct_blocks(elems[j],tree,1e-9) for j in 2:4])
            future=record(rep,pos,mode,"true_future_discovery",@timed (rep=block_payload(terminal.V,tree,1e-9,Any[],:terminal),v=copy(terminal.v)))
            v=record(rep,pos,mode,"compact_updates",@timed cc_build((p=ps,terminal=future)))
            record(rep,pos,mode,"condition_probe",@timed cc_condition(v.value))
            pre=record(rep,pos,mode,"closed1_actions",@timed di_prefix(e1,v,P))
            observed=record(rep,pos,mode,"directions",@timed di_consume(v,e1,packet,saved.mu))
        end
        g=cc_compare_dirs(dirs,"repeat$(rep)_$mode",observed,packet,rule,refpath)
        pe=max(scaled_error(pre.Y,full[1].V*P),scaled_error(pre.v,full[1].v))
        emit(checks,rep,mode,g.scaled,g.relative,pe,g.pass&&pe<1e-7);@assert g.pass&&pe<1e-7
        record(rep,pos,mode,"cleanup_gc",@timed GC.gc())
    end
    for io in (dirs,errio,timings,checks);close(io);end
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,Dict("accurate"=>true,"nB"=>n,"N"=>24,"window"=>4,"available_stage_records"=>available,
          "direct_generated_stages"=>[1],"unavailable_stage_records"=>[2,3,4],"true_future_direct"=>false,
          "poisoned_Sh0_independence_pass"=>true,"ordinary_error"=>ordinary_error,"saved_probe_error"=>probe_error,
          "consumer_reference_scaled"=>consumer_gate.scaled,"direction_scaled"=>gate.scaled,"direction_relative"=>gate.relative,
          "negative_control_rejected"=>!negative.pass,"negative_scaled"=>negative.scaled,"value_max_scaled"=>maxerr,
          "generator_stats"=>stats,"chain_stats"=>buildstats,"consumer_stats"=>runstats,"prefix_stats"=>prefixstats,
          "direct_leaf1_coefficient_bytes"=>p1.bytes,"direct_leaf1_object_graph_bytes_not_rss"=>Base.summarysize(p1),
          "direct_cold_s_not_benchmark"=>direct.time,"retained_interior_dimensions"=>[n,2n,3n],
          "chain_coefficient_bytes"=>8 .* [c.first.K.entries,c.chain.H.entries,c.value.H.entries],
          "maxrss_all_paths_bytes"=>Sys.maxrss(),"initial_full_boundary_matrices_constructed_by_direct_path"=>0,
          "initial_dense_SVD_for_direct_stage"=>false,"algorithm_dense_exports"=>0,
          "stage_network_factorization_is_given_input"=>true,"optimizer_run"=>false,"medium_extended"=>false,
          "new_capture"=>false,"reuse_across_iterates"=>false,"full_direct_pipeline_tested"=>false,
          "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input)))))
    end
    println("DIRECT_INITIAL_COMPLETE accurate=true full_direct_pipeline=false optimizer=false");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @assert length(ARGS)==3;di_run(ARGS...)
end
