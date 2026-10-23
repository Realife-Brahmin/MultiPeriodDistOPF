include(joinpath(@__DIR__,"compact_probe.jl"))
include(joinpath(@__DIR__,"compact_chain.jl"))

function cc_dense(elems,terminal)
    n=length(first(elems).ma);w=PairWork(n)
    a=_pp_combine!(deepcopy(elems[3]),elems[4],w)
    b=_pp_combine!(deepcopy(elems[2]),a,w)
    v=_pp_combine!(deepcopy(b),terminal,w)
    (first=a,chain=b,value=v)
end
function cc_inputs(elems,terminal,tree)
    p=[construct_blocks(elems[j],tree,1e-9) for j in 2:4]
    reports=Any[];rep=block_payload(terminal.V,tree,1e-9,reports,:terminal)
    (p=p,terminal=(rep=rep,v=copy(terminal.v)),terminal_reports=reports)
end
function cc_build(inputs;stats=Dict{String,Int}())
    p=inputs.p;a=ck_merge(p[2],p[3];stats=stats)
    b=cc_prepend(p[1],a;stats=stats)
    v=cc_close(b,inputs.terminal;stats=stats)
    (first=a,chain=b,value=v)
end
function cc_condition(c)
    n=length(c.H.tree.ids);x=fill(1/n,n);invnorm=0.0
    for _ in 1:4
        y=ck_solve(c.factor,x);invnorm=max(invnorm,norm(y,1))
        z=ck_solve(c.factor,sign.(y));x=zeros(n);x[argmax(abs.(z))]=1
    end
    norm(ck_apply(c.H,ones(n)),Inf)*invnorm # diagnostic estimate, NOT certificate
end
function cc_consume(c,packet,mu,elems;stats=Dict{String,Int}(),corrupt=false)
    n=length(c.value.v);future=c.value
    v=copy(future.v);corrupt && (v .+= 1e6sin.(collect(1:n)))
    action=x->cc_value(future,x)
    ds=[cc_stage(packet,v,action,mu,pr.dx;stats=stats) for pr in packet.probes]
    # Battery trajectory: first state fixed at zero, then solve x1 and recover
    # x2:x4 from the retained interior. No dense boundary matrix is exported.
    e=elems[1];k=x->MatrixUpperAction(e.Mbb,x)+action(x)
    sol=cc_cg(k,-(e.mb+v),diag(e.Mbb))
    stats["battery_cg_iterations"]=sol.iterations
    x1=sol.x
    rhs=copy(future.gradient);rhs[1:n] .+= ck_apply(future.chain.left.reps[2],x1;trans=true)
    rest=-ck_solve(future.factor,rhs)
    states=[zeros(n),x1,rest[1:n],rest[n+1:2n],rest[2n+1:3n]]
    path=[(dx=states[t],dp=(states[t+1]-states[t])./packet.f) for t in 1:4]
    (directions=ds,path=path)
end
MatrixUpperAction(A,x)=Symmetric(A,:U)*x
function cc_compare_dirs(io,label,a,packet,rule,referencepath)
    worst=0.0;relative=0.0
    function record(probe,field,A,B)
        d=norm(A-B);r=norm(B);s=d/max(1.0,r);q=d/max(1e-10,r)
        @assert isfinite(s)&&isfinite(q)
        emit(io,label,probe,field,d,r,s,q);worst=max(worst,s);relative=max(relative,q)
    end
    for (j,pr) in enumerate(packet.probes)
        d=directions(rule,pr.dx)
        for f in (:du,:dphi,:dzl,:dzu);record("stage1_probe$j",f,getproperty(a.directions[j],f),getproperty(d,f));end
    end
    for t in 1:4,f in (:dx,:dp);record("battery$t",f,getproperty(a.path[t],f),getproperty(referencepath[t],f));end
    (scaled=worst,relative=relative,pass=worst<DIR_SCALED&&relative<DIR_RELATIVE)
end
function cc_dense_consume(d,packet,mu,elems,terminal)
    rule=policy(packet,d.value,mu)
    ds=[directions(rule,pr.dx) for pr in packet.probes]
    w=PairWork(length(d.value.v))
    v3=_pp_combine!(deepcopy(d.first),terminal,w)
    v4=_pp_combine!(deepcopy(elems[4]),terminal,w)
    vals=Any[d.value,d.value,v3,v4,terminal] # suffix1 is unused by battery_path
    (directions=ds,path=battery_path(elems,vals,packet.f))
end
function cc_run(cap,input,out)
    @assert !ispath(out);mkpath(out)
    meta=TOML.parsefile(joinpath(cap,"metadata.toml"));leaf=joinpath(cap,"leaves.jls")
    @assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
    for r in meta["stage_files"]
        path=joinpath(cap,r["file"]);@assert filesize(path)==r["bytes"]&&bytes2hex(sha256(read(path)))==r["sha256"]
    end
    saved=deserialize(leaf);raw=saved.elems;n=length(first(raw).ma)
    @assert n==1020&&length(raw)==24&&"late" in meta["tags"]
    data=deserialize(input);@assert length(data[:Bset])==n
    full=production_values(raw,PintVal(zeros(n,n),zeros(n)))
    ov=ordinary(raw,PintVal(zeros(n,n),zeros(n)))
    ordinary_error=maximum(max(scaled_error(full[j].V,ov[j].V),scaled_error(full[j].v,ov[j].v)) for j in 1:24)
    @assert ordinary_error<1e-7
    packet=only([deserialize(joinpath(cap,r["file"])) for r in meta["stage_files"] if deserialize(joinpath(cap,r["file"])).stage==1])
    rule=policy(packet,full[2],saved.mu);probe_error=0.0
    for pr in packet.probes,f in (:du,:dphi)
        probe_error=max(probe_error,scaled_error(getproperty(directions(rule,pr.dx),f),getproperty(pr,f)))
    end
    @assert probe_error<1e-7
    elems=raw[1:4];terminal=full[5];referencepath=battery_path(elems,full,packet.f)
    dirs=open(joinpath(out,"directions.csv"),"w");emit(dirs,"mode,probe,field,absolute,reference,scaled,relative_floor_1e_10")
    # Independent consumer check with the unchanged production future, before
    # introducing compression. A failure here is a consumer/reference issue.
    consumerstats=Dict{String,Int}()
    exactds=[cc_stage(packet,full[2].v,x->full[2].V*x,saved.mu,pr.dx;stats=consumerstats) for pr in packet.probes]
    consumer_gate=cc_compare_dirs(dirs,"dense_operator_consumer",(directions=exactds,path=referencepath),packet,rule,referencepath)
    println("CHAIN_CONSUMER_REFERENCE pass=$(consumer_gate.pass) scaled=$(consumer_gate.scaled) relative=$(consumer_gate.relative)");flush(stdout)
    @assert consumer_gate.pass "independent operator direction consumer gate"
    topology=@timed topology_hierarchy(data;maxleaf=8);tree=topology.value
    inputs=@timed cc_inputs(elems,terminal,tree)
    stats=Dict{String,Int}()
    result=try
        cc_build(inputs.value;stats=stats)
    catch ex
        # Representation rejection is a scientific result, not a timing win.
        open(joinpath(out,"rejection.toml"),"w") do io
            TOML.print(io,Dict("reason"=>sprint(showerror,ex),"stats"=>stats,"optimizer"=>false))
        end
        rethrow()
    end
    open(joinpath(out,"ranks.csv"),"w") do io
        emit(io,"representation,depth,node_rows,left_rows,right_rows,lr_rank,rl_rank,lr_dense,rl_dense,entries")
        function walk(a,label,depth=0)
            if isnothing(a.localblock)
                emit(io,label,depth,length(a.tree.ids),length(a.left.tree.ids),length(a.right.tree.ids),size(ck_factors(a.lr)[1],2),size(ck_factors(a.rl)[1],2),a.lr.kind==:dense,a.rl.kind==:dense,a.entries)
                walk(a.left,label,depth+1);walk(a.right,label,depth+1)
            else
                emit(io,label,depth,length(a.tree.ids),0,0,0,0,false,false,a.entries)
            end
        end
        for (i,p) in enumerate(inputs.value.p),(j,a) in enumerate(p.reps);walk(a,"leaf$(i+1)_block$j");end
        walk(inputs.value.terminal.rep,"true_future")
        walk(result.first.K,"merge34");walk(result.chain.H,"merge234");walk(result.value.H,"closed234")
    end
    P=hcat(ones(n),sin.(collect(1:n)),cos.(collect(1:n)))
    ref=cc_dense(elems,terminal)
    errio=open(joinpath(out,"values.csv"),"w");emit(errio,"depth,field,absolute,reference,scaled")
    maxerr=0.0
    # Final validation-only reconstruction is outside EVERY execution timer.
    for (depth,c,d) in ((1,result.first,ref.first),(2,result.chain,ref.chain))
        exported=ck_export(c)
        for f in (:Maa,:Mab,:Mbb,:ma,:mb)
            A=getproperty(exported,f);B=getproperty(d,f);err=scaled_error(A,B)
            emit(errio,depth,f,norm(A-B),norm(B),err);maxerr=max(maxerr,err)
        end
    end
    V=zeros(n,n)
    for j in 1:16:n
        ix=j:min(j+15,n);E=zeros(n,length(ix));for (k,i) in enumerate(ix);E[i,k]=1;end
        V[:,ix]=cc_value(result.value,E)
    end
    for (f,A,B) in ((:V,V,ref.value.V),(:v,result.value.v,ref.value.v),(:production_V,V,full[2].V),(:production_v,result.value.v,full[2].v))
        err=scaled_error(A,B);emit(errio,3,f,norm(A-B),norm(B),err);maxerr=max(maxerr,err)
    end
    skew=norm(V-V')/max(1.0,norm(V))
    close(errio)
    runstats=Dict{String,Int}();consumed=cc_consume(result,packet,saved.mu,elems;stats=runstats)
    gate=cc_compare_dirs(dirs,"compact_composed",consumed,packet,rule,referencepath)
    negative=cc_consume(result,packet,saved.mu,elems;corrupt=true)
    negativegate=cc_compare_dirs(dirs,"corrupt_composed_gradient",negative,packet,rule,referencepath)
    @assert !negativegate.pass "negative control must depend on composed operations"
    accurate=gate.pass&&maxerr<1e-7
    println("CHAIN_ACCURACY pass=$accurate values=$maxerr scaled=$(gate.scaled) relative=$(gate.relative) stats=$stats consumer=$runstats");flush(stdout)
    timing=open(joinpath(out,"timings.csv"),"w");emit(timing,"repetition,order,mode,phase,seconds,alloc_bytes,gc_s")
    checks=open(joinpath(out,"repeat_checks.csv"),"w");emit(checks,"repetition,mode,scaled,relative,pass")
    function record(rep,pos,mode,phase,z)
        emit(timing,rep,pos,mode,phase,z.time,z.bytes,z.gctime);z.value
    end
    if accurate
        # Warm EXACT timed paths, including discovery and both consumers.
        cc_inputs(elems,terminal,tree);cc_build(inputs.value);cc_condition(result.value)
        cc_consume(result,packet,saved.mu,elems);cc_dense_consume(ref,packet,saved.mu,elems,terminal)
        for rep in 1:3,(pos,mode) in enumerate(isodd(rep) ? ("dense","compact") : ("compact","dense"))
            GC.gc()
            if mode=="dense"
                d=record(rep,pos,mode,"three_merges",@timed cc_dense(elems,terminal))
                z=record(rep,pos,mode,"downstream_directions",@timed cc_dense_consume(d,packet,saved.mu,elems,terminal))
            else
                fresh=record(rep,pos,mode,"fresh_input_discovery",@timed cc_inputs(elems,terminal,tree))
                st=Dict{String,Int}()
                r=record(rep,pos,mode,"compact_merges_and_close",@timed cc_build(fresh;stats=st))
                record(rep,pos,mode,"condition_probe",@timed cc_condition(r.value))
                z=record(rep,pos,mode,"downstream_directions",@timed cc_consume(r,packet,saved.mu,elems))
            end
            g=cc_compare_dirs(dirs,"repeat$(rep)_$mode",z,packet,rule,referencepath)
            emit(checks,rep,mode,g.scaled,g.relative,g.pass);@assert g.pass
            record(rep,pos,mode,"cleanup_gc",@timed GC.gc())
        end
    end
    close(dirs);close(timing);close(checks)
    # Actual upper-triangle inertia is checked outside timers on the explicitly
    # materialized VALIDATION matrix only, not fed back into execution.
    ev=eigvals(Symmetric(elems[1].Mbb+V,:U));rv=eigvals(Symmetric(elems[1].Mbb+full[2].V,:U))
    threshold=n*eps(Float64)*maximum(abs,rv)
    inertia_ok=count(<(-threshold),ev)==count(<(-threshold),rv)&&count(x->abs(x)<=threshold,ev)==count(x->abs(x)<=threshold,rv)
    @assert inertia_ok
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,Dict("accurate"=>accurate,"nB"=>n,"window"=>4,"N"=>24,"stats"=>stats,"consumer_stats"=>runstats,
          "ordinary_error"=>ordinary_error,"saved_probe_error"=>probe_error,"consumer_reference_scaled"=>consumer_gate.scaled,
          "max_value_scaled"=>maxerr,"direction_scaled"=>gate.scaled,"direction_relative"=>gate.relative,
          "negative_control_rejected"=>!negativegate.pass,"negative_scaled"=>negativegate.scaled,"inertia_ok"=>inertia_ok,
          "symmetry_defect"=>skew,"interior_dimensions"=>[n,2n,3n],
          "augmented_coefficient_bytes"=>8 .* [result.first.K.entries,result.chain.H.entries,result.value.H.entries],
          "input_coefficient_bytes"=>sum(p.bytes for p in inputs.value.p)+8inputs.value.terminal.rep.entries,
          "object_graph_bytes_not_rss"=>Base.summarysize(result),"maxrss_all_paths_bytes"=>Sys.maxrss(),
          "topology_cold_s_not_benchmark"=>topology.time,"factor_discovery_cold_s_not_benchmark"=>inputs.time,
          "closed_future_dense_base_bytes"=>8n*n,"dense_stage_base_bytes"=>8length(packet.S.b0.Sh0),
          "algorithm_dense_exports"=>0,"local_dense_cap"=>CC_LOCAL_FALLBACK_CAP,
          "optimizer_run"=>false,"medium_extended"=>false,"initial_factors_from_network_elimination"=>false,
          "reuse"=>"only identical frozen inputs within this iterate; exclude fresh_input_discovery only for this explicitly conditional reuse case",
          "closure"=>"retained interior dimension grows n to 2n to 3n; no fixed-rank or bounded-size boundary closure claim",
          "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input)))))
    end
    println("COMPACT_CHAIN_COMPLETE accurate=$accurate optimizer=false medium=false");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @assert length(ARGS)==3;cc_run(ARGS...)
end
