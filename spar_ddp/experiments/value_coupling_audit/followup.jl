# No optimizer invocation. Existing production reference and policies are reused.
include(joinpath(ENV["VALUE_FOLLOWUP_BASE"],"replay.jl"))
include(joinpath(@__DIR__,"followup_blocks.jl"))
const DIR_SCALED=1e-7
const DIR_RELATIVE=1e-5

function measure_dirs(io,label,elems,values,reference,packets,mu,ref_rules)
    worst=0.0;relative=0.0
    function record(t,probe,f,A,B)
        absolute_error=norm(A-B);reference_norm=norm(B)
        scaled_value=absolute_error/max(1.0,reference_norm)
        relative_value=absolute_error/max(reference_norm,1e-10)
        @assert isfinite(scaled_value) && isfinite(relative_value)
        emit(io,label,t,probe,f,absolute_error,reference_norm,scaled_value,relative_value)
        worst=max(worst,scaled_value);relative=max(relative,relative_value)
    end
    a=battery_path(elems,values,packets[1].f);b=battery_path(elems,reference,packets[1].f)
    for t in eachindex(elems),f in (:dx,:dp);record(t,"path",f,getproperty(a[t],f),getproperty(b[t],f));end
    for (t,p) in sort(collect(packets);by=first)
        rule=policy(p,values[t+1],mu)
        for (i,pr) in enumerate(p.probes)
            d=directions(rule,pr.dx);r=directions(ref_rules[t],pr.dx)
            for f in (:du,:dphi,:dzl,:dzu);record(t,i,f,getproperty(d,f),getproperty(r,f));end
        end
    end
    (scaled=worst,relative=relative,pass=worst<DIR_SCALED && relative<DIR_RELATIVE)
end
function followup_scan(elems,terminal,tree,base_tol,conditioned,io,label;detail=nothing,inertia_io=nothing)
    n=length(first(elems).ma);payloads=Int[];ranks=Int[];fallbacks=0
    # Whole-call allocation/GC additionally records wrappers and deepcopy costs.
    function transform(e,depth,t,kappa)
        e isa PintVal && return e
        tau=conditioned ? controlled_tol(base_tol,kappa) : base_tol
        a=@timed construct_blocks(e,tree,tau)
        b=@timed materialize_blocks(a.value,e,tau)
        push!(io,(label,depth,t,"construction",a.time,a.bytes,a.gctime,b.value.bytes,8*(3n*n+2n),tau,kappa,b.value.fallback))
        push!(io,(label,depth,t,"materialization",b.time,b.bytes,b.gctime,b.value.bytes,8*(3n*n+2n),tau,kappa,b.value.fallback))
        push!(payloads,b.value.bytes);fallbacks+=Int(b.value.fallback)
        for r in a.value.reports
            push!(ranks,r.rank)
            isnothing(detail) || emit(detail,label,depth,t,r.field,r.depth,r.side,r.rows,r.cols,r.rank,r.dense,r.tail,r.scale)
        end
        b.value.elem
    end
    function condition(K,depth,t)
        !conditioned && return 1.0
        z=@timed elimination_condition(K)
        push!(io,(label,depth,t,"condition_control",z.time,z.bytes,z.gctime,0,0,base_tol,z.value,false))
        z.value
    end
    xs=Any[]
    for (t,e) in enumerate(elems)
        k=condition(e.Maa+e.Mbb,0,t)
        push!(xs,transform(deepcopy(e),0,t,k))
    end
    push!(xs,deepcopy(terminal));h=1;depth=0;work=PairWork(n)
    refs=isnothing(inertia_io) ? nothing : Any[deepcopy(elems)...,deepcopy(terminal)]
    refwork=isnothing(inertia_io) ? nothing : PairWork(n);inertia_ok=true
    while h<length(xs)
        depth+=1
        for i in 1:2h:length(xs)
            mid=i+h;mid>length(xs) && continue
            for t in i:mid-1
                a=xs[t];b=xs[mid]
                K=a.Mbb+(b isa PintElem ? b.Maa : b.V)
                k=condition(K,depth,t)
                if !isnothing(inertia_io)
                    ra=refs[t];rb=refs[mid];RK=ra.Mbb+(rb isa PintElem ? rb.Maa : rb.V)
                    ev=eigvals(Symmetric(K,:U));rv=eigvals(Symmetric(RK,:U))
                    threshold=n*eps(Float64)*maximum(abs,rv)
                    neg=count(<(-threshold),ev);rneg=count(<(-threshold),rv)
                    near=count(x->abs(x)<=threshold,ev);rnear=count(x->abs(x)<=threshold,rv)
                    same=neg==rneg && near==rnear;inertia_ok &= same
                    emit(inertia_io,label,depth,t,minimum(ev),maximum(ev),threshold,neg,near,rneg,rnear,same)
                    refs[t]=_pp_combine!(ra,rb,refwork)
                end
                z=@timed _pp_combine!(a,b,work)
                push!(io,(label,depth,t,"combination",z.time,z.bytes,z.gctime,0,0,base_tol,k,false))
                xs[t]=transform(z.value,depth,t,k)
            end
        end
        h*=2
    end
    (values=xs,payload_ratio=isempty(payloads) ? 1.0 : sum(payloads)/(length(payloads)*8*(3n*n+2n)),
     rank_max=isempty(ranks) ? 0 : maximum(ranks),fallbacks=fallbacks,inertia_ok=inertia_ok)
end
function run_followup(cap,input,out,window)
    @assert !ispath(out);mkpath(out)
    meta=TOML.parsefile(joinpath(cap,"metadata.toml"))
    leaf=joinpath(cap,"leaves.jls");@assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
    for r in meta["stage_files"]
        p=joinpath(cap,r["file"]);@assert filesize(p)==r["bytes"] && bytes2hex(sha256(read(p)))==r["sha256"]
    end
    payload=deserialize(leaf);raw=payload.elems;n=length(first(raw).ma);N=length(raw)
    data=deserialize(input);@assert length(data[:Bset])==n
    tree=topology_hierarchy(data;maxleaf=8)
    full=production_values(raw,PintVal(zeros(n,n),zeros(n)))
    ord=ordinary(raw,PintVal(zeros(n,n),zeros(n)))
    ordinary_error=maximum(max(scaled_error(full[t].V,ord[t].V),scaled_error(full[t].v,ord[t].v)) for t in 1:N)
    @assert ordinary_error<1e-7
    last=window==0 ? N : window;elems=raw[1:last];terminal=full[last+1]
    packets=Dict{Int,Any}()
    for r in meta["stage_files"]
        p=deserialize(joinpath(cap,r["file"]));p.stage<=last && (packets[p.stage]=p)
    end
    rules=Dict{Int,Any}();probe_error=0.0
    for (t,p) in packets
        rule=policy(p,full[t+1],payload.mu);rules[t]=rule
        for pr in p.probes
            d=directions(rule,pr.dx)
            for f in (:du,:dphi)
                err=scaled_error(getproperty(d,f),getproperty(pr,f));probe_error=max(probe_error,err)
                @assert err<1e-7 "unchanged dense saved-probe gate"
            end
        end
    end
    println("FOLLOWUP_REFERENCE_PASS n=$n ordinary=$ordinary_error probes=$probe_error");flush(stdout)
    dirs=open(joinpath(out,"directions.csv"),"w");emit(dirs,"mode,stage,probe,field,absolute,reference,scaled,relative_floor_1e_10")
    costs=open(joinpath(out,"costs.csv"),"w");emit(costs,"mode,depth,stage,component,seconds,alloc_bytes,gc_s,payload_bytes,dense_bytes,tol,kappa,raw_fallback")
    costrows=Any[]
    rankio=open(joinpath(out,"ranks.csv"),"w");emit(rankio,"mode,scan_depth,stage,block,hierarchy_depth,side,rows,cols,rank,dense,tail,scale")
    inertio=open(joinpath(out,"inertia.csv"),"w");emit(inertio,"mode,depth,stage,min_eig,max_eig,threshold,negative,near_zero,reference_negative,reference_near_zero,same")
    trials=Any[];chosen=nothing
    # Baseline is evaluated first using the same window association and future.
    baseline=production_values(elems,terminal)
    gate=measure_dirs(dirs,"dense_window",elems,baseline,full,packets,payload.mu,rules)
    @assert gate.pass "dense window directions must satisfy unchanged accuracy targets"
    for (label,tol,controlled) in (("block_1e-6",1e-6,false),("condition_1e-6",1e-6,true),("condition_1e-9",1e-9,true),("condition_1e-12",1e-12,true))
        z=@timed followup_scan(elems,terminal,tree,tol,controlled,costrows,label;detail=rankio,inertia_io=inertio)
        d=measure_dirs(dirs,label,elems,z.value.values,full,packets,payload.mu,rules)
        push!(trials,Dict("mode"=>label,"tol"=>tol,"conditioned"=>controlled,"scaled"=>d.scaled,"relative"=>d.relative,"pass"=>d.pass,
            "payload_ratio"=>z.value.payload_ratio,"rank_max"=>z.value.rank_max,"raw_fallbacks"=>z.value.fallbacks,"inertia_ok"=>z.value.inertia_ok,"diagnostic_s"=>z.time))
        println("FOLLOWUP_TRIAL mode=$label scaled=$(d.scaled) relative=$(d.relative) pass=$(d.pass) payload=$(z.value.payload_ratio)");flush(stdout)
        if controlled && d.pass && z.value.inertia_ok && z.value.payload_ratio<0.99
            chosen=(label=label,tol=tol,conditioned=controlled);break
        end
        GC.gc()
    end
    # If no candidate passes, raw production fallback is explicit, not a
    # compression success. Do not time a fake compact algorithm or run solves.
    bench=open(joinpath(out,"warmed.csv"),"w")
    emit(bench,"repetition,order,mode,total_s,alloc_bytes,gc_s,payload_ratio,rank_max,fallbacks,scaled,relative")
    if !isnothing(chosen)
        # Warm the exact chosen code paths and dense baseline completely.
        warm=followup_scan(elems,terminal,tree,chosen.tol,true,costrows,"warmup")
        @assert measure_dirs(dirs,"warmup",elems,warm.values,full,packets,payload.mu,rules).pass
        production_values(elems,terminal);GC.gc()
        for rep in 1:3
            order=isodd(rep) ? ("dense","candidate") : ("candidate","dense")
            for (pos,mode) in enumerate(order)
                GC.gc() # Same boundary policy; explicit GC is outside each timer.
                if mode=="dense"
                    z=@timed production_values(elems,terminal)
                    d=measure_dirs(dirs,"repeat$(rep)_dense",elems,z.value,full,packets,payload.mu,rules)
                    emit(bench,rep,pos,mode,z.time,z.bytes,z.gctime,1.0,0,0,d.scaled,d.relative)
                else
                    z=@timed followup_scan(elems,terminal,tree,chosen.tol,true,costrows,"repeat$rep")
                    d=measure_dirs(dirs,"repeat$(rep)_candidate",elems,z.value.values,full,packets,payload.mu,rules)
                    @assert d.pass
                    emit(bench,rep,pos,mode,z.time,z.bytes,z.gctime,z.value.payload_ratio,z.value.rank_max,z.value.fallbacks,d.scaled,d.relative)
                end
                @assert d.pass
            end
        end
    else
        println("FOLLOWUP_NO_ACCURATE_COMPRESSED_CANDIDATE raw_dense_fallback=true");flush(stdout)
    end
    for row in costrows;emit(costs,row...);end
    for io in (dirs,costs,rankio,inertio,bench);close(io);end
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,Dict("trials"=>trials,"chosen"=>isnothing(chosen) ? "none_raw_fallback" : chosen.label,
          "nB"=>n,"N"=>N,"window"=>last,"tags"=>meta["tags"],"probe_error"=>probe_error,"ordinary_error"=>ordinary_error,
          "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input))),
          "direction_scaled_gate"=>DIR_SCALED,"direction_relative_gate"=>DIR_RELATIVE,"relative_floor"=>1e-10,
          "compact_kernel"=>false,"resident_memory_saving_demonstrated"=>false,"optimizer_run"=>false,
          "condition_control_charged"=>true,"offline_direction_oracle_used_for_selection"=>true,
          "component_io_in_timing"=>false,"in_memory_timing_records_charged"=>true,"maxrss_bytes"=>Sys.maxrss(),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads()))
    end
    println("FOLLOWUP_COMPLETE out=$out");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @assert length(ARGS)==4
    run_followup(ARGS[1],ARGS[2],ARGS[3],parse(Int,ARGS[4]))
end
