# Isolated frozen-iterate accuracy/storage diagnostic. NEVER an optimizer timer.
using Serialization, LinearAlgebra, SparseArrays, SHA, TOML
BLAS.set_num_threads(1)
ENV["FILTERDDP_TREE_KKT"]="1"
ENV["FILTERDDP_FACTOR_BACKED_POLICY"]="1"
const ROOT=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(ROOT,"spar_ddp/examples/power_system/ieee123c_filterddp.jl"))
include(joinpath(ROOT,"spar_ddp/examples/power_system/procs_pass.jl"))
include(joinpath(@__DIR__,"hierarchy.jl"))
include(joinpath(@__DIR__,"production_reference.jl"))

scaled_error(a,b)=norm(a-b)/max(1.0,norm(b))
function emit(io,x...)
    println(io,join(x,','));flush(io)
end
function compare!(io,mode,depth,lo,hi,a,b)
    fields=a isa PintElem ? (:Maa,:Mab,:Mbb,:ma,:mb) : (:V,:v)
    for f in fields
        A=getproperty(a,f);B=getproperty(b,f);d=norm(A-B);s=norm(B)
        @assert isfinite(d) && isfinite(s)
        emit(io,mode,depth,lo,hi,f,d,s,d/max(s,floatmin(Float64)),d/max(1.0,s),maximum(abs,A-B))
    end
end
function inertia!(io,mode,depth,lo,hi,K)
    asym=norm(K-K')/max(1.0,norm(K));ts=time_ns()
    ev=eigvals(Symmetric((K+K')/2));threshold=length(ev)*eps(Float64)*maximum(abs,ev)
    amin=minimum(abs,ev);amax=maximum(abs,ev)
    emit(io,mode,depth,lo,hi,minimum(ev),maximum(ev),threshold,
         count(<(-threshold),ev),count(x->abs(x)<=threshold,ev),
         amin==0 ? Inf : amax/amin,asym,(time_ns()-ts)/1e9)
end
function project(e,groups)
    mats=[zeros(size(e.Maa)) for _ in 1:3]
    for g in groups,(j,f) in enumerate((:Maa,:Mab,:Mbb));mats[j][g,g]=getproperty(e,f)[g,g];end
    PintElem(mats...,copy(e.ma),copy(e.mb))
end
function ordinary(elems,terminal)
    out=Vector{PintVal}(undef,length(elems)+1);out[end]=terminal
    for t in reverse(eachindex(elems));out[t]=pint_combine(elems[t],out[t+1]);end
    out
end
function policy(packet,future,mu)
    S=deepcopy(packet.S);r=_pp_solve(S,future.v,future.V);@assert !isnothing(r)
    rule=DDP4OPF.UpdateRule(Float64,packet.nx,packet.nu,packet.nc)
    # Dummy lambda affects _tp_post's reported dual scalar, NOT the policy.
    # We intentionally never report that scalar as an original dual residual.
    post=_tp_post((update=[rule],),S,r,1,mu,zeros(packet.nx),packet.nx,packet.nu,packet.nc,false,true,packet.cl)
    @assert !isnothing(post)
    rule
end
function directions(rule,dx)
    beta,omega=DDP4OPF.policy_actions!(rule,dx)
    (du=rule.α+beta,dphi=rule.ψ+omega,
     dzl=rule.χl-rule.Σ_L.*beta,dzu=rule.χu+rule.Σ_U.*beta)
end
function battery_path(elems,values,f)
    n=length(elems[1].ma);a=zeros(n);out=Any[]
    for t in eachindex(elems)
        e=elems[t];v=values[t+1]
        b=-(e.Mbb+v.V)\(e.Mab'*a+e.mb+v.v)
        push!(out,(dx=copy(a),dp=(b-a)./f));a=b
    end
    out
end

function run_replay(capturedir,input,outdir;window=0)
    @assert !ispath(outdir);mkpath(outdir)
    meta=TOML.parsefile(joinpath(capturedir,"metadata.toml"))
    leafpath=joinpath(capturedir,"leaves.jls")
    @assert bytes2hex(sha256(read(leafpath)))==meta["leaves_sha256"]
    for r in meta["stage_files"]
        p=joinpath(capturedir,r["file"])
        @assert filesize(p)==r["bytes"] && bytes2hex(sha256(read(p)))==r["sha256"]
    end
    payload=deserialize(leafpath);allleaves=payload.elems;n=length(first(allleaves).ma)
    N=length(allleaves);@assert N==meta["stages"]
    data=deserialize(input);@assert length(data[:Bset])==n
    tree=topology_hierarchy(data;maxleaf=8);groups=coupling_groups(tree)
    @assert sort(vcat(groups...))==collect(1:n)
    # True dense future at the window exit; never an artificial zero boundary.
    ts=time_ns();terminal_zero=PintVal(zeros(n,n),zeros(n))
    full=production_values(allleaves,terminal_zero)
    ordinary_values=ordinary(allleaves,terminal_zero)
    ordinary_error=maximum(max(scaled_error(full[t].V,ordinary_values[t].V),scaled_error(full[t].v,ordinary_values[t].v)) for t in 1:N)
    @assert ordinary_error<1e-7 "ordinary versus production value gate"
    reference_s=(time_ns()-ts)/1e9
    last=window==0 ? N : min(window,N);elems=allleaves[1:last];terminal=full[last+1]
    packets=Dict{Int,Any}()
    for r in meta["stage_files"]
        p=deserialize(joinpath(capturedir,r["file"]));p.stage<=last && (packets[p.stage]=p)
    end
    @assert haskey(packets,1)
    # First gate: reproduced FULL control/equality policies at identical saved
    # zero/sine probes, not merely the battery slice of the Newton direction.
    probes=open(joinpath(outdir,"probes.csv"),"w")
    emit(probes,"stage,probe,field,absolute_fro,reference_fro,scaled_error")
    rules=Dict{Int,Any}()
    for (t,p) in sort(collect(packets);by=first)
        rule=policy(p,full[t+1],payload.mu);rules[t]=rule
        for (i,pr) in enumerate(p.probes)
            d=directions(rule,pr.dx)
            for f in (:du,:dphi)
                A=getproperty(d,f);B=getproperty(pr,f);err=scaled_error(A,B)
                emit(probes,t,i,f,norm(A-B),norm(B),err)
                @assert err<1e-7 "dense full-direction probe gate t=$t probe=$i field=$f error=$err"
            end
        end
    end
    close(probes)
    println("DENSE_DIRECTION_PROBES_PASS tags=$(join(meta["tags"],';')) nB=$n N=$N window=$last");flush(stdout)
    # Record actual stage tree interfaces. These are NOT a temporal rank bound.
    open(joinpath(outdir,"separators.csv"),"w") do io
        emit(io,"stage,bus,own_width,parent_interface_width,retained_battery_energy_rows")
        for (t,p) in sort(collect(packets);by=first),i in eachindex(p.S.lay.own)
            emit(io,t,data[:Nset][i],length(p.S.lay.own[i]),length(p.S.lay.iface[i]),length(p.S.lay.E))
        end
    end
    errors=open(joinpath(outdir,"errors.csv"),"w")
    emit(errors,"mode,depth,lo,hi,block,absolute_fro,reference_fro,relative_fro,unit_floor_scaled,max_abs")
    for t in 1:N;compare!(errors,"dense",-1,t-1,N,full[t],ordinary_values[t]);end
    inert=open(joinpath(outdir,"inertia.csv"),"w")
    emit(inert,"mode,depth,lo,hi,min_eig,max_eig,zero_threshold,negative,near_zero,abs_condition,symmetry_error,eigen_s")
    costs=open(joinpath(outdir,"costs.csv"),"w")
    emit(costs,"mode,depth,lo,hi,combine_s,combine_alloc_bytes,combine_gc_s,construction_s,materialization_s,payload_bytes,dense_three_blocks_bytes,reconstruction_relative_fro")
    ranks=open(joinpath(outdir,"ranks.csv"),"w")
    emit(ranks,"mode,scan_depth,lo,hi,hierarchy_depth,rows,cols,rank,dense_fallback,cross_fro,tail_fro,svd_s,separator_bus")
    spectra=open(joinpath(outdir,"spectra.csv"),"w")
    emit(spectra,"mode,scan_depth,lo,hi,hierarchy_depth,node,singular_index,singular_value")
    function transform(e,mode,depth,lo,hi,tm=(time=0.0,bytes=0,gctime=0.0))
        if e isa PintVal || mode=="dense"
            bytes=e isa PintVal ? 8*(n*n+n) : 8*(3n*n+2n)
            emit(costs,mode,depth,lo,hi,tm.time,tm.bytes,tm.gctime,0.0,0.0,bytes,bytes,0.0)
            return e
        end
        if mode=="drop"
            z=@timed project(e,groups)
            emit(costs,mode,depth,lo,hi,tm.time,tm.bytes,tm.gctime,z.time,0.0,8*(3sum(length(g)^2 for g in groups)+2n),8*(3n*n+2n),norm(joint_value(z.value)-joint_value(e))/max(norm(joint_value(e)),floatmin(Float64)))
            return z.value
        end
        tol=mode=="exact" ? 0.0 : parse(Float64,mode)
        z=retain_coupling(e,tree;tol=tol)
        emit(costs,mode,depth,lo,hi,tm.time,tm.bytes,tm.gctime,z.construction_s,z.materialization_s,z.numeric_payload_bytes,8*(3n*n+2n),z.relative_frobenius_error)
        for (j,r) in enumerate(z.reports)
            emit(ranks,mode,depth,lo,hi,r.depth,r.rows,r.cols,r.rank,r.dense_fallback,r.cross_frobenius,r.tail_frobenius,r.svd_s,join(r.separator_buses,';'))
            # Full spectra recorded once in zero-truncation mode, not duplicated
            # across approximate budgets at the same evolving representation.
            if mode=="exact"
                for (i,s) in enumerate(r.singular_values);emit(spectra,mode,depth,lo,hi,r.depth,j,i,s);end
            end
        end
        z.elem
    end
    values_by_mode=Dict{String,Any}()
    for mode in ("dense","exact","drop","1e-8","1e-5")
        xs=Any[transform(deepcopy(e),mode,0,i-1,i) for (i,e) in enumerate(elems)];push!(xs,deepcopy(terminal))
        for t in 1:last;compare!(errors,mode,0,t-1,t,xs[t],elems[t]);end
        refs=Any[deepcopy(e) for e in elems];push!(refs,deepcopy(terminal))
        ends=vcat(collect(1:last),N);h=1;depth=0
        work=PairWork(n);reference_work=PairWork(n)
        while h<length(xs)
            depth+=1
            for i in 1:2h:length(xs)
                mid=i+h;mid>length(xs) && continue
                # The right operand is unchanged throughout this block/round.
                for t in i:mid-1
                    hi=ends[mid];a=xs[t];b=xs[mid]
                    K=a.Mbb+(b isa PintElem ? b.Maa : b.V)
                    inertia!(inert,mode,depth,t-1,hi,K)
                    timed=@timed _pp_combine!(a,b,work)
                    exact=_pp_combine!(refs[t],refs[mid],reference_work)
                    xs[t]=transform(timed.value,mode,depth,t-1,hi,timed)
                    refs[t]=exact;ends[t]=hi
                    compare!(errors,mode,depth,t-1,hi,xs[t],exact)
                    if mode=="dense" || mode=="exact"
                        fields=exact isa PintElem ? (:Maa,:Mab,:Mbb,:ma,:mb) : (:V,:v)
                        @assert maximum(scaled_error(getproperty(xs[t],f),getproperty(exact,f)) for f in fields)<1e-7
                    end
                end
            end
            h*=2
        end
        for t in 1:last
            @assert xs[t] isa PintVal
            compare!(errors,mode,depth+1,t-1,N,xs[t],full[t])
            if mode=="dense"
                @assert max(scaled_error(xs[t].V,full[t].V),scaled_error(xs[t].v,full[t].v))<1e-7 "ordinary versus scan"
            end
        end
        values_by_mode[mode]=xs
        println("FROZEN_SCAN_COMPLETE mode=$mode tags=$(join(meta["tags"],';')) nB=$n window=$last");flush(stdout)
        GC.gc()
    end
    dirs=open(joinpath(outdir,"directions.csv"),"w")
    emit(dirs,"mode,stage,probe,field,absolute_fro,reference_fro,relative_fro,unit_floor_scaled")
    refpath=battery_path(elems,full,packets[1].f)
    for (mode,values) in values_by_mode
        path=battery_path(elems,values,packets[1].f)
        for t in 1:last,f in (:dx,:dp)
            A=getproperty(path[t],f);B=getproperty(refpath[t],f)
            emit(dirs,mode,t,"path",f,norm(A-B),norm(B),norm(A-B)/max(norm(B),floatmin(Float64)),scaled_error(A,B))
        end
        for (t,p) in sort(collect(packets);by=first)
            rule=policy(p,values[t+1],payload.mu)
            for (i,pr) in enumerate(p.probes)
                d=directions(rule,pr.dx);ref=directions(rules[t],pr.dx)
                for f in (:du,:dphi,:dzl,:dzu)
                    A=getproperty(d,f);B=getproperty(ref,f)
                    emit(dirs,mode,t,i,f,norm(A-B),norm(B),norm(A-B)/max(norm(B),floatmin(Float64)),scaled_error(A,B))
                end
            end
        end
    end
    for io in (errors,inert,costs,ranks,spectra,dirs);close(io);end
    open(joinpath(outdir,"metadata.toml"),"w") do io
        TOML.print(io,Dict("tags"=>meta["tags"],"iteration"=>payload.k,"barrier_iteration"=>payload.j,"mu"=>payload.mu,
            "nB"=>n,"N"=>N,"window"=>last,"groups"=>groups,"reference_s"=>reference_s,
            "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input))),
            "julia"=>string(VERSION),"maxrss_bytes"=>Sys.maxrss(),"compact_kernel"=>false,
            "timing_is_benchmark"=>false,"full_direction_probe_gate"=>true,
            "reference_path"=>"production_inplace","ordinary_value_max_error"=>ordinary_error))
    end
    println("FROZEN_REPLAY_COMPLETE out=$outdir");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==4 || error("replay.jl CAPTURE_DIR INPUT OUTPUT WINDOW(0=full)")
    run_replay(ARGS[1],ARGS[2],ARGS[3];window=parse(Int,ARGS[4]))
end
