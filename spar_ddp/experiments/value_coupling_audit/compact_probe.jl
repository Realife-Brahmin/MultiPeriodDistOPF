# One actual leaf/leaf component, plus charged legacy-consumer hybrid replay.
# This file never invokes an optimizer or expands medium cases automatically.
include(joinpath(@__DIR__,"followup.jl"))
include(joinpath(@__DIR__,"compact_kernel.jl"))
function ck_continue(elems,terminal,firstmerge)
    xs=Any[deepcopy(e) for e in elems];push!(xs,deepcopy(terminal))
    work=PairWork(length(first(elems).ma));h=1
    while h<length(xs)
        for i in 1:2h:length(xs)
            mid=i+h;mid>length(xs) && continue
            for t in i:mid-1
                xs[t]=(h==1 && t==3 && mid==4) ? deepcopy(firstmerge) : _pp_combine!(xs[t],xs[mid],work)
            end
        end
        h*=2
    end
    xs
end
function ck_probe_actions(c,P)
    [ck_action(c,f,P) for f in (:Maa,:Mab,:Mbb)]
end
ck_dense_actions(d,P)=[getproperty(d,f)*P for f in (:Maa,:Mab,:Mbb)]
function ck_run(capturedir,input,out)
    @assert !ispath(out);mkpath(out)
    meta=TOML.parsefile(joinpath(capturedir,"metadata.toml"))
    leaf=joinpath(capturedir,"leaves.jls")
    @assert bytes2hex(sha256(read(leaf)))==meta["leaves_sha256"]
    for r in meta["stage_files"]
        f=joinpath(capturedir,r["file"])
        @assert filesize(f)==r["bytes"] && bytes2hex(sha256(read(f)))==r["sha256"]
    end
    payload=deserialize(leaf);raw=payload.elems;n=length(first(raw).ma)
    @assert n==1020 && payload.j>0 && length(raw)==24 && "late" in meta["tags"]
    data=deserialize(input);@assert length(data[:Bset])==n
    hierarchy=@timed topology_hierarchy(data;maxleaf=8);tree=hierarchy.value
    full=production_values(raw,PintVal(zeros(n,n),zeros(n)))
    ord=ordinary(raw,PintVal(zeros(n,n),zeros(n)))
    @assert maximum(max(scaled_error(full[t].V,ord[t].V),scaled_error(full[t].v,ord[t].v)) for t in eachindex(raw))<1e-7
    elems=raw[1:4];terminal=full[5];x=elems[3];y=elems[4]
    packets=Dict{Int,Any}();rules=Dict{Int,Any}();probe_error=0.0
    for r in meta["stage_files"]
        packet=deserialize(joinpath(capturedir,r["file"]))
        packet.stage>4 && continue
        packets[packet.stage]=packet;rule=policy(packet,full[packet.stage+1],payload.mu);rules[packet.stage]=rule
        for pr in packet.probes,f in (:du,:dphi)
            err=scaled_error(getproperty(directions(rule,pr.dx),f),getproperty(pr,f))
            probe_error=max(probe_error,err);@assert err<1e-7
        end
    end
    dirs=open(joinpath(out,"directions.csv"),"w");emit(dirs,"mode,stage,probe,field,absolute,reference,scaled,relative_floor_1e_10")
    @assert measure_dirs(dirs,"dense",elems,production_values(elems,terminal),full,packets,payload.mu,rules).pass
    println("COMPACT_REFERENCE_PASS n=$n saved_probe_error=$probe_error");flush(stdout)
    input_tol=1e-9;update_tol=1e-12
    setup=@timed (construct_blocks(x,tree,input_tol),construct_blocks(y,tree,input_tol))
    p,q=setup.value
    stats=Dict{String,Int}()
    attempt=try
        @timed ck_merge(p,q;tol=update_tol,stats=stats)
    catch exception
        if exception isa SingularException || occursin("COMPACT_RANK_GROWTH",sprint(showerror,exception))
            # Explicit scientific/representation rejection, not a weakened gate.
            fallback=@timed _pp_combine!(deepcopy(x),y,PairWork(n))
            open(joinpath(out,"metadata.toml"),"w") do io
                TOML.print(io,Dict("outcome"=>"compact_rejected_dense_fallback","reason"=>sprint(showerror,exception),
                  "stats"=>stats,"fallback_s"=>fallback.time,"fallback_alloc_bytes"=>fallback.bytes,
                  "input_setup_s_cold_not_benchmark"=>setup.time,"nB"=>n,"probe_error"=>probe_error,
                  "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input))),
                  "optimizer_run"=>false,"medium_extended"=>false))
            end
            close(dirs);println("COMPACT_COMPLETE rejected=true reason=$(sprint(showerror,exception))");return
        end
        rethrow()
    end
    c=attempt.value
    exported=ck_export(c) # Validation-only full export; timed separately below.
    reference=_pp_combine!(deepcopy(x),y,PairWork(n))
    errors=open(joinpath(out,"value_errors.csv"),"w")
    emit(errors,"phase,field,absolute,reference,scaled")
    value_error=0.0
    for f in (:Maa,:Mab,:Mbb,:ma,:mb)
        A=getproperty(exported,f);B=getproperty(reference,f);err=scaled_error(A,B)
        emit(errors,"merge",f,norm(A-B),norm(B),err);value_error=max(value_error,err)
    end
    hybrid=ck_continue(elems,terminal,exported)
    for t in 1:4,f in (:V,:v)
        A=getproperty(hybrid[t],f);B=getproperty(full[t],f);err=scaled_error(A,B)
        emit(errors,"suffix$t",f,norm(A-B),norm(B),err);value_error=max(value_error,err)
    end
    gate=measure_dirs(dirs,"compact_hybrid",elems,hybrid,full,packets,payload.mu,rules)
    # Prove that the saved stage1 future/policy actually depends on this merge.
    # The old first merge changes only suffix1, unused by these direction
    # probes. Merge3/4 feeds suffix2 and thus the selected stage1 policy.
    perturbed=deepcopy(exported)
    perturbed.ma .+= 1e6 .* sin.(collect(1:n))
    perturbed.mb .+= 1e6 .* cos.(collect(1:n))
    negative_values=ck_continue(elems,terminal,perturbed)
    negative_gate=measure_dirs(dirs,"sensitivity_negative_control",elems,negative_values,full,packets,payload.mu,rules)
    negative_full=0.0
    for (t,packet) in packets
        nrule=policy(packet,negative_values[t+1],payload.mu)
        for pr in packet.probes,f in (:du,:dphi,:dzl,:dzu)
            delta=scaled_error(getproperty(directions(nrule,pr.dx),f),getproperty(directions(rules[t],pr.dx),f))
            negative_full=max(negative_full,delta)
        end
    end
    @assert !negative_gate.pass && negative_full>1e-7 "direction audit must detect a corrupted tested merge"
    P=hcat(ones(n),sin.(collect(1:n)),cos.(collect(1:n)))
    actions=ck_probe_actions(c,P);action_error=maximum(scaled_error(actions[i],getproperty(reference,f)*P) for (i,f) in enumerate((:Maa,:Mab,:Mbb)))
    solve_error=scaled_error(ck_apply(c.K,ck_solve(c.factor,P)),P)
    ce=ck_condition_estimate(c)
    println("COMPACT_ACCURACY values=$value_error directions=$(gate.scaled) relative=$(gate.relative) actions=$action_error solve=$solve_error pass=$(gate.pass)");flush(stdout)
    timings=open(joinpath(out,"timings.csv"),"w")
    emit(timings,"repetition,order,mode,phase,seconds,alloc_bytes,gc_s")
    checks=open(joinpath(out,"repeat_checks.csv"),"w");emit(checks,"repetition,values_scaled,directions_scaled,directions_relative,pass")
    function timed_row(rep,pos,mode,phase,z)
        emit(timings,rep,pos,mode,phase,z.time,z.bytes,z.gctime)
        z.value
    end
    accurate=gate.pass && value_error<1e-7 && action_error<1e-7 && solve_error<1e-7
    if accurate
        # All exact paths warmed; no promise of zero OS/noise, retain all pairs.
        ck_condition_estimate(c);ck_probe_actions(c,P);ck_export(c);ck_continue(elems,terminal,exported)
        warm_dense=_pp_combine!(deepcopy(x),y,PairWork(n));ck_dense_actions(warm_dense,P);GC.gc()
        for rep in 1:3
            for (pos,mode) in enumerate(isodd(rep) ? ("dense","compact") : ("compact","dense"))
                GC.gc()
                if mode=="dense"
                    w=timed_row(rep,pos,mode,"workspace",@timed PairWork(n))
                    d=timed_row(rep,pos,mode,"merge",@timed _pp_combine!(deepcopy(x),y,w))
                    timed_row(rep,pos,mode,"actions",@timed ck_dense_actions(d,P))
                    timed_row(rep,pos,mode,"legacy_downstream",@timed ck_continue(elems,terminal,d))
                else
                    pq=timed_row(rep,pos,mode,"fresh_input_factor_discovery",@timed (construct_blocks(x,tree,input_tol),construct_blocks(y,tree,input_tol)))
                    repeatstats=Dict{String,Int}()
                    compact=timed_row(rep,pos,mode,"merge_factor_updates",@timed ck_merge(pq[1],pq[2];tol=update_tol,stats=repeatstats))
                    timed_row(rep,pos,mode,"condition_probe",@timed ck_condition_estimate(compact))
                    timed_row(rep,pos,mode,"actions",@timed ck_probe_actions(compact,P))
                    d=timed_row(rep,pos,mode,"legacy_dense_export",@timed ck_export(compact))
                    future=timed_row(rep,pos,mode,"legacy_downstream",@timed ck_continue(elems,terminal,d))
                    # All accuracy verification remains outside phase timers.
                    g=measure_dirs(dirs,"repeat$rep",elems,future,full,packets,payload.mu,rules)
                    ve=maximum(scaled_error(getproperty(d,f),getproperty(reference,f)) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
                    emit(checks,rep,ve,g.scaled,g.relative,g.pass && ve<1e-7);@assert g.pass && ve<1e-7
                end
                timed_row(rep,pos,mode,"cleanup_gc",@timed GC.gc())
            end
        end
    end
    close(errors);close(dirs);close(timings);close(checks)
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,Dict("outcome"=>accurate ? "compact_accuracy_pass" : "compact_accuracy_rejected",
          "nB"=>n,"window"=>4,"N"=>24,"compact_component"=>true,"compact_full_scan"=>false,
          "tested_merge_left_period"=>3,"direction_audit_negative_control_rejected"=>!negative_gate.pass,
          "negative_control_full_direction_scaled"=>negative_full,"dense_actions_exact_path_warmed"=>true,
          "legacy_export_required"=>true,"full_optimizer_run"=>false,"medium_extended"=>false,
          "input_tol"=>input_tol,"update_tol"=>update_tol,"stats"=>stats,
          "value_max_scaled"=>value_error,"action_max_scaled"=>action_error,"solve_residual_scaled"=>solve_error,
          "direction_max_scaled"=>gate.scaled,"direction_relative_floor_1e_10"=>gate.relative,
          "saved_probe_error"=>probe_error,"condition_probe"=>ce.probe_condition,"condition_probe_certified"=>false,
          "input_coefficient_bytes"=>p.bytes+q.bytes,"raw_six_blocks_bytes"=>8*(6n*n+4n),
          "raw_input_Maa_skew"=>norm(x.Maa-x.Maa')/max(1.0,norm(x.Maa)),
          "raw_input_Mbb_skew"=>norm(x.Mbb-x.Mbb')/max(1.0,norm(x.Mbb)),
          "compact_object_graph_bytes_not_rss"=>Base.summarysize(c),"maxrss_all_phases_bytes"=>Sys.maxrss(),
          "topology_setup_s"=>hierarchy.time,"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
          "leaves_sha256"=>meta["leaves_sha256"],"input_sha256"=>bytes2hex(sha256(read(input))),
          "reuse_assumption"=>"inputs may be reused only while this frozen iterate is unchanged; fresh factor discovery always charged separately",
          "operator_deferred_cost"=>"Schur products remain factored operators; actions and dense legacy export charged; no free output-materialization win",
          "initial_cold_clocks_are_benchmarks"=>false))
    end
    println("COMPACT_COMPLETE accurate=$accurate medium_extended=false optimizer=false");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @assert length(ARGS)==3
    ck_run(ARGS...)
end
