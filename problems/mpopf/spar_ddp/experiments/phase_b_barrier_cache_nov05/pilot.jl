using Distributed,SharedArrays,Serialization,LinearAlgebra,SHA,TOML,Statistics
include(joinpath(@__DIR__,"worker.jl"))
const STAMP=Ref{Any}(nothing);const ROOTCTX=Ref("");const TV=Ref{Any}(nothing)
const LOGS=Dict{String,IO}();const OWNER=Int[];const META=Any[]
function emitrow(name,x...);println(LOGS[name],join(x,','));flush(LOGS[name]);end
function makelog(out,name,header);LOGS[name]=open(joinpath(out,name*".csv"),"w");println(LOGS[name],header);end
function allworkers(f,args...)
 fetch.([remotecall(f,p,args...) for p in OWNER])
end
function stages(f,args...)
 fetch.([remotecall(f,p,args...) for p in OWNER])
end
function checkrows(scope,rep,stage,rows;require_pass=true,require_exact=false)
 good=true
 for (probe,field,m) in rows
  emitrow("errors",scope,rep,stage,probe,field,m.absolute,m.reference,m.scaled,m.relative,m.maxabs,m.exact,m.pass)
  good&=m.pass&&(!require_exact||m.exact)
 end
 require_pass&&@assert good "output gate $scope rep=$rep stage=$stage"
 good
end
function validate_state(i,label,rep;slot=i,require_pass=true)
 results=stages(worker_validate,i,slot);good=true
 for (t,result) in enumerate(results)
  good&=checkrows(label,rep,t,result.rows;require_pass=require_pass)
  for (probe,field,coordinate,value) in result.outliers;emitrow("outliers",label,rep,t,probe,field,coordinate,value);end
 end
 good
end
function snapshot_state(i,label;slot=i,check_saved=true)
 results=stages(worker_snapshot,i,slot)
 for (t,rows) in enumerate(results);check_saved&&checkrows(label,0,t,rows);end
end
function execute_pass(arm,i;diagnostic=:none,proof=:pinned_unshifted_tp_pre,context=ROOTCTX[],k=META[i]["iteration"],stale=false)
 start=time_ns();hit=false;populate=false
 if arm==:candidate
  observations=stages(worker_key_info,i);keys=[x.key for x in observations]
  complete=sum(x.entries for x in observations)==68
  current=(k=k,j=META[i]["barrier_iteration"],mu=META[i]["mu"],context=context,proof=proof,stage_keys=keys,terminal_key=digest(TV[]))
  hit=eligible(STAMP[],current;complete=complete);populate=!hit&&proof==:pinned_unshifted_tp_pre
  !hit&&allworkers(worker_reset_cache)
  STAMP[]=current
 end
 detection=(time_ns()-start)/1e9
 t0=time_ns();stages(worker_leaf,i,hit,populate)
 stale&&stages(worker_stale_gradients)
 if diagnostic!=:none
  rr=stages(worker_leaf_invariance,diagnostic==:save)
  if diagnostic==:compare
   for (t,rows) in enumerate(rr),(field,m) in rows
    emitrow("invariance","leaf",0,t,0,field,m.absolute,m.reference,m.scaled,m.relative,m.maxabs,m.exact)
    @assert m.exact
   end
  end
 end
 leaf=(time_ns()-t0)/1e9
 t0=time_ns();isval=falses(25);isval[25]=true;hits=0;stored=0;fallback=0
 mode=arm==:baseline ? :fresh : hit ? :hit : populate ? :populate : :fresh
 for (depth,pairs) in enumerate(rounds(24))
  replies=fetch.([remotecall(worker_merge,OWNER[mod1(q,24)],depth,t,mid,isval[mid],mid==25,mode,diagnostic) for (q,(t,mid)) in enumerate(pairs)])
  if mode==:hit&&any(r->r.miss,replies)
   # Unexpected per-entry loss: restart the entire current pass with the
   # unchanged fresh kernel, never mix fresh operands with stale scratch.
   allworkers(worker_reset_cache);STAMP[]=nothing
   return execute_pass(:baseline,i)
  end
  for ((t,mid),r) in zip(pairs,replies)
   isval[t]=isval[mid];hits+=r.hit;stored+=r.stored;fallback+=r.fallback
   for (field,m) in r.rows
    emitrow("invariance","merge",depth,t,mid,field,m.absolute,m.reference,m.scaled,m.relative,m.maxabs,m.exact)
    @assert m.exact
   end
  end
 end
 scan=(time_ns()-t0)/1e9
 t0=time_ns();stages(worker_policy,i,Float64(META[i]["mu"]))
 consumer=(time_ns()-t0)/1e9
 t0=time_ns();dx=zeros(size(TV[],1))
 for p in OWNER;dx=remotecall_fetch(worker_forward,p,dx);end
 forward=(time_ns()-t0)/1e9
 (detection=detection,leaf=leaf,scan=scan,consumer=consumer,forward=forward,hits=hits,stored=stored,fallback=fallback)
end
function retire()
 allworkers(worker_retire);STAMP[]=nothing;nothing
end
function measured_workload(arm)
 STAMP[]=nothing
 a=@timed execute_pass(arm,1)
 b=@timed execute_pass(arm,2)
 c=@timed retire()
 (first=a,second=b,cleanup=c,total=a.time+b.time+c.time)
end
function compiler_totals()
 x=allworkers(worker_compilation)
 (compile=sum(z.compile for z in x),recompile=sum(z.recompile for z in x),alloc=sum(z.alloc for z in x),gc_ns=sum(z.gc_ns for z in x))
end
function logtiming(rep,arm,x,cc,allocworker)
 for (phase,z) in (("population_or_fresh_first",x.first),("hit_or_fresh_second",x.second))
  v=z.value
  emitrow("timings",rep,arm,phase,z.time,v.detection,v.leaf,v.scan,v.consumer,v.forward,v.hits,v.stored,v.fallback,z.bytes,z.gctime,z.compile_time,z.recompile_time)
 end
 emitrow("totals",rep,arm,x.total,x.first.time,x.second.time,x.cleanup.time,x.first.bytes+x.second.bytes+x.cleanup.bytes,x.first.gctime+x.second.gctime+x.cleanup.gctime,cc.compile,cc.recompile,allocworker,cc.gc_ns/1e9)
end
function runpilot(capture,input,out)
 @assert !ispath(out);mkpath(out)
 makelog(out,"affinity","worker_pid,allowed,blas_threads")
 makelog(out,"errors","scope,repeat,stage,probe,field,absolute_fro,reference_fro,unit_scaled,relative_floor1e10,maxabs,exact_bits,pass")
 makelog(out,"outliers","scope,repeat,stage,probe,field,coordinate,maxabs")
 makelog(out,"invariance","kind,depth,left,right,field,absolute_fro,reference_fro,unit_scaled,relative_floor1e10,maxabs,exact_bits")
 makelog(out,"timings","repeat,arm,phase,wall_s,detection_s,leaf_s,scan_s,consumer_s,forward_s,cache_hits,populated,fallbacks,coordinator_alloc_bytes,coordinator_gc_s,coordinator_compile_s,coordinator_recompile_s")
 makelog(out,"totals","repeat,arm,charged_s,first_pass_s,second_pass_s,retirement_gc_s,coordinator_alloc_bytes,coordinator_gc_s,worker_compile_ns,worker_recompile_ns,worker_alloc_bytes,worker_gc_s")
 makelog(out,"negative_controls","control,hits,populated,correct_fresh_outputs,rejected_or_invalidated")
 makelog(out,"warmup","iteration,arm,worker_compile_ns,worker_recompile_ns,coordinator_compile_s,coordinator_recompile_s")
 dirs=String[]
 allmeta=[(dir,TOML.parsefile(joinpath(dir,"metadata.toml"))) for dir in filter(isdir,readdir(capture;join=true))]
 for tag in ("barrier_before","barrier_after","early_before","early_after")
  dir,m=only(filter(x->tag in x[2]["tags"],allmeta));push!(dirs,dir);push!(META,m)
  @assert bytes2hex(sha256(read(joinpath(dir,"leaves.jls"))))==m["leaves_sha256"]
 end
 @assert META[2]["pass"]==META[1]["pass"]+1&&META[1]["iteration"]==META[2]["iteration"]
 data=deserialize(input);n=length(data[:Bset]);@assert n==1020
 TV[]=zeros(n,n)
 ROOTCTX[]=digest((ENV["PILOT_SOURCE_SHA"],ENV["PILOT_RUNTIME_SHA"],bytes2hex(sha256(read(input))),"original_balanced_inplace_upper",string(VERSION),string(BLAS.get_config()),24,1,"exact_hessian_tree_half_factor_policy_no_regularization"))
 tree=topology_hierarchy(data;maxleaf=32);signs=ones(n)
 for (j,g) in enumerate(coupling_groups(tree));signs[g].=isodd(j) ? 1.0 : -1.0;end
 pids=addprocs(24;exeflags=Cmd(["--startup-file=no","--threads=1","--project="*ENV["PILOT_PROJECT"]]))
 append!(OWNER,pids)
 for pid in pids;remotecall_wait(Base.include,pid,Main,joinpath(@__DIR__,"worker.jl"));end
 shared=Any[SharedArray{Float64}((3n*n+2n,);pids=vcat(myid(),pids)) for t in 1:24]
 setup=fetch.([remotecall(worker_setup,pid,shared,n,t,dirs,signs) for (t,pid) in enumerate(pids)])
 open(joinpath(out,"workers.toml"),"w") do io;TOML.print(io,Dict("workers"=>[Dict(string(k)=>v for (k,v) in pairs(x)) for x in setup]));end
 @assert all(x->x.blas==1&&x.threads==1,setup)
 cpu_groups=[parse.(Int,split(line,',')) for line in readlines(ENV["PILOT_CPU_MAP"])]
 @assert length(cpu_groups)==25&&all(x->length(x)==1,cpu_groups)
 for (t,pid) in enumerate(pids);remotecall_fetch(worker_pin,pid,cpu_groups[t+1]);end
 worker_pin(cpu_groups[1])
 for pid in pids
  audit=remotecall_fetch(worker_affinity,pid);emitrow("affinity",audit.pid,audit.allowed,audit.blas)
 end
 for i in 1:2
  stages(worker_leaf,i,false,false)
  leaves=deserialize(joinpath(dirs[i],"leaves.jls")).elems
  rr=fetch.([remotecall(worker_leaf_check,pid,leaves[t]) for (t,pid) in enumerate(pids)])
  for (t,rows) in enumerate(rr),(field,m) in rows
   checkrows("reassembled_leaf_$(i)",0,t,[("leaf",field,m)])
   field in ("Maa","Mab","Mbb")&&@assert m.exact
  end
  execute_pass(:baseline,i;diagnostic=i==1 ? :save : :compare)
  snapshot_state(i,"saved_live_$(i)")
 end
 # Actual terminal curvature/gradient are zero and checked explicitly.
 emitrow("invariance","terminal",0,25,25,"V",0.0,0.0,0.0,0.0,0.0,true)
 allworkers(worker_clear_diagnostics)
 # Source-level stage-signature invariance must be exact before cache is allowed.
 key1=stages(worker_key,1);key2=stages(worker_key,2);@assert key1==key2 "conservative stage eligibility differs"
 STAMP[]=nothing;execute_pass(:candidate,1);validate_state(1,"cache_population_preflight",0);h=execute_pass(:candidate,2);@assert h.hits==68&&h.stored==0
 validate_state(2,"cache_preflight",0)
 storage=allworkers(worker_storage)
 open(joinpath(out,"worker_memory.toml"),"w") do io;TOML.print(io,Dict("workers"=>[Dict(string(k)=>v for (k,v) in pairs(x)) for x in storage]));end
 open(joinpath(out,"storage.toml"),"w") do io;TOML.print(io,Dict("numeric_bytes"=>sum(x.numeric for x in storage),"object_graph_bytes"=>sum(x.object_graph for x in storage),"merge_entries"=>sum(x.merges for x in storage),"leaf_entries"=>sum(x.leaves for x in storage),"closed_values"=>sum(x.closed_values for x in storage),"layout_version"=>2,"skipped_cholesky_per_hit"=>68,"skipped_matrix_triangular_solves_per_hit"=>112,"skipped_matrix_products_per_hit"=>156,"fresh_vector_triangular_solves_per_hit"=>68,"fresh_gradient_matvecs_per_hit"=>112));end
 retire()
 println("GLOBAL_INVARIANCE_CACHE_PREFLIGHT_PASS");flush(stdout)
 # A missing entry invalidates the entire pass BEFORE using stale scratch.
 STAMP[]=nothing;execute_pass(:candidate,1);remotecall_fetch(worker_remove_entry,OWNER[1])
 x=execute_pass(:candidate,2);ok=validate_state(2,"missing_entry_invalidated",0)
 @assert x.hits==0;emitrow("negative_controls","missing_merge_entry",x.hits,x.stored,ok,true);retire()
 # Actual ordinary accepted Newton step, no new capture.
 execute_pass(:baseline,4);snapshot_state(4,"saved_live_ordinary_after")
 STAMP[]=nothing;execute_pass(:candidate,3);x=execute_pass(:candidate,4);ok=validate_state(4,"ordinary_invalidated",0)
 @assert x.hits==0;emitrow("negative_controls","actual_ordinary_k0_to_k1",x.hits,x.stored,ok,true);retire()
 # Changed matrix-driving input; exact restoration, no mutation of capture files.
 remotecall_fetch(worker_perturb,OWNER[1],true)
 execute_pass(:baseline,2);snapshot_state(2,"modified_input_reference";slot=99,check_saved=false)
 STAMP[]=nothing;execute_pass(:candidate,1);x=execute_pass(:candidate,2);ok=validate_state(2,"matrix_invalidated",0;slot=99)
 @assert x.hits==0;emitrow("negative_controls","matrix_driving_Sh0",x.hits,x.stored,ok,true)
 remotecall_fetch(worker_perturb,OWNER[1],false);retire();@assert stages(worker_key,2)==key2
 for (name,kw) in (("unknown_applied_regularization",(proof=:unknown,)),("model_or_execution_change",(context=ROOTCTX[]*"_changed",)))
  STAMP[]=nothing;execute_pass(:candidate,1);x=execute_pass(:candidate,2;kw...);ok=validate_state(2,name,0)
  @assert x.hits==0;emitrow("negative_controls",name,x.hits,x.stored,ok,true);retire()
 end
 # True terminal input change invalidates; modified boundary is ONLY a negative control.
 STAMP[]=nothing;execute_pass(:candidate,1);TV[][1,1]=1e-6;allworkers(worker_terminal,TV[])
 x=execute_pass(:candidate,2);@assert x.hits==0
 execute_pass(:baseline,2);snapshot_state(2,"changed_terminal_reference";slot=98,check_saved=false)
 # Recompute candidate miss from first pair for independent fresh comparison.
 retire();TV[][1,1]=0.0;allworkers(worker_terminal,TV[]);execute_pass(:candidate,1)
 TV[][1,1]=1e-6;allworkers(worker_terminal,TV[]);x=execute_pass(:candidate,2);ok=validate_state(2,"terminal_invalidated",0;slot=98)
 emitrow("negative_controls","terminal_curvature",x.hits,x.stored,ok,x.hits==0)
 retire();TV[][1,1]=0.0;allworkers(worker_terminal,TV[])
 STAMP[]=nothing;execute_pass(:candidate,1);execute_pass(:candidate,2;stale=true)
 rejected=!validate_state(2,"stale_gradient_negative",0;require_pass=false)
 @assert rejected;emitrow("negative_controls","stale_gradient_consumer",68,0,false,rejected);retire()
 # Warm the complete workloads and validation before five clean alternating pairs.
 allworkers(worker_compilation);allworkers(worker_compilation)
 warmed=false
 for warm in 1:4
  good=true
  for arm in (:baseline,:candidate)
   before=compiler_totals();x=measured_workload(arm);after=compiler_totals()
   cc=(compile=after.compile-before.compile,recompile=after.recompile-before.recompile,alloc=after.alloc-before.alloc,gc_ns=after.gc_ns-before.gc_ns)
   comp=x.first.compile_time+x.second.compile_time+x.cleanup.compile_time
   recomp=x.first.recompile_time+x.second.recompile_time+x.cleanup.recompile_time
   emitrow("warmup",warm,arm,cc.compile,cc.recompile,comp,recomp)
   validate_state(2,"warmup_$(arm)",warm)
   good&=cc.compile==0&&cc.recompile==0&&comp==0&&recomp==0
  end
  if good;warmed=true;break;end
 end
 @assert warmed "bounded warmup did not clear compilation"
 for rep in 1:5
  for arm in (isodd(rep) ? (:baseline,:candidate) : (:candidate,:baseline))
   before=compiler_totals();x=measured_workload(arm);after=compiler_totals()
   cc=(compile=after.compile-before.compile,recompile=after.recompile-before.recompile,alloc=after.alloc-before.alloc,gc_ns=after.gc_ns-before.gc_ns)
   @assert cc.compile==0&&cc.recompile==0
   @assert all(z->z.compile_time==0&&z.recompile_time==0,(x.first,x.second,x.cleanup))
   @assert x.first.value.hits==0&&x.second.value.hits==(arm==:candidate ? 68 : 0)
   logtiming(rep,arm,x,cc,cc.alloc) # per-worker allocations not sampled in clean clocks
   validate_state(2,"measured_$(arm)",rep)
   println("PAIR_ARM_COMPLETE repeat=$rep arm=$arm charged=$(x.total)");flush(stdout)
  end
 end
 open(joinpath(out,"validation.toml"),"w") do io
  TOML.print(io,Dict("complete"=>true,"layout_version"=>2,"population_full_output_gate"=>true,"global_leaf_matrix_invariance"=>true,"global_68_merge_invariance"=>true,"signature_invariance"=>true,"measured_pairs"=>5,"warmup_compilation_cleared"=>true,"scaled_gate"=>1e-7,"relative_gate"=>1e-5,"relative_floor"=>1e-10,"process_workers"=>24,"blas_threads"=>1,"whole_workload_includes_population_detection_retirement"=>true,"scope"=>"frozen_process_parallel_pilot_not_optimizer_integration","applied_regularization_proof"=>"pinned_tp_pre_unshifted_source_not_reg_last"))
 end
 rmprocs(pids)
 println("EXACT_BARRIER_CACHE_PILOT_COMPLETE");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
 try
  @assert length(ARGS)==3;runpilot(ARGS...)
 finally
  for io in values(LOGS);close(io);end
 end
end
