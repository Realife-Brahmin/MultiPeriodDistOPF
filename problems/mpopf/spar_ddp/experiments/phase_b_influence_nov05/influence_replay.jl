# Frozen full-horizon sensitivity audit; never runs an approximate optimizer.
using Serialization,LinearAlgebra,SparseArrays,SHA,TOML,Statistics
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"lib","replay.jl"))
include(joinpath(@__DIR__,"influence_core.jl"))
const HANDLES=Dict{String,IO}()
function table(name,header)
 io=open(joinpath(OUTPUT,name*".csv"),"w");println(io,header);HANDLES[name]=io;io
end
row(name,args...)=emit(HANDLES[name],args...)
function cost(label,context,f)
 z=@timed f();row("costs",context,label,z.time,z.bytes,z.gctime);z.value
end
# Julia do-block convention.
cost(f::Function,label,context)=cost(label,context,f)
function err(name,context,field,a,b;stage=0,packet=nothing)
 m=metric(a,b)
 row(name,context,stage,field,m.absolute,m.reference,m.scaled,m.relative,m.maxabs,m.reference_near_zero,passes_metric(m))
 D=abs.(a-b);order=partialsortperm(vec(D),1:min(5,length(D));rev=true)
 for idx in order
  coord=Tuple(CartesianIndices(D)[idx]);bus="unmapped"
  if packet!==nothing&&length(coord)==1
   g=field in (:dphi,) ? packet.nu+coord[1] : coord[1]
   for (j,own) in enumerate(packet.S.lay.own)
    g in own&&(bus=string(DATA[:Nset][j]);break)
   end
   for j in 1:packet.nx
    g==packet.S.active[j]&&(bus=string(DATA[:Bset][j]);break)
   end
  end
  row("outliers",name,context,stage,field,join(coord,';'),bus,D[idx])
 end
 m
end
function scan(leaves;context="baseline",previous=nothing,mode=:current,trace=false)
 n=length(first(leaves).ma);N=length(leaves)
 xs=Any[deepcopy(e) for e in leaves];push!(xs,PintVal(zeros(n,n),zeros(n)))
 refs=trace ? Any[deepcopy(e) for e in leaves] : Any[]
 trace&&push!(refs,PintVal(zeros(n,n),zeros(n)))
 w=PairWork(n);rw=PairWork(n);points=Dict{Tuple{Int,Int,Int},Any}();affected=falses(N+1)
 spans=[(i,i) for i in 1:N];push!(spans,(N+1,N))
 for op in scan_ops(N)
  depth,t,mid=op;lo=spans[t][1];hi=max(spans[t][2],spans[mid][2])
  if op in LOCATIONS
   @assert xs[t] isa PintElem&&xs[mid] isa PintElem
   points[op]=merge_inputs(xs[t],xs[mid])
   row("locations",context,depth,t,mid,spans[t][1],spans[t][2],spans[mid][1],spans[mid][2],lo,hi)
  end
  hit=previous!==nothing&&op==(2,5,7);wasleft=affected[t];wasright=affected[mid]
  changed=hit||wasleft||wasright
  if hit
   z=retained_inputs(merge_inputs(xs[t],xs[mid]),previous,mode)
   xs[t]=cost("injected_combine",context) do;input_combine!(xs[t],xs[mid],w,z);end
  else
   if trace&&changed
    H=xs[t].Mbb+(xs[mid] isa PintElem ? xs[mid].Maa : xs[mid].V)
    ok=cost("fallback_detection",context) do;isposdef(Symmetric(H,:U));end
    !ok&&row("failures",context,"dense_fallback",depth,t,mid,"original_production_fallback_retained")
   end
   xs[t]=cost("production_combine",context) do;_pp_combine!(xs[t],xs[mid],w);end
  end
  affected[t]=changed;spans[t]=(lo,hi)
  if trace
   refs[t]=_pp_combine!(refs[t],refs[mid],rw)
   row("dependencies",context,depth,t,mid,hit,wasleft,wasright,changed,lo,hi)
   for f in fields_of(xs[t])
    err("merge_errors",context*"_d$(depth)_l$(t)_r$(mid)",f,getproperty(xs[t],f),getproperty(refs[t],f);stage=t)
   end
  end
 end
 @assert all(x->x isa PintVal,xs)
 (values=xs,points=points,affected=affected)
end
function packet(dir,t);deserialize(joinpath(dir,"stage_$(t).jls"));end
function extras(n)
 signs=ones(n)
 for (j,g) in enumerate(GROUPS);signs[g].=isodd(j) ? 1.0 : -1.0;end
 [(kind="cosine_small",dx=1e-3.*cos.(collect(1:n))),
  (kind="cosine_large",dx=cos.(collect(1:n))),
  (kind="topology_small",dx=1e-3.*signs),(kind="topology_large",dx=signs)]
end
function independent_next(e,v,dx)
 -(Symmetric(e.Mbb+v.V,:U)\(e.Mab'*dx+e.mb+v.v))
end
function check_baseline(state,dir,payload,values,points)
 good=true;negative_rejected=false
 # No-change substitution kernel must reproduce the original selected operation.
 # Recreate selected target input operands by a fresh production scan.
 xs=Any[deepcopy(e) for e in payload.elems];n=length(first(payload.elems).ma)
 push!(xs,PintVal(zeros(n,n),zeros(n)));w=PairWork(n)
 for op in scan_ops(length(payload.elems))
  _,t,mid=op
  if op==(2,5,7)
   x=deepcopy(xs[t]);input_combine!(x,xs[mid],PairWork(n),points[op])
   exact=_pp_combine!(deepcopy(xs[t]),xs[mid],PairWork(n))
   for f in fields_of(x)
    a=getproperty(x,f);b=getproperty(exact,f)
    good&=a==b;err("baseline",state*"_current_kernel",f,a,b)
   end
  end
  xs[t]=_pp_combine!(xs[t],xs[mid],w)
 end
 for t in 1:24
  p=packet(dir,t);@assert p.stage==t&&p.nx==n
  rule=cost("reference_consumer",state*"_stage$(t)") do;policy(p,values[t+1],payload.mu);end
  for pr in p.probes
   d=directions(rule,pr.dx)
   for f in (:du,:dphi,:dzl,:dzu)
    m=err("baseline",state*"_"*pr.kind,f,getproperty(d,f),getproperty(pr,f);stage=t,packet=p);good&=passes_metric(m)
   end
   live_next=p.S.fx*pr.dx+p.S.fu*d.du
   battery_next=cost("independent_battery_consumer",state*"_stage$(t)") do;independent_next(payload.elems[t],values[t+1],pr.dx);end
   err("independent",state*"_"*pr.kind,:next_dx,battery_next,live_next;stage=t)
   pr.kind=="actual_soc"&&(good&=passes_metric(err("baseline",state*"_actual_forward",:next_dx,live_next,p.actual_next_dx;stage=t)))
  end
  if t==1
   future=values[2];bad=PintVal(copy(future.V),future.v.+1e6/sqrt(n))
   badrule=policy(p,bad,payload.mu);d=directions(badrule,zeros(n));r=directions(rule,zeros(n))
   m=err("negative",state*"_composed_gradient",:du,d.du,r.du;stage=t,packet=p)
   negbattery=independent_next(payload.elems[t],bad,zeros(n))
   refbattery=independent_next(payload.elems[t],future,zeros(n))
   mb=err("negative",state*"_independent_composed_gradient",:next_dx,negbattery,refbattery;stage=t)
   negative_rejected=!passes_metric(m)&&!passes_metric(mb)
  end
 end
 row("baseline_status",state,good,negative_rejected)
 good&&negative_rejected
end
function distribution(context,field,A,B)
 ts=time_ns();D=A-B;n=size(A,1);m=metric(A,B);sq=sum(abs2,D)
 diagenergy=ndims(D)==2 ? sum(abs2,diag(D)) : sq
 row("changes",context,field,m.absolute,m.reference,m.scaled,m.relative,m.maxabs,m.reference_near_zero,sq==0 ? 0.0 : diagenergy/sq,(time_ns()-ts)/1e9)
 row("costs",context,"change_detection",(time_ns()-ts)/1e9,0,0)
 for (axis,vals) in (("coordinates",vec(abs.(D))),("rows",ndims(D)==2 ? vec(sqrt.(sum(abs2,D;dims=2))) : abs.(D)),("columns",ndims(D)==2 ? vec(sqrt.(sum(abs2,D;dims=1))) : abs.(D)))
  for q in (0.0,0.5,0.9,0.99,1.0);row("quantiles",context,field,axis,q,quantile(vals,q));end
  if axis!="coordinates"
   for (i,v) in enumerate(vals);row("row_column",context,field,axis,i,DATA[:Bset][i],v);end
  end
 end
 err("input_errors",context,field,A,B)
 for (i,g) in enumerate(GROUPS)
  if ndims(D)==2
   for (j,h) in enumerate(GROUPS)
    v=@view D[g,h];r=@view B[g,h];row("blocks",context,field,i,j,length(g),length(h),norm(v),norm(v)/max(1e-10,norm(r)),maximum(abs,v))
   end
  else
   row("blocks",context,field,i,0,length(g),1,norm(D[g]),norm(D[g])/max(1e-10,norm(B[g])),maximum(abs,D[g]))
  end
 end
 if ndims(D)==2
  ss=cost("delta_svd",context*"_"*string(field)) do;svdvals(D);end
  fullss=cost("current_svd",context*"_"*string(field)) do;svdvals(A);end
  for (kind,s) in (("delta",ss),("current",fullss))
   total=sum(abs2,s)
   for (i,v) in enumerate(s);row("spectra",context,field,kind,i,v);end
   for rel in (1e-3,1e-6,1e-9,1e-12)
    rank=count(>(rel*first(s)),s)
    row("ranks",context,field,kind,rel,rank,total==0 ? 0.0 : sqrt(sum(abs2,s[rank+1:end])/total))
   end
  end
 end
end
function condition(context,H)
 ev=cost("actual_upper_eigen",context) do;eigvals(Symmetric(H,:U));end
 threshold=length(ev)*eps()*maximum(abs,ev);lo=minimum(abs,ev);hi=maximum(abs,ev)
 row("condition",context,minimum(ev),maximum(ev),count(<(-threshold),ev),count(x->abs(x)<=threshold,ev),count(>(threshold),ev),lo==0 ? Inf : hi/lo,norm(H-H')/max(1.0,norm(H)),"original_coordinate_upper")
end
function influence_pair(name,before,after)
 a=STATES[before];b=STATES[after];prev=a.scan.points[(2,5,7)]
 for loc in LOCATIONS
  for f in (:H,:BL,:BR,:h)
   distribution(name*"_d$(loc[1])_l$(loc[2])_r$(loc[3])",f,getproperty(b.scan.points[loc],f),getproperty(a.scan.points[loc],f))
  end
  condition(name*"_before_d$(loc[1])",a.scan.points[loc].H);condition(name*"_after_d$(loc[1])",b.scan.points[loc].H)
 end
 candidates=Dict{Symbol,Any}()
 for mode in INTERVENTIONS
  context=name*"_"*string(mode)
  try
   s=scan(b.payload.elems;context=context,previous=prev,mode=mode,trace=true)
   @assert findall(s.affected)==collect(1:5)
   candidates[mode]=s.values
   for t in 1:25,f in (:V,:v)
    err("suffix_errors",context,f,getproperty(s.values[t],f),getproperty(b.scan.values[t],f);stage=t)
   end
  catch e
   row("failures",context,string(typeof(e)),2,5,7,replace(sprint(showerror,e),','=>';','\n'=>' '))
   println("INTERVENTION_FAILURE $context $(typeof(e))");flush(stdout)
  end
 end
 n=length(first(b.payload.elems).ma);paths=Dict(mode=>zeros(n) for mode in keys(candidates));refdx=zeros(n)
 for t in 1:24
  p=packet(b.dir,t)
  ref=cost("reference_consumer",name*"_stage$(t)") do;policy(p,b.scan.values[t+1],b.payload.mu);end
  reftrajectory=directions(ref,refdx);refnext=p.S.fx*refdx+p.S.fu*reftrajectory.du
  probes=vcat(collect(p.probes),extras(n))
  for (mode,values) in candidates
   context=name*"_"*string(mode)
   try
    r=cost("modified_consumer",context*"_stage$(t)") do;policy(p,values[t+1],b.payload.mu);end
    for pr in probes
     d=directions(r,pr.dx);base=directions(ref,pr.dx)
     for f in (:du,:dphi,:dzl,:dzu)
      err("directions",context*"_"*pr.kind,f,getproperty(d,f),getproperty(base,f);stage=t,packet=p)
     end
     for field in (:V,:v)
      # Hessian action and complete gradient; preserve full norms separately.
      av=field==:V ? values[t+1].V*pr.dx : values[t+1].v
      bv=field==:V ? b.scan.values[t+1].V*pr.dx : b.scan.values[t+1].v
      err("actions",context*"_"*pr.kind,field,av,bv;stage=t)
     end
    end
    dx=paths[mode];d=directions(r,dx);nextdx=p.S.fx*dx+p.S.fu*d.du
    for f in (:du,:dphi,:dzl,:dzu)
     err("trajectory",context,f,getproperty(d,f),getproperty(reftrajectory,f);stage=t,packet=p)
    end
    err("trajectory",context,:dx,dx,refdx;stage=t)
    err("trajectory",context,:next_dx,nextdx,refnext;stage=t)
    independent=independent_next(b.payload.elems[t],values[t+1],dx)
    err("independent",context,:next_dx,independent,nextdx;stage=t)
    paths[mode]=nextdx
   catch e
    row("failures",context,string(typeof(e)),0,t,t+1,replace(sprint(showerror,e),','=>';','\n'=>' '))
    # Failed path stays failed, never silently replaced with a reference path.
    paths[mode].=NaN
   end
  end
  refdx=refnext
 end
 println("INFLUENCE_PAIR_COMPLETE pair=$name candidates=$(length(candidates))");flush(stdout)
end
function main(capture,input,output)
 global OUTPUT=output;@assert !ispath(output);mkpath(output)
 global DATA=deserialize(input)
 global GROUPS=coupling_groups(topology_hierarchy(DATA;maxleaf=32))
 global STATES=Dict{String,Any}()
 headers=Dict("costs"=>"context,operation,seconds,allocation_bytes,gc_s",
 "locations"=>"context,depth,left,right,left_lo,left_hi,right_lo,right_hi,out_lo,out_hi",
 "dependencies"=>"context,depth,left,right,injected,left_dependency,right_dependency,affected,out_lo,out_hi",
 "outliers"=>"table,context,stage,field,coordinate,bus_or_unmapped,absolute_error",
 "baseline_status"=>"state,production_pass,negative_rejected",
 "failures"=>"context,outcome,depth,left,right,detail",
 "changes"=>"context,field,absolute_fro,previous_fro,unit_scaled,relative_floor1e10,maxabs,previous_near_zero,delta_diagonal_energy_fraction,change_detection_s",
 "quantiles"=>"context,field,distribution,quantile,value",
 "row_column"=>"context,field,axis,index,battery_bus,delta_norm",
 "blocks"=>"context,field,row_group,col_group,rows,cols,delta_fro,relative_floor1e10,maxabs",
 "spectra"=>"context,field,kind,index,singular_value",
 "ranks"=>"context,field,kind,relative_cutoff,rank,relative_tail_fro",
 "condition"=>"context,min_eig,max_eig,negative,near_zero,positive,abs_condition,raw_asymmetry,completion")
 for name in ("baseline","independent","negative","merge_errors","input_errors","suffix_errors","directions","actions","trajectory")
  headers[name]="context,stage,field,absolute_fro,reference_fro,unit_scaled,relative_floor1e10,maxabs,reference_near_zero,unchanged_error_gate"
 end
 for (name,header) in headers;table(name,header);end
 open(joinpath(output,"topology_groups.toml"),"w") do io;TOML.print(io,Dict("groups"=>GROUPS,"buses"=>DATA[:Bset]));end
 allgood=true
 for dir in sort(filter(isdir,readdir(capture;join=true)))
  meta=TOML.parsefile(joinpath(dir,"metadata.toml"));state=only(meta["tags"])
  @assert bytes2hex(sha256(read(joinpath(dir,"leaves.jls"))))==meta["leaves_sha256"]
  for r in meta["stage_files"]
   file=joinpath(dir,r["file"]);@assert bytes2hex(sha256(read(file)))==r["sha256"]&&filesize(file)==r["bytes"]
  end
  payload=deserialize(joinpath(dir,"leaves.jls"));@assert length(payload.elems)==24
  s=scan(payload.elems;context=state)
  good=check_baseline(state,dir,payload,s.values,s.points);allgood&=good
  STATES[state]=(dir=dir,payload=payload,scan=s,meta=meta)
  println("BASELINE_STATE_COMPLETE state=$state pass=$good");flush(stdout);GC.gc()
 end
 @assert length(STATES)==8
 open(joinpath(output,"baseline_gate.toml"),"w") do io;TOML.print(io,Dict("pass"=>allgood,"states"=>8,"stages_per_state"=>24,"unit_scaled_gate"=>1e-7,"relative_gate"=>1e-5,"relative_floor"=>1e-10));end
 if !allgood
  println("BASELINE_REFERENCE_GATE_FAILED no_interventions_run");return
 end
 for name in ("early","middle","late","barrier")
  influence_pair(name,name*"_before",name*"_after");GC.gc()
 end
 open(joinpath(output,"audit_metadata.toml"),"w") do io
  TOML.print(io,Dict("completed"=>true,"baseline_gate"=>true,"interventions_planned"=>24,"timing_is_benchmark"=>false,"julia"=>string(VERSION),"blas"=>string(BLAS.get_config()),"blas_threads"=>BLAS.get_num_threads(),"input_sha256"=>bytes2hex(sha256(read(input))),"capture_root"=>capture,"maxrss_bytes"=>Sys.maxrss(),"interpretation"=>"mixed_input_sensitivity_not_deployable_reuse"))
 end
 println("INFLUENCE_AUDIT_COMPLETE");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
 try
  @assert length(ARGS)==3
  main(ARGS...)
 finally
  for io in values(HANDLES);close(io);end
 end
end
