using Distributed,SharedArrays,Serialization,LinearAlgebra,SHA,TOML,SparseArrays
include(joinpath(@__DIR__,"lib","replay.jl"))
include(joinpath(@__DIR__,"cache_core.jl"))
const WPACK=Dict{Int,Any}();const WREF=Dict{Int,Any}()
const ENTRIES=Dict{Tuple{Int,Int,Int},MergeMatrixCache}()
const DIAG=Dict{Any,Any}()
const WSHARED=Ref{Any}(nothing);const WRULE=Ref{Any}(nothing);const WDIR=Ref{Any}(nothing)
const WSTAGE=Ref(0);const WWHICH=Ref(0);const WSIGNS=Ref(Float64[])
function worker_setup(shared,n,t,dirs,signs)
 BLAS.set_num_threads(1);ENV["FILTERDDP_PAIRED_INPLACE"]="1"
 Base.cumulative_compile_timing(true)
 WSHARED[]=shared;WSTAGE[]=t;WSIGNS[]=signs
 empty!(_PP_SLOTS);empty!(_PP_WORK)
 for a in shared
  push!(_PP_SLOTS,shared_slot_view(sdata(a),n))
 end
 push!(_PP_WORK,PairWork(n));_PP_END[]=PintVal(zeros(n,n),zeros(n))
 for (i,dir) in enumerate(dirs)
  meta=TOML.parsefile(joinpath(dir,"metadata.toml"));r=only(filter(r->r["file"]=="stage_$(t).jls",meta["stage_files"]))
  file=joinpath(dir,r["file"]);@assert bytes2hex(sha256(read(file)))==r["sha256"]
  WPACK[i]=deserialize(file)
 end
 (pid=getpid(),blas=BLAS.get_num_threads(),backend=string(BLAS.get_config()),threads=Threads.nthreads(),rss=Sys.maxrss())
end
worker_key(i)=stage_signature(WPACK[i])
worker_key_info(i)=(key=worker_key(i),entries=length(ENTRIES))
function worker_reset_cache()
 empty!(ENTRIES);nothing
end
function worker_leaf(i,hit,populate)
 t=WSTAGE[];p=WPACK[i];WWHICH[]=i
 fill!(_PP_END[].v,0.0) # current actual terminal gradient, never cached
 if hit
  leaf_vectors!(_PP_SLOTS[t],p.S)
 else
  @assert _pp_elem(p.S,t,p.nx)!==nothing
 end
 nothing
end
function worker_leaf_check(expected)
 e=_PP_SLOTS[WSTAGE[]]
 [(string(f),errors(getproperty(e,f),getproperty(expected,f))) for f in (:Maa,:Mab,:Mbb,:ma,:mb)]
end
function worker_merge(depth,t,mid,yval,last,mode,diagnostic)
 x=_PP_SLOTS[t];y=last ? _PP_END[] : yval ? PintVal(_PP_SLOTS[mid].Maa,_PP_SLOTS[mid].ma) : _PP_SLOTS[mid]
 key=(depth,t,mid);val=y isa PintVal;rows=Any[]
 if diagnostic!=:none
  inputs=(H=x.Mbb+(val ? y.V : y.Maa),BL=copy(x.Mab),BR=val ? nothing : copy(y.Mab),leftA=copy(x.Maa),rightC=val ? nothing : copy(y.Mbb))
  if diagnostic==:save;DIAG[(key,:in)]=inputs
  else
   old=DIAG[(key,:in)]
   for f in keys(inputs)
    aa=getproperty(inputs,f);bb=getproperty(old,f)
    isnothing(aa)||push!(rows,(string(f),errors(aa,bb)))
   end
  end
 end
 hit=false;stored=false;fallback=false;miss=false
 if mode==:hit&&haskey(ENTRIES,key)
  hit_merge!(x,y,_PP_WORK[1],ENTRIES[key]);hit=true
 elseif mode==:populate
  _,entry=populate_merge!(x,y,_PP_WORK[1]);fallback=entry===nothing
  if entry!==nothing;ENTRIES[key]=entry;stored=true;end
 elseif mode==:hit
  # Do NOT use an unmaterialized intermediate as a fresh operand.
  miss=true
 else
  _pp_combine!(x,y,_PP_WORK[1])
 end
 if diagnostic!=:none
  out=val ? (V=copy(x.Maa),) : (Maa=copy(x.Maa),Mab=copy(x.Mab),Mbb=copy(x.Mbb))
  if diagnostic==:save;DIAG[(key,:out)]=out
  else
   old=DIAG[(key,:out)]
   for f in keys(out);push!(rows,("output_"*string(f),errors(getproperty(out,f),getproperty(old,f))));end
  end
 end
 (hit=hit,stored=stored,fallback=fallback,miss=miss,rows=rows)
end
function worker_policy(i,mu)
 t=WSTAGE[];p=WPACK[i];v=t==24 ? _PP_END[] : PintVal(_PP_SLOTS[t+1].Maa,_PP_SLOTS[t+1].ma)
 WRULE[]=policy(p,v,mu);nothing
end
function worker_forward(dx)
 p=WPACK[WWHICH[]];d=directions(WRULE[],dx)
 WDIR[]=(dx=copy(dx),du=copy(d.du),dphi=copy(d.dphi),dzl=copy(d.dzl),dzu=copy(d.dzu))
 p.S.fx*dx+p.S.fu*d.du
end
function worker_probes()
 p=WPACK[WWHICH[]];n=p.nx
 vcat(collect(p.probes),[(kind="cosine_small",dx=1e-3.*cos.(collect(1:n))),(kind="cosine_large",dx=cos.(collect(1:n))),(kind="topology_small",dx=1e-3.*WSIGNS[]),(kind="topology_large",dx=WSIGNS[])])
end
function worker_snapshot(i,slot=i)
 t=WSTAGE[];p=WPACK[i]
 v=t==24 ? _PP_END[] : PintVal(_PP_SLOTS[t+1].Maa,_PP_SLOTS[t+1].ma)
 probes=Any[]
 for pr in worker_probes()
  d=directions(WRULE[],pr.dx)
  push!(probes,(kind=pr.kind,dirs=deepcopy(d),action=copy(v.V*pr.dx),gradient=copy(v.v)))
 end
 WREF[slot]=(V=copy(_PP_SLOTS[t].Maa),v=copy(_PP_SLOTS[t].ma),trajectory=deepcopy(WDIR[]),probes=probes)
 # Also verify exact original saved live full controls/bounds, not just new replay.
 rows=Any[]
 for pr in p.probes
  d=directions(WRULE[],pr.dx)
  for f in (:du,:dphi,:dzl,:dzu);push!(rows,(pr.kind,string(f),errors(getproperty(d,f),getproperty(pr,f))));end
 end
 rows
end
function worker_validate(i,slot=i)
 t=WSTAGE[];p=WPACK[i];r=WREF[slot];rows=Any[];outliers=Any[]
 function add(probe,field,aa,bb)
  m=errors(aa,bb);push!(rows,(probe,string(field),m))
  delta=abs.(aa-bb)
  for idx in partialsortperm(vec(delta),1:min(3,length(delta));rev=true)
   push!(outliers,(probe,string(field),idx,delta[idx]))
  end
 end
 add("suffix",:V,_PP_SLOTS[t].Maa,r.V);add("suffix",:v,_PP_SLOTS[t].ma,r.v)
 for f in (:dx,:du,:dphi,:dzl,:dzu);add("trajectory",f,getproperty(WDIR[],f),getproperty(r.trajectory,f));end
 v=t==24 ? _PP_END[] : PintVal(_PP_SLOTS[t+1].Maa,_PP_SLOTS[t+1].ma)
 for (pr,ref) in zip(worker_probes(),r.probes)
  d=directions(WRULE[],pr.dx)
  for f in (:du,:dphi,:dzl,:dzu);add(pr.kind,f,getproperty(d,f),getproperty(ref.dirs,f));end
  add(pr.kind,:V_action,v.V*pr.dx,ref.action);add(pr.kind,:gradient,v.v,ref.gradient)
 end
 (rows=rows,outliers=outliers)
end
function worker_storage()
 (numeric=sum(cache_numeric_bytes,values(ENTRIES);init=0),object_graph=Base.summarysize(ENTRIES),merges=length(ENTRIES),leaves=0,closed_values=count(x->x.V!==nothing,values(ENTRIES)),rss=Sys.maxrss())
end
function worker_retire()
 empty!(ENTRIES);GC.gc(false);nothing
end
function worker_compilation()
 c,r=Base.cumulative_compile_time_ns();g=Base.gc_num()
 (compile=c,recompile=r,alloc=Base.gc_bytes(),gc_ns=g.total_time)
end
worker_clear_diagnostics()=(empty!(DIAG);GC.gc(false);nothing)
const OLD_COEFFICIENT=Ref(0.0)
function worker_perturb(on)
 p=WPACK[2]
 if on;OLD_COEFFICIENT[]=p.S.b0.Sh0[1,1];p.S.b0.Sh0[1,1]+=1e-3
 else;p.S.b0.Sh0[1,1]=OLD_COEFFICIENT[];end
 nothing
end
function worker_stale_gradients()
 t=WSTAGE[];p=WPACK[1];e=_PP_SLOTS[t]
 leaf_vectors!(e,p.S);nothing
end

function worker_leaf_invariance(save)
 e=_PP_SLOTS[WSTAGE[]]
 if save
  DIAG[:leaf]=(Maa=copy(e.Maa),Mab=copy(e.Mab),Mbb=copy(e.Mbb));return Any[]
 end
 [(string(f),errors(getproperty(e,f),getproperty(DIAG[:leaf],f))) for f in (:Maa,:Mab,:Mbb)]
end
function worker_terminal(value)
 copyto!(_PP_END[].V,value);fill!(_PP_END[].v,0.0);nothing
end

function worker_pin(cpus)
 mask=zeros(UInt64,cld(maximum(cpus)+1,64))
 for cpu in cpus;mask[fld(cpu,64)+1]|=UInt64(1)<<mod(cpu,64);end
 for tid in readdir("/proc/self/task")
  rc=ccall(:sched_setaffinity,Cint,(Cint,Csize_t,Ptr{UInt64}),parse(Int,tid),sizeof(mask),mask)
  @assert rc==0
 end
 nothing
end
function worker_affinity()
 masks=[strip(only(filter(s->startswith(s,"Cpus_allowed_list:"),split(read("/proc/self/task/$tid/status",String),'\n')))[19:end]) for tid in readdir("/proc/self/task")]
 @assert length(unique(masks))==1
 (pid=getpid(),allowed=first(masks),blas=BLAS.get_num_threads())
end

worker_remove_entry()=(isempty(ENTRIES)||delete!(ENTRIES,first(keys(ENTRIES)));nothing)
