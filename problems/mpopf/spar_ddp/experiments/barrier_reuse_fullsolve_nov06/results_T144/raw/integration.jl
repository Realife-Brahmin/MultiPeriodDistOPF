# Loaded only in the immutable opt-in experiment snapshot; production defaults unchanged.
include(joinpath(@__DIR__,"br_cache_kernel.jl"))
using SHA,Serialization
const BR_LOCAL_INVALID=Ref(false)
const BR_ON=get(ENV,"STUDY_BARRIER_REUSE","0")=="1"
const BR_DEBUG=get(ENV,"STUDY_BARRIER_DEBUG","0")=="1"
const BR_GATE=BRGate()
const BR_ENTRIES=Dict{Tuple{Int,Int,Int},BRMatrixEntry}()
const BR_MODE=Ref(:fresh)
const BR_CONTEXT=Ref{Any}(nothing)
const BR_EXPECTED=Ref(0)
const BR_VALID=Ref(false)
const BR_FAILED=Ref(false)
const BR_NEGATIVE_DONE=Set{Int}()
const BR_OLD_VECTORS=Dict{Any,Any}()
const BR_STAGE_KEYS=Dict{Int,String}()
const BR_REF_RULES=Dict{Int,Any}()
const BR_SHADOW=Ref{Any}(nothing)
const BR_DEBUG_ROWS=Ref(0)
# Explicit invalidation API for future mutation sites. The pinned driver has no
# model/data/terminal writes between install and exit. Direct external mutation
# is unsupported; debug mode fingerprints every actual stage to audit that proof.
barrier_cache_invalidate!(reason::Symbol=:unknown)=(BR_LOCAL_INVALID[]=true;br_event!(BR_GATE,reason))
function br_context()
 Tuple(sort!([(k,v) for (k,v) in ENV if startswith(k,"FILTERDDP_") || startswith(k,"REDUCED_") || startswith(k,"TERMINAL_SOC_") || k in
  ("STUDY_SCAN","STUDY_KERNEL","STUDY_WORKERS","STUDY_BLAS_THREADS","REDUCED_CB","TERMINAL_SOC_GAMMA","REDUCED_PROFILE")]))
end
function br_supported()
 get(ENV,"STUDY_SCAN","")=="grouped" && get(ENV,"STUDY_KERNEL","")=="right" &&
 get(ENV,"FILTERDDP_PROCS_BALANCED","")=="1" &&
 isempty(get(ENV,"FILTERDDP_HYBRID_MU","")) && isempty(get(ENV,"FILTERDDP_VALUE_GROUPING","")) &&
 isempty(get(ENV,"FILTERDDP_FROZEN_CAPTURE_DIR","")) &&
 get(ENV,"FILTERDDP_LEAN_EXACT","")=="1" && get(ENV,"FILTERDDP_NETWORK_BACKEND","")=="tree"
end
function br_begin(d)
 br_request!(BR_GATE,d.k,d.j,Float64(d.μ),br_context(),BR_ON&&br_supported())
end
function br_release!()
 empty!(BR_ENTRIES);empty!(BR_STAGE_KEYS);empty!(BR_REF_RULES)
 empty!(BR_NEGATIVE_DONE);empty!(BR_OLD_VECTORS)
 BR_EXPECTED[]=0;BR_VALID[]=false;BR_CONTEXT[]=nothing;BR_MODE[]=:fresh
 nothing
end
function br_worker_start(request,context,expected,reset)
 reset && br_release!() # Drop warm-up cache payload BEFORE the worker's full GC.
 ready=!BR_LOCAL_INVALID[] && BR_VALID[] && length(BR_ENTRIES)==expected && BR_EXPECTED[]==expected && BR_CONTEXT[]==context
 BR_MODE[]=request&&ready ? :hit : BR_ON ? :populate : :fresh
 BR_LOCAL_INVALID[]=false
 BR_CONTEXT[]=context;BR_EXPECTED[]=expected;BR_FAILED[]=false
 empty!(BR_REF_RULES)
 ready
end
function br_force_population()
 BR_MODE[]=BR_ON ? :populate : :fresh
 nothing
end
function br_digest(x)
 io=IOBuffer();serialize(io,x);bytes2hex(sha256(take!(io)))
end
function br_stage_key(S,W)
 fields=(:x,:u,Symbol("\u03d5"),:zl,:zu,:lxx,:lux,:cxx,:cx,:cu,:fx,:fu,:fxx,:active,:fu_active,:lux_in_B,:K,:lay,:stat,:fac,:Bd,:Up,:b0,:nzc,:cxE,:inv_ul,:inv_uu,Symbol("\u03a3_L"),Symbol("\u03a3_U"))
 br_digest(([(f,getproperty(S,f)) for f in fields],W.ocp.control_limits,W.nx,W.nu,W.nc))
end
function br_check(a,b,label;exact=false)
 delta=norm(a-b);ref=norm(b);scaled=delta/max(1.0,ref);rel=delta/max(1e-10,ref)
 bits=size(a)==size(b)&&reinterpret(UInt64,vec(a))==reinterpret(UInt64,vec(b))
 @assert isfinite(delta)&&scaled<1e-7&&rel<1e-5&&(!exact||bits) "BR_DEBUG_FAIL $label scaled=$scaled relative=$rel bits=$bits"
 BR_DEBUG_ROWS[]+=1
 println("BR_DEBUG_CHECK label=$label absolute=$delta scaled=$scaled relative=$rel maxabs=$(maximum(abs,a-b;init=0.0)) exact=$bits")
 nothing
end
function br_leaf!(S,t,n)
 @assert !(BR_MODE[]===:hit && BR_LOCAL_INVALID[]) "Invalidation during hit pass"
 if BR_DEBUG
  key=br_stage_key(S,_PW[])
  if !(t in BR_NEGATIVE_DONE)
   for f in (:x,:u,Symbol("\u03d5"),:zl,:zu,:lxx,:K)
    changed=deepcopy(S);a=getproperty(changed,f)
    if a isa SparseMatrixCSC
     a[1,1]+=0.01
    else
     a[1]+=0.01
    end
    @assert br_stage_key(changed,_PW[])!=key
    println("BR_NEGATIVE stage=$t kind=$f rejected=true")
   end
   push!(BR_NEGATIVE_DONE,t)
  end
  if BR_MODE[]===:hit;@assert get(BR_STAGE_KEYS,t,"")==key "Stage curvature/state fingerprint changed at eligible hit";end
  BR_STAGE_KEYS[t]=key
 end
 if BR_MODE[]===:hit && !BR_DEBUG
  br_leaf_vectors!(_PP_SLOTS[t],S)
 else
  isnothing(_pp_elem(S,t,n)) && return nothing
  if BR_DEBUG
   vectors=deepcopy(_PP_SLOTS[t]);br_leaf_vectors!(vectors,S)
   br_check(vectors.ma,_PP_SLOTS[t].ma,"leaf_$(t)_ma";exact=true)
   br_check(vectors.mb,_PP_SLOTS[t].mb,"leaf_$(t)_mb";exact=true)
  end
 end
 true
end
function br_shadow_worker(N,n)
 BR_DEBUG||return nothing
 GC.gc(false)
 # Phase A has completed on every owner; shared leaves are read-only here.
 # Run the diagnostic reference on an existing worker with two physical cores,
 # matching the real kernel thread count without oversubscribing the coordinator.
 study_thread_check()
 xs=Any[deepcopy(_PP_SLOTS[t]) for t in 1:N]
 push!(xs,PintVal(zeros(n,n),zeros(n)))
 nodes=Dict{Any,Any}();w=PairWork(n)
 for (depth,ops) in enumerate(study_scan_rounds(N,true)),(t,mid) in ops
  xx=deepcopy(xs[t]);yy=deepcopy(xs[mid])
  result,entry=br_populate!(xx,yy,PairWork(n))
  xs[t]=right_combine!(xs[t],xs[mid],w)
  for f in (xs[t] isa PintVal ? (:V,:v) : (:Maa,:Mab,:Mbb,:ma,:mb))
   br_check(getproperty(result,f),getproperty(xs[t],f),"shadow_$(depth)_$(t)_$(f)";exact=true)
  end
  nodes[(depth,t,mid)]=(entry=entry===nothing ? nothing : br_digest(entry),
    ma=copy(xs[t] isa PintVal ? xs[t].v : xs[t].ma),mb=xs[t] isa PintVal ? nothing : copy(xs[t].mb),
    V=xs[t] isa PintVal ? br_digest(xs[t].V) : nothing)
 end
 (nodes=nodes,values=xs)
end
function br_debug_node(depth,t,mid)
 BR_DEBUG ? BR_SHADOW[].nodes[(depth,t,mid)] : nothing
end
function br_merge(depth,t,mid,yval,last,reference)
 @assert !(BR_MODE[]===:hit && BR_LOCAL_INVALID[]) "Invalidation during hit scan"
 study_thread_check();t0=time_ns();g0=_pw_gc()
 x=_PP_SLOTS[t];y=last ? _PP_END[] : yval ? PintVal(_PP_SLOTS[mid].Maa,_PP_SLOTS[mid].ma) : _PP_SLOTS[mid]
 key=(depth,t,mid);entry=nothing;hit=0;population=0;fallback=0
 if BR_MODE[]===:hit
  # A missing entry cannot be combined with unmaterialized hit intermediates.
  # Preflight checks completeness; unexpected loss fails, never silently mixes paths.
  @assert haskey(BR_ENTRIES,key) "Unexpected mid-scan cache loss"
  entry=BR_ENTRIES[key];br_hit!(x,y,_PP_WORK[1],entry);hit=1
 elseif BR_MODE[]===:populate
  _,entry=br_populate!(x,y,_PP_WORK[1],get(BR_ENTRIES,key,nothing))
  if entry===nothing
   delete!(BR_ENTRIES,key);BR_FAILED[]=true;fallback=1
  else
   BR_ENTRIES[key]=entry;population=1
  end
 else
  STUDY_MERGE(x,y,_PP_WORK[1])
 end
 if BR_DEBUG
  if BR_MODE[]!==:fresh
   @assert (entry===nothing ? nothing : br_digest(entry))==reference.entry "Cache factor/cross/output matrix mismatch"
  end
  if BR_MODE[]===:hit && haskey(BR_OLD_VECTORS,key)
   old=BR_OLD_VECTORS[key]
   delta=norm(old-reference.ma);scale=norm(reference.ma)
   rejected=delta/max(1.,scale)>=1e-7 || delta/max(1e-10,scale)>=1e-5
   println("BR_STALE_VECTOR depth=$depth stage=$t rejected=$rejected scaled=$(delta/max(1.,scale)) relative=$(delta/max(1e-10,scale))")
  end
  BR_OLD_VECTORS[key]=copy(x.ma)
  br_check(x.ma,reference.ma,"merge_$(depth)_$(t)_ma";exact=true)
  if reference.mb!==nothing;br_check(x.mb,reference.mb,"merge_$(depth)_$(t)_mb";exact=true);end
  if reference.V!==nothing;@assert br_digest(x.Maa)==reference.V "Closed value curvature mismatch";end
 end
 (el=(time_ns()-t0)/1e9,gc=_pw_gc()-g0,hit=hit,population=population,fallback=fallback)
end
function br_worker_finish()
 BR_VALID[]=BR_ON && !BR_FAILED[] && length(BR_ENTRIES)==BR_EXPECTED[]
 (entries=length(BR_ENTRIES),bytes=sum(br_bytes,values(BR_ENTRIES);init=0),valid=BR_VALID[],mode=BR_MODE[],debug_rows=BR_DEBUG_ROWS[])
end
function br_finish(M,d,requested,hit,stats)
 info=[remotecall_fetch(br_worker_finish,p) for p in M.pids]
 println("BR_PASS solve=$(M.exits) iteration=$(d.k) barrier_iteration=$(d.j) requested=$requested hit=$hit populated=$(stats[1]) reused=$(stats[2]) fallback=$(stats[3]) entries=$(sum(x.entries for x in info)) bytes=$(sum(x.bytes for x in info)) valid=$(all(x.valid for x in info)) debug=$BR_DEBUG epoch=$(BR_GATE.epoch)")
 nothing
end
function br_ref_value(t,last)
 BR_DEBUG ? (last ? BR_SHADOW[].values[end] : BR_SHADOW[].values[t+1]) : nothing
end
function br_debug_recover(S,t,mu,lambda,Vn,R)
 BR_DEBUG||return nothing
 W=_PW[];s=W.solver;SS=deepcopy(S)
 RR=_pp_solve(SS,Vn.v,Vn.V);@assert RR!==nothing
 for f in (:xE,Symbol("Q\u00fb"),:B_active,:XB)
  # Current production field Q-u uses a combining circumflex in some snapshots.
  hasproperty(R,f)&&hasproperty(RR,f)&&br_check(getproperty(R,f),getproperty(RR,f),"recovery_$(t)_$(f)")
 end
 old=s.update[t];s.update[t]=deepcopy(old)
 try
  P=_tp_post(s,SS,RR,t,mu,lambda,W.nx,W.nu,W.nc,false,true,W.ocp.control_limits)
  @assert P!==nothing
  BR_REF_RULES[t]=s.update[t]
 finally
  s.update[t]=old
 end
 nothing
end
function br_debug_rule(t)
 BR_DEBUG||return nothing
 W=_PW[];rule=W.solver.update[t];ref=BR_REF_RULES[t]
 for f in (Symbol("\u03b1"),Symbol("\u03c8"),Symbol("\u03c7l"),Symbol("\u03c7u"),Symbol("\u03a3_L"),Symbol("\u03a3_U"))
  br_check(getproperty(rule,f),getproperty(ref,f),"rule_$(t)_$(f)")
 end
 br_check(rule.factor_policy.βB,ref.factor_policy.βB,"feedback_$(t)")
 for scale in (1e-4,1.0,10.0),kind in (1,2,3)
  dx=scale.*[kind==1 ? cos(i) : kind==2 ? sin(i) : (isodd(i÷32) ? 1. : -1.) for i in 1:W.nx]
  a,b=DDP4OPF.policy_actions!(rule,dx);c,e=DDP4OPF.policy_actions!(ref,dx)
  br_check(a,c,"control_action_$(t)_$(kind)_$(scale)")
  br_check(b,e,"equality_action_$(t)_$(kind)_$(scale)")
  br_check(rule.χl.-rule.Σ_L.*a,ref.χl.-ref.Σ_L.*c,"lower_action_$(t)_$(kind)_$(scale)")
  br_check(rule.χu.+rule.Σ_U.*a,ref.χu.+ref.Σ_U.*c,"upper_action_$(t)_$(kind)_$(scale)")
 end
 nothing
end
function br_debug_state(t,x,actual)
 BR_DEBUG||return nothing
 W=_PW[];nom=W.solver.nominal[t];ref=BR_REF_RULES[t];rows=ref.factor_policy.active_rows
 ur=copy(nom.u);ur[rows].+=ref.α[rows].+ref.factor_policy.βB*(x.-nom.x)
 br_check(actual,Vector{Float64}(W.ocp.dynamics.f(x,ur)),"actual_forward_$(t)")
 nothing
end
function br_debug_directions(t,dx,du,dphi,dzl,dzu)
 BR_DEBUG||return nothing
 ref=BR_REF_RULES[t];a,b=DDP4OPF.policy_actions!(ref,dx)
 for (name,actual,expected) in (("control",du,ref.α.+a),("equality",dphi,ref.ψ.+b),("lower",dzl,ref.χl.-ref.Σ_L.*a),("upper",dzu,ref.χu.+ref.Σ_U.*a))
  br_check(actual,expected,"actual_$(t)_$(name)")
 end
 nothing
end

function br_retire!()
 bytes=sum(br_bytes,values(BR_ENTRIES);init=0)
 br_release!()
 if BR_ON
  old=GC.enable(true);GC.gc(true);GC.enable(old)
 end
 println("BR_RETIRE bytes=$bytes collected=$BR_ON")
 nothing
end

function br_debug_before(M,N,n)
 BR_DEBUG||return nothing
 BR_SHADOW[]=remotecall_fetch(br_shadow_worker,M.pids[1],N,n)
 nothing
end
