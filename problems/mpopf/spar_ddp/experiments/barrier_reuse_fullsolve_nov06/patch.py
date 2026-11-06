# Apply ONLY to a new immutable copy of the validated grouped/right runtime.
from pathlib import Path
import sys
root=Path(sys.argv[1]);ps=root/'spar_ddp/examples/power_system'
p=ps/'procs_pass.jl';s=p.read_text(encoding='utf-8')
module='DDP4OPF' if 'DDP4OPF.Solver' in s else 'SPARDDP'
def once(a,b):
 global s
 assert s.count(a)==1,(a,s.count(a))
 s=s.replace(a,b)
once('function pw_full_gc()\n','function pw_full_gc()\n    br_release!()\n')
once('function pw_pass_start(advance::Bool)','function pw_pass_start(advance::Bool, request, context, expected, reset)')
a=s.index('function pw_pass_start(');b=s.index('\nend',a)
part=s[a:b];part=part.replace('    if advance\n','    ready=br_worker_start(request,context,expected,reset)\n    if advance\n')
part=part.replace('return (el=','return (ready=ready, el=');s=s[:a]+part+s[b:]
once('isnothing(_pp_elem(S, t, W.nx)) && return nothing','isnothing(br_leaf!(S, t, W.nx)) && return nothing')
a=s.index('function pw_combine(');b=s.index('\nend',a)+len('\nend')
s=s[:a]+s[b:]
once('function pw_step3(t::Int, μ::Float64, λ::Vector{Float64}, last::Bool)','function pw_step3(t::Int, μ::Float64, λ::Vector{Float64}, last::Bool, refV)')
once('    isnothing(R) && return nothing\n    P = _tp_post','    isnothing(R) && return nothing\n    br_debug_recover(S,t,μ,λ,refV,R)\n    P = _tp_post')
once('    isnothing(P) && return nothing\n    study_thread_check()','    isnothing(P) && return nothing\n    br_debug_rule(t)\n    study_thread_check()')
once('    return Vector{Float64}(W.ocp.dynamics.f(x, ur))','    result=Vector{Float64}(W.ocp.dynamics.f(x, ur))\n    br_debug_state(t,x,result)\n    return result')
once('    W.dir[t] = (dx=δx, du=du, dϕ=dϕ, dzl=dzl, dzu=dzu)','    br_debug_directions(t,δx,du,dϕ,dzl,dzu)\n    W.dir[t] = (dx=δx, du=du, dϕ=dϕ, dzl=dzl, dzu=dzu)')
once('    t_start = time_ns()\n    fresh =','    t_start = time_ns()\n    requested=br_begin(d);context=br_context()\n    expected=Dict(p=>0 for p in M.pids)\n    for ops in study_scan_rounds(N,true), (q,(t,mid)) in enumerate(ops)\n        expected[M.pids[mod1(q,length(M.pids))]]+=1\n    end\n    br_stats=zeros(Int,3)\n    fresh =')
once('remotecall_fetch(pw_pass_start, p, advance)','remotecall_fetch(pw_pass_start, p, advance, requested, context, expected[p], fresh)')
once('    t_0 = time_ns()\n','''    hit=requested && all(r.ready for r in r0)
    if requested && !hit
        @sync for p in M.pids;@async remotecall_fetch(br_force_population,p);end
    end
    t_0 = time_ns()
''')
once('    # 2. the rounds\n','    br_debug_before(M,N,nx)\n    # 2. the rounds\n')
once('for pairs in study_scan_rounds(N, get(ENV,"STUDY_SCAN","balanced")=="grouped")','for (depth,pairs) in enumerate(study_scan_rounds(N, get(ENV,"STUDY_SCAN","balanced")=="grouped"))')
once('remotecall_fetch(pw_combine,','remotecall_fetch(br_merge,')
once('                t, mid, isval[mid], mid == N + 1)','                depth, t, mid, isval[mid], mid == N + 1, br_debug_node(depth,t,mid))')
once('        for (t, mid) in pairs; isval[t] = isval[mid]; end','        br_stats .+= [sum(r.population for r in res),sum(r.hit for r in res),sum(r.fallback for r in res)]\n        for (t, mid) in pairs; isval[t] = isval[mid]; end')
once('remotecall_fetch(pw_step3, w, t, μ, λs[t], t == N)','remotecall_fetch(pw_step3, w, t, μ, λs[t], t == N, br_ref_value(t,t==N))')
once('    t_3 = time_ns()\n','    t_3 = time_ns()\n    br_finish(M,d,requested,hit,br_stats)\n')
once('    M.exits += 1','    M.exits += 1\n    @sync for p in M.pids;@async remotecall_fetch(br_retire!,p);end')
once('    _PM[] = ProcsMain(pids, owner, shared, -1, Float64[], 0)','    _PM[] = ProcsMain(pids, owner, shared, -1, Float64[], 0)\n    '+module+'.STUDY_STATE_EVENT[]=event->br_event!(BR_GATE,event)')
s+='\ninclude(joinpath(@__DIR__,"br_integration.jl"))\n'
p.write_text(s,encoding='utf-8')
paths=list(root.glob('spar_ddp/*/src/solve.jl'));assert len(paths)==1,paths
p=paths[0];s=p.read_text(encoding='utf-8')
s='const STUDY_STATE_EVENT=Ref{Any}(nothing)\nstudy_state_event(e)=(STUDY_STATE_EVENT[]===nothing ? nothing : STUDY_STATE_EVENT[](e))\n'+s
once('    reuse_reset!()','    reuse_reset!()\n    study_state_event(:solve_start)')
once('            data.j += 1','            study_state_event(:barrier)\n            data.j += 1')
once('        update_nominal_trajectory!(solver)','        study_state_event(:accepted)\n        update_nominal_trajectory!(solver)')
p.write_text(s,encoding='utf-8')
p=ps/'br_integration.jl';text=p.read_text(encoding='utf-8').replace('DDP4OPF.',module+'.');p.write_text(text,encoding='utf-8')

# Preserve the separate diagnostic warmup transcript; clean warmup unchanged.
p=ps/'ieee123c_filterddp.jl';s=p.read_text(encoding='utf-8')
once('    rm(scratch; force=true)','''    if get(ENV,"STUDY_BARRIER_DEBUG","0")=="1"
        cp(scratch,joinpath(ENV["STUDY_OUT"],"warmup-debug.log");force=false)
    end
    rm(scratch; force=true)''')
p.write_text(s,encoding='utf-8')
