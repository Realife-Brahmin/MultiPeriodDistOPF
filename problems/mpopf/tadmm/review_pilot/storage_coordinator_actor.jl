# Worker-side functions of the parallel storage-coordinated tADMM. Included on
# the coordinator and on every worker. Each worker owns a set of periods (their
# OPF and repair models) and a set of batteries (their whole-horizon step).
isdefined(Main,:run_storage_coordinator) || include("storage_coordinator.jl")
const SC=Dict{Symbol,Any}()
include("storage_coordinator_conic_stage.jl")
sc_conic()=SC_STAGE_SOLVER!="ipopt"

function sc_affinity()
    Sys.islinux() || return "not-linux"
    return String(strip(match(r"Cpus_allowed_list:\s*([^\n]+)",read("/proc/self/status",String)).captures[1]))
end

function sc_prepare!(datafile,system,periods,chunk,rho)
    BLAS.set_num_threads(1)
    empty!(SC); GC.gc()
    d=deserialize(datafile); apply_review_cb!(d,system)
    SC[:d]=d; SC[:periods]=sort(collect(periods)); SC[:chunk]=collect(chunk)
    SC[:stages]=Dict(t=>make_power_stage(d,t;rho) for t in periods)
    SC[:repair]=Dict(t=>make_power_stage(d,t;rho=0.0) for t in periods)
    SC[:coord]=isempty(chunk) ? nothing : make_coordinator(d;idx=chunk)
    return (;worker=myid(),pid=getpid(),host=gethostname(),affinity=sc_affinity(),periods=SC[:periods],batteries=length(chunk),
        julia_threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads())
end

# Update I for this worker's periods. zc and uc map period => vector over batteries.
function sc_xstep!(zc,uc,rho,k)
    d=SC[:d]
    return map(SC[:periods]) do t
        tick=time_ns()
        x,solver_s,iterations=power_stage_step!(SC[:stages][t],d,zc[t],uc[t],rho,k)
        (;t,x,solver_s,iterations,elapsed_s=(time_ns()-tick)/1e9)
    end
end

# Update II for this worker's batteries. target has one row per battery of the chunk.
function sc_zstep!(target,k)
    tick=time_ns()
    if sc_conic()
        row=coordinator_step_conic!(SC[:coord],SC[:d],SC[:chunk],target)
        row.ok || error("battery step failed")
        return (;z=row.z,solver_s=row.solver_s,elapsed_s=(time_ns()-tick)/1e9)
    end
    z,solver_s=coordinator_step!(SC[:coord],target,k;start=pop!(SC,:coord_resume,nothing))
    return (;z,solver_s,elapsed_s=(time_ns()-tick)/1e9)
end

sc_dispatch()=[(t,stage_dispatch(SC[:stages][t],SC[:d])) for t in SC[:periods]]

# Fix the batteries of each owned period at the given schedule and re-solve the network.
function sc_repair!(zc)
    d=SC[:d]
    return map(SC[:periods]) do t
        tick=time_ns(); s=SC[:repair][t]
        if has_values(s.m)
            vals=value.(s.vars); set_start_value.(s.vars,vals)
        end
        for (n,j) in enumerate(d[:Bset]); fix(s.pb[j],zc[t][n];force=true); end
        optimize!(s.m)
        # A tightly constrained period can stop at Ipopt's acceptable level instead
        # of the strict tolerance. Return any finite point with its status; the
        # caller judges it by the measured violation of the original model.
        strict=termination_status(s.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
        ok=has_values(s.m) && all(isfinite,value.(s.vars))
        (;t,ok,strict,status=string(raw_status(s.m)),u=ok ? stage_dispatch(s,d) : Float64[],solver_s=solve_time(s.m),elapsed_s=(time_ns()-tick)/1e9)
    end
end

# ---- Lean variant: each worker receives only its own periods' data and keeps
# one model per period, used for the iterations and, when asked, for the repair.

"""Columns of a per-period data matrix for a subset of periods, indexed like the full matrix."""
struct PeriodColumns <: AbstractMatrix{Float64}
    rows::Int
    periods::Int
    position::Dict{Int,Int}
    columns::Matrix{Float64}
end
Base.size(a::PeriodColumns)=(a.rows,a.periods)
Base.getindex(a::PeriodColumns,j::Int,t::Int)=a.columns[j,a.position[t]]

"""Copy of the instance holding per-period matrices for `periods` only. Call after
the screening and layout caches exist, so no worker needs the other periods."""
function slim_instance(d,periods)
    T=d[:T]; position=Dict(t=>n for (n,t) in enumerate(periods)); out=Dict{Symbol,Any}()
    for (key,value) in d
        out[key]=value isa Matrix{Float64} && size(value,2)==T && size(value,1)>1 ?
            PeriodColumns(size(value,1),T,position,value[:,collect(periods)]) : value
    end
    return out
end

function sc_prepare_lean!(d,periods,chunk,rho)
    BLAS.set_num_threads(1)
    empty!(SC); GC.gc()
    SC[:d]=d; SC[:periods]=sort(collect(periods)); SC[:chunk]=collect(chunk)
    SC[:stages]=Dict(t=>(sc_conic() ? make_power_stage_conic(d,t;rho) : make_power_stage(d,t;rho)) for t in periods)
    SC[:coord]=isempty(chunk) ? nothing : (sc_conic() ? make_coordinator_conic(d;idx=chunk) : make_coordinator(d;idx=chunk))
    SC[:resume]=Dict{Int,Vector{Float64}}()
    return (;worker=myid(),pid=getpid(),host=gethostname(),affinity=sc_affinity(),periods=SC[:periods],batteries=length(chunk),
        julia_threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads())
end

function sc_xstep_lean!(zc,uc,rho,k)
    d=SC[:d]
    return map(SC[:periods]) do t
        tick=time_ns()
        x,solver_s,iterations=sc_conic() ? power_stage_step_conic!(SC[:stages][t],d,zc[t],uc[t],rho) :
            power_stage_step!(SC[:stages][t],d,zc[t],uc[t],rho,k;start=pop!(SC[:resume],t,nothing))
        (;t,x,solver_s,iterations,elapsed_s=(time_ns()-tick)/1e9)
    end
end

# Repair in the iteration model itself: drop the penalty, fix the batteries at the
# schedule, solve, then restore the bounds and remember the iterate to resume from.
function sc_repair_shared!(zc;max_iter=3000)
    d=SC[:d]
    return map(SC[:periods]) do t
        tick=time_ns(); s=SC[:stages][t]
        if sc_conic()
            row=repair_stage_conic!(s,d,zc[t])
            return (;t,ok=row.ok,strict=row.strict,status=row.status,u=row.u,solver_s=row.solver_s,elapsed_s=(time_ns()-tick)/1e9)
        end
        # The iterate to resume from: saved by an earlier repair of this iteration, else the model's own.
        iterate=haskey(SC[:resume],t) ? SC[:resume][t] : value.(s.vars)
        set_start_value.(s.vars,iterate)
        for (n,j) in enumerate(d[:Bset])
            set_objective_coefficient(s.m,s.pb[j],0.0)
            set_objective_coefficient(s.m,s.pb[j],s.pb[j],s.cb)
            fix(s.pb[j],zc[t][n];force=true)
        end
        # A period whose network cannot absorb the fixed batteries runs long before
        # Ipopt gives up; the cap makes such an attempt cheap. It is then not strict.
        set_optimizer_attribute(s.m,"max_iter",max_iter)
        optimize!(s.m)
        # SC_REPAIR_ACCEPT=1 also takes Ipopt's acceptable level: the dispatch is scored
        # against the original model afterwards (violation <= no_tol), which decides.
        accepted=get(ENV,"SC_REPAIR_ACCEPT","0")=="1" ? (MOI.LOCALLY_SOLVED,MOI.OPTIMAL,MOI.ALMOST_LOCALLY_SOLVED) : (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
        strict=termination_status(s.m) in accepted
        ok=has_values(s.m) && all(isfinite,value.(s.vars))
        u=ok ? stage_dispatch(s,d) : Float64[]
        status=string(raw_status(s.m)); solver_s=solve_time(s.m)   # read before the model changes again
        set_optimizer_attribute(s.m,"max_iter",3000)
        for j in d[:Bset]
            unfix(s.pb[j]); set_lower_bound(s.pb[j],-d[:P_B_R_pu][j]); set_upper_bound(s.pb[j],d[:P_B_R_pu][j])
        end
        SC[:resume][t]=iterate
        (;t,ok,strict,status,u,solver_s,elapsed_s=(time_ns()-tick)/1e9)
    end
end

# Battery step with some periods pinned: the schedule closest to `target` that
# the batteries can follow when their powers in the pinned periods are given.
# pinned maps period => powers of this worker's batteries.
function sc_zstep_pinned!(target,pinned)
    c=SC[:coord]; d=SC[:d]; tick=time_ns(); limits=[d[:P_B_R_pu][j] for j in d[:Bset][SC[:chunk]]]
    if sc_conic()
        row=coordinator_step_conic!(c,d,SC[:chunk],target;pinned)
        return (;ok=row.ok,z=row.z,elapsed_s=(time_ns()-tick)/1e9)
    end
    iterate=haskey(SC,:coord_resume) ? SC[:coord_resume] : value.(c.vars)
    set_start_value.(c.vars,iterate)
    for i in eachindex(target); set_objective_coefficient(c.m,c.z[i],-2*target[i]); end
    for (t,p) in pinned, n in eachindex(p); fix(c.z[n,t],p[n];force=true); end
    optimize!(c.m)
    ok=termination_status(c.m) in (MOI.LOCALLY_SOLVED,MOI.OPTIMAL)
    z=ok ? value.(c.z) : zeros(0,0)
    for (t,p) in pinned, n in eachindex(p)
        unfix(c.z[n,t]); set_lower_bound(c.z[n,t],-limits[n]); set_upper_bound(c.z[n,t],limits[n])
    end
    SC[:coord_resume]=iterate
    return (;ok,z,elapsed_s=(time_ns()-tick)/1e9)
end
