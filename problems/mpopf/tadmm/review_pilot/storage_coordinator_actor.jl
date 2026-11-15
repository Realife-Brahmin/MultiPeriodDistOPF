# Worker-side functions of the parallel storage-coordinated tADMM. Included on
# the coordinator and on every worker. Each worker owns a set of periods (their
# OPF and repair models) and a set of batteries (their whole-horizon step).
isdefined(Main,:run_storage_coordinator) || include("storage_coordinator.jl")
const SC=Dict{Symbol,Any}()

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
    z,solver_s=coordinator_step!(SC[:coord],target,k)
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
