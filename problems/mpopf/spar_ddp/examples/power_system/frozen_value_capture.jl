# Diagnostic capture only. No value approximation and no solver decisions.
using SHA, TOML
const _FC_STAGE=Dict{Int,Any}()
const _FC_STATE=Dict{Symbol,Any}()
fc_active(exits)=exits>0 || get(ENV,"FILTERDDP_WARMUP","0")=="0"
function fc_cache_stage(t,S,W)
    # Capture all stages' quadratic leaves at coordinator; only three full
    # local KKT/derivative records are needed for recovered-direction probes.
    t in unique((1,cld(W.ocp.N,2),W.ocp.N)) || return
    _FC_STAGE[t]=(S=deepcopy(S),cl=deepcopy(W.ocp.control_limits),nx=W.nx,nu=W.nu,nc=W.nc)
end
function fc_before(M,n,k,j,mu)
    fc_active(M.exits) || return nothing
    if k==0 && j==0
        empty!(_FC_STATE);_FC_STATE[:index]=0;_FC_STATE[:tags]=Set{String}()
    end
    _FC_STATE[:index]+=1
    slots=PintElem[]
    for a in M.shared
        p=pointer(sdata(a));blk=i->copy(unsafe_wrap(Array,p+(i-1)*n*n*sizeof(Float64),(n,n)))
        vec=i->copy(unsafe_wrap(Array,p+(3n*n+(i-1)*n)*sizeof(Float64),(n,)))
        push!(slots,PintElem(blk(1),blk(2),blk(3),vec(1),vec(2)))
    end
    (elems=slots,k=k,j=j,mu=mu,index=_FC_STATE[:index])
end
function fc_write_stage(t,dir)
    W=_PW[];rec=_FC_STAGE[t];rule=W.solver.update[t]
    probes=Vector{Any}()
    for dx in (zeros(W.nx),[1e-3*sin(i) for i in 1:W.nx])
        du,dp=SPARDDP.policy_actions!(rule,dx)
        push!(probes,(dx=copy(dx),du=rule.α+du,dphi=rule.ψ+dp))
    end
    # S was copied BEFORE phase C mutates scratch right-hand sides.
    f=[nonzeros(rec.S.fu_active)[rec.S.fu_active.colptr[i]] for i in 1:W.nx]
    packet=(stage=t,S=rec.S,cl=rec.cl,nx=rec.nx,nu=rec.nu,nc=rec.nc,
            nominal=(x=copy(rec.S.x),u=copy(rec.S.u),phi=copy(rec.S.ϕ),zl=copy(rec.S.zl),zu=copy(rec.S.zu)),
            f=f,probes=probes)
    file=joinpath(dir,"stage_$(t).jls");serialize(file*".partial",packet);mv(file*".partial",file)
    (file=basename(file),sha=bytes2hex(sha256(read(file))),bytes=filesize(file))
end
function fc_select!(state,d,options,ocp)
    tags=state[:tags];want=String[]
    !("early" in tags) && push!(want,"early")
    middle=parse(Int,get(ENV,"FILTERDDP_FROZEN_MIDDLE_K","5"))
    d.k>=middle && !("middle" in tags) && push!(want,"middle")
    if get(state,:await_after,false) && d.j>state[:before_j]
        @assert d.k==state[:before_k] && d.μ<state[:before_mu]
        push!(want,"barrier_after");state[:await_after]=false
    end
    err0=max(d.dual_inf,d.cs_inf_0,d.primal_inf)
    errmu=max(d.dual_inf,d.cs_inf_μ,d.primal_inf)
    ni=(ocp.control_limits.nl+ocp.control_limits.nu)*ocp.N
    changes=err0>=options.optimality_tolerance && errmu<=options.κ_ϵ*d.μ && ni>0 && d.μ>options.optimality_tolerance/10
    if changes && !("barrier_before" in tags)
        push!(want,"barrier_before");state[:await_after]=true
        state[:before_j]=d.j;state[:before_k]=d.k;state[:before_mu]=d.μ
    end
    err0<options.optimality_tolerance && !("late" in tags) && push!(want,"late")
    want
end
function fc_after(payload,M,d,options,ocp)
    isnothing(payload) && return
    state=_FC_STATE;tags=state[:tags];want=fc_select!(state,d,options,ocp)
    isempty(want) && return
    root=ENV["FILTERDDP_FROZEN_CAPTURE_DIR"];mkpath(root)
    name=@sprintf("pass%04d_k%04d_j%03d",payload.index,d.k,d.j)
    dir=joinpath(root,name);@assert !ispath(dir);mkpath(dir)
    file=joinpath(dir,"leaves.jls");serialize(file*".partial",payload);mv(file*".partial",file)
    selected=unique((1,cld(ocp.N,2),ocp.N));records=Any[]
    for t in selected;push!(records,remotecall_fetch(fc_write_stage,M.owner[t],t,dir));end
    meta=Dict("tags"=>want,"pass"=>payload.index,"iteration"=>Int(d.k),"barrier_iteration"=>Int(d.j),"mu"=>Float64(d.μ),
              "objective"=>Float64(d.objective),"primal"=>Float64(d.primal_inf),"dual"=>Float64(d.dual_inf),"complementarity"=>Float64(d.cs_inf_0),
              "stages"=>ocp.N,"leaves_sha256"=>bytes2hex(sha256(read(file))),"stage_files"=>[Dict("file"=>r.file,"sha256"=>r.sha,"bytes"=>r.bytes) for r in records])
    open(joinpath(dir,"metadata.toml"),"w") do io;TOML.print(io,meta);end
    coverage=joinpath(root,"coverage.csv");fresh=!isfile(coverage)
    open(coverage,"a") do io
        fresh && println(io,"tag,snapshot,pass,iteration,barrier_iteration,mu,objective,primal,dual,complementarity")
        for tag in want
            println(io,join((tag,name,payload.index,d.k,d.j,d.μ,d.objective,d.primal_inf,d.dual_inf,d.cs_inf_0),','))
            push!(tags,tag)
        end
    end
    println("FROZEN_VALUE_CAPTURE tags=$(join(want,';')) snapshot=$name iteration=$(d.k) barrier_iteration=$(d.j) mu=$(d.μ) exact_dense=true")
    flush(stdout)
end
