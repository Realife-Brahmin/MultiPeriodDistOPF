# Opt-in diagnostic capture only; baseline arithmetic/decisions unchanged.
using SHA,TOML
const _FC_STAGE=Dict{Int,Any}()
const _FC_STATE=Dict{Symbol,Any}()
const AUDIT_KEYS=Dict((0,0)=>"early_before",(1,0)=>"early_after",
 (4,0)=>"barrier_before",(4,1)=>"barrier_after",
 (7,3)=>"middle_before",(8,3)=>"middle_after",
 (12,7)=>"late_before",(13,7)=>"late_after")
audit_tag(k,j)=get(AUDIT_KEYS,(Int(k),Int(j)),"")
fc_active(exits)=exits>0||get(ENV,"FILTERDDP_WARMUP","0")=="0"
function fc_cache_stage(t,S,W)
    _FC_STAGE[t]=(S=deepcopy(S),cl=deepcopy(W.ocp.control_limits),nx=W.nx,nu=W.nu,nc=W.nc)
end
function fc_before(M,n,k,j,mu,reg_before=NaN)
    fc_active(M.exits)||return nothing
    if k==0&&j==0
        empty!(_FC_STATE);_FC_STATE[:index]=0;_FC_STATE[:tags]=Set{String}()
    end
    _FC_STATE[:index]+=1
    tag=audit_tag(k,j);isempty(tag)&&return nothing
    slots=PintElem[]
    for a in M.shared
        p=pointer(sdata(a));blk=i->copy(unsafe_wrap(Array,p+(i-1)*n*n*sizeof(Float64),(n,n)))
        vec=i->copy(unsafe_wrap(Array,p+(3n*n+(i-1)*n)*sizeof(Float64),(n,)))
        push!(slots,PintElem(blk(1),blk(2),blk(3),vec(1),vec(2)))
    end
    (elems=slots,k=k,j=j,mu=mu,index=_FC_STATE[:index],tag=tag,reg_before=reg_before)
end
function fc_write_stage(t,dir,actual_dx)
    W=_PW[];rec=_FC_STAGE[t];rule=W.solver.update[t];probes=Any[]
    actual_du=nothing
    for (kind,dx) in (("zero",zeros(W.nx)),("sine",[1e-3*sin(i) for i in 1:W.nx]),("actual_soc",actual_dx))
        beta,omega=DDP4OPF.policy_actions!(rule,dx)
        du=copy(rule.α+beta);dphi=copy(rule.ψ+omega)
        dzl=copy(rule.χl-rule.Σ_L.*beta);dzu=copy(rule.χu+rule.Σ_U.*beta)
        push!(probes,(kind=kind,dx=copy(dx),du=du,dphi=dphi,dzl=dzl,dzu=dzu))
        kind=="actual_soc"&&(actual_du=du)
    end
    next_dx=rec.S.fx*actual_dx+rec.S.fu*actual_du
    @assert all(isfinite,next_dx)
    f=[nonzeros(rec.S.fu_active)[rec.S.fu_active.colptr[i]] for i in 1:W.nx]
    packet=(stage=t,S=rec.S,cl=rec.cl,nx=rec.nx,nu=rec.nu,nc=rec.nc,
      nominal=(x=copy(rec.S.x),u=copy(rec.S.u),phi=copy(rec.S.ϕ),zl=copy(rec.S.zl),zu=copy(rec.S.zu)),
      f=f,probes=probes,actual_next_dx=next_dx)
    file=joinpath(dir,"stage_$(t).jls");serialize(file*".partial",packet);mv(file*".partial",file)
    (file=basename(file),sha=bytes2hex(sha256(read(file))),bytes=filesize(file),next_dx=next_dx)
end
function fc_after(payload,M,d,options,ocp)
    isnothing(payload)&&return
    state=_FC_STATE;tag=payload.tag;@assert !(tag in state[:tags])
    root=ENV["FILTERDDP_FROZEN_CAPTURE_DIR"];mkpath(root)
    name=@sprintf("pass%04d_k%04d_j%03d",payload.index,d.k,d.j)
    dir=joinpath(root,name);@assert !ispath(dir);mkpath(dir)
    file=joinpath(dir,"leaves.jls");serialize(file*".partial",payload);mv(file*".partial",file)
    records=Any[];actual_dx=zeros(length(first(payload.elems).ma))
    for t in 1:ocp.N
        r=remotecall_fetch(fc_write_stage,M.owner[t],t,dir,actual_dx)
        push!(records,r);actual_dx=r.next_dx
    end
    meta=Dict("tags"=>[tag],"pass"=>payload.index,"iteration"=>Int(d.k),"barrier_iteration"=>Int(d.j),
      "mu"=>Float64(d.μ),"objective"=>Float64(d.objective),"primal"=>Float64(d.primal_inf),
      "dual"=>Float64(d.dual_inf),"complementarity"=>Float64(d.cs_inf_0),"stages"=>ocp.N,
      "regularization_before_pass"=>Float64(payload.reg_before),"regularization_reported_after_pass"=>Float64(d.reg_last),
      "regularization_semantics"=>"procs hook records before-pass reg_last; original successful procs pass sets reported reg_last=0",
      "schedule"=>"original_balanced","kernel"=>"original_inplace","full_stage_records"=>true,
      "leaves_sha256"=>bytes2hex(sha256(read(file))),
      "stage_files"=>[Dict("file"=>r.file,"sha256"=>r.sha,"bytes"=>r.bytes) for r in records])
    open(joinpath(dir,"metadata.toml"),"w") do io;TOML.print(io,meta);end
    coverage=joinpath(root,"coverage.csv");fresh=!isfile(coverage)
    open(coverage,"a") do io
        fresh&&println(io,"tag,snapshot,pass,iteration,barrier_iteration,mu,objective,primal,dual,complementarity")
        println(io,join((tag,name,payload.index,d.k,d.j,d.μ,d.objective,d.primal_inf,d.dual_inf,d.cs_inf_0),','))
    end
    push!(state[:tags],tag)
    println("INFLUENCE_CAPTURE tag=$tag snapshot=$name k=$(d.k) j=$(d.j) mu=$(d.μ) reg_before=$(payload.reg_before) stages=$(ocp.N)")
    flush(stdout)
end

