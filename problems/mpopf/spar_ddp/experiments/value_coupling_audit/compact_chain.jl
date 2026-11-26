# Isolated retained-interior composition. No production hook. A chain keeps
# the preceding Schur object and grows its interior: NOT bounded-rank closure.
const CC_LOCAL_FALLBACK_CAP=192
function cc_factor(a,stats)
    n=length(a.tree.ids);empty=zeros(0,0)
    u,v=(!isnothing(a.localblock) ? (empty,empty) : ck_factors(a.lr))
    w,z=(!isnothing(a.localblock) ? (empty,empty) : ck_factors(a.rl))
    rankgrowth=n>CK_LOCAL_DENSE_CAP&&size(u,2)+size(w,2)>=n
    if n<=CK_LOCAL_DENSE_CAP || rankgrowth&&n<=CC_LOCAL_FALLBACK_CAP
        A=ck_apply_ordered(a,Matrix{Float64}(I,n,n))
        rankgrowth&&ck_bump!(stats,"rank_growth_local_dense_fallbacks")
        stats["max_local_dense_dim"]=max(get(stats,"max_local_dense_dim",0),n)
        return CKFactor(a.tree.ids,nothing,nothing,lu(A),nothing,empty,empty,empty,empty)
    end
    rankgrowth&&error("CC_RANK_GROWTH node=$n core=$(size(u,2)+size(w,2)) exceeds bounded local fallback cap")
    lf=cc_factor(a.left,stats);rf=cc_factor(a.right,stats);r=size(u,2);s=size(w,2)
    stats["max_woodbury_core_dim"]=max(get(stats,"max_woodbury_core_dim",0),r+s)
    au=ck_solve_ordered(lf,u);bw=ck_solve_ordered(rf,w)
    core=[Matrix{Float64}(I,r,r) v'*bw;z'*au Matrix{Float64}(I,s,s)]
    CKFactor(a.tree.ids,lf,rf,nothing,lu(core),au,bw,v,z)
end
function cc_lift(ds,bs,n,tol,stats)
    m=length(ds);@assert length(bs)==m-1
    tree=ds[1].tree
    function walk(d,b)
        t=d[1].tree
        ids=vcat([t.ids .+ (j-1)*n for j in 1:m]...)
        if !isnothing(t.left)
            l=walk([a.left for a in d],[a.left for a in b])
            r=walk([a.right for a in d],[a.right for a in b])
            tr=CouplingTree(vcat(l.tree.ids,r.tree.ids),l.tree,r.tree,Int[])
            us=Matrix{Float64}[];vs=Matrix{Float64}[]
            function addoff(o,i,j)
                u,v=ck_factors(o)
                U=zeros(length(l.tree.ids),size(u,2));V=zeros(length(r.tree.ids),size(v,2))
                li=indexin(t.left.ids .+ (i-1)*n,l.tree.ids)
                ri=indexin(t.right.ids .+ (j-1)*n,r.tree.ids)
                U[Int.(li),:]=u;V[Int.(ri),:]=v;push!(us,U);push!(vs,V)
            end
            for j in 1:m;addoff(d[j].lr,j,j);end
            for j in 1:m-1
                addoff(b[j].lr,j,j+1);addoff(ck_transpose(b[j].rl),j+1,j)
            end
            off=ck_off(hcat(us...),hcat(vs...),tol,stats)
            return FollowupBlockNode(tr,nothing,l,r,off,ck_transpose(off),l.entries+r.entries+2off.entries)
        end
        tr=coupling_leaf(ids);A=zeros(length(ids),length(ids))
        for j in 1:m
            ii=Int.(indexin(t.ids .+ (j-1)*n,tr.ids));A[ii,ii]=d[j].localblock
            if j<m
                jj=Int.(indexin(t.ids .+ j*n,tr.ids));A[ii,jj]=b[j].localblock;A[jj,ii]=b[j].localblock'
            end
        end
        FollowupBlockNode(tr,A,nothing,nothing,nothing,nothing,length(A))
    end
    walk(ds,bs)
end
struct CCChain
    previous::CKSchur
    left::Any
    H::FollowupBlockNode
    factor::CKFactor
    aa::FollowupBlockNode
    bb::FollowupBlockNode
    ma::Vector{Float64}
    mb::Vector{Float64}
end
function cc_prepend(p,c::CKSchur;tol=1e-12,stats=Dict{String,Int}())
    n=length(p.ma)
    a=ck_upper(ck_add(p.reps[3],c.x.reps[1],tol,stats),tol,stats)
    H=cc_lift([a,c.K],[c.x.reps[2]],n,tol,stats);F=cc_factor(H,stats)
    z=ck_solve(F,vcat(p.mb+c.x.ma,c.x.mb+c.y.ma))
    aa=ck_upper(p.reps[1],tol,stats)
    CCChain(c,p,H,F,aa,c.bb,p.ma-ck_apply(p.reps[2],z[1:n]),c.y.mb-ck_apply(c.y.reps[2],z[n+1:end];trans=true))
end
function ck_action(c::CCChain,block,X)
    n=length(c.ma);x=reshape(X,n,:);Z=zeros(size(x));B=c.left.reps[2];D=c.previous.y.reps[2]
    if block==:Maa
        z=ck_solve(c.factor,vcat(ck_apply(B,x;trans=true),Z))
        y=ck_apply(c.aa,x)-ck_apply(B,z[1:n,:])
    elseif block==:Mab
        z=ck_solve(c.factor,vcat(Z,ck_apply(D,x)))
        y=-ck_apply(B,z[1:n,:])
    elseif block==:Mbb
        z=ck_solve(c.factor,vcat(Z,ck_apply(D,x)))
        y=ck_apply(c.bb,x)-ck_apply(D,z[n+1:end,:];trans=true)
    else;error("unknown chain block");end
    X isa AbstractVector ? vec(y) : y
end
struct CCClosed
    chain::CCChain
    H::FollowupBlockNode
    factor::CKFactor
    gradient::Vector{Float64}
    v::Vector{Float64}
end
function cc_close(c::CCChain,terminal;tol=1e-12,stats=Dict{String,Int}())
    p=c.left;q=c.previous.x;r=c.previous.y;n=length(p.ma)
    a=ck_upper(ck_add(p.reps[3],q.reps[1],tol,stats),tol,stats)
    d=ck_upper(ck_add(r.reps[3],terminal.rep,tol,stats),tol,stats)
    H=cc_lift([a,c.previous.K,d],[q.reps[2],r.reps[2]],n,tol,stats)
    F=cc_factor(H,stats);g=vcat(p.mb+q.ma,q.mb+r.ma,r.mb+terminal.v)
    z=ck_solve(F,g)
    CCClosed(c,H,F,g,p.ma-ck_apply(p.reps[2],z[1:n]))
end
function cc_value(c::CCClosed,X)
    n=length(c.v);x=reshape(X,n,:);z=ck_solve(c.factor,vcat(ck_apply(c.chain.left.reps[2],x;trans=true),zeros(2n,size(x,2))))
    y=ck_apply(c.chain.aa,x)-ck_apply(c.chain.left.reps[2],z[1:n,:])
    X isa AbstractVector ? vec(y) : y
end
# Bounded matrix-free consumer. Non-SPD breakdown, iteration limit or failed
# true residual is rejection, never silently accepted or a relaxed gate.
function cc_cg(action,b,diagonal;tol=2e-13,maxiter=300)
    x=zeros(length(b));r=copy(b);d=max.(abs.(diagonal),eps());z=r./d;p=copy(z);rz=dot(r,z)
    target=tol*max(norm(b),1e-12)
    norm(r)<=target && return (x=x,iterations=0,residual=norm(r))
    for k in 1:maxiter
        q=action(p);pq=dot(p,q);pq>0 || error("CC_CONSUMER_NOT_SPD iteration=$k")
        alpha=rz/pq;x .+= alpha.*p;r .-= alpha.*q
        if norm(r)<=target
            true_r=norm(action(x)-b)
            true_r<=10target || error("CC_CONSUMER_TRUE_RESIDUAL residual=$true_r target=$target")
            return (x=x,iterations=k,residual=true_r)
        end
        z=r./d;newrz=dot(r,z);p=z+(newrz/rz).*p;rz=newrz
    end
    error("CC_CONSUMER_ITERATION_CAP maxiter=$maxiter residual=$(norm(r)) target=$target")
end
function cc_stage(packet,v,action,mu,dx;stats=Dict{String,Int}())
    S=deepcopy(packet.S);lay=S.lay;stat=S.stat;b0=S.b0;n=length(v)
    f=[nonzeros(S.fu_active)[S.fu_active.colptr[k]] for k in 1:n]
    @assert all(k->rowvals(S.fu_active)[S.fu_active.colptr[k]]==k,1:n)
    # Raw Sh0 is retained. Production's grouped-order upper convention is
    # explicit here; future action symmetry is independently audited.
    Sh=Symmetric(b0.Sh0,:U);bat=b0.bat
    function schur(x)
        a=zeros(n);a[bat]=x
        Sh*x+(f.*action(f.*a))[bat]
    end
    b=copy(S.b);b[1:packet.nu]=-(S.lu_+S.fu'*v+S.mubar+S.cuϕ)
    b[S.active] .-= S.fu_active'*action(dx)
    S.lux_in_B && (b[S.active] .-= S.lux[S.active,:]*dx)
    b[packet.nu+1:end]=-S.c-S.cx*dx
    y=zeros(length(b));y2=similar(y);t2=zeros(length(b))
    _tp_network_solve!(y,lay,S.fac,stat.ooff,S.pre,S.z,b)
    r=reshape(copy(b[lay.E]),:,1)
    for (a,t,p) in stat.kef;r[a]-=S.nzc[p]*y[lay.Fc[t]];end
    for (si,i) in enumerate(lay.iso),(a,p) in stat.iso_e[si];r[a]-=S.nzc[p]*b[i]/S.nzc[stat.iso_d[si]];end
    rp=r[stat.perm,:]
    for g in eachindex(stat.granges)
        rg=stat.granges[g];np=b0.nP[g];rN=view(rp,rg[np+1:end],:)
        rN .*= b0.Ninv[g];rp[rg[1:np],:] .-= b0.PN[g]*rN
    end
    sol=cc_cg(schur,vec(rp[b0.prow,:]),diag(Sh))
    stats["stage_cg_maxiter"]=max(get(stats,"stage_cg_maxiter",0),sol.iterations)
    rp[b0.prow,:]=reshape(sol.x,:,1)
    for g in eachindex(stat.granges)
        rg=stat.granges[g];np=b0.nP[g]
        rp[rg[np+1:end],:] .-= b0.Ninv[g].*(b0.NP[g]*rp[rg[1:np],:])
    end
    xe=similar(rp);xe[stat.perm,:]=rp
    for (a,j,p) in stat.kef;t2[lay.Fc[j]]+=S.nzc[p]*xe[a];end
    _tp_network_solve!(y2,lay,S.fac,stat.ooff,S.pre,S.z,t2)
    for (si,i) in enumerate(lay.iso)
        z=b[i];for (a,p) in stat.iso_e[si];z-=S.nzc[p]*xe[a];end
        b[i]=z/S.nzc[stat.iso_d[si]]
    end
    for own in lay.own,i in own;b[i]=y[i]-y2[i];end
    b[lay.E]=vec(xe);du=b[1:packet.nu]
    (du=du,dphi=b[packet.nu+1:end],dzl=S.inv_ul.*mu-S.zl-S.Σ_L.*du,dzu=S.inv_uu.*mu-S.zu+S.Σ_U.*du)
end
