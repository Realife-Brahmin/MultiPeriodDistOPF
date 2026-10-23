# Diagnostic only. Direct sparse-local + separator-factor generation. Never
# reads S.b0.Sh0 or an already assembled boundary quadratic for generation.
function di_record_value(io,phase,field,A,B)
    absolute_error=norm(A-B);reference_norm=norm(B)
    scaled=absolute_error/max(1.0,reference_norm);@assert isfinite(scaled)
    println(io,join((phase,field,absolute_error,reference_norm,scaled),','));flush(io)
    scaled
end
function di_payload(P,U,V,tree,stats;tol=1e-12)
    ids=tree.ids
    if isnothing(tree.left)
        B=Matrix(P[ids,ids])+U[ids,:]*V[ids,:]'
        return FollowupBlockNode(tree,B,nothing,nothing,nothing,nothing,length(B))
    end
    l=di_payload(P,U,V,tree.left,stats;tol=tol);r=di_payload(P,U,V,tree.right,stats;tol=tol)
    function cross(I,J)
        F=P[I,J];ii,jj,_=findnz(F);ri=unique(ii);ci=unique(jj)
        if length(ri)<=length(ci)
            a=zeros(length(I),length(ri));for (k,i) in enumerate(ri);a[i,k]=1;end
            b=Matrix(F[ri,:]')
        else
            a=Matrix(F[:,ci]);b=zeros(length(J),length(ci));for (k,j) in enumerate(ci);b[j,k]=1;end
        end
        stats["max_sparse_cross_interface"]=max(get(stats,"max_sparse_cross_interface",0),size(a,2))
        ck_off(hcat(a,U[I,:]),hcat(b,V[J,:]),tol,stats)
    end
    lr=cross(tree.left.ids,tree.right.ids);rl=cross(tree.right.ids,tree.left.ids)
    FollowupBlockNode(tree,nothing,l,r,lr,rl,l.entries+r.entries+lr.entries+rl.entries)
end

function di_generate(S,tree;stats=Dict{String,Int}())
    S.fac isa TreeKKTFactor || error("DIRECT_UNAVAILABLE non-tree network factor")
    lay=S.lay;stat=S.stat;b0=S.b0;n=length(S.active);nz=S.nzc
    @assert length(tree.ids)==n
    fa=S.fu_active
    @assert all(k->fa.colptr[k+1]==fa.colptr[k]+1&&rowvals(fa)[fa.colptr[k]]==k,1:n)
    f=[nonzeros(fa)[fa.colptr[k]] for k in 1:n];p=1 ./ f
    # Recreate only LOCAL elimination blocks from the already available network
    # factors. Never allocate the global dense Sh0 or Maa/Mab/Mbb.
    I_=Int[];J_=Int[];V_=Float64[];offset=0
    for gi in eachindex(stat.gkeys)
        rg=stat.granges[gi];ng=length(rg);Dg=zeros(ng,ng)
        stats["max_local_elimination_dim"]=max(get(stats,"max_local_elimination_dim",0),ng)
        for (b,a,pos) in stat.kee[gi];Dg[b,a]=nz[pos];end
        if stat.gkeys[gi] != -1
            Bf=S.Bd[stat.gkeys[gi]];ge=stat.gkef[gi]
            for (a,aa,pa) in ge,(b,bb,pb) in ge;Dg[a,b]-=nz[pa]*Bf[aa,bb]*nz[pb];end
        end
        for (a,b,pa,pb,si) in stat.giso[gi];Dg[a,b]-=nz[pa]*nz[pb]/nz[stat.iso_d[si]];end
        np=b0.nP[gi];blk=Dg[1:np,1:np]-Matrix(b0.PN[gi]*(Diagonal(b0.Ninv[gi])*b0.NP[gi]))
        ids=b0.bat[offset+1:offset+np]
        for j in 1:np,i in 1:np
            z=p[ids[i]]*blk[i,j]*p[ids[j]]
            if z!=0;push!(I_,ids[i]);push!(J_,ids[j]);push!(V_,z);end
        end
        offset+=np
    end
    @assert offset==n
    localbb=sparse(I_,J_,V_,n,n)
    up=S.Up[b0.prow,:];response=S.fac.Mroot\up'
    U=zeros(n,size(up,2));V=zeros(n,size(up,2))
    U[b0.bat,:]=-p[b0.bat].*up
    V[b0.bat,:]=p[b0.bat].*response'
    stats["separator_width"]=size(up,2)
    # Same original-coordinate gradient and sparse coupling algebra as _pp_elem.
    nE=length(lay.E);nN=nE-n;rE=copy(S.b[lay.E])
    for (a,t,pos) in stat.kef;rE[a]-=nz[pos]*S.y1[lay.Fc[t]];end
    for (si,i) in enumerate(lay.iso),(a,pos) in stat.iso_e[si]
        rE[a]-=nz[pos]*(S.b[i]/nz[stat.iso_d[si]])
    end
    ninv=zeros(nN);ii=Int[];jj=Int[];vv=Float64[]
    for (g,rg) in enumerate(stat.granges)
        np=b0.nP[g];PN=b0.PN[g];rv=rowvals(PN);nzv=nonzeros(PN)
        for l in 1:length(rg)-np
            e=stat.perm[rg[np+l]]-n;ninv[e]=b0.Ninv[g][l]
            for q in nzrange(PN,l);push!(ii,stat.perm[rg[rv[q]]]);push!(jj,e);push!(vv,nzv[q]);end
        end
    end
    SPN=sparse(ii,jj,vv,n,nN);Nd=Diagonal(ninv)
    rhat=rE[1:n]-SPN*(ninv.*rE[n+1:end])
    Rhat=sparse(-sparse(S.lux)[S.active,:]+SPN*(Nd*S.cxE))
    D=sparse(sparse(S.lxx)-S.cxE'*(Nd*S.cxE))
    dvec=S.lx+S.cxϕ+S.cxE'*(ninv.*rE[n+1:end])
    R=sparse(Rhat'*Diagonal(p));RT=sparse(Diagonal(p)*Rhat)
    ps=[sparse(localbb+R+RT+D),sparse(-localbb-R),localbb]
    reps=[di_payload(ps[1],U,V,tree,stats),di_payload(ps[2],-U,V,tree,stats),di_payload(ps[3],U,V,tree,stats)]
    ma=p.*rhat+dvec;mb=-p.*rhat
    (reps=reps,ma=ma,mb=mb,reports=Any[],bytes=8*(sum(x.entries for x in reps)+2n))
end

function di_consume(c,p1,packet,mu;stats=Dict{String,Int}(),corrupt=false)
    n=length(c.value.v);v=copy(c.value.v)
    corrupt&&(v .+=1e6sin.(collect(1:n)))
    action=x->cc_value(c.value,x)
    ds=[cc_stage(packet,v,action,mu,pr.dx;stats=stats) for pr in packet.probes]
    B=ck_upper(p1.reps[3],1e-12,stats)
    diagonal_probe=ck_apply(B,ones(n))
    sol=cc_cg(x->ck_apply(B,x)+action(x),-(p1.mb+v),diagonal_probe)
    # The diagonal approximation above is a row-sum probe, not a certified
    # preconditioner/condition bound. Actual residual and direction gates apply.
    stats["battery_cg_iterations"]=sol.iterations;x1=sol.x
    rhs=copy(c.value.gradient);rhs[1:n] .+=ck_apply(c.value.chain.left.reps[2],x1;trans=true)
    rest=-ck_solve(c.value.factor,rhs)
    xs=[zeros(n),x1,rest[1:n],rest[n+1:2n],rest[2n+1:3n]]
    path=[(dx=xs[t],dp=(xs[t+1]-xs[t])./packet.f) for t in 1:4]
    (directions=ds,path=path)
end

function di_prefix(p1,c,P;stats=Dict{String,Int}())
    n=length(p1.ma);aa=ck_upper(p1.reps[1],1e-12,stats);bb=ck_upper(p1.reps[3],1e-12,stats)
    B=p1.reps[2];k=x->ck_apply(bb,x)+cc_value(c.value,x);d=ck_apply(bb,ones(n))
    solve=b->begin
        z=cc_cg(k,b,d);stats["prefix_cg_maxiter"]=max(get(stats,"prefix_cg_maxiter",0),z.iterations);z.x
    end
    v=p1.ma-ck_apply(B,solve(p1.mb+c.value.v))
    Y=hcat([ck_apply(aa,P[:,j])-ck_apply(B,solve(ck_apply(B,P[:,j];trans=true))) for j in axes(P,2)]...)
    (v=v,Y=Y)
end
function di_dense_prefix(p1,c,P)
    K=cholesky(Symmetric(p1.Mbb+c.value.V,:U))
    (v=p1.ma-p1.Mab*(K\(p1.mb+c.value.v)),Y=Symmetric(p1.Maa,:U)*P-p1.Mab*(K\(p1.Mab'*P)))
end
