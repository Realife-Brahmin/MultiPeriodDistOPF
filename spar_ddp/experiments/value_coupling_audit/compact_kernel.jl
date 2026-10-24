# Isolated HODLR/recursive Woodbury component prototype. No production hook.
# Input matrices remain raw; upper completion is explicit only where production
# uses Symmetric(K,:U) or copytri! on the Schur diagonal blocks.
const CK_LOCAL_DENSE_CAP=64
ck_bump!(s,k,v=1)=(s[k]=get(s,k,0)+v)
function ck_factors(o::FollowupOffBlock)
    o.kind==:lowrank && return (o.U,o.V)
    m,n=size(o.A)
    m<=n ? (Matrix{Float64}(I,m,m),copy(o.A')) : (copy(o.A),Matrix{Float64}(I,n,n))
end
function ck_off(U,V,tol,stats)
    m,r=size(U);n=size(V,1);@assert size(V,2)==r
    r==0 && return FollowupOffBlock(:lowrank,nothing,U,V,0)
    # Recompression uses factors' QR and a small core, never a reconstructed
    # full input matrix. Rank can grow and is recorded, not fixed by horizon.
    a=qr(U);b=qr(V);ra=Matrix(a.R);rb=Matrix(b.R)
    core=ra*rb';ck_bump!(stats,"core_svd_calls")
    stats["max_core_svd_dim"]=max(get(stats,"max_core_svd_dim",0),maximum(size(core)))
    f=svd(core;full=false);k=length(f.S);tail=0.0;budget=tol*norm(f.S)
    while k>0 && tail+f.S[k]^2<=budget^2;tail+=f.S[k]^2;k-=1;end
    stats["max_relative_recompression_tail_times_1e15"]=max(get(stats,"max_relative_recompression_tail_times_1e15",0),round(Int,1e15sqrt(tail)/max(norm(f.S),floatmin(Float64))))
    cu=Matrix(a.Q)*(f.U[:,1:k].*f.S[1:k]');cv=Matrix(b.Q)*Matrix(f.Vt[1:k,:]')
    stats["max_retained_rank"]=max(get(stats,"max_retained_rank",0),k)
    if k*(m+n)>=m*n
        ck_bump!(stats,"dense_cross_fallbacks")
        return FollowupOffBlock(:dense,cu*cv',nothing,nothing,m*n)
    end
    FollowupOffBlock(:lowrank,nothing,cu,cv,k*(m+n))
end
function ck_add(a::FollowupBlockNode,b::FollowupBlockNode,tol,stats)
    @assert a.tree.ids==b.tree.ids
    if !isnothing(a.localblock)
        c=a.localblock+b.localblock
        return FollowupBlockNode(a.tree,c,nothing,nothing,nothing,nothing,length(c))
    end
    l=ck_add(a.left,b.left,tol,stats);r=ck_add(a.right,b.right,tol,stats)
    au,av=ck_factors(a.lr);bu,bv=ck_factors(b.lr)
    lr=ck_off(hcat(au,bu),hcat(av,bv),tol,stats)
    au,av=ck_factors(a.rl);bu,bv=ck_factors(b.rl)
    rl=ck_off(hcat(au,bu),hcat(av,bv),tol,stats)
    FollowupBlockNode(a.tree,nothing,l,r,lr,rl,l.entries+r.entries+lr.entries+rl.entries)
end
ck_transpose(o::FollowupOffBlock)=o.kind==:dense ? FollowupOffBlock(:dense,copy(o.A'),nothing,nothing,o.entries) : FollowupOffBlock(:lowrank,nothing,o.V,o.U,o.entries)
function ck_upper_cross(a,L,R,tol,stats)
    maximum(L)<minimum(R) && return a.lr
    maximum(R)<minimum(L) && return ck_transpose(a.rl)
    # Interleaved original indices: mask factors by runs, not a full dense
    # matrix. This preserves the ORIGINAL upper-triangle convention under the
    # topology permutation. Separator/interface rank may consequently grow.
    ck_bump!(stats,"interleaved_upper_nodes")
    u,v=ck_factors(a.lr);w,z=ck_factors(a.rl)
    q=[count(<(i),R) for i in L];us=Matrix{Float64}[];vs=Matrix{Float64}[]
    for bucket in unique(q)
        rows=findall(==(bucket),q);pivot=L[first(rows)]
        high=findall(>(pivot),R);low=findall(<(pivot),R)
        if !isempty(high)
            x=zeros(length(L),size(u,2));y=zeros(length(R),size(v,2))
            x[rows,:]=u[rows,:];y[high,:]=v[high,:];push!(us,x);push!(vs,y)
        end
        if !isempty(low)
            x=zeros(length(L),size(z,2));y=zeros(length(R),size(w,2))
            x[rows,:]=z[rows,:];y[low,:]=w[low,:];push!(us,x);push!(vs,y)
        end
    end
    ck_off(hcat(us...),hcat(vs...),tol,stats)
end
function ck_upper(a::FollowupBlockNode,tol,stats)
    if !isnothing(a.localblock)
        @assert issorted(a.tree.ids)
        c=Matrix(Symmetric(a.localblock,:U))
        return FollowupBlockNode(a.tree,c,nothing,nothing,nothing,nothing,length(c))
    end
    l=ck_upper(a.left,tol,stats);r=ck_upper(a.right,tol,stats)
    off=ck_upper_cross(a,a.left.tree.ids,a.right.tree.ids,tol,stats)
    FollowupBlockNode(a.tree,nothing,l,r,off,ck_transpose(off),l.entries+r.entries+2off.entries)
end
function ck_apply_ordered(a::FollowupBlockNode,x;trans=false)
    !isnothing(a.localblock) && return (trans ? a.localblock' : a.localblock)*x
    k=length(a.left.tree.ids);xl=view(x,1:k,:);xr=view(x,k+1:size(x,1),:)
    top=ck_apply_ordered(a.left,xl;trans=trans);bot=ck_apply_ordered(a.right,xr;trans=trans)
    if trans
        u,v=ck_factors(a.rl);top+=v*(u'*xr)
        u,v=ck_factors(a.lr);bot+=v*(u'*xl)
    else
        u,v=ck_factors(a.lr);top+=u*(v'*xr)
        u,v=ck_factors(a.rl);bot+=u*(v'*xl)
    end
    vcat(top,bot)
end
function ck_apply(a,x;trans=false)
    xx=reshape(x,size(x,1),:);p=a.tree.ids
    y=ck_apply_ordered(a,xx[p,:];trans=trans)[invperm(p),:]
    x isa AbstractVector ? vec(y) : y
end
mutable struct CKFactor
    ids::Vector{Int}
    left::Any
    right::Any
    localfactor::Any
    corefactor::Any
    AinvU::Matrix{Float64}
    BinvW::Matrix{Float64}
    V::Matrix{Float64}
    Z::Matrix{Float64}
end
function ck_solve_ordered(f::CKFactor,b)
    !isnothing(f.localfactor) && return f.localfactor\b
    k=length(f.left.ids);a=ck_solve_ordered(f.left,view(b,1:k,:));c=ck_solve_ordered(f.right,view(b,k+1:size(b,1),:))
    r=size(f.AinvU,2);z=f.corefactor\vcat(f.V'*c,f.Z'*a)
    vcat(a-f.AinvU*view(z,1:r,:),c-f.BinvW*view(z,r+1:size(z,1),:))
end
function ck_solve(f,b)
    xx=reshape(b,size(b,1),:);p=f.ids
    y=ck_solve_ordered(f,xx[p,:])[invperm(p),:]
    b isa AbstractVector ? vec(y) : y
end
function ck_factor(a::FollowupBlockNode,stats)
    n=length(a.tree.ids);empty=zeros(0,0)
    if n<=CK_LOCAL_DENSE_CAP
        # Only bounded local dense fallbacks, never a full nB-by-nB matrix.
        A=ck_apply_ordered(a,Matrix{Float64}(I,n,n))
        ck_bump!(stats,"local_dense_factors");stats["max_local_dense_dim"]=max(get(stats,"max_local_dense_dim",0),n)
        return CKFactor(a.tree.ids,nothing,nothing,lu(A),nothing,empty,empty,empty,empty)
    end
    lf=ck_factor(a.left,stats);rf=ck_factor(a.right,stats)
    u,v=ck_factors(a.lr);w,z=ck_factors(a.rl);r=size(u,2);s=size(w,2)
    stats["max_woodbury_core_dim"]=max(get(stats,"max_woodbury_core_dim",0),r+s)
    r+s>=n && error("COMPACT_RANK_GROWTH node=$n core=$(r+s); dense fallback required")
    au=ck_solve_ordered(lf,u);bw=ck_solve_ordered(rf,w)
    core=[Matrix{Float64}(I,r,r) v'*bw;z'*au Matrix{Float64}(I,s,s)]
    CKFactor(a.tree.ids,lf,rf,nothing,lu(core),au,bw,v,z)
end
struct CKSchur
    x::Any
    y::Any
    aa::FollowupBlockNode
    bb::FollowupBlockNode
    K::FollowupBlockNode
    factor::CKFactor
    ma::Vector{Float64}
    mb::Vector{Float64}
end
function ck_merge(x,y;tol=1e-12,stats=Dict{String,Int}())
    K=ck_upper(ck_add(x.reps[3],y.reps[1],tol,stats),tol,stats)
    fac=ck_factor(K,stats)
    aa=ck_upper(x.reps[1],tol,stats);bb=ck_upper(y.reps[3],tol,stats)
    solved=ck_solve(fac,x.mb+y.ma)
    CKSchur(x,y,aa,bb,K,fac,x.ma-ck_apply(x.reps[2],solved),y.mb-ck_apply(y.reps[2],solved;trans=true))
end
function ck_action(c::CKSchur,block,x)
    B=c.x.reps[2];C=c.y.reps[2]
    if block==:Maa
        ck_apply(c.aa,x)-ck_apply(B,ck_solve(c.factor,ck_apply(B,x;trans=true)))
    elseif block==:Mab
        -ck_apply(B,ck_solve(c.factor,ck_apply(C,x)))
    elseif block==:Mbb
        ck_apply(c.bb,x)-ck_apply(C,ck_solve(c.factor,ck_apply(C,x));trans=true)
    else
        error("unknown Schur block")
    end
end
function ck_condition_estimate(c;steps=4)
    # Hager-style probe estimate, not a certified condition bound. Reuses the
    # same compact factor; no dense inverse or dense K construction.
    n=length(c.ma);x=fill(1/n,n);invnorm=0.0
    for _ in 1:steps
        y=ck_solve(c.factor,x);invnorm=max(invnorm,norm(y,1))
        z=ck_solve(c.factor,sign.(y));x=zeros(n);x[argmax(abs.(z))]=1
    end
    knorm=norm(ck_apply(c.K,ones(n)),Inf)
    (probe_condition=knorm*invnorm,inverse_one_norm_estimate=invnorm,
     norm_probe=knorm,certified=false)
end
function ck_export(c;batch=16)
    # EXPLICIT legacy-consumer/validation boundary, NOT part of compact kernel.
    # Never send a full identity through the operator inside compact execution.
    n=length(c.ma);mats=Matrix{Float64}[]
    for block in (:Maa,:Mab,:Mbb)
        A=zeros(n,n)
        for j in 1:batch:n
            ix=j:min(j+batch-1,n);E=zeros(n,length(ix))
            for (k,i) in enumerate(ix);E[i,k]=1;end
            A[:,ix]=ck_action(c,block,E)
        end
        push!(mats,A)
    end
    PintElem(mats...,copy(c.ma),copy(c.mb))
end
