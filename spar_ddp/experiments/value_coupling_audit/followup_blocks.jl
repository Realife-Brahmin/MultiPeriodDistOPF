# Opt-in frozen-data experiment only. Raw local blocks and vectors are exact.
# Unlike the joint symmetric prototype, no raw skew is projected away.
function block_payload(A,tree,tol,reports,field,depth=0)
    I=tree.ids
    if isnothing(tree.left)
        B=copy(A[I,I])
        return (tree=tree,localblock=B,left=nothing,right=nothing,lr=nothing,rl=nothing,entries=length(B))
    end
    a=block_payload(A,tree.left,tol,reports,field,depth+1)
    b=block_payload(A,tree.right,tol,reports,field,depth+1)
    function cross(L,R,side)
        B=copy(A[L,R]);scale=norm(B)
        if tol==0 || scale==0
            push!(reports,(field=field,depth=depth,side=side,rows=size(B,1),cols=size(B,2),rank=min(size(B)...),dense=true,tail=0.0,scale=scale))
            return (kind=:dense,A=B,U=nothing,V=nothing,entries=length(B))
        end
        F=svd(B;full=false);r=length(F.S);tail=0.0
        while r>0 && tail+F.S[r]^2<=(tol*scale)^2
            tail+=F.S[r]^2;r-=1
        end
        dense=r*sum(size(B))>=length(B)
        push!(reports,(field=field,depth=depth,side=side,rows=size(B,1),cols=size(B,2),rank=r,dense=dense,tail=dense ? 0.0 : sqrt(tail),scale=scale))
        dense && return (kind=:dense,A=B,U=nothing,V=nothing,entries=length(B))
        (kind=:lowrank,A=nothing,U=F.U[:,1:r].*F.S[1:r]',V=Matrix(F.Vt[1:r,:]'),entries=r*sum(size(B)))
    end
    lr=cross(tree.left.ids,tree.right.ids,"LR");rl=cross(tree.right.ids,tree.left.ids,"RL")
    (tree=tree,localblock=nothing,left=a,right=b,lr=lr,rl=rl,entries=a.entries+b.entries+lr.entries+rl.entries)
end
function expand_block(p,n)
    A=zeros(n,n)
    function visit(x)
        if !isnothing(x.localblock);A[x.tree.ids,x.tree.ids]=x.localblock;return;end
        visit(x.left);visit(x.right);L=x.tree.left.ids;R=x.tree.right.ids
        A[L,R]=x.lr.kind==:dense ? x.lr.A : x.lr.U*x.lr.V'
        A[R,L]=x.rl.kind==:dense ? x.rl.A : x.rl.U*x.rl.V'
    end
    visit(p);A
end
function construct_blocks(e,tree,tol)
    @assert 0<=tol<1 && isfinite(tol)
    reports=Any[]
    reps=[block_payload(getproperty(e,f),tree,tol,reports,string(f)) for f in (:Maa,:Mab,:Mbb)]
    (reps=reps,ma=copy(e.ma),mb=copy(e.mb),reports=reports,
     bytes=8*(sum(p.entries for p in reps)+2length(e.ma)))
end
function materialize_blocks(p,raw,tol)
    n=length(raw.ma);mats=[expand_block(q,n) for q in p.reps]
    err=maximum(norm(mats[i]-getproperty(raw,f))/max(norm(getproperty(raw,f)),floatmin(Float64)) for (i,f) in enumerate((:Maa,:Mab,:Mbb)))
    @assert isfinite(err)
    # Rounding and raw asymmetry are never hidden by symmetrizing the input.
    fallback=err>tol+32eps(Float64)
    out=fallback ? deepcopy(raw) : PintElem(mats...,copy(p.ma),copy(p.mb))
    (elem=out,fallback=fallback,error=fallback ? 0.0 : err,
     bytes=fallback ? 8*(3n*n+2n) : p.bytes)
end
function elimination_condition(K)
    # Explicit Inf-norm condition of the ACTUAL upper-triangle completion.
    # Expensive prototype control, measured and charged, not a free oracle.
    c=cond(Matrix(Symmetric(K,:U)),Inf)
    isfinite(c) ? max(1.0,c) : Inf
end
controlled_tol(tol,kappa)=isfinite(kappa) ? tol/max(1.0,sqrt(kappa)) : 0.0
