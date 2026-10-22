# Accuracy/storage prototype only; combination still materializes dense blocks.
# No production hook uses this file. Low rank is applied to cross-cluster
# blocks, never to the whole value Hessian or to individual battery variables.
struct CouplingTree
    ids::Vector{Int}
    left::Any
    right::Any
    separators::Vector{Int}
end
coupling_leaf(ids)=CouplingTree(sort(ids),nothing,nothing,Int[])
function topology_hierarchy(data;maxleaf=8)
    @assert maxleaf>=1
    rows=Dict(b=>i for (i,b) in enumerate(data[:Bset]))
    adj=Dict(v=>Int[] for v in data[:Nset])
    for (a,b) in data[:Lset];push!(adj[a],b);push!(adj[b],a);end
    for a in values(adj);sort!(a);end
    function joinchildren(children,separator)
        length(children)==1 && return only(children)
        total=sum(length(x.ids) for x in children);acc=0;split=1
        for i in 1:length(children)-1
            acc+=length(children[i].ids);split=i;acc>=total/2 && break
        end
        a=joinchildren(children[1:split],separator);b=joinchildren(children[split+1:end],separator)
        CouplingTree(vcat(a.ids,b.ids),a,b,[separator])
    end
    function build(nodes)
        ids=[rows[v] for v in nodes if haskey(rows,v)]
        length(ids)<=maxleaf && return coupling_leaf(ids)
        allowed=Set(nodes);start=minimum(nodes);parent=Dict(start=>0);order=[start];q=1
        while q<=length(order)
            v=order[q];q+=1
            for w in adj[v]
                if w in allowed && !haskey(parent,w);parent[w]=v;push!(order,w);end
            end
        end
        @assert length(order)==length(nodes) "component must be connected"
        mass=Dict(v=>Int(haskey(rows,v)) for v in order)
        for v in reverse(order);parent[v]!=0 && (mass[parent[v]]+=mass[v]);end
        total=length(ids)
        score(v)=maximum(vcat(total-mass[v],[mass[w] for w in adj[v] if get(parent,w,-1)==v]))
        sep=first(sort(order;by=v->(score(v),v)))
        visited=Set([sep]);children=CouplingTree[]
        haskey(rows,sep) && push!(children,coupling_leaf([rows[sep]]))
        for w in adj[sep]
            w in allowed && !(w in visited) || continue
            component=[w];push!(visited,w);i=1
            while i<=length(component)
                v=component[i];i+=1
                for u in adj[v]
                    if u in allowed && !(u in visited);push!(visited,u);push!(component,u);end
                end
            end
            any(v->haskey(rows,v),component) && push!(children,build(component))
        end
        @assert length(children)>1
        joinchildren(children,sep)
    end
    tree=build(collect(data[:Nset]));@assert sort(tree.ids)==collect(1:length(rows))
    tree
end
coupling_internals(t)=isnothing(t.left) ? 0 : 1+coupling_internals(t.left)+coupling_internals(t.right)
coupling_groups(t)=isnothing(t.left) ? [t.ids] : vcat(coupling_groups(t.left),coupling_groups(t.right))
joint_value(e)=[e.Maa e.Mab;e.Mab' e.Mbb]
function retain_coupling(e,tree;tol=1e-8)
    @assert 0<=tol<1
    start=time_ns();n=length(e.ma);original=joint_value(e)
    # Diagonal quadratic blocks are symmetric mathematically; include any
    # symmetrization discrepancy in the final reconstruction error, not Mab.
    H=(original+original')/2
    budget=tol*norm(H)/sqrt(max(1,2coupling_internals(tree)))
    reports=Any[]
    idx(t)=vcat(t.ids,n.+t.ids)
    function compress(t,depth)
        I=idx(t)
        isnothing(t.left) && return (tree=t,leaf=H[I,I],left=nothing,right=nothing,off=nothing,entries=length(I)^2)
        a=compress(t.left,depth+1);b=compress(t.right,depth+1)
        L=idx(t.left);R=idx(t.right);A=H[L,R];ts=time_ns();S=svd(A;full=false)
        r=length(S.S);tail2=0.0
        while r>0 && tail2+S.S[r]^2<=budget^2;tail2+=S.S[r]^2;r-=1;end
        dense=r*(length(L)+length(R))>=length(A)
        off=dense ? (kind=:dense,A=A) : (kind=:lowrank,U=S.U[:,1:r].*S.S[1:r]',V=Matrix(S.Vt[1:r,:]'))
        entries=dense ? length(A) : r*(length(L)+length(R))
        push!(reports,(depth=depth,rows=length(L),cols=length(R),rank=r,dense_fallback=dense,
                       cross_frobenius=norm(A),tail_frobenius=dense ? 0.0 : sqrt(tail2),svd_s=(time_ns()-ts)/1e9,
                       separator_buses=t.separators,singular_values=copy(S.S)))
        (tree=t,leaf=nothing,left=a,right=b,off=off,entries=a.entries+b.entries+entries)
    end
    representation=compress(tree,0);construction_s=(time_ns()-start)/1e9
    start=time_ns();out=zeros(2n,2n)
    function expand(x)
        I=idx(x.tree)
        if !isnothing(x.leaf);out[I,I]=x.leaf;return;end
        expand(x.left);expand(x.right);L=idx(x.tree.left);R=idx(x.tree.right)
        A=x.off.kind==:dense ? x.off.A : x.off.U*x.off.V'
        out[L,R]=A;out[R,L]=A'
    end
    expand(representation);materialization_s=(time_ns()-start)/1e9
    err=norm(out-original)/max(norm(original),floatmin(Float64))
    @assert isfinite(err) && err<=tol+1e-12 "global quadratic error budget"
    elem=PintElem(out[1:n,1:n],out[1:n,n+1:end],out[n+1:end,n+1:end],copy(e.ma),copy(e.mb))
    (elem=elem,representation=representation,reports=reports,relative_frobenius_error=err,
     numeric_payload_bytes=8*(representation.entries+2n),construction_s=construction_s,materialization_s=materialization_s)
end
