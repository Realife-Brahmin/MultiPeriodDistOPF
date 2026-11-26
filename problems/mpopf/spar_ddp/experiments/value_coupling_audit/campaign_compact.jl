# Complete four-period frozen workload. Reuses known stage factors; no optimizer
# hook. A growing retained interior is explicit, not fixed-rank Schur closure.
function campaign_tail(raw,first=5)
    # Keep GLOBAL production round association while skipping unneeded prefix
    # left operands. An isolated reindexed tail can round differently.
    n=length(raw[1].ma);N=length(raw)
    xs=Any[t<first ? nothing : deepcopy(raw[t]) for t in 1:N]
    push!(xs,PintVal(zeros(n,n),zeros(n)));w=PairWork(n);h=1
    while h<length(xs)
        for i in 1:2h:length(xs)
            mid=i+h;mid>length(xs)&&continue
            for t in max(i,first):mid-1;xs[t]=_pp_combine!(xs[t],xs[mid],w);end
        end
        h*=2
    end
    @assert xs[first] isa PintVal;xs[first]
end
struct CampaignValue
    aa::FollowupBlockNode
    B::FollowupBlockNode
    H::FollowupBlockNode
    factor::CKFactor
    g::Vector{Float64}
    v::Vector{Float64}
    n::Int
end
function campaign_suffix(ps,terminal,stats)
    m=length(ps);n=length(first(ps).ma)
    ds=FollowupBlockNode[]
    for j in 1:m
        next=j==m ? terminal.rep : ps[j+1].reps[1]
        push!(ds,ck_upper(ck_add(ps[j].reps[3],next,1e-12,stats),1e-12,stats))
    end
    bs=[ps[j+1].reps[2] for j in 1:m-1]
    H=cc_lift(ds,bs,n,1e-12,stats);F=cc_factor(H,stats)
    g=vcat([ps[j].mb+(j==m ? terminal.v : ps[j+1].ma) for j in 1:m]...)
    z=ck_solve(F,g);aa=ck_upper(ps[1].reps[1],1e-12,stats);B=ps[1].reps[2]
    CampaignValue(aa,B,H,F,g,ps[1].ma-ck_apply(B,z[1:n]),n)
end
function campaign_action(c::CampaignValue,X)
    x=reshape(X,c.n,:);h=length(c.g)-c.n
    z=ck_solve(c.factor,vcat(ck_apply(c.B,x;trans=true),zeros(h,size(x,2))))
    y=ck_apply(c.aa,x)-ck_apply(c.B,z[1:c.n,:])
    X isa AbstractVector ? vec(y) : y
end
struct CampaignPrefix
    aa::FollowupBlockNode
    B::FollowupBlockNode
    bb::FollowupBlockNode
    future::CampaignValue
    diagonal_probe::Vector{Float64}
    first_exit::Vector{Float64}
    v::Vector{Float64}
end
campaign_prefix_K(p,x)=ck_apply(p.bb,x)+campaign_action(p.future,x)
function campaign_prefix(p,future,stats)
    aa=ck_upper(p.reps[1],1e-12,stats);bb=ck_upper(p.reps[3],1e-12,stats)
    d=ck_apply(bb,ones(length(p.ma)))
    z=cc_cg(x->ck_apply(bb,x)+campaign_action(future,x),p.mb+future.v,d)
    stats["prefix_gradient_cg_iterations"]=z.iterations
    CampaignPrefix(aa,p.reps[2],bb,future,d,-z.x,p.ma-ck_apply(p.reps[2],z.x))
end
function campaign_action(p::CampaignPrefix,X)
    n=length(p.v);xx=reshape(X,n,:)
    yy=hcat([begin
        sol=cc_cg(x->campaign_prefix_K(p,x),ck_apply(p.B,xx[:,j];trans=true),p.diagonal_probe)
        ck_apply(p.aa,xx[:,j])-ck_apply(p.B,sol.x)
    end for j in axes(xx,2)]...)
    X isa AbstractVector ? vec(yy) : yy
end
function campaign_build(ps,terminal;stats=Dict{String,Int}())
    @assert length(ps)==4
    vals=Any[nothing]
    for j in 2:4;push!(vals,campaign_suffix(ps[j:4],terminal,stats));end
    vals[1]=campaign_prefix(ps[1],vals[2],stats);vals
end
function campaign_conditions(values)
    # Charged uncertified solve/norm probes for every factored suffix.
    [begin
        n=length(c.g);x=fill(1/n,n);iv=0.0
        for _ in 1:4
            y=ck_solve(c.factor,x);iv=max(iv,norm(y,1))
            z=ck_solve(c.factor,sign.(y));x=zeros(n);x[argmax(abs.(z))]=1
        end
        iv*norm(ck_apply(c.H,ones(n)),Inf)
    end for c in values[2:4]]
end
campaign_direct_inputs(packets,tree)=[di_generate(p.S,tree) for p in packets]
campaign_dense_inputs(packets,n)=[_pp_elem(p.S,p.stage,n) for p in packets]
campaign_boundary(terminal,tree)=(rep=block_payload(terminal.V,tree,1e-9,Any[],:terminal),v=copy(terminal.v))
campaign_value_actions(values,P)=[campaign_action(c,P) for c in values]
campaign_dense_actions(values,P)=[c.V*P for c in values[1:4]]
function campaign_consume(packets,values,terminal,mu;corrupt=false)
    all_dirs=Any[]
    for (t,packet) in enumerate(packets)
        future=t==4 ? terminal : values[t+1]
        vx=copy(future.v);corrupt&&t==1&&(vx .+=1e6sin.(collect(1:length(vx))))
        action=t==4 ? (x->ck_apply(future.rep,x)) : (x->campaign_action(future,x))
        push!(all_dirs,[cc_stage(packet,vx,action,mu,pr.dx) for pr in packet.probes])
    end
    n=length(values[1].v);x1=values[1].first_exit;future=values[2]
    rhs=copy(future.g);rhs[1:n] .+=ck_apply(future.B,x1;trans=true)
    rest=-ck_solve(future.factor,rhs)
    xs=[zeros(n),x1,rest[1:n],rest[n+1:2n],rest[2n+1:3n]]
    (directions=all_dirs,path=[(dx=xs[t],dp=(xs[t+1]-xs[t])./packets[t].f) for t in 1:4])
end
function campaign_dense_consume(packets,values,elems,mu)
    ds=Any[]
    for (t,packet) in enumerate(packets)
        rule=policy(packet,values[t+1],mu)
        push!(ds,[directions(rule,pr.dx) for pr in packet.probes])
    end
    (directions=ds,path=battery_path(elems,values,packets[1].f))
end
function campaign_compare(io,label,observed,packets,rules,pathref)
    worst=0.0;relative=0.0
    function record(t,probe,field,A,B)
        local ae=norm(A-B);local rn=norm(B)
        local se=ae/max(1.0,rn);local re=ae/max(1e-10,rn)
        @assert isfinite(se)&&isfinite(re)
        emit(io,label,t,probe,field,ae,rn,se,re);worst=max(worst,se);relative=max(relative,re)
    end
    for (t,packet) in enumerate(packets),(j,pr) in enumerate(packet.probes)
        d=directions(rules[t],pr.dx)
        for f in (:du,:dphi,:dzl,:dzu);record(t,j,f,getproperty(observed.directions[t][j],f),getproperty(d,f));end
    end
    for t in 1:4,f in (:dx,:dp);record(t,"path",f,getproperty(observed.path[t],f),getproperty(pathref[t],f));end
    (scaled=worst,relative=relative,pass=worst<DIR_SCALED&&relative<DIR_RELATIVE)
end
