# Opt-in research only. Never loaded without FILTERDDP_VALUE_GROUPING.
const _VG_GROUPS=Ref(Vector{Vector{Int}}())
const _VG_LABELS=Ref(Int[])
const _VG_WORK=Dict{Int,PairWork}()
function vg_install(data)
    mode=ENV["FILTERDDP_VALUE_GROUPING"];n=length(data[:Bset])
    groups=if mode=="chunks16"
        p=sortperm(data[:Bset]);[p[i:min(i+15,n)] for i in 1:16:n]
    elseif mode in ("root","dense")
        function top(b)
            b==data[:substationBus] && return b
            while data[:parent][b]!=data[:substationBus];b=data[:parent][b];end
            b
        end
        labels=top.(data[:Bset]);[findall(==(k),labels) for k in sort(unique(labels))]
    else
        error("Unknown diagnostic value grouping: $mode")
    end
    @assert sort(vcat(groups...))==collect(1:n)
    _VG_GROUPS[]=groups;labs=zeros(Int,n)
    empty!(_VG_WORK)
    for (i,g) in enumerate(groups);labs[g].=i;get!(_VG_WORK,length(g)) do;PairWork(length(g));end;end
    _VG_LABELS[]=labs
    println("VALUE_GROUP_DIAGNOSTIC mode=$mode states=$n groups=$(length(groups)) max_group=$(maximum(length,groups)) approximate=$(mode!="dense") dense_storage_retained=true")
end
function vg_zero_cross!(A)
    labs=_VG_LABELS[]
    for j in axes(A,2),i in axes(A,1);labs[i]==labs[j] || (A[i,j]=0.0);end
end
function vg_combine!(x::PintElem,y)
    for g in _VG_GROUPS[]
        a=PintElem(x.Maa[g,g],x.Mab[g,g],x.Mbb[g,g],x.ma[g],x.mb[g])
        b=y isa PintVal ? PintVal(y.V[g,g],y.v[g]) : PintElem(y.Maa[g,g],y.Mab[g,g],y.Mbb[g,g],y.ma[g],y.mb[g])
        z=_pp_combine!(a,b,_VG_WORK[length(g)])
        if z isa PintVal
            x.Maa[g,g]=z.V;x.ma[g]=z.v
        else
            for f in (:Maa,:Mab,:Mbb,:ma,:mb)
                A=getproperty(x,f);B=getproperty(z,f)
                ndims(A)==2 ? (A[g,g]=B) : (A[g]=B)
            end
        end
    end
    vg_zero_cross!(x.Maa)
    if y isa PintVal;return PintVal(x.Maa,x.ma);end
    vg_zero_cross!(x.Mab);vg_zero_cross!(x.Mbb)
    x
end
function vg_slot(a,n)
    p=pointer(sdata(a));blk=i->unsafe_wrap(Array,p+(i-1)*n*n*sizeof(Float64),(n,n))
    vec=i->unsafe_wrap(Array,p+(3n*n+(i-1)*n)*sizeof(Float64),(n,))
    PintElem(blk(1),blk(2),blk(3),vec(1),vec(2))
end
function vg_reference_before(shared,n,k,j,exits)
    k in (0,5,10,15) && j==0 || return nothing
    start=time_ns();N=length(shared);ref=Vector{PintVal}(undef,N);v=PintVal(zeros(n,n),zeros(n));w=PairWork(n)
    for t in N:-1:1
        e=vg_slot(shared[t],n)
        dest=get(ENV,"FILTERDDP_VALUE_REFERENCE_DIR","")
        if exits>0 && !isempty(dest) && t in unique((1,cld(N,2),N))
            mkpath(dest);serialize(joinpath(dest,"leaf_t$(t)_k$(k)_j$(j).jls"),e)
        end
        x=PintElem(copy(e.Maa),copy(e.Mab),copy(e.Mbb),copy(e.ma),copy(e.mb))
        v=_pp_combine!(x,v,w);ref[t]=v
    end
    println("VALUE_REFERENCE_BUILD iteration=$k barrier_iteration=$j phase=$(exits>0 ? "timed" : "warm") wall_s=$((time_ns()-start)/1e9) diagnostic_overhead_not_benchmark=true")
    ref
end
function vg_reference_after(shared,n,ref,k,j)
    isnothing(ref) && return
    eV=0.0;ev=0.0
    for t in eachindex(shared)
        x=vg_slot(shared[t],n);r=ref[t]
        eV=max(eV,norm(x.Maa-r.V,Inf)/max(1,norm(r.V,Inf)))
        ev=max(ev,norm(x.ma-r.v,Inf)/max(1,norm(r.v,Inf)))
    end
    @assert isfinite(eV) && isfinite(ev)
    ENV["FILTERDDP_VALUE_GROUPING"]=="dense" && @assert max(eV,ev)<1e-7 "dense-reference gate"
    println("VALUE_REFERENCE_CHECK iteration=$k barrier_iteration=$j mode=$(ENV["FILTERDDP_VALUE_GROUPING"]) relative_Vxx=$eV relative_Vx=$ev")
end
