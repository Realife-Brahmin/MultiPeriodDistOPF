# First controlled representation test. Standard libraries; no optimizer/model.
using Serialization, LinearAlgebra, SparseArrays, SHA, Statistics, Printf
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"../../../.."))
include(joinpath(ROOT,"spar_ddp/examples/power_system/pint_scan.jl"))
length(ARGS)==4 || error("Usage: replay.jl CAPTURE INPUT_JLS OUTPUT N_PERIODS")
const CAP=abspath(ARGS[1]); const OUT=abspath(ARGS[3]); mkpath(OUT)
const N=parse(Int,ARGS[4]); @assert N in (4,8)
const data=deserialize(ARGS[2]); const n=length(data[:Bset]); @assert n==1020
function rootbranch(b)
    seen=Set{Int}()
    while data[:parent][b]!=data[:substationBus]
        b in seen && error("cyclic parent map");push!(seen,b);b=data[:parent][b]
    end
    b
end
const gids=rootbranch.(data[:Bset]); const gkeys=sort(unique(gids))
const groups=[findall(==(g),gids) for g in gkeys]
@assert sort(vcat(groups...))==collect(1:n)
const left=vcat(groups[1:length(groups)÷2]...)
const right=setdiff(1:n,left)
const manifest_path=isfile(joinpath(CAP,"capture.sha256")) ? joinpath(CAP,"capture.sha256") : joinpath(dirname(CAP),"capture.sha256")
const manifest=Dict(basename(split(l)[end])=>split(l)[1] for l in readlines(manifest_path))
for t in 1:N
    f="stage_$t.jls";@assert bytes2hex(sha256(read(joinpath(CAP,f))))==manifest[f]
end
const elems=[deserialize(joinpath(CAP,"stage_$t.jls")) for t in 1:N]
@assert all(e->size(e.Maa)==size(e.Mab)==size(e.Mbb)==(n,n) && length(e.ma)==length(e.mb)==n,elems)
println("REPLAY_METADATA nB=$n periods=$N groups=$(length(groups)) Julia=$(VERSION) BLAS=$(BLAS.get_num_threads()) first_optimization_pass_only=true original_input_sha=$(bytes2hex(sha256(read(ARGS[2]))))")

# This path REMOVES cross-group quadratic terms. It is an approximation, not
# physical aggregation and not an exact structured solver. Every state remains.
function block_project(e)
    mats=[zeros(n,n) for _ in 1:3]
    for g in groups, (i,f) in enumerate((:Maa,:Mab,:Mbb));mats[i][g,g]=getproperty(e,f)[g,g];end
    PintElem(mats...,copy(e.ma),copy(e.mb))
end
function block_combine(x,y)
    mats=[zeros(n,n) for _ in 1:3];a=zeros(n);b=zeros(n)
    for g in groups
        sub(e)=PintElem(e.Maa[g,g],e.Mab[g,g],e.Mbb[g,g],e.ma[g],e.mb[g])
        z=pint_combine(sub(x),sub(y))
        for (i,f) in enumerate((:Maa,:Mab,:Mbb));mats[i][g,g]=getproperty(z,f);end
        a[g]=z.ma;b[g]=z.mb
    end
    PintElem(mats...,a,b)
end
function errors(x,y)
    maximum(norm(getproperty(x,f)-getproperty(y,f),Inf)/max(1,norm(getproperty(y,f),Inf)) for f in (:Maa,:Mab,:Mbb,:ma,:mb))
end

open(joinpath(OUT,"spectra.csv"),"w") do io
    println(io,"depth,interval,block,singular_index,singular_value,relative_to_cross_sigma1,relative_to_full_frobenius")
    open(joinpath(OUT,"rounds.csv"),"w") do rr
        println(rr,"depth,interval,block_diagonal_approx_relative_error,dense_s_median,block_approx_s_median,dense_alloc_bytes,block_alloc_bytes,dense_gc_s,block_gc_s,diagonal_block_payload_bytes,dense_payload_bytes")
        exact=deepcopy(elems); approx=block_project.(elems); intervals=[(i-1,i) for i in 1:N]
        depth=0
        while true
            for i in eachindex(exact)
                x=exact[i];label="$(intervals[i][1]):$(intervals[i][2])"
                for f in (:Maa,:Mab,:Mbb)
                    A=getproperty(x,f);s=svdvals(A[left,right]);scale=max(first(s),floatmin(Float64));whole=max(norm(A),floatmin(Float64))
                    for (j,v) in enumerate(s);println(io,join((depth,label,f,j,v,v/scale,v/whole),','));end
                end
            end
            length(exact)==1 && break
            nextx=PintElem[];nexta=PintElem[];nextinterval=Tuple{Int,Int}[]
            for i in 1:2:length(exact)
                x,y=exact[i:i+1];a,b=approx[i:i+1]
                # Both paths warmed; no model/JIT setup counted as optimizer time.
                ref=pint_combine(x,y);test=block_combine(a,b)
                # Two different association orders below also checked separately.
                dense=Float64[];block=Float64[];db=Int[];bb=Int[];dg=Float64[];bg=Float64[]
                for rep in 1:3
                    GC.gc()
                    if isodd(rep)
                        u=@timed pint_combine(x,y);v=@timed block_combine(a,b)
                    else
                        v=@timed block_combine(a,b);u=@timed pint_combine(x,y)
                    end
                    push!(dense,u.time);push!(block,v.time);push!(db,u.bytes);push!(bb,v.bytes);push!(dg,u.gctime);push!(bg,v.gctime)
                end
                label="$(intervals[i][1]):$(intervals[i+1][2])"
                err=errors(test,ref);@assert isfinite(err)
                stored=8*(3sum(length(g)^2 for g in groups)+2n)
                println(rr,join((depth+1,label,err,median(dense),median(block),median(db),median(bb),median(dg),median(bg),stored,8*(3n*n+2n)),','))
                println("REPLAY_COMPOSITION depth=$(depth+1) interval=$label approx_error=$err dense_s=$(median(dense)) block_approx_s=$(median(block))")
                push!(nextx,ref);push!(nexta,test);push!(nextinterval,(intervals[i][1],intervals[i+1][2]))
            end
            exact=nextx;approx=nexta;intervals=nextinterval;depth+=1
        end
        # Full dense left fold versus balanced association is the exactness gate.
        fold=reduce(pint_combine,elems);err=errors(exact[1],fold)
        @assert err<1e-7 "dense association regression: $err"
        # Check closing to a zero terminal quadratic, retaining all nB states.
        zero=PintVal(zeros(n,n),zeros(n));a=pint_combine(approx[1],zero);d=pint_combine(exact[1],zero)
        ev=norm(a.V-d.V,Inf)/max(1,norm(d.V,Inf));eg=norm(a.v-d.v,Inf)/max(1,norm(d.v,Inf))
        @assert isfinite(ev) && isfinite(eg)
        println("REPLAY_CHECK dense_association_error=$err approximate_closed_hessian_error=$ev approximate_closed_gradient_error=$eg chol_fallbacks=$(_PINT_CHOL_FAIL[]) maxrss_bytes=$(Sys.maxrss())")
    end
end
println("REPLAY_COMPLETE diagnostic_only=true no_feasibility_or_full_solve_claim=true")
