# Opt-in integrated grouped/right scan cache; gradients are never retained.
using LinearAlgebra,SparseArrays
mutable struct BRMatrixEntry
 U::Matrix{Float64}
 Y1::Matrix{Float64}
 Y2::Union{Nothing,Matrix{Float64}}
 V::Union{Nothing,Matrix{Float64}}
end
function br_store(old,A)
 old===nothing ? copy(A) : (copyto!(old,A);old)
end
# Same arithmetic, order, raw triangle, and fallback as study_right_kernel.jl.
function br_populate!(x,y,w,old=nothing)
 val=y isa PintVal
 w.K .= x.Mbb .+ (val ? y.V : y.Maa)
 w.g .= x.mb .+ (val ? y.v : y.ma)
 C=cholesky!(Symmetric(w.K);check=false)
 if !issuccess(C)
  return right_combine!(x,y,w),nothing
 end
 U=C.U
 copyto!(w.Y1,x.Mab)
 BLAS.trsm!('R','U','N','N',1.0,U.data,w.Y1)
 if !val;copyto!(w.Y2,y.Mab);ldiv!(U',w.Y2);end
 ldiv!(U',w.g)
 BLAS.syrk!('U','N',-1.0,w.Y1,1.0,x.Maa);LinearAlgebra.copytri!(x.Maa,'U')
 if !val
  copyto!(x.Mbb,y.Mbb)
  BLAS.syrk!('U','T',-1.0,w.Y2,1.0,x.Mbb);LinearAlgebra.copytri!(x.Mbb,'U')
  mul!(x.Mab,w.Y1,w.Y2,-1.0,0.0)
 end
 mul!(x.ma,w.Y1,w.g,-1.0,1.0)
 if !val;copyto!(x.mb,y.mb);mul!(x.mb,w.Y2',w.g,-1.0,1.0);end
 compatible=old!==nothing && (old.V!==nothing)==val && size(old.U)==size(x.Maa)
 c=compatible ? old : BRMatrixEntry(similar(x.Maa),similar(x.Maa),val ? nothing : similar(x.Maa),val ? similar(x.Maa) : nothing)
 copyto!(c.U,U);copyto!(c.Y1,w.Y1)
 val ? copyto!(c.V,x.Maa) : copyto!(c.Y2,w.Y2)
 return val ? PintVal(x.Maa,x.ma) : x,c
end
function br_hit!(x,y,w,c::BRMatrixEntry)
 val=y isa PintVal
 @assert val==(c.V!==nothing)
 w.g .= x.mb .+ (val ? y.v : y.ma)
 ldiv!(UpperTriangular(c.U)',w.g)
 mul!(x.ma,c.Y1,w.g,-1.0,1.0)
 if val
  copyto!(x.Maa,c.V)
 else
  copyto!(x.mb,y.mb);mul!(x.mb,c.Y2',w.g,-1.0,1.0)
 end
 val ? PintVal(x.Maa,x.ma) : x
end
br_bytes(c)=sum(x===nothing ? 0 : sizeof(x) for x in (c.U,c.Y1,c.Y2,c.V))
mutable struct BRGate
 epoch::UInt
 event::Symbol
 previous::Any
end
BRGate()=BRGate(0,:unknown,nothing)
function br_event!(g::BRGate,event::Symbol)
 g.event=event
 if event!==:barrier;g.epoch+=1;end
 event===:solve_start && (g.previous=nothing)
 nothing
end
function br_request!(g::BRGate,k,j,mu,context,supported)
 p=g.previous
 hit=supported && g.event===:barrier && p!==nothing && p.epoch==g.epoch &&
  p.k==k && p.j+1==j && mu<p.mu && p.context==context
 g.previous=(epoch=g.epoch,k=k,j=j,mu=mu,context=context)
 g.event=:unknown # authorization consumed once; repeated calls cannot hit
 hit
end
function br_leaf_vectors!(el,S)
 lay,stat,b0=S.lay,S.stat,S.b0;nB=length(S.active);nE=length(lay.E);nN=nE-nB
 f=[nonzeros(S.fu_active)[S.fu_active.colptr[k]] for k in 1:nB]
 rE=Vector{Float64}(undef,nE)
 for a in 1:nE;rE[a]=S.b[lay.E[a]];end
 for (a,t,p) in stat.kef;rE[a]-=S.nzc[p]*S.y1[lay.Fc[t]];end
 for (si,i) in enumerate(lay.iso)
  z=S.b[i]/S.nzc[stat.iso_d[si]]
  for (a,p) in stat.iso_e[si];rE[a]-=S.nzc[p]*z;end
 end
 ninv=zeros(nN);I_=Int[];J_=Int[];V_=Float64[]
 for (g,rg) in enumerate(stat.granges)
  npg=b0.nP[g];PNg=b0.PN[g];rv=rowvals(PNg);nzv=nonzeros(PNg)
  for l in 1:length(rg)-npg
   e=stat.perm[rg[npg+l]]-nB;ninv[e]=b0.Ninv[g][l]
   for q in nzrange(PNg,l);push!(I_,stat.perm[rg[rv[q]]]);push!(J_,e);push!(V_,nzv[q]);end
  end
 end
 SPN=sparse(I_,J_,V_,nB,nN);r0P=rE[1:nB];r0N=rE[nB+1:end]
 rhat=r0P-SPN*(ninv.*r0N);dvec=S.lx+getproperty(S,Symbol("cx\u03d5"))+S.cxE'*(ninv.*r0N);invf=1 ./f
 for i in 1:nB
  pr=invf[i]*rhat[i];el.ma[i]=pr+dvec[i];el.mb[i]=-pr
 end
 el
end
