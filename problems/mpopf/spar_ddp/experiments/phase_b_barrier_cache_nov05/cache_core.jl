# Opt-in frozen scan cache only; no production hooks or defaults.
using LinearAlgebra,Serialization,SHA,SparseArrays
struct MergeMatrixCache
 U::Matrix{Float64}
 Y1::Matrix{Float64}
 Y2::Union{Nothing,Matrix{Float64}}
 A::Matrix{Float64}
 B::Union{Nothing,Matrix{Float64}}
 C::Union{Nothing,Matrix{Float64}}
end
curvfields(x)=hasproperty(x,:V) ? (:V,) : (:Maa,:Mab,:Mbb)
function populate_merge!(x,y,w)
 val=hasproperty(y,:V)
 w.K .= x.Mbb .+ (val ? y.V : y.Maa)
 w.g .= x.mb .+ (val ? y.v : y.ma)
 C=cholesky!(Symmetric(w.K,:U);check=false)
 if !issuccess(C)
  r=pint_combine(x,y) # EXACT original fallback, never a near-equality shortcut.
  if val
   copyto!(x.Maa,r.V);copyto!(x.ma,r.v);return PintVal(x.Maa,x.ma),nothing
  end
  for f in (:Maa,:Mab,:Mbb,:ma,:mb);copyto!(getproperty(x,f),getproperty(r,f));end
  return x,nothing
 end
 U=C.U
 transpose!(w.Y1,x.Mab);ldiv!(U',w.Y1)
 if !val;copyto!(w.Y2,y.Mab);ldiv!(U',w.Y2);end
 ldiv!(U',w.g)
 BLAS.syrk!('U','T',-1.0,w.Y1,1.0,x.Maa);LinearAlgebra.copytri!(x.Maa,'U')
 if !val
  copyto!(x.Mbb,y.Mbb);BLAS.syrk!('U','T',-1.0,w.Y2,1.0,x.Mbb);LinearAlgebra.copytri!(x.Mbb,'U')
  mul!(x.Mab,w.Y1',w.Y2,-1.0,0.0)
 end
 mul!(x.ma,w.Y1',w.g,-1.0,1.0)
 if !val;copyto!(x.mb,y.mb);mul!(x.mb,w.Y2',w.g,-1.0,1.0);end
 entry=MergeMatrixCache(Matrix(U),copy(w.Y1),val ? nothing : copy(w.Y2),copy(x.Maa),val ? nothing : copy(x.Mab),val ? nothing : copy(x.Mbb))
 return val ? PintVal(x.Maa,x.ma) : x,entry
end
function hit_merge!(x,y,w,c::MergeMatrixCache)
 val=hasproperty(y,:V);@assert val==isnothing(c.Y2)
 w.g .= x.mb .+ (val ? y.v : y.ma)
 ldiv!(UpperTriangular(c.U)',w.g)
 mul!(x.ma,c.Y1',w.g,-1.0,1.0)
 if !val;copyto!(x.mb,y.mb);mul!(x.mb,c.Y2',w.g,-1.0,1.0);end
 copyto!(x.Maa,c.A)
 if !val;copyto!(x.Mab,c.B);copyto!(x.Mbb,c.C);end
 val ? PintVal(x.Maa,x.ma) : x
end
function digest(x)
 io=IOBuffer();serialize(io,x);bytes2hex(sha256(take!(io)))
end
function stage_signature(p)
 S=p.S
 # Exclude current RHS/gradients/mu and mutable work buffers; include all
 # matrix-driving and nominal data, even factors which could coincide at leaves.
 fields=(:x,:u,Symbol("\u03d5"),:zl,:zu,:lxx,:lux,:cxx,:cx,:cu,:fx,:fu,:fxx,:active,:fu_active,:lux_in_B,:K,:lay,:stat,:fac,:Bd,:Up,:b0,:nzc,:cxE,:inv_ul,:inv_uu,Symbol("\u03a3_L"),Symbol("\u03a3_U"))
 digest(([(f,getproperty(S,f)) for f in fields],p.cl,p.nx,p.nu,p.nc))
end
function eligible(previous,current)
 previous!==nothing && previous.k==current.k && current.j==previous.j+1 &&
 current.mu<previous.mu && current.context==previous.context &&
 current.proof==:pinned_unshifted_tp_pre && previous.proof==current.proof &&
 current.stage_keys==previous.stage_keys && current.terminal_key==previous.terminal_key
end
# Current vector work from _pp_elem, same arithmetic and ordering. No ma/mb
# from the old iterate or the saved leaf is used on a hit.
function leaf_vectors!(el,S)
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
exactbits(a,b)=size(a)==size(b)&&reinterpret(UInt64,vec(a))==reinterpret(UInt64,vec(b))
function errors(a,b)
 d=norm(a-b);s=norm(b)
 (absolute=d,reference=s,scaled=d/max(1.0,s),relative=d/max(1e-10,s),maxabs=maximum(abs,a-b),exact=exactbits(a,b),pass=d/max(1.0,s)<1e-7&&d/max(1e-10,s)<1e-5)
end
function rounds(N)
 all=Vector{Vector{Tuple{Int,Int}}}();h=1
 while h<N+1
  pairs=Tuple{Int,Int}[]
  for i in 1:2h:N+1
   mid=i+h;mid>N+1&&continue
   for t in i:mid-1;push!(pairs,(t,mid));end
  end
  push!(all,pairs);h*=2
 end
 all
end
cache_numeric_bytes(c::MergeMatrixCache)=sum(isnothing(x) ? 0 : sizeof(x) for x in (c.U,c.Y1,c.Y2,c.A,c.B,c.C))

function shared_slot_view(a,n)
 ptr=pointer(a)
 mat(i)=unsafe_wrap(Array,ptr+(i-1)*n*n*sizeof(Float64),(n,n))
 vec(i)=unsafe_wrap(Array,ptr+(3n*n+(i-1)*n)*sizeof(Float64),(n,))
 PintElem(mat(1),mat(2),mat(3),vec(1),vec(2))
end
