# Offline sensitivity diagnostics only. Production source is never modified.
const LOCATIONS=((1,13,14),(2,5,7),(3,9,13))
const INTERVENTIONS=(:H,:BL,:BR,:h,:Hdiag,:BLoffdiag)
fields_of(x)=hasproperty(x,:V) ? (:V,:v) : (:Maa,:Mab,:Mbb,:ma,:mb)
merge_inputs(x,y)=(H=x.Mbb+y.Maa,BL=copy(x.Mab),BR=copy(y.Mab),h=x.mb+y.ma)
function retained_inputs(current,previous,mode)
    H=copy(current.H);BL=copy(current.BL);BR=copy(current.BR);h=copy(current.h)
    mode==:H && copyto!(H,previous.H)
    mode==:BL && copyto!(BL,previous.BL)
    mode==:BR && copyto!(BR,previous.BR)
    mode==:h && copyto!(h,previous.h)
    if mode==:Hdiag
        for i in axes(H,1);H[i,i]=previous.H[i,i];end
    elseif mode==:BLoffdiag
        for j in axes(BL,2),i in axes(BL,1);i!=j && (BL[i,j]=previous.BL[i,j]);end
    end
    @assert mode in INTERVENTIONS || mode==:current
    (H=H,BL=BL,BR=BR,h=h)
end
# Same original-coordinate UPPER completion and arithmetic as production.
# External Maa/Mbb/ma/mb remain current even for inconsistent mixed inputs.
# A failed modified factorization is an explicit outcome; no silent fallback.
function input_combine!(x,y,w,z)
    copyto!(w.K,z.H);copyto!(w.g,z.h)
    C=cholesky!(Symmetric(w.K,:U);check=false)
    issuccess(C)||throw(PosDefException(0))
    U=C.U
    transpose!(w.Y1,z.BL);ldiv!(U',w.Y1)
    copyto!(w.Y2,z.BR);ldiv!(U',w.Y2);ldiv!(U',w.g)
    BLAS.syrk!('U','T',-1.0,w.Y1,1.0,x.Maa);LinearAlgebra.copytri!(x.Maa,'U')
    copyto!(x.Mbb,y.Mbb);BLAS.syrk!('U','T',-1.0,w.Y2,1.0,x.Mbb);LinearAlgebra.copytri!(x.Mbb,'U')
    mul!(x.Mab,w.Y1',w.Y2,-1.0,0.0)
    mul!(x.ma,w.Y1',w.g,-1.0,1.0)
    copyto!(x.mb,y.mb);mul!(x.mb,w.Y2',w.g,-1.0,1.0)
    x
end
function scan_ops(N)
    out=Tuple{Int,Int,Int}[];h=1;depth=0
    while h<N+1
        depth+=1
        for i in 1:2h:N+1
            mid=i+h;mid>N+1&&continue
            for t in i:mid-1;push!(out,(depth,t,mid));end
        end
        h*=2
    end
    out
end
function dependency_trace(N,target=(2,5,7))
    changed=falses(N+1);out=Any[]
    for (depth,t,mid) in scan_ops(N)
        hit=(depth,t,mid)==target
        affected=hit||changed[t]||changed[mid]
        push!(out,(depth=depth,left=t,right=mid,injected=hit,left_dependency=changed[t],right_dependency=changed[mid],affected=affected))
        changed[t]=affected
    end
    out,changed
end
metric(a,b)=(absolute=norm(a-b),reference=norm(b),scaled=norm(a-b)/max(1.0,norm(b)),relative=norm(a-b)/max(1e-10,norm(b)),maxabs=maximum(abs,a-b),reference_near_zero=norm(b)<1e-10)
passes_metric(m)=isfinite(m.scaled)&&isfinite(m.relative)&&m.scaled<1e-7&&m.relative<1e-5
