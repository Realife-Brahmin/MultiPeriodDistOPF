# Same two-sided value algebra and fallback as the pinned production kernel.
# Y1 = x.Mab / U, while Y2 = U' \ y.Mab and g = U' \ g.
function right_combine!(x::PintElem, y::PintElem, w::PairWork)
    K = w.K; g = w.g
    K .= x.Mbb .+ y.Maa
    g .= x.mb .+ y.ma
    C = cholesky!(Symmetric(K); check=false)
    if !issuccess(C)
        r = pint_combine(x, y)
        copyto!(x.Maa, r.Maa); copyto!(x.Mab, r.Mab)
        copyto!(x.Mbb, r.Mbb); copyto!(x.ma, r.ma); copyto!(x.mb, r.mb)
        return x
    end
    U = C.U
    copyto!(w.Y1, x.Mab)
    BLAS.trsm!('R', 'U', 'N', 'N', 1.0, U.data, w.Y1)
    copyto!(w.Y2, y.Mab)
    ldiv!(U', w.Y2)
    ldiv!(U', g)
    BLAS.syrk!('U', 'N', -1.0, w.Y1, 1.0, x.Maa)
    LinearAlgebra.copytri!(x.Maa, 'U')
    copyto!(x.Mbb, y.Mbb)
    BLAS.syrk!('U', 'T', -1.0, w.Y2, 1.0, x.Mbb)
    LinearAlgebra.copytri!(x.Mbb, 'U')
    mul!(x.Mab, w.Y1, w.Y2, -1.0, 0.0)
    mul!(x.ma, w.Y1, g, -1.0, 1.0)
    copyto!(x.mb, y.mb)
    mul!(x.mb, w.Y2', g, -1.0, 1.0)
    x
end

function right_combine!(x::PintElem, y::PintVal, w::PairWork)
    K = w.K; g = w.g
    K .= x.Mbb .+ y.V
    g .= x.mb .+ y.v
    C = cholesky!(Symmetric(K); check=false)
    if !issuccess(C)
        r = pint_combine(x, y)
        copyto!(x.Maa, r.V); copyto!(x.ma, r.v)
        return PintVal(x.Maa, x.ma)
    end
    U = C.U
    copyto!(w.Y1, x.Mab)
    BLAS.trsm!('R', 'U', 'N', 'N', 1.0, U.data, w.Y1)
    ldiv!(U', g)
    BLAS.syrk!('U', 'N', -1.0, w.Y1, 1.0, x.Maa)
    LinearAlgebra.copytri!(x.Maa, 'U')
    mul!(x.ma, w.Y1, g, -1.0, 1.0)
    PintVal(x.Maa, x.ma)
end
