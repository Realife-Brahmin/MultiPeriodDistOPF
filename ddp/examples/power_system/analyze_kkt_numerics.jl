using LinearAlgebra
using Printf
using Serialization
using SparseArrays

input = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
payload = deserialize(input)
K = payload.K
rhs = hasproperty(payload, :rhs) ? payload.rhs : payload.rhs_multi

n = size(K, 1)
F = @timed lu(K)
fac = F.value
X = @timed (fac \ rhs)
relres = norm(K * X.value - rhs) / max(norm(rhs), eps())
fill_ratio = (nnz(fac.L) + nnz(fac.U)) / nnz(K)
diag_u = abs.(diag(fac.U))

println("input=$input")
@printf("K n=%d nnz=%d density=%.6e symmetry_error=%.3e\n",
    n, nnz(K), nnz(K) / n^2, norm(K - K') / max(norm(K), eps()))
@printf("LU factor_s=%.6f solve_s=%.6f fill_ratio=%.3f min_abs_Udiag=%.3e max_abs_Udiag=%.3e pivot_ratio=%.3e relative_residual=%.3e\n",
    F.time, X.time, fill_ratio, minimum(diag_u), maximum(diag_u),
    minimum(diag_u) / maximum(diag_u), relres)

# A complete inertia calculation is affordable only for the smallest captured
# system. It confirms whether the symmetric KKT matrix has positive, negative,
# or numerically zero eigenvalues. Larger systems use sparse-LU diagnostics.
if n <= 2500
    row_scale = vec(maximum(abs.(K), dims=2))
    positive_scale = row_scale[row_scale .> 0]
    @printf("coefficient_scale min_row_max=%.3e max_row_max=%.3e spread=%.3e\n",
        minimum(positive_scale), maximum(positive_scale),
        maximum(positive_scale) / minimum(positive_scale))
    d = 1 ./ sqrt.(max.(row_scale, eps()))
    K_equilibrated = Diagonal(d) * K * Diagonal(d)
    vals = eigvals(Symmetric(Matrix(K_equilibrated)))
    scale = maximum(abs, vals)
    tol = n * eps(scale)
    npos = count(>(tol), vals)
    nneg = count(<(-tol), vals)
    nzero = n - npos - nneg
    @printf("equilibrated_inertia positive=%d negative=%d near_zero=%d tolerance=%.3e min_abs_eigenvalue=%.3e max_abs_eigenvalue=%.3e spectral_condition=%.3e\n",
        npos, nneg, nzero, tol, minimum(abs, vals), scale,
        scale / minimum(abs, vals))
else
    println("inertia=skipped (dense eigendecomposition deliberately limited to n<=2500)")
end
