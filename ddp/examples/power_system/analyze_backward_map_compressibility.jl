# Test whether FilterDDP's state-to-control and state-to-dual sensitivity maps
# can be represented accurately by a low-rank basis or selected state columns.

using LinearAlgebra
using Printf
using Serialization
using SparseArrays

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
stage = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

capture = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp",
    "$(system)_T$(T)_stage$(stage)_oracle_capture.jls"))
beta = Matrix(capture.beta)
omega = Matrix(capture.omega)
Vxx = Matrix(Symmetric((capture.future_Vxx + capture.future_Vxx') / 2))
nx = size(beta, 2)

function relative_column_errors(M, approximation)
    errors = [norm(M[:, j] - approximation[:, j]) /
              max(norm(M[:, j]), eps()) for j in axes(M, 2)]
    return maximum(errors), sum(errors) / length(errors)
end

function svd_approximation(M, rank)
    F = svd(M)
    rank == 0 && return zeros(size(M))
    return F.U[:, 1:rank] * Diagonal(F.S[1:rank]) * F.Vt[1:rank, :]
end

# Normalize the two maps before selecting common columns so omega's numerical
# scale cannot overwhelm beta's contribution to the pivoting decision.
combined = [beta / norm(beta); omega / norm(omega)]
column_order = Vector(qr(combined, ColumnNorm()).p)
ranks = unique(sort([0, 1, 2, 4, 8, 12, 16, 24, 32, 40, 48, nx]))

rows = NamedTuple[]
for rank in ranks
    beta_svd = svd_approximation(beta, rank)
    omega_svd = svd_approximation(omega, rank)
    if rank == 0
        beta_selected = zeros(size(beta))
        omega_selected = zeros(size(omega))
    else
        selected = column_order[1:rank]
        C = combined[:, selected]
        coefficients = C \ combined
        beta_selected = beta[:, selected] * coefficients
        omega_selected = omega[:, selected] * coefficients
    end
    beta_selected_max, beta_selected_mean = relative_column_errors(beta, beta_selected)
    omega_selected_max, omega_selected_mean = relative_column_errors(omega, omega_selected)
    push!(rows, (
        rank=rank,
        beta_svd_fro=norm(beta-beta_svd)/norm(beta),
        omega_svd_fro=norm(omega-omega_svd)/norm(omega),
        beta_selected_fro=norm(beta-beta_selected)/norm(beta),
        omega_selected_fro=norm(omega-omega_selected)/norm(omega),
        beta_selected_col_max=beta_selected_max,
        beta_selected_col_mean=beta_selected_mean,
        omega_selected_col_max=omega_selected_max,
        omega_selected_col_mean=omega_selected_mean,
    ))
end

singular_beta = svdvals(beta)
singular_omega = svdvals(omega)
eigen_vxx = sort(abs.(eigvals(Symmetric(Vxx))), rev=true)

function optimal_rank_for_relative_frobenius(singular_values, tolerance)
    total = sum(abs2, singular_values)
    for rank in 0:length(singular_values)
        tail = rank == length(singular_values) ? 0.0 :
               sum(abs2, @view singular_values[rank+1:end])
        sqrt(tail / total) <= tolerance && return rank
    end
    return length(singular_values)
end

output = joinpath(REPO, "ddp", "results", "network_filterddp",
    "$(system)_T$(T)_stage$(stage)_backward_map_compressibility.csv")
open(output, "w") do io
    println(io, "rank,beta_svd_fro,omega_svd_fro,beta_selected_fro,omega_selected_fro,beta_selected_col_max,beta_selected_col_mean,omega_selected_col_max,omega_selected_col_mean")
    for row in rows
        @printf(io, "%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
            row.rank, row.beta_svd_fro, row.omega_svd_fro,
            row.beta_selected_fro, row.omega_selected_fro,
            row.beta_selected_col_max, row.beta_selected_col_mean,
            row.omega_selected_col_max, row.omega_selected_col_mean)
    end
end

spectrum_output = joinpath(REPO, "ddp", "results", "network_filterddp",
    "$(system)_T$(T)_stage$(stage)_backward_map_spectrum.csv")
open(spectrum_output, "w") do io
    println(io, "index,beta_singular,omega_singular,Vxx_abs_eigenvalue")
    for j in 1:nx
        @printf(io, "%d,%.12e,%.12e,%.12e\n",
            j, singular_beta[j], singular_omega[j], eigen_vxx[j])
    end
end

@printf("COMPRESSIBILITY system=%s T=%d stage=%d nx=%d\n", system, T, stage, nx)
for row in rows
    @printf("COMPRESSIBILITY rank=%d beta_opt=%.3e omega_opt=%.3e beta_selected=%.3e omega_selected=%.3e beta_col_max=%.3e omega_col_max=%.3e\n",
        row.rank, row.beta_svd_fro, row.omega_svd_fro,
        row.beta_selected_fro, row.omega_selected_fro,
        row.beta_selected_col_max, row.omega_selected_col_max)
end
@printf("COMPRESSIBILITY condition_beta=%.3e condition_omega=%.3e condition_Vxx_abs=%.3e\n",
    singular_beta[1]/singular_beta[end], singular_omega[1]/singular_omega[end],
    eigen_vxx[1]/eigen_vxx[end])
for tolerance in (0.10, 0.05, 0.01)
    @printf("COMPRESSIBILITY tolerance=%.2f optimal_rank_beta=%d optimal_rank_omega=%d optimal_rank_Vxx=%d\n",
        tolerance,
        optimal_rank_for_relative_frobenius(singular_beta, tolerance),
        optimal_rank_for_relative_frobenius(singular_omega, tolerance),
        optimal_rank_for_relative_frobenius(eigen_vxx, tolerance))
end
println("COMPRESSIBILITY wrote=$output")
println("COMPRESSIBILITY wrote=$spectrum_output")
