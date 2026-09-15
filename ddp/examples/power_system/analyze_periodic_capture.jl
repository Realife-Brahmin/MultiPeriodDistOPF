# Summarize shapes and compressibility (singular-value decay, effective rank)
# across a directory of periodic FilterDDP backward-pass captures produced by
# the FILTERDDP_PERIODIC_CAPTURE_DIR instrumentation in
# ddp/DDP4OPF.jl/src/backward_pass.jl (diagnostic-only, not part of
# the paper's optimization stack).
#
# The multi-RHS block captured is [-Qu | -B_active/-B ; -c | -cx], i.e. one
# feedforward column plus nx state-sensitivity columns. The feedforward column
# is typically far larger in magnitude than the state-sensitivity columns, so
# it is EXCLUDED from the compressibility analysis below (it would otherwise
# dominate any joint SVD and produce a misleadingly low effective rank). The
# question this script actually answers is: how compressible is the nx-column
# state-sensitivity block that FilterDDP solves for beta/omega/Vxx/Vx, since
# that block -- not the single feedforward solve -- is what scales with nx and
# dominates the per-stage KKT solve cost.
#
# Usage:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/analyze_periodic_capture.jl \
#         ddp/results/network_filterddp/periodic_capture_ieee123C_1ph_T3

using LinearAlgebra
using Printf
using Serialization
using SparseArrays

capture_dir = length(ARGS) >= 1 ? ARGS[1] : error("pass the periodic-capture directory")

function optimal_rank_for_relative_frobenius(singular_values, tolerance)
    isempty(singular_values) && return 0
    total = sum(abs2, singular_values)
    total == 0 && return 0
    for rank in 0:length(singular_values)
        tail = rank == length(singular_values) ? 0.0 :
               sum(abs2, @view singular_values[rank+1:end])
        sqrt(tail / total) <= tolerance && return rank
    end
    return length(singular_values)
end

files = sort(filter(f -> endswith(f, ".jls"), readdir(capture_dir)))
isempty(files) && error("no .jls captures found in $capture_dir")

rows = NamedTuple[]
for file in files
    cap = deserialize(joinpath(capture_dir, file))
    beta = Matrix(cap.beta)
    omega = Matrix(cap.omega)
    K_nnz = cap.K === nothing ? 0 : nnz(cap.K)
    K_n = cap.K === nothing ? 0 : size(cap.K, 1)
    K_density = K_n == 0 ? 0.0 : K_nnz / (K_n^2)
    rhs_full = cap.rhs_multi === nothing ? zeros(0, 0) : Matrix(cap.rhs_multi)
    rhs_state = size(rhs_full, 2) > 1 ? rhs_full[:, 2:end] : zeros(0, 0)

    sv_beta = isempty(beta) ? Float64[] : svdvals(beta)
    sv_omega = isempty(omega) ? Float64[] : svdvals(omega)
    sv_rhs = isempty(rhs_state) ? Float64[] : svdvals(rhs_state)

    push!(rows, (
        iteration=cap.iteration, stage=cap.stage, nx=cap.nx, nu=cap.nu, nc=cap.nc,
        K_shape="$(K_n)x$(K_n)", K_nnz=K_nnz, K_density=K_density,
        rhs_state_shape="$(size(rhs_state,1))x$(size(rhs_state,2))",
        beta_shape="$(size(beta,1))x$(size(beta,2))",
        omega_shape="$(size(omega,1))x$(size(omega,2))",
        alpha_len=length(cap.alpha), psi_len=length(cap.psi),
        SigmaL_len=length(cap.Sigma_L), SigmaU_len=length(cap.Sigma_U),
        Vxx_shape="$(size(cap.Vxx_incoming,1))x$(size(cap.Vxx_incoming,2))",
        Vx_len=length(cap.Vx_incoming),
        rhs_sv_max=isempty(sv_rhs) ? NaN : sv_rhs[1],
        rhs_sv_min=isempty(sv_rhs) ? NaN : sv_rhs[end],
        rhs_rank_1pct=optimal_rank_for_relative_frobenius(sv_rhs, 0.01),
        rhs_rank_5pct=optimal_rank_for_relative_frobenius(sv_rhs, 0.05),
        rhs_rank_10pct=optimal_rank_for_relative_frobenius(sv_rhs, 0.10),
        beta_rank_1pct=optimal_rank_for_relative_frobenius(sv_beta, 0.01),
        beta_rank_5pct=optimal_rank_for_relative_frobenius(sv_beta, 0.05),
        omega_rank_1pct=optimal_rank_for_relative_frobenius(sv_omega, 0.01),
        omega_rank_5pct=optimal_rank_for_relative_frobenius(sv_omega, 0.05),
        rhs_cond=isempty(sv_rhs) || sv_rhs[end] == 0 ? Inf : sv_rhs[1] / sv_rhs[end],
        beta_cond=isempty(sv_beta) || sv_beta[end] == 0 ? Inf : sv_beta[1] / sv_beta[end],
        barrier_mu=cap.barrier_mu, reg=cap.reg,
    ))
end
sort!(rows, by = r -> (r.iteration, r.stage))

output = joinpath(capture_dir, "periodic_capture_summary.csv")
open(output, "w") do io
    println(io, "iteration,stage,nx,nu,nc,K_shape,K_nnz,K_density,rhs_state_shape,beta_shape,omega_shape,alpha_len,psi_len,SigmaL_len,SigmaU_len,Vxx_shape,Vx_len,rhs_sv_max,rhs_sv_min,rhs_rank_1pct,rhs_rank_5pct,rhs_rank_10pct,beta_rank_1pct,beta_rank_5pct,omega_rank_1pct,omega_rank_5pct,rhs_cond,beta_cond,barrier_mu,reg")
    for r in rows
        @printf(io, "%d,%d,%d,%d,%d,%s,%d,%.6e,%s,%s,%s,%d,%d,%d,%d,%s,%d,%.6e,%.6e,%d,%d,%d,%d,%d,%d,%d,%.3e,%.3e,%.3e,%.3e\n",
            r.iteration, r.stage, r.nx, r.nu, r.nc, r.K_shape, r.K_nnz, r.K_density,
            r.rhs_state_shape, r.beta_shape, r.omega_shape, r.alpha_len, r.psi_len,
            r.SigmaL_len, r.SigmaU_len, r.Vxx_shape, r.Vx_len,
            r.rhs_sv_max, r.rhs_sv_min,
            r.rhs_rank_1pct, r.rhs_rank_5pct, r.rhs_rank_10pct,
            r.beta_rank_1pct, r.beta_rank_5pct, r.omega_rank_1pct, r.omega_rank_5pct,
            r.rhs_cond, r.beta_cond, r.barrier_mu, r.reg)
    end
end

@printf("PERIODIC_CAPTURE files=%d dir=%s\n", length(files), capture_dir)
@printf("Shapes (stage-invariant): K=%s rhs_state=%s beta=%s omega=%s Sigma_L/U=len %d Vxx=%s Vx=len %d\n",
    rows[1].K_shape, rows[1].rhs_state_shape, rows[1].beta_shape, rows[1].omega_shape,
    rows[1].SigmaL_len, rows[1].Vxx_shape, rows[1].Vx_len)
@printf("%-4s %-3s %-9s | rhs_state sv_max/sv_min (nx=%d cols) | rhs_rank(1%%/5%%/10%%) | beta_rank(1%%/5%%) | omega_rank(1%%/5%%) | rhs_cond | beta_cond\n",
    "iter", "stg", "K_dens%", rows[1].nx)
for r in rows
    @printf("%-4d %-3d %-9.4f | %10.3e / %10.3e         | %3d / %3d / %3d      | %3d / %3d       | %3d / %3d       | %.2e | %.2e\n",
        r.iteration, r.stage, 100r.K_density, r.rhs_sv_max, r.rhs_sv_min,
        r.rhs_rank_1pct, r.rhs_rank_5pct, r.rhs_rank_10pct,
        r.beta_rank_1pct, r.beta_rank_5pct, r.omega_rank_1pct, r.omega_rank_5pct,
        r.rhs_cond, r.beta_cond)
end
println("PERIODIC_CAPTURE wrote=$output")
