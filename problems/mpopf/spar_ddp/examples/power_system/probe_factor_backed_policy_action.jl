using LinearAlgebra
using Printf
using Serialization
using SparseArrays

capture_path = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
repeats = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20

capture = deserialize(capture_path)
K = capture.K
rhs_state = @view capture.rhs[:, 2:end]
beta = capture.beta
omega = capture.omega
nx = size(beta, 2)
nu = size(beta, 1)
nc = size(omega, 1)

size(K) == (nu + nc, nu + nc) || error("inconsistent captured KKT dimensions")
size(rhs_state) == (nu + nc, nx) || error("inconsistent captured state RHS")

# The state-column block already has FilterDDP's sign convention:
# K * [beta; omega] = rhs_state = -[B; cx].  A policy action can therefore
# be evaluated from one RHS without materializing either complete map.
F = lu(K)
directions = Vector{Vector{eltype(beta)}}()
for index in unique((1, cld(nx, 2), nx))
    direction = zeros(eltype(beta), nx)
    direction[index] = one(eltype(beta))
    push!(directions, direction)
end
push!(directions, collect(range(-one(eltype(beta)), one(eltype(beta)); length=nx)))

max_beta_relative_error = Ref(0.0)
max_omega_relative_error = Ref(0.0)
max_kkt_relative_residual = Ref(0.0)
for direction in directions
    action_rhs = rhs_state * direction
    action = F \ action_rhs
    beta_action = @view action[1:nu]
    omega_action = @view action[nu+1:end]
    beta_reference = beta * direction
    omega_reference = omega * direction
    beta_scale = max(norm(beta_reference), eps(real(eltype(beta))))
    omega_scale = max(norm(omega_reference), eps(real(eltype(beta))))
    max_beta_relative_error[] = max(max_beta_relative_error[],
        norm(beta_action - beta_reference) / beta_scale)
    max_omega_relative_error[] = max(max_omega_relative_error[],
        norm(omega_action - omega_reference) / omega_scale)
    max_kkt_relative_residual[] = max(max_kkt_relative_residual[],
        norm(K * action - action_rhs) / max(norm(action_rhs), eps(real(eltype(beta)))))
end

probe_direction = directions[end]
direct_times = Float64[]
factor_times = Float64[]
for _ in 1:repeats
    push!(direct_times, @elapsed begin
        beta * probe_direction
        omega * probe_direction
    end)
    action_rhs = rhs_state * probe_direction
    push!(factor_times, @elapsed F \ action_rhs)
end
median_value(values) = sort(values)[cld(length(values), 2)]

@printf("FACTOR_POLICY_ACTION nx=%d nu=%d nc=%d beta_rel=%.3e omega_rel=%.3e kkt_residual=%.3e direct_action_s=%.9f factor_action_s=%.9f dense_maps_MiB=%.3f factor_summary_MiB=%.3f\n",
    nx, nu, nc, max_beta_relative_error[], max_omega_relative_error[],
    max_kkt_relative_residual[], median_value(direct_times), median_value(factor_times),
    (Base.summarysize(beta) + Base.summarysize(omega)) / 2.0^20,
    Base.summarysize(F) / 2.0^20)
