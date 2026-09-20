_derivative_matrix(A, m, n) = issparse(A) ? sparse(A) : reshape(vec(A), m, n)

# ---------------------------------------------------------------- frozen KKT --
# Fast-decoupled-style reuse of the per-stage KKT factorisation. FDPF freezes the
# Jacobian because susceptance is genuinely constant; here K = [H cu'; cu 0] is
# NOT constant -- it carries the barrier terms Sigma_L + Sigma_U, which scale
# like mu/s^2 and move sharply as mu drops and iterates approach bounds. So this
# is a quasi-Newton approximation whose worth is an empirical question: it saves
# assembly+factorisation per stage per iteration, and pays for it in step
# quality. Opt in with FILTERDDP_FREEZE_KKT=N (refactor every N backward passes);
# N <= 1 is the default and reproduces exact behaviour bit for bit.
const _FROZEN_KKT = Dict{Int,Any}()
const _FROZEN_AT  = Dict{Int,Int}()
const _FROZEN_STATS = Dict{Symbol,Int}(:factorisations => 0, :reuses => 0)

_freeze_period() = parse(Int, get(ENV, "FILTERDDP_FREEZE_KKT", "1"))

function _frozen_reset!()
    empty!(_FROZEN_KKT); empty!(_FROZEN_AT)
    _FROZEN_STATS[:factorisations] = 0; _FROZEN_STATS[:reuses] = 0
    return nothing
end

function _frozen_lu(t::Int, K, iter::Int)
    period = _freeze_period()
    if period <= 1
        _FROZEN_STATS[:factorisations] += 1
        return lu(K)
    end
    if !haskey(_FROZEN_KKT, t) || (iter - _FROZEN_AT[t]) >= period
        _FROZEN_KKT[t] = lu(K)
        _FROZEN_AT[t] = iter
        _FROZEN_STATS[:factorisations] += 1
    else
        _FROZEN_STATS[:reuses] += 1
    end
    return _FROZEN_KKT[t]
end


# ----------------------------------------------------- diagonal Hessian --
# Agenda question (R. Gupta, 2026-09-18): "can we not just use a diagonalized
# matrix?". FILTERDDP_DIAG_HESSIAN=1 replaces the stage Hessian block of
# K = [H cu'; cu 0] with its diagonal, leaving the constraint Jacobian intact.
#
# Two things make this a real approximation rather than a cheap trick:
#   * On ieee123 T=3 about 83% of nnz(H) IS the dense nB x nB block
#     fu' * Vxx * fu -- the intertemporal curvature. Diagonalising deletes
#     exactly the information that makes this a second-order method.
#   * diag(H) is singular in general: at the terminal stage nnz(H) = 536
#     against nu = 791, because line flows, voltages and currents pick up a
#     diagonal entry only through cuu, and unbounded controls get no barrier
#     term at all. Hence FILTERDDP_DIAG_HESSIAN_FLOOR (default 1e-8), which
#     floors each diagonal entry from below.
_diag_hessian_enabled() = get(ENV, "FILTERDDP_DIAG_HESSIAN", "0") != "0"
_diag_hessian_floor() = parse(Float64, get(ENV, "FILTERDDP_DIAG_HESSIAN_FLOOR", "1e-8"))

function _diagonalise_hessian(H)
    floor_val = _diag_hessian_floor()
    d = diag(H)
    @inbounds for i in eachindex(d)
        d[i] = max(d[i], floor_val)
    end
    return issparse(H) ? spdiagm(0 => d) : diagm(d)
end


function backward_pass!(solver::Solver{T, nx, nu, nc, nux, ncx}, ocp::OCP{T, nx, nu, nc},
            traj::Vector{TrajectoryElement{T, nx, nu, nc}}, data::SolverData{T}, options::Options{T}; verbose::Bool=false
            ) where {T, nx, nu, nc, nux, ncx}
    # A stale factorisation must never leak between solves.
    data.k == 0 && _frozen_reset!()
    reg::T = 0.0
    μ = data.μ
    δ_c = 0.
    reg = 0.0

    dynamics = ocp.dynamics
    cl = ocp.control_limits
    ni = (cl.nl + cl.nu) * ocp.N
    rhs_workspace = nothing
    feasibility_diagnostic = get(ENV, "FILTERDDP_FEASIBILITY_DIAGNOSTIC", "0") == "1"
    equality_ssq = T(0); equality_count = 0; equality_max = T(0)
    equality_worst_stage = 0; equality_worst_index = 0
    dynamics_ssq = T(0); dynamics_count = 0; dynamics_max = T(0)
    dynamics_worst_stage = 0; dynamics_worst_index = 0
    bound_ssq = T(0); bound_count = 0; bound_max = T(0)
    bound_worst_stage = 0; bound_worst_index = 0; bound_worst_kind = "none"
    
    while reg <= options.reg_max
        equality_ssq = T(0); equality_count = 0; equality_max = T(0)
        equality_worst_stage = 0; equality_worst_index = 0
        dynamics_ssq = T(0); dynamics_count = 0; dynamics_max = T(0)
        dynamics_worst_stage = 0; dynamics_worst_index = 0
        bound_ssq = T(0); bound_count = 0; bound_max = T(0)
        bound_worst_stage = 0; bound_worst_index = 0; bound_worst_kind = "none"
        data.status = 0
        V̂x = zeros(T, nx)
        V̂xx = zeros(T, nx, nx)
        λ = zeros(T, nx)

        data.barrier_lagrangian_curr = T(0.0)
        data.primal_1_curr = T(0.0)
        data.primal_inf = T(0.0)
        data.cs_inf_μ = T(0.0)
        data.cs_inf_0 = T(0.0)
        data.objective = T(0.0)
        data.dual_inf = T(0.0)
        data.expected_change_L = T(0.0)
        ϕ_norm = T(0.0)
        z_norm = T(0.0) 
        
        for t = ocp.N:-1:1
            timing_diagnostic = get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"
            memory_diagnostic = get(ENV, "FILTERDDP_MEMORY_DIAGNOSTIC", "0") == "1"
            stage_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            stage_maxrss_start = memory_diagnostic ? Sys.maxrss() : 0
            derivative_start_ns = time_ns()
            x, u, ϕ, zl, zu = traj[t].x, traj[t].u, traj[t].ϕ, traj[t].zl, traj[t].zu
            # PATCH: per-stage data when supplied, else the shared function
            constraints = stage_con(ocp, t)
            objective_t = stage_obj(ocp, t)

            # evaluate derivatives

            first_order_ns = 0
            second_order_ns = 0
            dynamics_first_order_ns = 0
            if t == ocp.N
                callback_start_ns = time_ns()
                lx = Vector(ocp.term_objective.lx(x, u))
                lu_ = Vector(ocp.term_objective.lu(x, u))
                data.objective += ocp.term_objective.l(x, u)[1]
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                lxx = _derivative_matrix(ocp.term_objective.lxx(x, u), nx, nx)
                lux = _derivative_matrix(ocp.term_objective.lux(x, u), nu, nx)
                luu = _derivative_matrix(ocp.term_objective.luu(x, u), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
            else
                callback_start_ns = time_ns()
                lx = Vector(objective_t.lx(x, u))
                lu_ = Vector(objective_t.lu(x, u))
                data.objective += objective_t.l(x, u)[1]
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                lxx = _derivative_matrix(objective_t.lxx(x, u), nx, nx)
                lux = _derivative_matrix(objective_t.lux(x, u), nu, nx)
                luu = _derivative_matrix(objective_t.luu(x, u), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
            end
            
            if nc > 0
                callback_start_ns = time_ns()
                c = Vector{T}(constraints.c(x, u))
                cx = _derivative_matrix(constraints.cx(x, u), nc, nx)
                cu = _derivative_matrix(constraints.cu(x, u), nc, nu)
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                cxx = _derivative_matrix(constraints.cxx(x, u, ϕ), nx, nx)
                cux = _derivative_matrix(constraints.cux(x, u, ϕ), nu, nx)
                cuu = _derivative_matrix(constraints.cuu(x, u, ϕ), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
                
                # evaluate constraint violation norms
                data.primal_1_curr += norm(c, 1)
                data.primal_inf = max(data.primal_inf, norm(c, Inf))
                if feasibility_diagnostic
                    equality_ssq += sum(abs2, c)
                    equality_count += length(c)
                    if !isempty(c)
                        value, index = findmax(abs, c)
                        if value > equality_max
                            equality_max = value
                            equality_worst_stage = t
                            equality_worst_index = index
                        end
                    end
                end
            end

            if feasibility_diagnostic
                if t < ocp.N
                    dynamics_residual = traj[t+1].x - dynamics.f(x, u)
                    dynamics_ssq += sum(abs2, dynamics_residual)
                    dynamics_count += length(dynamics_residual)
                    if !isempty(dynamics_residual)
                        value, index = findmax(abs, dynamics_residual)
                        if value > dynamics_max
                            dynamics_max = value
                            dynamics_worst_stage = t
                            dynamics_worst_index = index
                        end
                    end
                end
                @inbounds for i in eachindex(u)
                    if cl.maskl[i]
                        violation = max(cl.l[i] - u[i], zero(T))
                        bound_ssq += violation^2; bound_count += 1
                        if violation > bound_max
                            bound_max = violation; bound_worst_stage = t
                            bound_worst_index = i; bound_worst_kind = "lower"
                        end
                    end
                    if cl.masku[i]
                        violation = max(u[i] - cl.u[i], zero(T))
                        bound_ssq += violation^2; bound_count += 1
                        if violation > bound_max
                            bound_max = violation; bound_worst_stage = t
                            bound_worst_index = i; bound_worst_kind = "upper"
                        end
                    end
                end
            end

            callback_start_ns = time_ns()
            fx = _derivative_matrix(dynamics.fx(x, u), nx, nx)
            fu = _derivative_matrix(dynamics.fu(x, u), nx, nu)
            first_order_ns += time_ns() - callback_start_ns
            callback_start_ns = time_ns()
            fxx = _derivative_matrix(nc == 0 ? dynamics.fxx(x, u, V̂x) : dynamics.fxx(x, u, λ), nx, nx)
            fux = _derivative_matrix(nc == 0 ? dynamics.fux(x, u, V̂x) : dynamics.fux(x, u, λ), nu, nx)
            fuu = _derivative_matrix(nc == 0 ? dynamics.fuu(x, u, V̂x) : dynamics.fuu(x, u, λ), nu, nu)
            second_order_ns += time_ns() - callback_start_ns
            derivative_s = (time_ns() - derivative_start_ns) / 1e9
            first_order_s = first_order_ns / 1e9
            second_order_s = second_order_ns / 1e9
            derivative_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - stage_alloc_start : 0
            algebra_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            algebra_start_ns = time_ns()

            # evaluate barrier Lagrangian

            ul = u - cl.l
            uu = cl.u - u
            data.barrier_lagrangian_curr -= μ* (dot(log.(ul), cl.maskl) + dot(log.(uu), cl.masku))

            # evaluate complementary slackness errors
            cs_l = ul .* zl
            cs_u = uu .* zu
            data.cs_inf_0 = max(data.cs_inf_0, norm(cs_l .* cl.maskl, Inf))
            data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_l .- μ).* cl.maskl, Inf))
            data.cs_inf_0 = max(data.cs_inf_0, norm(cs_u .* cl.masku, Inf))
            data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_u .- μ) .* cl.masku, Inf))

            inv_ul = inv.(ul) .* cl.maskl
            inv_uu = inv.(uu) .* cl.masku

            # Qû = Lu' -μŪ^{-1}e + fu' * V̂x
            Qû = lu_ + fu' * V̂x + μ .* (inv_uu - inv_ul)
            # C = Lxx + fx' * Vxx * fx + V̄x ⋅ fxx
            C = lxx + fx' * V̂xx * fx + fxx
    
            sparse_stage = issparse(luu) || issparse(fu) || (nc > 0 && (issparse(cu) || issparse(cuu)))
            structured_B = sparse_stage && issparse(fu) && nnz(lux) == 0 &&
                nnz(fux) == 0 && (nc == 0 || nnz(cux) == 0)
            active_B_rows = Int[]
            B_active = Matrix{T}(undef, 0, 0)
            ux_tmp = Matrix{T}(undef, 0, 0)
            # Ĥ = Luu + Σ + fu' * Vxx * fu + V̄x ⋅ fuu
            Σ_L = inv_ul .* zl
            Σ_U = inv_uu .* zu
            if sparse_stage
                Ĥ = sparse(luu) + spdiagm(0 => Σ_L + Σ_U) +
                     sparse(fu)' * sparse(V̂xx) * sparse(fu) + sparse(fuu)
            else
                ux_tmp = fu' * V̂xx
                Ĥ = luu + diagm(Σ_L) + diagm(Σ_U) + ux_tmp * fu + fuu
            end
            # B = Lux + fu' * Vxx * fx + V̄x ⋅ fux
            if structured_B
                active_B_rows = sort!(unique(findnz(fu)[2]))
                fu_active = Matrix(@view fu[:, active_B_rows])
                B_active = fu_active' * V̂xx * fx
            else
                isempty(ux_tmp) && (ux_tmp = fu' * V̂xx)
                B = ux_tmp * fx
                B .+= lux
                B .+= fux
            end

            if nc > 0
                data.barrier_lagrangian_curr += dot(c, ϕ)
                Qû = Qû + cu' * ϕ
                C = C + cxx
                Ĥ = Ĥ + cuu
                !structured_B && (B .+= cux)
            end
            
            # inertia correction / regularisation
            if !iszero(reg)
                @inbounds for i in axes(Ĥ, 1)
                    Ĥ[i, i] += reg
                end
            end

            # Opt-in diagonal curvature model; see _diagonalise_hessian above.
            # Applied after the inertia correction so that reg still reaches the
            # diagonal, and before sparse_kkt is decided so the branch is unchanged.
            _diag_hessian_enabled() && (Ĥ = _diagonalise_hessian(Ĥ))

            # Sparse network models use the full saddle-point system directly,
            # avoiding a dense QR basis and the explicit reduced Hessian Z'HZ.
            sparse_kkt = nc > 0 && (issparse(Ĥ) || issparse(cu))
            algebra_s = (time_ns() - algebra_start_ns) / 1e9
            algebra_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - algebra_alloc_start : 0
            kkt_assembly_s = 0.0
            factor_s = 0.0
            solve_s = 0.0
            kkt_alloc_bytes = 0
            factor_alloc_bytes = 0
            solve_alloc_bytes = 0
            blocked_value = false
            blocked_α = Vector{T}()
            blocked_ψ = Vector{T}()
            blocked_Vxx = Matrix{T}(undef, 0, 0)
            blocked_Vx = Vector{T}()
            if sparse_kkt
                kkt_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                kkt_start_ns = time_ns()
                Ĥ = sparse(Symmetric(Ĥ))
                K = [Ĥ sparse(cu'); sparse(cu) spzeros(T, nc, nc)]
                # Fill-in diagnostic: nnz of the coefficient before factorisation
                # and of the LU factors after, per stage. Opt-in; this is what
                # distinguishes "the matrix got denser" from "pivoting got worse"
                # as the explanation for factorisation cost rising ~6x on an
                # identically sized K when the batteries actually cycle.
                nnz_diagnostic = get(ENV, "FILTERDDP_NNZ_DIAGNOSTIC", "0") == "1"

                capture_this_kkt = haskey(ENV, "FILTERDDP_CAPTURE_KKT") &&
                    t == parse(Int, get(ENV, "FILTERDDP_CAPTURE_STAGE", "1"))
                blocked_value = get(ENV, "FILTERDDP_BLOCKED_VALUE_RHS", "0") == "1" &&
                    get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" &&
                    structured_B && !capture_this_kkt
                block_width = min(parse(Int, get(ENV, "FILTERDDP_VALUE_BLOCK_WIDTH", "128")), nx)
                rhs_width = blocked_value ? block_width : nx + 1
                if isnothing(rhs_workspace)
                    rhs_workspace = Matrix{T}(undef, nu + nc, rhs_width)
                end
                rhs = rhs_workspace
                @views @. rhs[1:nu, 1] = -Qû
                @views @. rhs[nu+1:end, 1] = -c
                if !blocked_value
                    @views begin
                    if structured_B
                        fill!(rhs[1:nu, 2:nx+1], zero(T))
                        @. rhs[active_B_rows, 2:nx+1] = -B_active
                    else
                        @. rhs[1:nu, 2:nx+1] = -B
                    end
                    @. rhs[nu+1:end, 2:nx+1] = -cx
                    end
                end
                captured_rhs = capture_this_kkt ? copy(rhs) : nothing
                # DIAGNOSTIC: snapshot the multi-RHS block *before* ldiv! overwrites
                # it in place, so the periodic capture below records the true input
                # to the KKT solve rather than its solution.
                periodic_capture_active = haskey(ENV, "FILTERDDP_PERIODIC_CAPTURE_DIR") &&
                    data.k % parse(Int, get(ENV, "FILTERDDP_PERIODIC_CAPTURE_STRIDE", "5")) == 0
                periodic_rhs_snapshot = periodic_capture_active ? copy(rhs) : nothing
                kkt_assembly_s = (time_ns() - kkt_start_ns) / 1e9
                kkt_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - kkt_alloc_start : 0
                kkt_solution = Matrix{T}(undef, 0, 0)
                F = nothing
                try
                    if timing_diagnostic || memory_diagnostic
                        factor_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                        factor_start_ns = time_ns()
                        F = _frozen_lu(t, K, data.k)
                        if nnz_diagnostic
                            # L and U extraction is costly, hence opt-in only
                            nL = nnz(F.L); nU = nnz(F.U)
                            @printf("FILTERDDP_NNZ iteration=%d stage=%d n=%d nnz_K=%d nnz_LU=%d fill_ratio=%.3f
",
                                    data.k, t, size(K, 1), nnz(K), nL + nU, (nL + nU) / max(nnz(K), 1))
                            flush(stdout)
                        end
                        factor_s = (time_ns() - factor_start_ns) / 1e9
                        factor_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - factor_alloc_start : 0
                        solve_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                        solve_start_ns = time_ns()
                        if blocked_value
                            ldiv!(F, @view(rhs[:, 1:1]))
                            blocked_α = copy(@view rhs[1:nu, 1])
                            blocked_ψ = copy(@view rhs[nu+1:end, 1])
                            blocked_Vxx = zeros(T, nx, nx)
                            blocked_Vx = zeros(T, nx)
                            for first_col in 1:block_width:nx
                                columns = first_col:min(first_col + block_width - 1, nx)
                                width = length(columns)
                                rhs_block = @view rhs[:, 1:width]
                                fill!(rhs_block, zero(T))
                                @views @. rhs[active_B_rows, 1:width] = -B_active[:, columns]
                                @views @. rhs[nu+1:end, 1:width] = -cx[:, columns]
                                ldiv!(F, rhs_block)
                                β_block = @view rhs[1:nu, 1:width]
                                ω_block = @view rhs[nu+1:end, 1:width]
                                @views blocked_Vxx[columns, :] .=
                                    β_block[active_B_rows, :]' * B_active + ω_block' * cx
                                @views blocked_Vx[columns] .= β_block' * Qû + ω_block' * c
                                all(isfinite, rhs_block) || (data.status = 1; break)
                            end
                        else
                            ldiv!(F, rhs)
                            kkt_solution = rhs
                        end
                        solve_s = (time_ns() - solve_start_ns) / 1e9
                        solve_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - solve_alloc_start : 0
                    else
                        F = _frozen_lu(t, K, data.k)
                        if nnz_diagnostic
                            # L and U extraction is costly, hence opt-in only
                            nL = nnz(F.L); nU = nnz(F.U)
                            @printf("FILTERDDP_NNZ iteration=%d stage=%d n=%d nnz_K=%d nnz_LU=%d fill_ratio=%.3f
",
                                    data.k, t, size(K, 1), nnz(K), nL + nU, (nL + nU) / max(nnz(K), 1))
                            flush(stdout)
                        end
                        if blocked_value
                            ldiv!(F, @view(rhs[:, 1:1]))
                            blocked_α = copy(@view rhs[1:nu, 1])
                            blocked_ψ = copy(@view rhs[nu+1:end, 1])
                            blocked_Vxx = zeros(T, nx, nx)
                            blocked_Vx = zeros(T, nx)
                            for first_col in 1:block_width:nx
                                columns = first_col:min(first_col + block_width - 1, nx)
                                width = length(columns)
                                rhs_block = @view rhs[:, 1:width]
                                fill!(rhs_block, zero(T))
                                @views @. rhs[active_B_rows, 1:width] = -B_active[:, columns]
                                @views @. rhs[nu+1:end, 1:width] = -cx[:, columns]
                                ldiv!(F, rhs_block)
                                β_block = @view rhs[1:nu, 1:width]
                                ω_block = @view rhs[nu+1:end, 1:width]
                                @views blocked_Vxx[columns, :] .=
                                    β_block[active_B_rows, :]' * B_active + ω_block' * cx
                                @views blocked_Vx[columns] .= β_block' * Qû + ω_block' * c
                                all(isfinite, rhs_block) || (data.status = 1; break)
                            end
                        else
                            ldiv!(F, rhs)
                            kkt_solution = rhs
                        end
                    end
                    !blocked_value && !all(isfinite, kkt_solution) && (data.status = 1)
                catch err
                    verbose && @warn "Sparse KKT factorization failed" exception=(err, catch_backtrace())
                    data.status = 1
                end
            elseif nc > 0
                A = cu'
                qrf = qr(A)
                Q = Matrix(qrf.Q * Matrix{T}(I, nu, nu))
                Y = Q[:, 1:nc]
                Z = Q[:, nc+1:nu]
                
                AY = A' * Y
                fk = lu(AY)
                α_β_y = fk \ [-c -cx]

                Ĥ = Symmetric(Ĥ)
                M = Symmetric(Z' * Ĥ * Z)
                ck = cholesky(M; check=false)
            else
                ck = cholesky(Symmetric(Ĥ); check=false)
            end
            !sparse_kkt && ck.info != 0 && (data.status = 1)
            
            if data.status == 1
                if iszero(reg) # initial setting of regularisation
                    reg = (data.reg_last == 0.0) ? options.reg_1 : max(options.reg_min, options.κ_w_m * data.reg_last)
                elseif get(ENV, "FILTERDDP_PERIODIC_CAPTURE_KKT_ONLY", "0") == "1"
                    capture_file = joinpath(capture_dir,
                        @sprintf("iter%04d_stage%03d.jls", data.k, t))
                    Serialization.serialize(capture_file, (
                        iteration=data.k, stage=t, nx=nx, nu=nu, nc=nc,
                        K=copy(K), rhs=copy(@view periodic_rhs_snapshot[:, 1:1]),
                        barrier_mu=μ, reg=reg))
                else
                    reg = (data.reg_last == 0.0) ? options.κ_̄w_p * reg : options.κ_w_p * reg
                end
                break
            end

            update_start_ns = time_ns()
            update_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            if sparse_kkt
                update_rule = solver.update[t]
                copyto!(update_rule.α, blocked_value ? blocked_α : @view(kkt_solution[1:nu, 1]))
                copyto!(update_rule.ψ, blocked_value ? blocked_ψ : @view(kkt_solution[nu+1:end, 1]))
                α = update_rule.α
                ψ = update_rule.ψ
                factor_backed_policy = get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" && structured_B
                if factor_backed_policy
                    β = blocked_value ? zeros(T, 0, 0) : @view(kkt_solution[1:nu, 2:nx+1])
                    ω = blocked_value ? zeros(T, 0, 0) : @view(kkt_solution[nu+1:end, 2:nx+1])
                    update_rule.β = zeros(T, 0, 0)
                    update_rule.ω = zeros(T, 0, 0)
                    update_rule.factor_policy = FactorBackedPolicy{T}(
                        F, copy(active_B_rows), copy(B_active), sparse(cx),
                        zeros(T, nu + nc), zeros(T, nx))
                else
                    size(update_rule.β) == (nu, nx) || (update_rule.β = zeros(T, nu, nx))
                    size(update_rule.ω) == (nc, nx) || (update_rule.ω = zeros(T, nc, nx))
                    copyto!(update_rule.β, @view kkt_solution[1:nu, 2:nx+1])
                    copyto!(update_rule.ω, @view kkt_solution[nu+1:end, 2:nx+1])
                    update_rule.factor_policy = nothing
                    β = update_rule.β
                    ω = update_rule.ω
                end
            elseif nc > 0
                α_β_z = ck \ (Z' * ([-Qû -B] - Ĥ * Y * α_β_y))
                α_β = Y * α_β_y + Z * α_β_z
                α = α_β[:, 1]
                β = α_β[:, 2:nx+1]

                ψ_ω = ((Y' * ([-Qû -B] - Ĥ * α_β))' / fk)'
                ψ = ψ_ω[:, 1]
                ω = ψ_ω[:, 2:nx+1]
            else
                α_β = ck \ [-Qû -B]
                α = α_β[:, 1]
                β = α_β[:, 2:nx+1]
                ψ = zeros(T, nc)
                ω = zeros(T, nc, nx)
            end

            if sparse_kkt && capture_this_kkt
                Serialization.serialize(ENV["FILTERDDP_CAPTURE_KKT"], (
                    K=K, rhs=captured_rhs, stage=t, nx=nx, nu=nu, nc=nc,
                    x=copy(x), u=copy(u), phi=copy(ϕ), zl=copy(zl), zu=copy(zu),
                    barrier_mu=μ, future_Vx=copy(V̂x), future_Vxx=copy(V̂xx),
                    next_state=copy(dynamics.f(x, u)), alpha=copy(α),
                    beta=copy(β), psi=copy(ψ), omega=copy(ω)))
            end

            # DIAGNOSTIC (not part of the paper's optimization stack): dump every
            # matrix/vector FilterDDP carries per stage, every FILTERDDP_PERIODIC_
            # CAPTURE_STRIDE-th outer iteration, for offline shape/compressibility
            # inspection. No effect unless FILTERDDP_PERIODIC_CAPTURE_DIR is set.
            #
            # Two modes, since raw dumps scale as O(nu*nx + nc*nx) per snapshot --
            # fine at ieee123/ieee2522 scale (KB-MB), multiple GB per snapshot at
            # large10k scale (nu~54665, nc~42303, nx~1020), which filled a disk
            # during the first large10k attempt. Set FILTERDDP_PERIODIC_CAPTURE_
            # INLINE_STATS=1 to compute the same shape/SVD-rank statistics that
            # analyze_periodic_capture.jl would have computed offline, right here
            # where the arrays already live in memory, and append one small CSV
            # row per snapshot instead of serializing the arrays themselves.
            if sparse_kkt && periodic_capture_active
                capture_dir = ENV["FILTERDDP_PERIODIC_CAPTURE_DIR"]
                mkpath(capture_dir)
                if get(ENV, "FILTERDDP_PERIODIC_CAPTURE_INLINE_STATS", "0") == "1"
                    function _optimal_rank(sv, tol)
                        isempty(sv) && return 0
                        total = sum(abs2, sv)
                        total == 0 && return 0
                        for r in 0:length(sv)
                            tail = r == length(sv) ? 0.0 : sum(abs2, @view sv[r+1:end])
                            sqrt(tail / total) <= tol && return r
                        end
                        return length(sv)
                    end
                    rhs_state = Matrix(@view periodic_rhs_snapshot[:, 2:end])
                    beta_d = Matrix(β)
                    omega_d = Matrix(ω)
                    sv_rhs = isempty(rhs_state) ? Float64[] : svdvals(rhs_state)
                    sv_beta = isempty(beta_d) ? Float64[] : svdvals(beta_d)
                    sv_omega = isempty(omega_d) ? Float64[] : svdvals(omega_d)
                    csv_path = joinpath(capture_dir, "periodic_capture_summary.csv")
                    header_needed = !isfile(csv_path)
                    open(csv_path, "a") do io
                        header_needed && println(io,
                            "iteration,stage,nx,nu,nc,K_shape,K_nnz,K_density,rhs_state_shape,beta_shape,omega_shape,alpha_len,psi_len,SigmaL_len,SigmaU_len,Vxx_shape,Vx_len,rhs_sv_max,rhs_sv_min,rhs_rank_1pct,rhs_rank_5pct,rhs_rank_10pct,beta_rank_1pct,beta_rank_5pct,omega_rank_1pct,omega_rank_5pct,rhs_cond,beta_cond,barrier_mu,reg")
                        @printf(io, "%d,%d,%d,%d,%d,%dx%d,%d,%.6e,%dx%d,%dx%d,%dx%d,%d,%d,%d,%d,%dx%d,%d,%.6e,%.6e,%d,%d,%d,%d,%d,%d,%d,%.3e,%.3e,%.3e,%.3e\n",
                            data.k, t, nx, nu, nc,
                            size(K,1), size(K,2), nnz(K), nnz(K) / (size(K,1)^2),
                            size(rhs_state,1), size(rhs_state,2),
                            size(beta_d,1), size(beta_d,2), size(omega_d,1), size(omega_d,2),
                            length(α), length(ψ), length(Σ_L), length(Σ_U),
                            size(V̂xx,1), size(V̂xx,2), length(V̂x),
                            isempty(sv_rhs) ? NaN : sv_rhs[1], isempty(sv_rhs) ? NaN : sv_rhs[end],
                            _optimal_rank(sv_rhs, 0.01), _optimal_rank(sv_rhs, 0.05), _optimal_rank(sv_rhs, 0.10),
                            _optimal_rank(sv_beta, 0.01), _optimal_rank(sv_beta, 0.05),
                            _optimal_rank(sv_omega, 0.01), _optimal_rank(sv_omega, 0.05),
                            (isempty(sv_rhs) || sv_rhs[end] == 0) ? Inf : sv_rhs[1]/sv_rhs[end],
                            (isempty(sv_beta) || sv_beta[end] == 0) ? Inf : sv_beta[1]/sv_beta[end],
                            μ, reg)
                    end
                else
                    capture_file = joinpath(capture_dir,
                        @sprintf("iter%04d_stage%03d.jls", data.k, t))
                    Serialization.serialize(capture_file, (
                        iteration=data.k, stage=t, nx=nx, nu=nu, nc=nc,
                        x=copy(x), u=copy(u),
                        K=copy(K),
                        rhs_multi=periodic_rhs_snapshot,
                        kkt_solution=copy(kkt_solution),
                        alpha=copy(α), beta=copy(Matrix(β)), psi=copy(ψ),
                        omega=copy(Matrix(ω)),
                        Sigma_L=copy(Σ_L), Sigma_U=copy(Σ_U),
                        Vxx_incoming=copy(V̂xx), Vx_incoming=copy(V̂x),
                        barrier_mu=μ, reg=reg))
                end
            end

            # update parameters of update rule for ineq. dual variables

            if sparse_kkt
                @. update_rule.χl = inv_ul * μ - zl - Σ_L * α
                @. update_rule.χu = inv_uu * μ - zu + Σ_U * α
                copyto!(update_rule.Σ_L, Σ_L)
                copyto!(update_rule.Σ_U, Σ_U)
                χl = update_rule.χl
                χu = update_rule.χu
            else
                χl = inv_ul .* μ - zl - Σ_L .* α
                χu = inv_uu .* μ - zu + Σ_U .* α
                solver.update[t] = UpdateRule{T, nx, nu, nc, nux, ncx}(α, ψ, β, ω, χl, χu, Σ_L, Σ_U)
            end

            # evaluate optimality dual error
            Lu = lu_ - zl + zu + fu' * λ
            if nc > 0
                Lu = Lu + cu' * ϕ
            end

            data.dual_inf = max(data.dual_inf, norm(Lu, Inf))
            z_norm += sum(zl)
            z_norm += sum(zu)
            ϕ_norm += norm(ϕ, 1)

            # Update return V derivatives for next timestep Vxx = C + β' * B + ω' cx
            if blocked_value
                V̂xx = C + blocked_Vxx
                V̂x = lx + blocked_Vx + fx' * V̂x
            elseif structured_B
                beta_active = Matrix(@view β[active_B_rows, :])
                beta_B = beta_active' * B_active
                V̂xx = C + beta_B
                V̂x = lx + β' * Qû + fx' * V̂x
            else
                beta_B = β' * B
                V̂xx = C + beta_B
                V̂x = lx + β' * Qû + fx' * V̂x
            end
            λ = lx + fx' * λ
            if nc > 0
                if !blocked_value
                    V̂xx = V̂xx + ω' * cx
                    V̂x = V̂x + ω' * c
                end
                V̂x = V̂x + cx' * ϕ
                λ = λ + cx' * ϕ
            end

            # evaluate sufficient decrease condition in forward pass
            data.expected_change_L += dot(Qû, α)
            nc > 0 && (data.expected_change_L += dot(c, ψ))
            update_s = (time_ns() - update_start_ns) / 1e9
            update_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - update_alloc_start : 0
            timing_diagnostic && @printf(
                "FILTERDDP_TIMING iteration=%d barrier_iteration=%d stage=%d derivative_s=%.9f first_order_s=%.9f second_order_s=%.9f algebra_s=%.9f kkt_assembly_s=%.9f factor_s=%.9f solve_s=%.9f update_s=%.9f K_nnz=%d rhs_cols=%d Vxx_nnz=%d\n",
                data.k, data.j, t, derivative_s, first_order_s, second_order_s, algebra_s, kkt_assembly_s, factor_s, solve_s,
                update_s, sparse_kkt ? nnz(K) : 0, sparse_kkt ? size(rhs, 2) : 0,
                count(!iszero, V̂xx))
            memory_diagnostic && @printf(
                "FILTERDDP_MEMORY stage=%d derivative_alloc_MiB=%.3f algebra_alloc_MiB=%.3f kkt_alloc_MiB=%.3f factor_alloc_MiB=%.3f solve_alloc_MiB=%.3f update_alloc_MiB=%.3f K_MiB=%.3f rhs_MiB=%.3f solution_MiB=%.3f beta_MiB=%.3f omega_MiB=%.3f bound_sens_MiB=%.3f Vxx_MiB=%.3f update_rule_MiB=%.3f maxrss_delta_MiB=%.3f\n",
                t, derivative_alloc_bytes / 2.0^20, algebra_alloc_bytes / 2.0^20,
                kkt_alloc_bytes / 2.0^20, factor_alloc_bytes / 2.0^20,
                solve_alloc_bytes / 2.0^20, update_alloc_bytes / 2.0^20,
                sparse_kkt ? Base.summarysize(K) / 2.0^20 : 0.0,
                sparse_kkt ? Base.summarysize(rhs) / 2.0^20 : 0.0,
                sparse_kkt && kkt_solution !== rhs ? Base.summarysize(kkt_solution) / 2.0^20 : 0.0,
                Base.summarysize(β) / 2.0^20, Base.summarysize(ω) / 2.0^20,
                (Base.summarysize(Σ_L) + Base.summarysize(Σ_U)) / 2.0^20,
                Base.summarysize(V̂xx) / 2.0^20,
                Base.summarysize(solver.update[t]) / 2.0^20,
                max(Sys.maxrss() - stage_maxrss_start, 0) / 2.0^20)
        end
        scaling_dual = max(options.s_max, (ϕ_norm + z_norm) / max(ni + nc * ocp.N, 1.0))  / options.s_max
        scaling_cs = max(options.s_max, z_norm / max(ni, 1.0))  / options.s_max
        data.dual_inf /= scaling_dual
        data.cs_inf_0 /= scaling_cs
        data.cs_inf_μ /= scaling_cs
        data.barrier_lagrangian_curr += data.objective
        data.status == 0 && break
    end
    data.reg_last = reg
    if feasibility_diagnostic && data.status == 0
        equality_rms = sqrt(equality_ssq / max(equality_count, 1))
        dynamics_rms = sqrt(dynamics_ssq / max(dynamics_count, 1))
        bound_rms = sqrt(bound_ssq / max(bound_count, 1))
        @printf("FILTERDDP_FEASIBILITY iteration=%d barrier_iteration=%d equality_count=%d equality_rms=%.12e equality_max=%.12e equality_worst_stage=%d equality_worst_index=%d dynamics_count=%d dynamics_rms=%.12e dynamics_max=%.12e dynamics_worst_stage=%d dynamics_worst_index=%d bound_count=%d bound_rms=%.12e bound_max=%.12e bound_worst_stage=%d bound_worst_index=%d bound_worst_kind=%s\n",
            data.k, data.j, equality_count, equality_rms, equality_max,
            equality_worst_stage, equality_worst_index,
            dynamics_count, dynamics_rms, dynamics_max,
            dynamics_worst_stage, dynamics_worst_index,
            bound_count, bound_rms, bound_max, bound_worst_stage,
            bound_worst_index, bound_worst_kind)
        flush(stdout)
    end
    data.status != 0 && (verbose && (@warn "Backward pass failure, unable to find an iteration matrix with correct inertia."))
end
