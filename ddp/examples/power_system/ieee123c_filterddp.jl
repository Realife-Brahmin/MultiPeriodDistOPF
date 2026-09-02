# Direct FilterDDP transcription of the tADMM workspace's IEEE123C BFM-SOCP.
#
# The network algebraic variables are FilterDDP controls. The only states are
# battery energies. Inequalities unsupported by FilterDDP are represented as
# box-bounded controls, or as equality constraints with nonnegative slacks.
#
# Run after export_ieee123c_data.jl:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/ieee123c_filterddp.jl \
#         [system] [T] [dimensions|build|solver|solve]

using FilterDDP
using JuMP
using LinearAlgebra
using Printf
using Serialization
using SparseArrays

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

function control_layout(data)
    N, L, B, D = length(data[:Nset]), length(data[:Lset]),
                  length(data[:Bset]), length(data[:Dset])
    k = 0
    take(n) = (r = (k+1):(k+n); k += n; r)
    idx = (ps=first(take(1)), qs=first(take(1)), P=take(L), Q=take(L),
           v=take(N), ell=take(L), pb=take(B), qnorm=take(D),
           soc_slack=take(L), energy_slack=take(B))
    return idx, k
end

function analytic_dynamics(nx, nu, pbidx, dt)
    f = (x,u) -> x .- dt .* u[pbidx]
    fx = (x,u) -> Matrix{Float64}(I, nx, nx)
    fu = function (x,u)
        J = spzeros(nx, nu)
        for b in 1:nx
            J[b, pbidx[b]] = -dt
        end
        J
    end
    zxx = (x,u,λ) -> zeros(nx, nx)
    zux = (x,u,λ) -> spzeros(nu, nx)
    zuu = (x,u,λ) -> spzeros(nu, nu)
    FilterDDP.Dynamics{nx,nu,typeof(f),typeof(fx),typeof(fu),typeof(zxx),typeof(zux),typeof(zuu)}(
        f, fx, fu, zxx, zux, zuu)
end

function analytic_objective(nx, nu, psidx, pbidx, price, pbase, dt, C_B)
    l = (x,u) -> [price*pbase*dt*u[psidx] + C_B*pbase^2*dt*sum(abs2, u[pbidx])]
    lx = (x,u) -> zeros(nx)
    lu = function (x,u)
        g = zeros(nu)
        g[psidx] = price*pbase*dt
        g[pbidx] .= 2C_B*pbase^2*dt .* u[pbidx]
        g
    end
    lxx = (x,u) -> zeros(nx, nx)
    lux = (x,u) -> spzeros(nu, nx)
    luu = function (x,u)
        H = spzeros(nu, nu)
        for k in pbidx
            H[k,k] = 2C_B*pbase^2*dt
        end
        H
    end
    FilterDDP.Objective{nx,nu,typeof(l),typeof(lx),typeof(lu),typeof(lxx),typeof(lux),typeof(luu)}(
        l, lx, lu, lxx, lux, luu)
end

function build_model(data)
    Nstage = data[:T]
    buses, lines = data[:Nset], data[:Lset]
    batteries, ders = data[:Bset], data[:Dset]
    root = data[:substationBus]
    nonroot = data[:Nm1set]
    idx, nu = control_layout(data)
    nx = length(batteries)

    buspos = Dict(j => k for (k,j) in enumerate(buses))
    linepos = Dict(e => k for (k,e) in enumerate(lines))
    batpos = Dict(j => k for (k,j) in enumerate(batteries))
    derpos = Dict(j => k for (k,j) in enumerate(ders))

    dt = data[:delta_t_h]
    pbase = data[:kVA_B]
    dyn = analytic_dynamics(nx, nu, idx.pb, dt)

    lower = fill(-Inf, nu); upper = fill(Inf, nu)
    lower[idx.ps] = 0.0
    for (k,j) in enumerate(buses)
        lower[idx.v[k]] = data[:Vminpu][j]^2
        upper[idx.v[k]] = data[:Vmaxpu][j]^2
    end
    lower[idx.ell] .= 0.0
    for (k,j) in enumerate(batteries)
        lower[idx.pb[k]] = -data[:P_B_R_pu][j]
        upper[idx.pb[k]] =  data[:P_B_R_pu][j]
        width = (data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j]
        lower[idx.energy_slack[k]] = 0.0
        upper[idx.energy_slack[k]] = width
    end
    lower[idx.qnorm] .= -1.0; upper[idx.qnorm] .= 1.0
    lower[idx.soc_slack] .= 0.0
    limits = ControlLimits(lower, upper)

    stage_objs = Any[]
    stage_cons = Any[]
    for t in 1:Nstage
        price = data[:LoadShapeCost][t]
        push!(stage_objs, analytic_objective(nx, nu, idx.ps, idx.pb,
            price, pbase, dt, data[:C_B]))

        qmax = Dict(j => sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2)) for j in ders)
        function equations(x,u)
            c = Any[]
            # Real-power balance: root, then every non-root bus.
            push!(c, u[idx.ps] - sum(u[idx.P[linepos[e]]] for e in data[:L1set]))
            for j in nonroot
                incoming = u[idx.P[linepos[(data[:parent][j],j)]]]
                outgoing = sum(u[idx.P[linepos[(j,k)]]] for k in data[:children][j]; init=0.0)
                pL = j in data[:NLset] ? data[:p_L_pu][j,t] : 0.0
                pD = j in ders ? data[:p_D_pu][j,t] : 0.0
                pb = j in batteries ? u[idx.pb[batpos[j]]] : 0.0
                incoming_line = (data[:parent][j],j)
                r = data[:rdict_pu][incoming_line]
                push!(c, outgoing - incoming + r*u[idx.ell[linepos[incoming_line]]] - pb - pD + pL)
            end
            # Reactive-power balance.
            push!(c, u[idx.qs] - sum(u[idx.Q[linepos[e]]] for e in data[:L1set]))
            for j in nonroot
                incoming = u[idx.Q[linepos[(data[:parent][j],j)]]]
                outgoing = sum(u[idx.Q[linepos[(j,k)]]] for k in data[:children][j]; init=0.0)
                qL = j in data[:NLset] ? data[:q_L_pu][j,t] : 0.0
                qD = j in ders ? qmax[j]*u[idx.qnorm[derpos[j]]] : 0.0
                incoming_line = (data[:parent][j],j)
                z = data[:xdict_pu][incoming_line]
                push!(c, outgoing - incoming + z*u[idx.ell[linepos[incoming_line]]] - qD + qL)
            end
            # BFM voltage drop.
            for (e,(i,j)) in enumerate(lines)
                r, z = data[:rdict_pu][(i,j)], data[:xdict_pu][(i,j)]
                push!(c, u[idx.v[buspos[j]]] - u[idx.v[buspos[i]]] +
                         2*(r*u[idx.P[e]] + z*u[idx.Q[e]]) - (r^2+z^2)*u[idx.ell[e]])
            end
            # SOCP inequality P^2+Q^2 <= v_i*ell, with nonnegative slack.
            for (e,(i,_)) in enumerate(lines)
                push!(c, u[idx.P[e]]^2 + u[idx.Q[e]]^2 -
                         u[idx.v[buspos[i]]]*u[idx.ell[e]] + u[idx.soc_slack[e]])
            end
            push!(c, u[idx.v[buspos[root]]] - 1.05^2)
            # Bound the post-action energy x - dt*P_B with a box slack.
            for (b,j) in enumerate(batteries)
                emin = data[:soc_min][j] * data[:B_R_pu][j]
                push!(c, x[b] - dt*u[idx.pb[b]] - emin - u[idx.energy_slack[b]])
            end
            c
        end
        nc = 2length(buses) + 2length(lines) + 1 + nx
        cx = function (x,u)
            J = spzeros(nc, nx)
            first_energy = 2length(buses) + 2length(lines) + 2
            for b in 1:nx
                J[first_energy+b-1, b] = 1.0
            end
            J
        end
        # The constraint Jacobian has a fixed sparsity pattern. Build its
        # constant entries once per stage and update only the four
        # iterate-dependent SOC entries per line on each callback.
        cu_J = spzeros(nc, nu)
        let J = cu_J
            row = 1
            J[row, idx.ps] = 1.0
            for e in data[:L1set]
                J[row, idx.P[linepos[e]]] -= 1.0
            end
            for j in nonroot
                row += 1
                incoming_line = (data[:parent][j],j)
                J[row, idx.P[linepos[incoming_line]]] -= 1.0
                J[row, idx.ell[linepos[incoming_line]]] += data[:rdict_pu][incoming_line]
                for k in data[:children][j]
                    J[row, idx.P[linepos[(j,k)]]] += 1.0
                end
                j in batteries && (J[row, idx.pb[batpos[j]]] -= 1.0)
            end
            row += 1
            J[row, idx.qs] = 1.0
            for e in data[:L1set]
                J[row, idx.Q[linepos[e]]] -= 1.0
            end
            for j in nonroot
                row += 1
                incoming_line = (data[:parent][j],j)
                J[row, idx.Q[linepos[incoming_line]]] -= 1.0
                J[row, idx.ell[linepos[incoming_line]]] += data[:xdict_pu][incoming_line]
                for k in data[:children][j]
                    J[row, idx.Q[linepos[(j,k)]]] += 1.0
                end
                j in ders && (J[row, idx.qnorm[derpos[j]]] -= qmax[j])
            end
            for (e,(i,j)) in enumerate(lines)
                row += 1
                r, z = data[:rdict_pu][(i,j)], data[:xdict_pu][(i,j)]
                J[row, idx.v[buspos[j]]] = 1.0
                J[row, idx.v[buspos[i]]] = -1.0
                J[row, idx.P[e]] = 2r
                J[row, idx.Q[e]] = 2z
                J[row, idx.ell[e]] = -(r^2+z^2)
            end
            for (e,(i,_)) in enumerate(lines)
                row += 1
                # Nonzero placeholders retain these positions in the CSC
                # structure until the first value update.
                J[row, idx.P[e]] = 1.0
                J[row, idx.Q[e]] = 1.0
                J[row, idx.v[buspos[i]]] = 1.0
                J[row, idx.ell[e]] = 1.0
                J[row, idx.soc_slack[e]] = 1.0
            end
            row += 1
            J[row, idx.v[buspos[root]]] = 1.0
            for b in 1:nx
                row += 1
                J[row, idx.pb[b]] = -dt
                J[row, idx.energy_slack[b]] = -1.0
            end
        end
        cu = let J = cu_J
            function (x,u)
                row = 2length(buses) + length(lines)
                for (e,(i,_)) in enumerate(lines)
                    row += 1
                    J[row, idx.P[e]] = 2u[idx.P[e]]
                    J[row, idx.Q[e]] = 2u[idx.Q[e]]
                    J[row, idx.v[buspos[i]]] = -u[idx.ell[e]]
                    J[row, idx.ell[e]] = -u[idx.v[buspos[i]]]
                end
                J
            end
        end
        cxx = (x,u,ϕ) -> zeros(nx, nx)
        cux = (x,u,ϕ) -> spzeros(nu, nx)
        cuu = function (x,u,ϕ)
            H = spzeros(nu, nu)
            first_soc = 2length(buses) + length(lines) + 1
            for (e,(i,_)) in enumerate(lines)
                weight = ϕ[first_soc+e-1]
                p, q = idx.P[e], idx.Q[e]
                v, ell = idx.v[buspos[i]], idx.ell[e]
                H[p,p] += 2weight
                H[q,q] += 2weight
                H[v,ell] -= weight
                H[ell,v] -= weight
            end
            H
        end
        con = FilterDDP.EqualityConstraints{nx,nu,nc,typeof(equations),typeof(cx),typeof(cu),typeof(cxx),typeof(cux),typeof(cuu)}(
            equations, cx, cu, cxx, cux, cuu)
        push!(stage_cons, con)
    end

    # FilterDDP has a control at its final stage; use that stage's ordinary cost.
    term = stage_objs[end]
    ocp = build_ocp(Nstage, stage_objs[1], term, dyn, stage_cons[1], limits;
                    stage_objectives=stage_objs, stage_constraints=stage_cons)
    return ocp, idx, nx, nu, length(stage_cons[1].c(zeros(nx), zeros(nu)))
end

function main(args=ARGS)
system = length(args) >= 1 ? args[1] : "ieee123C_1ph"
T = length(args) >= 2 ? parse(Int, args[2]) : 2
mode = length(args) >= 3 ? args[3] : "dimensions"
quiet = length(args) >= 4 && args[4] == "quiet"
datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                    "network_data_$(system)_T$(T).jls")
data = deserialize(datafile)
idx, nu = control_layout(data)
nx = length(data[:Bset])
nc = 2length(data[:Nset]) + 2length(data[:Lset]) + 1 + nx

@printf("%s T=%d: nx=%d, nu=%d, nc=%d, dense nu^2=%d, reduced=%d\n",
        system, T, nx, nu, nc, nu^2, nu-nc)
mode == "dimensions" && exit()

t0 = time()
ocp, idx, nx, nu, nc_actual = build_model(data)
@printf("build complete: %.3f s, nc=%d\n", time()-t0, nc_actual)
mode == "build" && exit()

println("constructing dynamic Solver storage...")
flush(stdout)
t_solver = time()
max_iterations = parse(Int, get(ENV, "FILTERDDP_MAX_ITERATIONS", "200"))
optimality_tolerance = parse(Float64, get(ENV, "FILTERDDP_OPTIMALITY_TOLERANCE", "1e-7"))
solver = Solver(ocp; options=Options{Float64}(verbose=!quiet,
    optimality_tolerance=optimality_tolerance, max_iterations=max_iterations))
@printf("optimality_tolerance=%.3e max_iterations=%d\n", optimality_tolerance, max_iterations)
@printf("Solver construction complete: %.3f s\n", time()-t_solver)
mode == "solver" && exit()
x0 = Float64[data[:B0_pu][j] for j in data[:Bset]]
u0 = zeros(nu)
u0[idx.ps] = 1.0
u0[idx.v] .= 1.0
u0[idx.ell] .= 1e-3
u0[idx.soc_slack] .= 1e-3
for (b,j) in enumerate(data[:Bset])
    emin = data[:soc_min][j] * data[:B_R_pu][j]
    u0[idx.energy_slack[b]] = x0[b] - emin
end
ubar = [copy(u0) for _ in 1:T]
t1 = time()
println("entering solve! ...")
flush(stdout)
status = solve!(solver, x0, ubar)
@printf("solve complete: %.3f s, iterations=%d, status=%s\n",
        time()-t1, solver.data.k, string(status))
@printf("final residuals: primal=%.12e dual=%.12e complementarity=%.12e\n",
        solver.data.primal_inf, solver.data.dual_inf, solver.data.cs_inf_0)

# Compare against the existing centralized JuMP/Ipopt result when available.
xddp, uddp = get_trajectory(solver)
solutionfile = joinpath(REPO, "ddp", "results", "network_filterddp",
                        "filterddp_solution_$(system)_T$(T).jls")
serialize(solutionfile, Dict(
    :system => system,
    :T => T,
    :status => status,
    :iterations => solver.data.k,
    :x => xddp,
    :u => uddp,
))
println("saved_solution=$solutionfile")
flush(stdout)
Jddp = 0.0
max_eq = 0.0
for t in 1:T
    Jddp += data[:LoadShapeCost][t]*data[:kVA_B]*data[:delta_t_h]*uddp[t][idx.ps]
    Jddp += data[:C_B]*data[:kVA_B]^2*data[:delta_t_h]*sum(uddp[t][k]^2 for k in idx.pb)
    max_eq = max(max_eq, norm(ocp.stage_constraints[t].c(xddp[t], uddp[t]), Inf))
end
@printf("FilterDDP objective=%.12f max_equality_residual=%.3e\n", Jddp, max_eq)

reffile = joinpath(REPO, "envs", "tadmm", "processedData", "$(system)_T$(T)", "sol_socp_bf.jls")
if isfile(reffile)
    ref = deserialize(reffile)
    maxdiff = Dict(:P_Subs=>0.0, :Q_Subs=>0.0, :P=>0.0, :Q=>0.0,
                   :v=>0.0, :ell=>0.0, :P_B=>0.0, :B=>0.0, :q_D=>0.0)
    for t in 1:T
        maxdiff[:P_Subs] = max(maxdiff[:P_Subs], abs(uddp[t][idx.ps] - ref[:P_Subs][t]))
        maxdiff[:Q_Subs] = max(maxdiff[:Q_Subs], abs(uddp[t][idx.qs] - ref[:Q_Subs][t]))
        for (e,line) in enumerate(data[:Lset])
            maxdiff[:P] = max(maxdiff[:P], abs(uddp[t][idx.P[e]] - ref[:P][line,t]))
            maxdiff[:Q] = max(maxdiff[:Q], abs(uddp[t][idx.Q[e]] - ref[:Q][line,t]))
            maxdiff[:ell] = max(maxdiff[:ell], abs(uddp[t][idx.ell[e]] - ref[:ℓ][line,t]))
        end
        for (k,j) in enumerate(data[:Nset])
            maxdiff[:v] = max(maxdiff[:v], abs(uddp[t][idx.v[k]] - ref[:v][j,t]))
        end
        for (b,j) in enumerate(data[:Bset])
            maxdiff[:P_B] = max(maxdiff[:P_B], abs(uddp[t][idx.pb[b]] - ref[:P_B][j,t]))
            Bpost = xddp[t][b] - data[:delta_t_h]*uddp[t][idx.pb[b]]
            maxdiff[:B] = max(maxdiff[:B], abs(Bpost - ref[:B][j,t]))
        end
        for (d,j) in enumerate(data[:Dset])
            qmax = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
            maxdiff[:q_D] = max(maxdiff[:q_D], abs(qmax*uddp[t][idx.qnorm[d]] - ref[:q_D][j,t]))
        end
    end
    @printf("Ipopt objective=%.12f objective_gap=%.3e\n", ref[:objective], abs(Jddp-ref[:objective]))
    for key in (:P_Subs,:Q_Subs,:P,:Q,:v,:ell,:P_B,:B,:q_D)
        @printf("max_abs_diff[%s]=%.3e\n", string(key), maxdiff[key])
    end
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
