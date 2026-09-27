# Can the stage KKT use a block-diagonal Jacobian? (Agenda, R. Gupta,
# 2026-10-02.) Splits the radial feeder into k areas (subtrees, cut greedily
# from the leaves so each holds about N/k buses), assigns every variable and
# every constraint row to an area, and drops the Jacobian entries that couple
# two areas. With the diagonal Hessian the KKT matrix then falls apart into k
# independent blocks (one factorization per area, in parallel). On a captured
# stage KKT, for each k it reports
#   dropped         Jacobian entries removed (both triangles)
#   dev_primal / dev_feedforward
#                   how far the (n_x+1)-column solution moves (control rows)
#   refine_steps    block-Jacobi iterative refinement, x += M \ (b - K x) with
#                   M the block-diagonal matrix: steps to a 1e-10 relative
#                   residual on the feedforward column (or "diverges"), i.e.
#                   whether the blocks work as a preconditioner even if not
#                   as a replacement
#   nnz_LU, factor_s, largest_block
#
#   julia --project=envs/ddp2026 ddp/examples/power_system/kkt_block_jacobian.jl \
#         <capture.jls> <network_data.jls> <system> <iter> <out.csv> [k list]

using Serialization, SparseArrays, LinearAlgebra, Printf, Statistics

capfile, datafile, system, iter, outcsv = ARGS[1:5]
klist = length(ARGS) >= 6 ? parse.(Int, split(ARGS[6], ',')) : [2, 4, 8, 16, 32, 64]
BLAS.set_num_threads(1)

c = deserialize(capfile)
K = SparseMatrixCSC{Float64, Int64}(c.K); B = Matrix{Float64}(c.rhs)
n = size(K, 1); nu, nc = c.nu, c.nc
d = deserialize(datafile)
buses, lines = d[:Nset], d[:Lset]; root = d[:substationBus]; nonroot = d[:Nm1set]
bats, ders = d[:Bset], d[:Dset]
N, L, nB, nD = length(buses), length(lines), length(bats), length(ders)
children = d[:children]

# Variable layout (control_layout in ieee123c_filterddp.jl) and row layout.
off = cumsum([0, 1, 1, L, L, N, L, nB, nD, L])     # ps qs P Q v ell pb qnorm soc_slack energy_slack
rng(k, m) = (off[k] + 1):(off[k] + m)
iP, iQ, iv, iell, ipb, iqn, isoc, ien = rng(3, L), rng(4, L), rng(5, N), rng(6, L), rng(7, nB), rng(8, nD), rng(9, L), rng(10, nB)

# Subtree sizes, then greedy cuts from the leaves.
order = Int[]                                     # post-order
stack = [(root, false)]
while !isempty(stack)
    j, done = pop!(stack)
    if done; push!(order, j); continue; end
    push!(stack, (j, true))
    for ch in get(children, j, Int[]); push!(stack, (ch, false)); end
end
function partition(k)
    target = N / k
    area = Dict{Int, Int}(); remaining = Dict{Int, Int}()
    nareas = 1
    for j in order
        remaining[j] = 1 + sum((remaining[ch] for ch in get(children, j, Int[]) if !haskey(area, ch)); init=0)
        if j != root && remaining[j] >= target && nareas < k
            nareas += 1
            # assign j's still-unassigned subtree to the new area
            st = [j]
            while !isempty(st)
                v = pop!(st)
                haskey(area, v) && continue
                area[v] = nareas
                append!(st, get(children, v, Int[]))
            end
            remaining[j] = 0
        end
    end
    for j in buses; haskey(area, j) || (area[j] = 1); end
    return area, nareas
end

function block_keep(area)
    busarea = [area[j] for j in buses]
    linearea = [area[j] for (_, j) in lines]
    colarea = zeros(Int, nu)
    colarea[1:2] .= area[root]
    colarea[iP] .= linearea; colarea[iQ] .= linearea; colarea[iell] .= linearea; colarea[isoc] .= linearea
    colarea[iv] .= busarea
    colarea[ipb] .= [area[j] for j in bats]; colarea[ien] .= [area[j] for j in bats]
    colarea[iqn] .= [area[j] for j in ders]
    rowarea = zeros(Int, nc)
    rowarea[1] = area[root]; rowarea[N+1] = area[root]
    for (k, j) in enumerate(nonroot); rowarea[1+k] = area[j]; rowarea[N+1+k] = area[j]; end
    rowarea[2N+1:2N+L] .= linearea; rowarea[2N+L+1:2N+2L] .= linearea
    rowarea[2N+2L+1] = area[root]
    rowarea[2N+2L+2:end] .= [area[j] for j in bats]
    a = vcat(colarea, rowarea)
    I, J, V = findnz(K)
    keep = [a[i] == a[j] for (i, j) in zip(I, J)]
    return sparse(I[keep], J[keep], V[keep], n, n), count(.!keep), a
end

lu(K); F0 = lu(K); X0 = F0 \ B
rel(A, R) = norm(A - R) / norm(R)
b1 = B[:, 1]; nb1 = norm(b1)
open(outcsv, "w") do io
    println(io, "system,iter,mu,k_requested,areas,dropped,dropped_frac,largest_block,nnz_LU,factor_s,dev_all,dev_primal,dev_feedforward_primal,refine_steps,refine_final_res")
    for k in klist
        k > N ÷ 4 && continue
        area, na = partition(k)
        Kb, dropped, a = block_keep(area)
        largest = maximum(count(==(g), a) for g in 1:na)
        line = try
            lu(Kb); tf = @elapsed Fb = lu(Kb)
            Xb = Fb \ B
            x = Fb \ b1; steps = 0; res = norm(b1 - K * x) / nb1
            while res > 1e-10 && steps < 200 && isfinite(res) && res < 1e6
                x .+= Fb \ (b1 - K * x); steps += 1
                res = norm(b1 - K * x) / nb1
            end
            @sprintf("%d,%.4f,%.3e,%.3e,%.3e,%s,%.3e", nnz(Fb.L) + nnz(Fb.U), tf, rel(Xb, X0),
                     rel(Xb[1:nu, :], X0[1:nu, :]), rel(Xb[1:nu, 1], X0[1:nu, 1]),
                     res <= 1e-10 ? string(steps) : "diverges", res)
        catch
            "singular,NaN,NaN,NaN,NaN,NaN,NaN"
        end
        @printf(io, "%s,%s,%.3e,%d,%d,%d,%.3e,%d,%s\n", system, iter, c.barrier_mu, k, na, dropped,
                dropped / nnz(K), largest, line)
        flush(io)
    end
end
println("KKT_BLOCK_JACOBIAN $system iter $iter done")
