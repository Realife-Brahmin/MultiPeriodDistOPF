# The battery value recursion by pairwise combination of the periods
# (parallel in time), after Sarkka and Garcia-Fernandez, "Temporal
# Parallelization of Dynamic Programming and Linear Quadratic Control", IEEE
# TAC 68(2), 2023 (arXiv 2104.03186). Diagnostic only (FILTERDDP_PINT_SCAN=1):
# the solve itself is the ordinary sequential one.
#
# Needs the exact battery-block Hessian (no FILTERDDP_DIAG_HESSIAN), the tree
# solver and FILTERDDP_HALF_BLOCK=1. Once the network of a stage is eliminated,
# what the backward sweep solves is a linear-quadratic problem in the battery
# energies. With a = dx (state before the stage), dp the battery powers and
# b = a + f .* dp the state after it, the stage leaves
#     q(a, dp) = 1/2 dp' S0 dp - (rhat + Rhat a)' dp + 1/2 a' D a + d' a,
# where none of S0, rhat, Rhat, D, d depends on the next stage's value function
# (S0 is the halved battery block without the curvature term). The sweep is
#     V_new(a) = min_dp q(a, dp) + V(a + f .* dp).
# Instead, each stage is written as a function of both of its end states,
#     G(a, b) = q(a, (b - a) ./ f),
# two neighbours are combined by minimising over the state between them, and
# the combinations are done pairwise: about log2(T) rounds instead of T steps.
#
# On the selected backward passes this file (1) redoes the sweep from the
# stage data and compares it with the solver's value functions, (2) does the
# same by pairwise combination, and (3) times both. The parallel time is a
# count: the combinations of a round are run one after another and the slowest
# is taken, as in parsim_from_log.py.

using LinearAlgebra, SparseArrays, Printf

struct PintLeaf
    t::Int
    S0::Matrix{Float64}
    rhat::Vector{Float64}
    Rhat::SparseMatrixCSC{Float64, Int}     # n_B x n_x
    D::SparseMatrixCSC{Float64, Int}        # n_x x n_x
    d::Vector{Float64}
    f::Vector{Float64}                      # b = a + f .* dp
    Vx_out::Vector{Float64}                 # the solver's value function after this stage
    Vxx_out::Matrix{Float64}
    alphaP::Vector{Float64}                 # the solver's battery-power step
    seq_s::Float64                          # the solver's sequential time for this stage
end

# G(a, b) = 1/2 a' Maa a + a' Mab b + 1/2 b' Mbb b + ma' a + mb' b
struct PintElem
    Maa::Matrix{Float64}; Mab::Matrix{Float64}; Mbb::Matrix{Float64}
    ma::Vector{Float64}; mb::Vector{Float64}
end
# V(a) = 1/2 a' V a + v' a
struct PintVal
    V::Matrix{Float64}; v::Vector{Float64}
end

const _PINT_LEAVES = PintLeaf[]
const _PINT_CHOL_FAIL = Ref(0)

function _pint_leaf(p)
    nB = length(p.active_B_rows); nu = p.nu
    fa = sparse(p.fu_active); nx = size(fa, 1)
    (nx == nB && all(k -> fa.colptr[k+1] == fa.colptr[k] + 1 && rowvals(fa)[fa.colptr[k]] == k, 1:nB)) ||
        error("pint_scan: the dynamics are not x + f .* p with one battery per state")
    f = [nonzeros(fa)[fa.colptr[k]] for k in 1:nB]
    cx_s = sparse(p.cx)
    energy_rows = sort!(unique(rowvals(cx_s)))
    p.F.lay.E == vcat(p.active_B_rows, nu .+ energy_rows) ||
        error("pint_scan: battery rows differ from the tree layout's")
    cxE = cx_s[energy_rows, :]
    S0, ninv, SPN = tree_kkt_block0(p.F, p.K, nB)
    rE = tree_kkt_reduced_rhs(p.F, p.rhs0)
    r0P = rE[1:nB] .+ f .* p.Vx_in               # the right-hand side without f' V_x
    r0N = rE[nB+1:end]
    Nd = Diagonal(ninv)
    rhat = r0P - SPN * (ninv .* r0N)
    L = sparse(p.lux)[p.active_B_rows, :]
    Rhat = sparse(-L + SPN * (Nd * cxE))
    D = sparse(sparse(p.lxx) - cxE' * (Nd * cxE))
    d = Vector(p.lx) + cx_s' * p.phi + cxE' * (ninv .* r0N)
    return PintLeaf(p.t, S0, rhat, Rhat, D, d, f, copy(p.Vx_out), p.Vxx_out,
                    p.alpha[p.active_B_rows], p.seq_s)
end

function pint_elem(lf::PintLeaf)
    p = 1 ./ lf.f; n = length(p)
    Np = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n; Np[i, j] = p[i] * lf.S0[i, j] * p[j]; end
    Rp = Matrix(lf.Rhat') .* p'                  # Rhat' * Diagonal(p)
    Maa = Np + Rp; Maa .+= Rp'; Maa .+= lf.D
    Mab = Np + Rp; Mab .*= -1
    pr = p .* lf.rhat
    return PintElem(Maa, Mab, Np, pr + lf.d, -pr)
end

# One step of the sweep from the stage data: (V, v) after the stage -> before it.
function pint_step(lf::PintLeaf, V::Matrix{Float64}, v::Vector{Float64})
    f = lf.f; n = length(f)
    S = copy(lf.S0)
    Rt = Matrix(lf.Rhat)
    @inbounds for j in 1:n, i in 1:n
        S[i, j] += f[i] * V[i, j] * f[j]
        Rt[i, j] -= f[i] * V[i, j]
    end
    rt = lf.rhat .- f .* v
    Vn = V + lf.D
    C = cholesky!(Symmetric(S); check=false)
    if issuccess(C)
        ldiv!(C.U', Rt)
        y = C.U' \ rt
        BLAS.syrk!('U', 'T', -1.0, Rt, 1.0, Vn); LinearAlgebra.copytri!(Vn, 'U')
        vn = lf.d + v - Rt' * y
        return PintVal(Vn, vn), C.U \ y
    end
    _PINT_CHOL_FAIL[] += 1
    S = copy(lf.S0)
    @inbounds for j in 1:n, i in 1:n; S[i, j] += f[i] * V[i, j] * f[j]; end
    Fl = lu!(S)
    X = Fl \ Rt; a = Fl \ rt
    mul!(Vn, Rt', X, -1.0, 1.0)
    return PintVal(Vn, lf.d + v - Rt' * a), a
end

# (a, z) combined with (z, b): minimise over z.
function pint_combine(x::PintElem, y::PintElem)
    K = x.Mbb + y.Maa
    g = x.mb + y.ma
    C = cholesky!(Symmetric(K); check=false)
    if issuccess(C)
        Y1 = copy(transpose(x.Mab)); ldiv!(C.U', Y1)
        Y2 = copy(y.Mab); ldiv!(C.U', Y2)
        w = C.U' \ g
        Maa = copy(x.Maa); BLAS.syrk!('U', 'T', -1.0, Y1, 1.0, Maa); LinearAlgebra.copytri!(Maa, 'U')
        Mbb = copy(y.Mbb); BLAS.syrk!('U', 'T', -1.0, Y2, 1.0, Mbb); LinearAlgebra.copytri!(Mbb, 'U')
        Mab = Matrix{Float64}(undef, size(K)); mul!(Mab, Y1', Y2, -1.0, 0.0)
        return PintElem(Maa, Mab, Mbb, x.ma - Y1' * w, y.mb - Y2' * w)
    end
    _PINT_CHOL_FAIL[] += 1
    Fl = lu!(x.Mbb + y.Maa)
    X1 = Fl \ Matrix(transpose(x.Mab)); X2 = Fl \ y.Mab; w = Fl \ g
    return PintElem(x.Maa - x.Mab * X1, -(x.Mab * X2), y.Mbb - y.Mab' * X2,
                    x.ma - x.Mab * w, y.mb - y.Mab' * w)
end
# (a, z) closed with a value function of z.
function pint_combine(x::PintElem, y::PintVal)
    K = x.Mbb + y.V
    g = x.mb + y.v
    C = cholesky!(Symmetric(K); check=false)
    if issuccess(C)
        Y1 = copy(transpose(x.Mab)); ldiv!(C.U', Y1)
        w = C.U' \ g
        Vn = copy(x.Maa); BLAS.syrk!('U', 'T', -1.0, Y1, 1.0, Vn); LinearAlgebra.copytri!(Vn, 'U')
        return PintVal(Vn, x.ma - Y1' * w)
    end
    _PINT_CHOL_FAIL[] += 1
    Fl = lu!(x.Mbb + y.V)
    X1 = Fl \ Matrix(transpose(x.Mab)); w = Fl \ g
    return PintVal(x.Maa - x.Mab * X1, x.ma - x.Mab * w)
end

# All the value functions by pairwise combination. xs[1:T] are the stages'
# two-ended functions and xs[T+1] the value function after the last stage. On
# return xs[t] is the value function before stage t. Level r combines, in every
# block of 2^(r+1) positions, each entry of the left half with the first entry
# of the right half; the combinations of a level are independent of each other.
# Returns the times of the full combinations and of the closing ones, by level.
function pint_scan!(xs::Vector{Any})
    n = length(xs); full = Vector{Float64}[]; closing = Vector{Float64}[]
    every = _pint_gc_every(size(xs[end].V, 1)); done = 0
    h = 1
    while h < n
        tf = Float64[]; tc = Float64[]
        for i in 1:2h:n
            mid = i + h
            mid > n && continue
            right = xs[mid]
            for t in i:mid-1
                dt, xs[t] = _pint_timed(() -> pint_combine(xs[t], right))
                push!(right isa PintVal ? tc : tf, dt)
                done += 1; done % every == 0 && GC.gc(false)
            end
        end
        push!(full, tf); push!(closing, tc)
        h *= 2
    end
    return full, closing
end

# Time one operation with the garbage collector off, so that a collection
# triggered by the harness's own allocations does not land in one combination
# (seen: 75-94 ms against a median of 5 ms on med2522 T=96). The collector is
# run between operations instead, every `_pint_gc_every` of them.
function _pint_timed(f)
    GC.enable(false)
    t0 = time_ns()
    r = f()
    dt = (time_ns() - t0) / 1e9
    GC.enable(true)
    return dt, r
end
_pint_gc_every(n::Int) = max(1, floor(Int, 1.5e9 / (48 * max(n, 1)^2)))

_pint_med(v) = isempty(v) ? 0.0 : sort(v)[(length(v) + 1) ÷ 2]
_pint_max(v) = isempty(v) ? 0.0 : maximum(v)
_pint_rel(A, B) = norm(A - B) / max(norm(B), floatmin())

function _pint_measure(leaves::Vector{PintLeaf})
    T = length(leaves); nB = length(leaves[1].f)
    _PINT_CHOL_FAIL[] = 0
    # (1) the sweep from the stage data
    GC.gc()
    V = PintVal(zeros(nB, nB), zeros(nB))
    ts = zeros(T); eV = 0.0; ev = 0.0; ea = 0.0
    every = _pint_gc_every(nB)
    for t in T:-1:1
        lf = leaves[t]
        Vp = V
        ts[t], (V, a) = _pint_timed(() -> pint_step(lf, Vp.V, Vp.v))
        eV = max(eV, _pint_rel(V.V, lf.Vxx_out)); ev = max(ev, _pint_rel(V.v, lf.Vx_out))
        ea = max(ea, _pint_rel(a, lf.alphaP))
        (T - t + 1) % every == 0 && GC.gc(false)
    end
    chol_sweep = _PINT_CHOL_FAIL[]; _PINT_CHOL_FAIL[] = 0
    # (2) the same by pairwise combination
    GC.gc()
    tl = zeros(T)
    xs = Vector{Any}(undef, T + 1)
    for t in 1:T
        tl[t], xs[t] = _pint_timed(() -> pint_elem(leaves[t]))
        t % every == 0 && GC.gc(false)
    end
    xs[T+1] = PintVal(zeros(nB, nB), zeros(nB))
    full, closing = pint_scan!(xs)
    sV = 0.0; sv = 0.0
    for t in 1:T
        sV = max(sV, _pint_rel(xs[t].V, leaves[t].Vxx_out)); sv = max(sv, _pint_rel(xs[t].v, leaves[t].Vx_out))
    end
    depth = count(i -> !isempty(full[i]) || !isempty(closing[i]), eachindex(full))
    lvl_max = sum(max(_pint_max(full[i]), _pint_max(closing[i])) for i in eachindex(full); init=0.0)
    lvl_med = sum(max(_pint_med(full[i]), _pint_med(closing[i])) for i in eachindex(full); init=0.0)
    allf = reduce(vcat, full; init=Float64[]); allc = reduce(vcat, closing; init=Float64[])
    # One worker per period: build the stage functions, the rounds, then every
    # stage's own battery solve with its value function, all stages at once.
    return (eV=eV, ev=ev, ea=ea, sV=sV, sv=sv, chol_sweep=chol_sweep, chol_scan=_PINT_CHOL_FAIL[],
            sweep_sum=sum(ts), step_med=_pint_med(ts), step_max=maximum(ts), depth=depth,
            merges=length(allf), closes=length(allc), merge_med=_pint_med(allf), merge_max=_pint_max(allf),
            close_med=_pint_med(allc), leaf_max=maximum(tl), rounds_max=lvl_max, rounds_med=lvl_med,
            work=sum(allf) + sum(allc),
            count_max=maximum(tl) + lvl_max + maximum(ts), count_med=_pint_med(tl) + lvl_med + _pint_med(ts))
end

function _pint_report(k::Int, mu, reg)
    leaves = sort(_PINT_LEAVES; by=lf -> lf.t)
    T = length(leaves)
    (T > 0 && [lf.t for lf in leaves] == collect(1:T)) || (println("PINT_SCAN_SKIP iteration=$k stages=$T"); return)
    nB = length(leaves[1].f)
    blas0 = BLAS.get_num_threads()
    solver_seq = sum(lf.seq_s for lf in leaves)
    for nth in unique((1, blas0))
        BLAS.set_num_threads(nth)
        m = _pint_measure(leaves)
        @printf("PINT_SCAN iteration=%d T=%d nB=%d blas_threads=%d mu=%.3e reg=%.3e sweep_err_Vxx=%.3e sweep_err_Vx=%.3e sweep_err_alpha=%.3e scan_err_Vxx=%.3e scan_err_Vx=%.3e chol_fail_sweep=%d chol_fail_scan=%d",
                k, T, nB, nth, mu, reg, m.eV, m.ev, m.ea, m.sV, m.sv, m.chol_sweep, m.chol_scan)
        @printf(" sweep_sum_s=%.6f sweep_step_med_s=%.6f sweep_step_max_s=%.6f depth=%d merges=%d closes=%d merge_med_s=%.6f merge_max_s=%.6f close_med_s=%.6f leaf_max_s=%.6f rounds_max_s=%.6f rounds_med_s=%.6f work_s=%.6f",
                m.sweep_sum, m.step_med, m.step_max, m.depth, m.merges, m.closes, m.merge_med, m.merge_max,
                m.close_med, m.leaf_max, m.rounds_max, m.rounds_med, m.work)
        @printf(" count_max_s=%.6f count_med_s=%.6f speedup_max=%.3f speedup_med=%.3f solver_seq_sum_s=%.6f\n",
                m.count_max, m.count_med, m.sweep_sum / m.count_max, m.sweep_sum / m.count_med, solver_seq)
        flush(stdout)
    end
    BLAS.set_num_threads(blas0)
    return
end

# Compile the timed functions on a small random problem, so that the first
# measured pass is not a compilation.
function _pint_warm()
    n = 4; T = 5
    lv = PintLeaf[]
    for t in 1:T
        A = randn(n, n)
        push!(lv, PintLeaf(t, A' * A + n * I, randn(n), sparse(0.1I, n, n), sparse(1.0I, n, n), randn(n),
                           fill(-0.5, n), zeros(n), zeros(n, n), zeros(n), 0.0))
    end
    _pint_measure(lv)
    return nothing
end

function install_pint_scan_hook()
    _pint_warm()
    DDP4OPF.PINT_SCAN_HOOK[] = function (kind::Symbol, p)
        if kind === :reset
            empty!(_PINT_LEAVES)
        elseif kind === :stage
            push!(_PINT_LEAVES, _pint_leaf(p))
        elseif kind === :finish
            _pint_report(p.k, p.mu, p.reg)
            empty!(_PINT_LEAVES)
        end
        return nothing
    end
    return nothing
end
