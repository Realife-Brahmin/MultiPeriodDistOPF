# Ordering and solver benchmark on one captured FilterDDP stage KKT system.
#
# Question (R. Gupta, agenda 2026-09-25): does a symmetric-indefinite LDL'
# factorization, or a different fill-reducing ordering, beat the UMFPACK LU that
# FilterDDP's backward pass uses -- on the real workload, which is one
# factorization followed by an (n_x+1)-column solve?
#
# Every configuration factorizes the SAME captured K and solves the SAME captured
# right-hand side. Per configuration, after one discarded warm-up (compilation),
# `repeats` timed repetitions give medians of:
#   ordering_s    UMFPACK symbolic analysis / MUMPS analysis (JOB=1), incl. ordering
#   factor_s      numerical factorization only
#   solve1_s      solve with the first RHS column (the feed-forward term)
#   solvewide_s   solve with all n_x+1 columns, exactly as FilterDDP does it
# plus factor size, fill, memory, pivoting statistics and relative residuals
# ||K X - B||_F / ||B||_F.
#
# Solvers actually available here (see ma57_availability_*.txt for MA57/HSL):
#   UMFPACK  unsymmetric LU; the baseline is Julia's default control, i.e. exactly
#            FilterDDP's lu(K). Orderings: AMD/COLAMD, METIS, CHOLMOD, BEST, each
#            under the unsymmetric or symmetric strategy.
#   MUMPS    sym = 2, general symmetric-indefinite LDL' with 1x1/2x2 pivoting (the
#            same algorithm family as MA57, and the solver inside this machine's
#            Ipopt). Orderings via ICNTL(7): AMD, AMF, SCOTCH, PORD, METIS, QAMD,
#            automatic; plus METIS on the compressed graph (ICNTL(12)=2) and
#            Ipopt's own MUMPS option set.
# Sequential: BLAS, OpenMP and OpenBLAS are pinned to one thread by the caller.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/kkt_ordering_benchmark.jl \
#         <capture.jls> <system> <arm> <out.csv> [repeats]

using LinearAlgebra, SparseArrays, Serialization, Printf, Statistics, Dates
using MUMPS

const U = SparseArrays.UMFPACK
const LSS = SparseArrays.LibSuiteSparse

input, system, arm, out_csv = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
repeats = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 5
families = Set(split(get(ENV, "KKT_ORDER_FAMILIES", "umfpack,mumps"), ','))
only = get(ENV, "KKT_ORDER_ONLY", "")            # optional comma list of config names

BLAS.set_num_threads(1)
threads = @sprintf("blas=%d omp=%s openblas=%s julia=%d", BLAS.get_num_threads(),
                   get(ENV, "OMP_NUM_THREADS", "unset"), get(ENV, "OPENBLAS_NUM_THREADS", "unset"),
                   Threads.nthreads())

cap = deserialize(input)
K = SparseMatrixCSC{Float64, Int64}(cap.K)
B = Matrix{Float64}(cap.rhs)
n = size(K, 1); nu, nc, nx = cap.nu, cap.nc, cap.nx
nnzK = nnz(K); nnz_triu = nnz(triu(K))
asym = norm(K - K', 1) / norm(K, 1)
nnzH = nnz(K[1:nu, 1:nu]); nnz22 = nnz(K[nu+1:end, nu+1:end])
b1 = B[:, 1:1]
@printf("KKT %s  system=%s arm=%s  n=%d (nu=%d nc=%d)  nnz=%d  nnz(H)=%d  nnz(22)=%d  rhs=%dx%d  asym=%.2e  repeats=%d  %s\n",
        basename(input), system, arm, n, nu, nc, nnzK, nnzH, nnz22, size(B)..., asym, repeats, threads)
flush(stdout)

relres(X, R) = norm(K * X - R) / norm(R)
bandwidth(p, q) = begin
    ip = invperm(p); iq = invperm(q); I, J, _ = findnz(K)
    maximum(abs(ip[I[k]] - iq[J[k]]) for k in eachindex(I))
end
millions(v) = v < 0 ? -Float64(v) * 1e6 : Float64(v)

# Ordering-intrinsic fill. UMFPACK counts nonzeros while MUMPS counts dense-front
# storage, so their factor sizes are not directly comparable. This applies a pivot
# order symmetrically to K's pattern and counts the exact simplicial Cholesky fill
# of an SPD surrogate (|K| plus a dominant diagonal): no pivoting, no supernode
# amalgamation -- it measures the ordering alone. Reported as 2 nnz(L) / nnz(K), the
# scale of UMFPACK's (nnz(L)+nnz(U)) / nnz(K).
const Ssurr = let A = abs.(K); A + spdiagm(0 => vec(sum(A; dims=2)) .+ 1.0) end
struct_fill(order) = 2 * nnz(sparse(ldlt(Symmetric(Ssurr); perm=order).LD)) / nnzK

const rows = Vector{Dict{String, Any}}()
function record!(r)
    push!(rows, r)
    @printf("  %-28s ord %8.4f  fact %8.4f  solve1 %8.5f  wide %8.4f  f+w %8.4f | fill %6.2f  res1 %.1e  resW %.1e  %s\n",
            r["config"], r["ordering_s"], r["factor_s"], r["solve1_s"], r["solvewide_s"],
            r["factor_s"] + r["solvewide_s"], r["fill_ratio"], r["relres_1"], r["relres_wide"], r["note"])
    flush(stdout)
end
base_row() = Dict{String, Any}("system" => system, "arm" => arm, "capture" => basename(input),
    "n" => n, "nu" => nu, "nc" => nc, "nx" => nx, "rhs_cols" => size(B, 2), "nnz_K" => nnzK,
    "nnz_triu_K" => nnz_triu, "nnz_H" => nnzH, "nnz_22" => nnz22, "K_asymmetry" => asym,
    "threads" => threads, "repeats" => repeats)

# ------------------------------------------------------------------ UMFPACK --
const ORDER_NAME = Dict(0 => "CHOLMOD", 1 => "AMD/COLAMD", 2 => "GIVEN", 3 => "METIS",
                        4 => "BEST", 5 => "NONE", 6 => "USER", 7 => "METIS_GUARD")
const STRAT_NAME = Dict(0 => "AUTO", 1 => "UNSYMMETRIC", 3 => "SYMMETRIC")
function umf_control(strategy, ordering)
    c = U.get_umfpack_control(Float64, Int64)
    strategy === nothing || (c[LSS.UMFPACK_STRATEGY + 1] = strategy)
    ordering === nothing || (c[LSS.UMFPACK_ORDERING + 1] = ordering)
    return c
end

function bench_umfpack(name, control, note)
    to = Float64[]; tf = Float64[]; t1 = Float64[]; tw = Float64[]
    local F, X1, XW
    for r in 0:repeats
        GC.gc()
        F = U.UmfpackLU(K; control=copy(control))
        a = @elapsed U.umfpack_symbolic!(F, nothing)
        f = @elapsed U.umfpack_numeric!(F)
        X1 = copy(b1); s1 = @elapsed ldiv!(F, X1)
        XW = copy(B);  sw = @elapsed ldiv!(F, XW)
        r == 0 && continue
        push!(to, a); push!(tf, f); push!(t1, s1); push!(tw, sw)
    end
    info = F.info
    unit = info[LSS.UMFPACK_SIZE_OF_UNIT + 1]
    Lf, Uf = F.L, F.U
    nLU = nnz(Lf) + nnz(Uf)
    r = base_row()
    merge!(r, Dict("solver" => "UMFPACK", "config" => name, "factorization" => "LU",
        "ordering_requested" => ORDER_NAME[Int(control[LSS.UMFPACK_ORDERING + 1])],
        "strategy_requested" => STRAT_NAME[Int(control[LSS.UMFPACK_STRATEGY + 1])],
        "ordering_used" => ORDER_NAME[Int(info[LSS.UMFPACK_ORDERING_USED + 1])],
        "strategy_used" => STRAT_NAME[Int(info[LSS.UMFPACK_STRATEGY_USED + 1])],
        "ordering_s" => median(to), "factor_s" => median(tf), "solve1_s" => median(t1),
        "solvewide_s" => median(tw), "ordering_s_min" => minimum(to), "factor_s_min" => minimum(tf),
        "solvewide_s_min" => minimum(tw),
        "factor_entries" => nLU, "factor_entries_def" => "nnz(L)+nnz(U) incl. diagonals",
        "fill_ratio" => nLU / nnzK,
        "factor_memory_MB" => info[LSS.UMFPACK_NUMERIC_SIZE + 1] * unit / 2^20,
        "peak_memory_MB" => info[LSS.UMFPACK_PEAK_MEMORY + 1] * unit / 2^20,
        "flops" => info[LSS.UMFPACK_FLOPS + 1],
        "pivot_stats" => @sprintf("offdiag_pivots=%d", Int(info[LSS.UMFPACK_NOFF_DIAG + 1])),
        "bandwidth_permuted" => bandwidth(F.p, F.q),
        "struct_fill_sym_order" => struct_fill(F.q),
        "relres_1" => relres(X1, b1), "relres_wide" => relres(XW, B),
        "status" => F.status == 0 ? "ok" : "umfpack_status_$(F.status)", "note" => note))
    record!(r)
end

if "umfpack" in families
    S_U, S_S = LSS.UMFPACK_STRATEGY_UNSYMMETRIC, LSS.UMFPACK_STRATEGY_SYMMETRIC
    for (name, strat, ord, note) in (
        ("UMFPACK_default",       nothing, nothing, "FilterDDP baseline: lu(K), Julia default control"),
        ("UMFPACK_unsym_colamd",  S_U, LSS.UMFPACK_ORDERING_AMD,     "COLAMD on A'A"),
        ("UMFPACK_sym_amd",       S_S, LSS.UMFPACK_ORDERING_AMD,     "AMD on A+A', diagonal preference"),
        ("UMFPACK_unsym_metis",   S_U, LSS.UMFPACK_ORDERING_METIS,   "METIS on A'A"),
        ("UMFPACK_sym_metis",     S_S, LSS.UMFPACK_ORDERING_METIS,   "METIS on A+A'"),
        ("UMFPACK_cholmod",       nothing, LSS.UMFPACK_ORDERING_CHOLMOD, "SuiteSparse C default: AMD/COLAMD, METIS if fill is high"),
        ("UMFPACK_best",          nothing, LSS.UMFPACK_ORDERING_BEST,    "tries AMD/COLAMD, METIS, NESDIS; keeps the least fill"))
        isempty(only) || name in split(only, ',') || continue
        try
            bench_umfpack(name, umf_control(strat, ord), note)
        catch err
            r = base_row(); merge!(r, Dict("solver" => "UMFPACK", "config" => name, "status" => "failed",
                "note" => first(split(sprint(showerror, err), '\n'))))
            push!(rows, r); @printf("  %-28s FAILED: %s\n", name, r["note"])
        end
    end
end

# -------------------------------------------------------------------- MUMPS --
const MUMPS_ORDER = Dict(0 => "AMD", 1 => "USER", 2 => "AMF", 3 => "SCOTCH", 4 => "PORD",
                         5 => "METIS", 6 => "QAMD", 7 => "AUTO")
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

function mumps_make(icntl_mods, cntl_mods)
    icntl = MUMPS.get_icntl()
    m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, icntl, MUMPS.default_cntl64)
    MUMPS.suppress_display!(m)
    for (i, v) in icntl_mods; MUMPS.set_icntl!(m, i, v; displaylevel=0); end
    for (i, v) in cntl_mods;  MUMPS.set_cntl!(m, i, v; displaylevel=0);  end
    MUMPS.associate_matrix!(m, K)
    return m
end
function mumps_phase!(m, job)
    MUMPS.set_job!(m, job); MUMPS.invoke_mumps!(m)
    m.infog[1] < 0 && error("MUMPS $(job) INFOG(1)=$(m.infog[1]) INFOG(2)=$(m.infog[2])")
    return m
end

function bench_mumps(name, icntl_mods, cntl_mods, note)
    to = Float64[]; tf = Float64[]; t1 = Float64[]; tw = Float64[]
    local m, X1, XW, perm, stats
    mem = 20
    for r in 0:repeats
        GC.gc()
        m = mumps_make(vcat(icntl_mods, [14 => mem]), cntl_mods)
        a = @elapsed mumps_phase!(m, MUMPS.ANALYZE)
        f = try
            @elapsed mumps_phase!(m, MUMPS.FACTOR)
        catch err
            # INFOG(1) = -8/-9: workspace estimate too small (delayed pivots).
            # Enlarge ICNTL(14) and restart the repetition; recorded in the note.
            (m.infog[1] in (-8, -9) && mem < 1000 && r == 0) || rethrow()
            finalize(m); mem *= 5; m = mumps_make(vcat(icntl_mods, [14 => mem]), cntl_mods)
            a = @elapsed mumps_phase!(m, MUMPS.ANALYZE)
            @elapsed mumps_phase!(m, MUMPS.FACTOR)
        end
        X1 = copy(b1); MUMPS.associate_rhs!(m, X1; unsafe=true)
        s1 = @elapsed mumps_phase!(m, MUMPS.SOLVE)
        XW = copy(B);  MUMPS.associate_rhs!(m, XW; unsafe=true)
        sw = @elapsed mumps_phase!(m, MUMPS.SOLVE)
        if r == repeats
            perm = copy(unsafe_wrap(Array, m.sym_perm, n))
            stats = (ordering_used = MUMPS_ORDER[Int(m.infog[7])], neg = m.infog[12],
                     delayed = m.infog[13], null = m.infog[28], entries = millions(m.infog[29]),
                     est_entries = millions(m.infog[20]), memMB = m.infog[22],
                     flops = m.rinfog[3])
        end
        finalize(m)
        r == 0 && continue
        push!(to, a); push!(tf, f); push!(t1, s1); push!(tw, sw)
    end
    r = base_row()
    merge!(r, Dict("solver" => "MUMPS", "config" => name, "factorization" => "LDLT",
        "ordering_requested" => join(["ICNTL($i)=$v" for (i, v) in icntl_mods], " "),
        "strategy_requested" => "sym=2", "ordering_used" => stats.ordering_used, "strategy_used" => "sym=2",
        "ordering_s" => median(to), "factor_s" => median(tf), "solve1_s" => median(t1),
        "solvewide_s" => median(tw), "ordering_s_min" => minimum(to), "factor_s_min" => minimum(tf),
        "solvewide_s_min" => minimum(tw),
        "factor_entries" => stats.entries,
        "factor_entries_def" => "INFOG(29): entries in L and D (one triangle)",
        # An LU with the same pattern stores L and U = D L', i.e. twice the one-sided
        # count; this makes fill_ratio comparable with UMFPACK's (nnz(L)+nnz(U))/nnz(K).
        "fill_ratio" => 2 * stats.entries / nnzK,
        "fill_ratio_one_sided" => stats.entries / nnz_triu,
        "factor_memory_MB" => stats.memMB, "peak_memory_MB" => stats.memMB, "flops" => stats.flops,
        "pivot_stats" => @sprintf("negative=%d delayed=%d null=%d est_entries=%.0f", stats.neg,
                                  stats.delayed, stats.null, stats.est_entries),
        "inertia_ok" => stats.neg == nc && stats.null == 0,
        # SYM_PERM(i) is the position of variable i in the pivot order.
        "bandwidth_permuted" => bandwidth(invperm(perm), invperm(perm)),
        "struct_fill_sym_order" => struct_fill(invperm(perm)),
        "relres_1" => relres(X1, b1), "relres_wide" => relres(XW, B),
        "status" => "ok", "note" => note * (mem > 20 ? " (ICNTL(14)=$mem)" : "")))
    record!(r)
    return perm
end

if "mumps" in families
    ipopt_like = [7 => 7, 8 => 77, 6 => 7]              # Ipopt: pivot_order, scaling, permuting_scaling
    for (name, icntl, cntl, note) in (
        ("MUMPS_sym_auto",    [7 => 7], [], "MUMPS default ordering choice"),
        ("MUMPS_sym_amd",     [7 => 0], [], "AMD"),
        ("MUMPS_sym_amf",     [7 => 2], [], "approximate minimum fill"),
        ("MUMPS_sym_qamd",    [7 => 6], [], "AMD with quasi-dense row detection"),
        ("MUMPS_sym_pord",    [7 => 4], [], "PORD"),
        ("MUMPS_sym_scotch",  [7 => 3], [], "SCOTCH nested dissection"),
        ("MUMPS_sym_metis",   [7 => 5], [], "METIS nested dissection"),
        ("MUMPS_sym_metis_compressed", [7 => 5, 12 => 2], [], "METIS on the 2x2-compressed graph"),
        # ICNTL(12)=0 (automatic) compresses the graph into 2x2 pivot candidates on
        # these zero-diagonal KKT matrices; =1 orders the plain graph instead.
        ("MUMPS_sym_amd_plain",   [7 => 0, 12 => 1], [], "AMD on the uncompressed graph"),
        ("MUMPS_sym_amf_plain",   [7 => 2, 12 => 1], [], "AMF on the uncompressed graph"),
        ("MUMPS_sym_metis_plain", [7 => 5, 12 => 1], [], "METIS on the uncompressed graph"),
        ("MUMPS_sym_ipopt",   ipopt_like, [1 => 1e-6], "Ipopt's MUMPS options (pivtol 1e-6)"))
        isempty(only) || name in split(only, ',') || continue
        try
            bench_mumps(name, icntl, cntl, note)
        catch err
            r = base_row(); merge!(r, Dict("solver" => "MUMPS", "config" => name, "status" => "failed",
                "note" => first(split(sprint(showerror, err), '\n'))))
            push!(rows, r); @printf("  %-28s FAILED: %s\n", name, r["note"])
        end
    end
end

cols = ["system", "arm", "capture", "n", "nu", "nc", "nx", "rhs_cols", "nnz_K", "nnz_triu_K", "nnz_H",
        "nnz_22", "K_asymmetry", "solver", "config", "factorization", "strategy_requested",
        "ordering_requested", "strategy_used", "ordering_used", "threads", "repeats",
        "ordering_s", "factor_s", "solve1_s", "solvewide_s", "ordering_s_min", "factor_s_min",
        "solvewide_s_min", "factor_plus_wide_s", "total_s", "factor_entries", "factor_entries_def",
        "fill_ratio", "fill_ratio_one_sided", "factor_memory_MB", "peak_memory_MB", "flops",
        "pivot_stats", "inertia_ok", "bandwidth_permuted", "struct_fill_sym_order",
        "relres_1", "relres_wide", "status", "note"]
csvcell(v) = v isa AbstractString ? (occursin(r"[,\"]", v) ? "\"" * replace(v, "\"" => "\"\"") * "\"" : v) :
             v isa AbstractFloat ? @sprintf("%.6g", v) : string(v)
mkpath(dirname(out_csv))
open(out_csv, "w") do io
    println(io, join(cols, ","))
    for r in rows
        if haskey(r, "factor_s")
            r["factor_plus_wide_s"] = r["factor_s"] + r["solvewide_s"]
            r["total_s"] = r["ordering_s"] + r["factor_s"] + r["solvewide_s"]
        end
        println(io, join((csvcell(get(r, c, "")) for c in cols), ","))
    end
end
println("wrote $out_csv  ($(length(rows)) rows)  finished=$(now())")
MUMPS.MPI.Finalized() || MUMPS.MPI.Finalize()
