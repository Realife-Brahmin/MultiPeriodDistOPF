# Linear-solver sweep on a captured per-stage KKT system.
#
# Supersedes benchmark_captured_kkt.jl, which compared UMFPACK against MUMPS but
# never recorded or controlled a thread count -- so its numbers cannot answer
# "what about MUMPS on multiple threads?" (R. Gupta, agenda 2026-09-18).  Every
# row here logs OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and BLAS.get_num_threads().
#
# Four families are swept, because the KKT K = [H cu'; cu 0] has properties the
# default configuration does not exploit:
#
#   umfpack        default: unsymmetric LU, COLAMD, pivot tolerance 0.1
#   umfpack_sym    UMFPACK_STRATEGY_SYMMETRIC -- K is structurally symmetric, so
#                  AMD on A+A' with diagonal preference should cut fill.  Exact.
#   umfpack_pivot  the pivot tolerance swept down.  Factorisation cost on an
#                  identically sparse K rises 6.4x once batteries cycle while
#                  nnz(K) is unchanged and nnz(LU) rises only 36%; delayed pivots
#                  are the standing hypothesis and this is the knob for them.
#   mumps          unsymmetric and symmetric, with ICNTL(16) OpenMP threads.
#
# Usage:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/benchmark_kkt_solvers.jl <capture.jls> [repeats] [out.csv]
# Set KKT_BENCH_FAMILIES to a comma list to restrict (default all).

using LinearAlgebra
using MUMPS
using Printf
using Serialization
using SparseArrays

const UMF = SparseArrays.UMFPACK
# 1-based Julia indices into UMFPACK's control vector.
const CTRL_STRATEGY = 6          # C UMFPACK_STRATEGY (0-based 5)
const STRATEGY_SYMMETRIC = 3.0   # C UMFPACK_STRATEGY_SYMMETRIC

input   = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
repeats = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
out_csv = length(ARGS) >= 3 ? ARGS[3] : replace(input, r"\.jls$" => "_solver_sweep.csv")
families = Set(split(get(ENV, "KKT_BENCH_FAMILIES",
                         "umfpack,umfpack_sym,umfpack_pivot,mumps"), ','))

payload = deserialize(input)
K, rhs = payload.K, payload.rhs
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

omp   = get(ENV, "OMP_NUM_THREADS", "unset")
obl   = get(ENV, "OPENBLAS_NUM_THREADS", "unset")
blas  = BLAS.get_num_threads()
@printf("KKT %s  size=%s  nnz=%d  rhs=%s  repeats=%d\n",
        basename(input), string(size(K)), nnz(K), string(size(rhs)), repeats)
@printf("threads: OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s BLAS.get_num_threads=%d CPU_THREADS=%d\n",
        omp, obl, blas, Sys.CPU_THREADS)

middle(xs) = sort(xs)[cld(length(xs), 2)]
relres(X) = norm(K * X - rhs) / norm(rhs)

rows = NamedTuple[]
function record!(name, factor_s, solve_s, X, nlu, note)
    r = (solver=name, omp=omp, openblas=obl, blas_threads=blas,
         factor_s=factor_s, solve_s=solve_s, total_s=factor_s + solve_s,
         nnz_LU=nlu, relative_residual=relres(X), note=note)
    push!(rows, r)
    @printf("  %-22s factor %9.4f s  solve %9.4f s  total %9.4f s  nnz_LU %10d  res %.3e  %s\n",
            name, r.factor_s, r.solve_s, r.total_s, nlu, r.relative_residual, note)
    return r
end

function bench_umfpack(name, control, note)
    F = lu(K; control=control); X = F \ rhs            # warm / compile
    ft = Float64[]; st = Float64[]
    local Ft, Xt
    for _ in 1:repeats
        GC.gc()
        push!(ft, @elapsed Ft = lu(K; control=control))
        push!(st, @elapsed Xt = Ft \ rhs)
        X = Xt
    end
    record!(name, middle(ft), middle(st), X, nnz(Ft.L) + nnz(Ft.U), note)
end

if "umfpack" in families
    bench_umfpack("UMFPACK", UMF.get_umfpack_control(Float64, Int64), "default strategy, pivtol 0.1")
end

if "umfpack_sym" in families
    local c = UMF.get_umfpack_control(Float64, Int64)
    c[CTRL_STRATEGY] = STRATEGY_SYMMETRIC
    bench_umfpack("UMFPACK_symmetric", c, "STRATEGY_SYMMETRIC (AMD on A+A')")
end

if "umfpack_pivot" in families
    for pt in (1e-2, 1e-3, 1e-4, 1e-6)
        local c = UMF.get_umfpack_control(Float64, Int64)
        c[UMF.JL_UMFPACK_PIVOT_TOLERANCE] = pt
        bench_umfpack(@sprintf("UMFPACK_pivtol_%.0e", pt), c, "relaxed partial pivoting")
    end
end

if "mumps" in families
    icntl_threads = parse(Int, get(ENV, "KKT_BENCH_MUMPS_ICNTL16", "0"))
    for (sym, label) in ((MUMPS.mumps_unsymmetric, "MUMPS_unsymmetric"),
                         (MUMPS.mumps_symmetric,   "MUMPS_symmetric"))
        try
            make() = begin
                icntl = MUMPS.get_icntl()
                icntl_threads > 0 && (icntl[16] = icntl_threads)
                F = MUMPS.Mumps{eltype(K)}(sym, icntl, MUMPS.default_cntl64)
                MUMPS.associate_matrix!(F, K); MUMPS.suppress_display!(F)
                MUMPS.mumps_factorize!(F); F
            end
            F = make(); X = F \ rhs; finalize(F)       # warm
            ft = Float64[]; st = Float64[]
            local Ft, Xt
            for _ in 1:repeats
                GC.gc()
                push!(ft, @elapsed Ft = make())
                push!(st, @elapsed Xt = Ft \ rhs)
                X = Xt; finalize(Ft)
            end
            record!(label, middle(ft), middle(st), X, -1,
                    icntl_threads > 0 ? "ICNTL(16)=$icntl_threads" : "ICNTL(16) default")
        catch err
            @warn "$label failed" exception=(err, catch_backtrace())
        end
    end
end

open(out_csv, "w") do io
    println(io, "solver,omp_num_threads,openblas_num_threads,blas_threads,factor_s,solve_s,total_s,nnz_LU,relative_residual,note")
    for r in rows
        @printf(io, "%s,%s,%s,%d,%.9f,%.9f,%.9f,%d,%.6e,%s\n",
                r.solver, r.omp, r.openblas, r.blas_threads,
                r.factor_s, r.solve_s, r.total_s, r.nnz_LU, r.relative_residual, r.note)
    end
end
println("wrote $out_csv")
MUMPS.MPI.Finalized() || MUMPS.MPI.Finalize()
