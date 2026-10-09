# How much of MUMPS's Schur-complement call is the (pattern-only) analysis?
# Analyse once on the first captured stage, then refactor every capture with
# the same pattern, against a fresh analysis+factorization per call. Checks
# that the reused-analysis Schur complement equals the fresh one.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/battery_schur_reuse_bench.jl <captures...>

using LinearAlgebra, Printf, Serialization, SparseArrays
import MUMPS
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

function new_instance(relax)
    m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
    MUMPS.suppress_display!(m)
    MUMPS.set_icntl!(m, 8, 0; displaylevel=0)
    MUMPS.set_icntl!(m, 14, relax; displaylevel=0)
    return m
end

function fresh_schur(K, E; relax=400)
    m = new_instance(relax)
    MUMPS.associate_matrix!(m, K)
    MUMPS.mumps_schur_complement!(m, E)            # analysis + factorization
    m.infog[1] < 0 && error("fresh: INFOG(1) = $(m.infog[1])")
    S = MUMPS.get_schur_complement(m); MUMPS.finalize!(m)
    return S
end

caps = [deserialize(f) for f in ARGS]
E = findall(i -> any(!iszero, @view caps[1].rhs_multi[i, 2:end]), 1:size(caps[1].K, 1))
fresh_schur(caps[1].K, E)                            # compile

m = new_instance(400)
MUMPS.associate_matrix!(m, caps[1].K)
MUMPS.set_schur_centralized_by_column!(m, E)
m.job = MUMPS.ANALYZE                                 # analysis only
t_an = @elapsed MUMPS.invoke_mumps!(m)
m.infog[1] < 0 && error("analysis: INFOG(1) = $(m.infog[1])")
@printf("ANALYSIS_ONCE n=%d nE=%d analysis_s=%.3f\n", size(caps[1].K, 1), length(E), t_an)
for (k, d) in enumerate(caps)
    t_fresh = @elapsed S_fresh = fresh_schur(d.K, E)
    MUMPS.associate_matrix!(m, d.K)                  # same pattern, new values
    m.job = MUMPS.FACTOR                             # factorization only
    t_fact = @elapsed MUMPS.invoke_mumps!(m)
    m.infog[1] < 0 && error("refactor: INFOG(1) = $(m.infog[1])")
    S = MUMPS.get_schur_complement(m)
    @printf("REUSE capture=%d fresh_s=%.3f factor_only_s=%.3f rel_diff=%.2e\n",
            k, t_fresh, t_fact, norm(S - S_fresh) / norm(S_fresh))
end
MUMPS.finalize!(m)
