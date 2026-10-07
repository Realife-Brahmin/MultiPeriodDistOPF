# Schur complement of the stage KKT onto the battery rows, for
# FILTERDDP_BATTERY_SCHUR=1 (DDP4OPF.BATTERY_SCHUR_HOOK). MUMPS returns it
# directly (ICNTL(19)); UMFPACK has no such interface. Symmetric mode, so
# MUMPS reads the upper triangle of K. Checked against full solves on captured
# stages by battery_block_check.jl.
#
# The stage KKT has the same sparsity pattern at every stage and iteration
# (the KKT pattern cache holds it fixed), and MUMPS's analysis depends only on
# the pattern. So the analysis runs once and every later call only refactors
# (battery_schur_reuse_bench.jl: 0.87 -> 0.33 s per call at large10k); it
# re-runs if the pattern or the battery rows ever change. A too-small
# workspace estimate (INFOG(1) = -8 or -9) is retried with a larger relaxation
# ICNTL(14); MUMPS otherwise returns garbage silently.

import MUMPS   # not `using`: MUMPS exports solve!, which would shadow DDP4OPF.solve!
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

const _SCHUR_RELAX = Ref(400)
const _SCHUR_STATE = Ref{Any}(nothing)   # (m, colptr, rowval, E) of the analysed pattern

function _schur_analyse(K, E)
    st = _SCHUR_STATE[]
    isnothing(st) || MUMPS.finalize!(st.m)
    m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
    MUMPS.suppress_display!(m)
    MUMPS.set_icntl!(m, 8, 0; displaylevel=0)    # no scaling with a Schur complement
    MUMPS.set_icntl!(m, 14, _SCHUR_RELAX[]; displaylevel=0)
    MUMPS.associate_matrix!(m, K)
    MUMPS.set_schur_centralized_by_column!(m, E)
    m.job = MUMPS.ANALYZE
    MUMPS.invoke_mumps!(m)
    m.infog[1] < 0 && error("MUMPS analysis failed, INFOG(1) = $(m.infog[1])")
    st = (m=m, colptr=copy(K.colptr), rowval=copy(K.rowval), E=copy(E))
    _SCHUR_STATE[] = st
    println("BATTERY_SCHUR analysis n=$(size(K, 1)) nE=$(length(E))")
    return st
end

function mumps_battery_schur(K, E)
    E = collect(E)
    st = _SCHUR_STATE[]
    if isnothing(st) || st.E != E || st.colptr != K.colptr || st.rowval != K.rowval
        st = _schur_analyse(K, E)
    end
    MUMPS.associate_matrix!(st.m, K)
    while true
        st.m.job = MUMPS.FACTOR
        MUMPS.invoke_mumps!(st.m)
        status = st.m.infog[1]
        status >= 0 && return copy(MUMPS.get_schur_complement(st.m))
        status in (-8, -9) || error("MUMPS Schur complement failed, INFOG(1) = $status")
        _SCHUR_RELAX[] *= 2
        MUMPS.set_icntl!(st.m, 14, _SCHUR_RELAX[]; displaylevel=0)
        @printf("MUMPS workspace retry: ICNTL(14) = %d\n", _SCHUR_RELAX[])
    end
end

DDP4OPF.BATTERY_SCHUR_HOOK[] = mumps_battery_schur
println("BATTERY_SCHUR hook: MUMPS, analysis reused")
