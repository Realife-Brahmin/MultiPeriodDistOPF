# Schur complement of the stage KKT onto the battery rows, for
# FILTERDDP_BATTERY_SCHUR=1 (DDP4OPF.BATTERY_SCHUR_HOOK). MUMPS returns it
# directly (ICNTL(19)); UMFPACK has no such interface. Symmetric mode, so
# MUMPS reads the upper triangle of K. A too-small workspace estimate
# (INFOG(1) = -8 or -9) is retried with a larger relaxation ICNTL(14); it
# otherwise returns garbage silently. Checked against full solves on captured
# stages by battery_block_check.jl.

import MUMPS   # not `using`: MUMPS exports solve!, which would shadow DDP4OPF.solve!
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

const _SCHUR_RELAX = Ref(50)

function mumps_battery_schur(K, E)
    while true
        m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
        MUMPS.suppress_display!(m)
        MUMPS.set_icntl!(m, 8, 0; displaylevel=0)    # no scaling with a Schur complement
        MUMPS.set_icntl!(m, 14, _SCHUR_RELAX[]; displaylevel=0)
        MUMPS.associate_matrix!(m, K)
        MUMPS.mumps_schur_complement!(m, collect(E))
        status = m.infog[1]
        if status >= 0
            S = MUMPS.get_schur_complement(m)
            MUMPS.finalize!(m)
            return S
        end
        MUMPS.finalize!(m)
        status in (-8, -9) || error("MUMPS Schur complement failed, INFOG(1) = $status")
        _SCHUR_RELAX[] *= 2
        @printf("MUMPS workspace retry: ICNTL(14) = %d\n", _SCHUR_RELAX[])
    end
end

DDP4OPF.BATTERY_SCHUR_HOOK[] = mumps_battery_schur
println("BATTERY_SCHUR hook: MUMPS")
