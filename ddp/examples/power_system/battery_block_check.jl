# Does the value update need only the battery block of the stage KKT inverse?
#
# With the factor-backed policy, the backward pass uses the n_x feedback
# columns of K \ rhs only in the rows E = {battery powers, energy constraints}
# (the only nonzero rows of those right-hand-side columns), and uses
# beta' Qu + omega' c, which by symmetry of K equals R' [alpha; psi].  If so,
# the feedback rows needed are X_E = S^{-1} R_E with S the Schur complement of
# K onto E, and no n_x-column solve is required.
#
# On captured stage systems (FILTERDDP_PERIODIC_CAPTURE_DIR files) this checks
#   1. E has 2 n_x entries;
#   2. X_E from MUMPS's Schur complement against the full solution;
#   3. the symmetry identity used for V_x;
# and times, informally, UMFPACK factor + full solve against MUMPS factor
# with the Schur complement + dense solve.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/battery_block_check.jl <capture files...>

using LinearAlgebra, Printf, Serialization, SparseArrays
using MUMPS
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

# MUMPS reports a too-small workspace estimate as INFO(1) = -9 (or -8) and
# returns garbage; retry with a larger relaxation ICNTL(14) until it succeeds.
function schur_onto(K, E; relax=50)
    while true
        m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.default_icntl, MUMPS.default_cntl64)
        MUMPS.suppress_display!(m)
        MUMPS.set_icntl!(m, 8, 0; displaylevel=0)   # no scaling with a Schur complement
        MUMPS.set_icntl!(m, 14, relax; displaylevel=0)
        MUMPS.associate_matrix!(m, K)
        MUMPS.mumps_schur_complement!(m, E)
        status = m.infog[1]
        if status >= 0
            S = MUMPS.get_schur_complement(m)
            MUMPS.finalize!(m)
            return S, relax
        end
        MUMPS.finalize!(m)
        status in (-8, -9) || error("MUMPS failed with INFOG(1) = $status")
        relax *= 2
    end
end

for file in ARGS
    d = deserialize(file)
    K, R, X = d.K, d.rhs_multi, d.kkt_solution
    nx = d.nx
    E = findall(i -> any(!iszero, @view R[i, 2:end]), 1:size(R, 1))
    t_lu = @elapsed F = lu(K)
    t_full = @elapsed Xfull = F \ R
    schur_onto(K, E)                                  # warm-up / compile
    t_schur = @elapsed ((S, relax) = schur_onto(K, E))
    t_dense = @elapsed XE = lu(S) \ R[E, 2:end]
    ref = @view X[E, 2:end]
    err_block = norm(XE - ref) / norm(ref)
    lhs = X[:, 2:end]' * R[:, 1]
    rhs = R[:, 2:end]' * X[:, 1]
    err_vx = norm(lhs - rhs) / max(norm(lhs), eps())
    asym = norm(S - S') / norm(S)
    @printf("BATTERY_BLOCK file=%s relax=%d n=%d nx=%d nE=%d share=%.4f err_block=%.3e err_vx_identity=%.3e schur_asym=%.1e umf_lu_s=%.3f umf_full_solve_s=%.3f mumps_schur_s=%.3f dense_s=%.3f\n",
            basename(dirname(file)) * "/" * basename(file), relax, size(K, 1), nx, length(E),
            length(E) / size(K, 1), err_block, err_vx, asym, t_lu, t_full, t_schur, t_dense)
    flush(stdout)
end
MUMPS.MPI.Finalized() || MUMPS.MPI.Finalize()
