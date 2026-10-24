# Independent, solver-free check of an exported relaxation certificate (2026-10-22).
#   julia problems/msopf/root_level/verify_certificate.jl <certificate.txt> [<more files>]
# Needs only Base and LinearAlgebra (no JuMP, no solver). Reads the file written by export_certificate
# (feasibility_certificate.jl) and recomputes the Lagrangian lower bound
#     min over the box of  c'x - sum_i lam_i (a_i'x + a0_i - b_i) - sum_k <Z_k, M_k(x)>
# in EXACT rational arithmetic (every Float64 in the file is a dyadic rational, read exactly):
#   * scalar multipliers projected to their sign (GE: >= 0, LE: <= 0, EQ: free), as in dual_bound;
#   * PSD multipliers: the raw dual triangle Z is used if an exact LDL' factorization shows Z PSD; otherwise the
#     Float64 eigenvalue clip of Z (as dual_bound does) is shifted by tau*I, tau increased until the exact LDL'
#     of the shifted matrix succeeds (all pivots > 0) -- an exact proof that the multiplier used is PSD;
#   * second-order-cone multipliers projected onto the cone as dual_bound does, then their head raised (if needed)
#     until exact arithmetic confirms the multiplier lies in the cone;
#   * the Frobenius pairing (off-diagonal entries count twice), the residual coefficients and the box minimum
#     are all exact; nothing is rounded until the final result is printed.
# The result is an exact lower bound on min c'x over {constraints, box} GIVEN the data in the file: the model
# coefficients and boxes as written. Whether those boxes and the model's derived bounds are themselves valid is
# a separate question (the derivations), and the effect of widening every box by a relative 1e-12 (to cover
# Float64 rounding in the box construction) is reported too.
using LinearAlgebra, Printf

const Rat = Rational{BigInt}
ratof(s) = (x = parse(Float64, s); isfinite(x) ? Rat(x) : x)

"Exact LDL' PSD test of a symmetric rational matrix: true if PSD, false otherwise."
function exact_psd(A::Matrix{Rat})
    n = size(A, 1); A = copy(A)
    for k in 1:n
        d = A[k, k]
        if d < 0
            return false
        elseif d == 0
            all(A[i, k] == 0 for i in k+1:n) || return false
            continue
        end
        for i in k+1:n
            l = A[i, k] / d
            l == 0 && continue
            for j in k+1:i
                A[i, j] -= l * A[j, k]
            end
        end
        for i in k+1:n, j in k+1:i; A[j, i] = A[i, j]; end
    end
    return true
end

function verify_file(path; widen = 1e-12)
    t0 = time()
    r = Dict{Int,Rat}(); c0 = Rat(0)
    lo = Dict{Int,Any}(); hi = Dict{Int,Any}()
    nlin, npsd, nshift, nraw, maxtau = 0, 0, 0, 0, 0.0; nsoc, nsocfix = 0, 0
    sign_viol = 0.0; psd_min_eig = Inf
    lines = readlines(path)
    i = 1
    while i <= length(lines)
        ln = lines[i]; i += 1
        (isempty(ln) || ln[1] == '#') && continue
        f = split(ln)
        if f[1] == "OBJ"
            c0 += ratof(f[2]); nt = parse(Int, f[3])
            for k in 1:nt; v = parse(Int, f[2+2k]); r[v] = get(r, v, Rat(0)) + ratof(f[3+2k]); end
        elseif f[1] == "VAR"
            v = parse(Int, f[2]); lo[v] = ratof(f[3]); hi[v] = ratof(f[4])
        elseif f[1] == "LIN"
            nlin += 1
            sense = f[2]; b = ratof(f[3]); lam = parse(Float64, f[4]); a0 = ratof(f[5]); nt = parse(Int, f[6])
            lp = sense == "GE" ? max(lam, 0.0) : sense == "LE" ? min(lam, 0.0) : lam
            sign_viol = max(sign_viol, abs(lam - lp))
            lq = Rat(lp)
            lq == 0 && continue
            for k in 1:nt; v = parse(Int, f[5+2k]); r[v] = get(r, v, Rat(0)) - lq * ratof(f[6+2k]); end
            c0 -= lq * (a0 - b)
        elseif f[1] == "SOC"
            nsoc += 1
            y = parse.(Float64, f[3:end]); n = length(y)
            # project as dual_bound does, then make the head exactly large enough: y1^2 >= sum yi^2, y1 >= 0, in exact arithmetic
            tt, uu = y[1], y[2:end]; nu = norm(uu)
            yp = nu <= tt ? copy(y) : nu <= -tt ? zero(y) : vcat((tt + nu) / 2, (tt + nu) / 2 .* uu ./ nu)
            yq = Rat.(yp); tail = sum(yq[2:end] .^ 2)
            if !(yq[1] >= 0 && yq[1]^2 >= tail)
                h = yp[1]; k = 0
                while !(Rat(h) >= 0 && Rat(h)^2 >= tail)
                    h = nextfloat(max(h, sqrt(Float64(tail))) * (1 + 1e-15)); k += 1
                end
                yq[1] = Rat(h); nsocfix += 1
            end
            for _ in 1:n
                g = split(lines[i]); i += 1
                q = parse(Int, g[2]); w = yq[q]
                w == 0 && continue
                c0 -= w * ratof(g[3]); nt = parse(Int, g[4])
                for k in 1:nt; v = parse(Int, g[3+2k]); r[v] = get(r, v, Rat(0)) - w * ratof(g[4+2k]); end
            end
        elseif f[1] == "PSD"
            npsd += 1
            n = parse(Int, f[2]); tri = parse.(Float64, f[3:end]); ntri = length(tri)
            L = zeros(n, n); idx = Dict{Int,Tuple{Int,Int}}()
            for j in 1:n, ii in 1:j; q = (j * (j - 1)) ÷ 2 + ii; L[ii, j] = L[j, ii] = tri[q]; idx[q] = (ii, j); end
            Zq = Rat.(L)
            if exact_psd(Zq)
                nraw += 1
            else
                ev = eigen(Symmetric(L)); psd_min_eig = min(psd_min_eig, minimum(ev.values))
                Lp = ev.vectors * Diagonal(max.(ev.values, 0.0)) * ev.vectors'; Lp = (Lp + Lp') / 2
                Zq = Rat.(Lp); tau = 0.0; nrm = max(opnorm(Lp), 1e-300)
                while !exact_psd(Zq + Rat(tau) * I)
                    tau = tau == 0 ? 1e-15 * nrm : 10tau
                end
                Zq = Zq + Rat(tau) * I; nshift += 1; maxtau = max(maxtau, tau / nrm)
            end
            for _ in 1:ntri
                g = split(lines[i]); i += 1
                q = parse(Int, g[2]); (ii, j) = idx[q]; w = ii == j ? Zq[ii, j] : 2 * Zq[ii, j]
                w == 0 && continue
                c0 -= w * ratof(g[3]); nt = parse(Int, g[4])
                for k in 1:nt; v = parse(Int, g[3+2k]); r[v] = get(r, v, Rat(0)) - w * ratof(g[4+2k]); end
            end
        end
    end
    # box minimum, exact; and with every box widened by a relative `widen` (plus 1e-300 absolute)
    lb = c0; lbw = c0; boxterm = Rat(0); unb = ""
    for (v, a) in r
        a == 0 && continue
        l, u = lo[v], hi[v]
        x = a > 0 ? l : u
        if !(x isa Rat)
            unb = string(v); lb = -Inf; break
        end
        lb += a * x; boxterm += abs(a * x)
        xw = x + (a > 0 ? -1 : 1) * (abs(x) * Rat(widen) + Rat(1e-300))
        lbw += a * xw
    end
    return (; bound = lb, bound_widened = lbw, boxterm, nlin, npsd, nraw, nshift, maxtau, sign_viol, psd_min_eig, unbounded = unb, nsoc, nsocfix, seconds = time() - t0)
end

function report(path)
    v = verify_file(path)
    println(@sprintf("%s\n  exact lower bound %s pu (= %.10f; x1000 for kW)\n  boxes widened by 1e-12 relative: %.10f\n  box term %.4e; %d scalar rows, %d PSD blocks (%d raw duals PSD as given, %d clipped and shifted, largest shift %.1e of the block norm); largest sign violation %.2e; most negative raw PSD-dual eigenvalue %.2e; %.1f s",
                     path, isfinite(v.bound) ? string(Float64(BigFloat(v.bound))) : "-Inf", isfinite(v.bound) ? Float64(BigFloat(v.bound)) : -Inf,
                     isfinite(v.bound) ? Float64(BigFloat(v.bound_widened)) : -Inf, Float64(v.boxterm), v.nlin, v.npsd, v.nraw, v.nshift, v.maxtau, v.sign_viol,
                     v.psd_min_eig, v.seconds))
    return v
end

if abspath(PROGRAM_FILE) == (@__FILE__); foreach(report, ARGS); end
