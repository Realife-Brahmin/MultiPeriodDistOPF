# current_bounds.jl -- source-current bounds derived from the original model (2026-10-22, extended for batteries).
#
# Derivation: results/multiperiod/source_current_bounds.md (idle batteries) and daily_pilot.md (dispatched
# batteries). For a fixed angle vector and one period t:
#   * every terminal voltage is nonzero (lemma), so power-form KCL is current-form KCL with device currents
#     J_b = conj(S_b / V_b), and Y V = J + sum_s (E_s / z_s) e_t(s) with Y nonsingular (Re Y positive definite);
#   * I_s = a_s + sum_b c_sb J_b,   a_s = (E_s - V0_t(s)) / z_s,   c_sb = -Z_t(s),b / z_s,   V = V0 + Z J;
#   * the device power S_b of bus b lies in the rectangle [Plo, Phi] x [Qlo, Qhi] (loads fixed, PV real power
#     fixed, PV reactive power within its rating, battery power within its rating); if the rectangle contains
#     the origin, arg S_b is unrestricted, otherwise it lies in the arc the rectangle subtends (< 180 deg);
#   * |V_b - V0_b| <= rho_b, so |V_b| >= max(lo, |V0_b| - rho_b) and arg V_b is within asin(rho_b / |V0_b|) of
#     arg V0_b; hence J_b lies in a circular sector (arg J = arg V - arg S, |J| <= Smax / Vmin), and rho_b is
#     recomputed from those sectors (support functions) until it stops decreasing;
#   * |I_s| <= max over a direction grid of [Re(conj(u) a_s) + sum_b h_b(u)] / cos(pi / ndir): every vector
#     has a grid direction within pi/ndir of its own, so this bounds |I_s| for every allowed state.
# Generic in the number type: Float64 (with an outward allowance `safe`) or BigFloat (an independent
# high-precision recomputation of the same bound).

wrapang(x) = x - 2 * oftype(x, pi) * round(x / (2 * oftype(x, pi)))
asin_or_pi(r, v) = r < v ? asin(r / v) : oftype(r, pi)

"Add the support function of the hull of the sector {r e^{j phi}: 0 <= r <= R, |phi - phic| <= w} over the grid."
function sector_support!(acc, psis, phic, w, R)
    T = eltype(acc); pT = T(pi)
    w >= pT && (acc .+= R; return acc)
    @inbounds for k in eachindex(psis)
        d = abs(wrapang(psis[k] - phic))
        acc[k] += d <= w ? R : (d - w < pT / 2 ? R * cos(d - w) : zero(T))
    end
    return acc
end

"Value of the support function of one sector hull at angular distance d from its centre."
sector_value(d, w, R) = (pT = oftype(d, pi); w >= pT ? R : d <= w ? R : (d - w < pT / 2 ? R * cos(d - w) : zero(d)))

"Add the support function of the union of two opposite sectors (centres phic and phic + pi, half-width w): the max of the two."
function double_support!(acc, psis, phic, w, R)
    pT = eltype(acc)(pi)
    @inbounds for k in eachindex(psis)
        d1 = abs(wrapang(psis[k] - phic)); d2 = abs(wrapang(psis[k] - phic - pT))
        acc[k] += max(sector_value(d1, w, R), sector_value(d2, w, R))
    end
    return acc
end

"""
Corner enclosure. For a fixed V, J = conj(S) V / |V|^2 maps the power rectangle to a rotated, scaled rectangle whose
support function in direction psi is max_i |c_i| cos(psi - arg V + arg c_i) / |V| over its corners c_i. Over
arg V within w of phi0 and |V| >= Vmin, each corner term is at most sector_value(|psi - phi0 + arg c_i|, w, |c_i| / Vmin)
(the floor at zero only over-estimates). Adds the max over corners, times `scale` (|coefficient| / Vmin).
`cs` holds (|c_i|, arg c_i).
"""
function corner_support!(acc, psis, phi0, w, scale, cs)
    T = eltype(acc)
    @inbounds for k in eachindex(psis)
        m = zero(T)
        for (R, a) in cs
            m = max(m, sector_value(abs(wrapang(psis[k] - phi0 + a)), w, scale * R))
        end
        acc[k] += m
    end
    return acc
end

"""
    device_rectangles(net, kw, t; soc_tighten = false)

Device power set of every bus with devices at period t, as (bus, Plo, Phi, Qlo, Qhi) in pu: loads at the period's
load factor (fixed), PV real power at the period's PV factor (fixed), PV reactive power in [-q, q] with
q = sqrt(S_rated^2 - P^2) when it is dispatched, battery power in [-P_rated, P_rated] when batteries are dispatched
(zero when idle). The rectangle encloses every allowed device power; nothing is conditioned on a schedule.
"""
function device_rectangles(net, kw, t)
    S = net.S_base_kVA; lam = kw.load_factor[t]; pvf = get(kw, :pv_factor, zeros(kw.T))[t]
    battery = get(kw, :battery, :optimize); pv_q = get(kw, :pv_q, :optimize)
    P = Dict(b => 0.0 for b in net.buses); Q0 = Dict(b => 0.0 for b in net.buses); qd = Dict(b => 0.0 for b in net.buses); pb = Dict(b => 0.0 for b in net.buses)
    for d in net.loads; P[d.bus] -= lam * d.kW / S; Q0[d.bus] -= lam * d.kvar / S; end
    for p in get(kw, :pvs, [])
        P[p.bus] += pvf * p.kW / S
        pv_q === :optimize && (qd[p.bus] += sqrt(max(p.kVA^2 - (pvf * p.kW)^2, 0.0)) / S)
    end
    battery === :optimize && for bt in net.batteries; pb[bt.bus] += bt.kW_rated / S; end
    src = Set(s.bus for s in net.sources)
    out = [(b, P[b] - pb[b], P[b] + pb[b], Q0[b] - qd[b], Q0[b] + qd[b]) for b in net.buses if P[b] != 0 || Q0[b] != 0 || qd[b] != 0 || pb[b] != 0]
    any(r[1] in src for r in out) && error("device at a source terminal: the band (Vmin) argument does not cover it")
    return out
end

"Arc of arg S over the rectangle: (centre, half-width); half-width pi (unrestricted) when it contains the origin."
function rect_arc(::Type{T}, Plo, Phi, Qlo, Qhi) where {T}
    (Plo <= 0 <= Phi && Qlo <= 0 <= Qhi) && return (zero(T), T(pi))
    mid = complex(T(Plo + Phi) / 2, T(Qlo + Qhi) / 2); phc = angle(mid)
    rel = [wrapang(angle(complex(T(p), T(q))) - phc) for p in (Plo, Phi), q in (Qlo, Qhi)]
    lo_, hi_ = minimum(rel), maximum(rel)
    return (phc + (lo_ + hi_) / 2, (hi_ - lo_) / 2)
end

"""
    current_bounds(RM, ang, rects, t; T = Float64, ndir = 3600, ndir_rho = ndir, lo = 0.95, safe = 1 + 1e-9)

Bounds on |I_s| (pu, type T) for every source at angle vector `ang`, period t, device rectangles `rects`, using
the network response `RM` (Z in BigFloat). Returns the sector bound, the termwise bound, the no-load term and
the disc radii history; also the sectors (for diagnostics).
"""
# segments = true: a device set that is a segment through zero power along one axis (a battery bus without PV or
# reactive load: S = P_B real, |P_B| <= rating) has arg S in {0, pi} (or {-pi/2, pi/2}), so J = conj(S) V / |V|^2
# lies in two opposite thin sectors around arg V, not in a full disc. segments = false: the full disc (as before).
function current_bounds(RM, ang, rects, t; T = Float64, ndir = 3600, ndir_rho = ndir, lo = 0.95, safe = T == Float64 ? 1 + 1e-9 : one(T), iters = 8, segments = false, corners = false)
    id = RM.id; Yb = big(net.Y_base); n = length(net.buses)
    b0 = zeros(Complex{BigFloat}, n); Ed = Dict{String,Complex{BigFloat}}(); ysd = Dict{String,Complex{BigFloat}}()
    for s in net.sources
        E = big(s.pu) * cis(deg2rad(big(angle_at(ang, s.name, t)))); y = 1 / (big(s.z_series) * Yb)
        b0[id[s.bus]] += y * E; Ed[s.name] = E; ysd[s.name] = y
    end
    V0b = RM.Z * b0
    CT = Complex{T}
    V0 = CT.(V0b); Zt = CT.(RM.Z)
    D = [r[1] for r in rects]; k_ = [id[b] for b in D]; nd = length(D)
    Smax = [maximum(abs(complex(T(p), T(q))) for p in (r[2], r[3]), q in (r[4], r[5])) for r in rects]
    sarc = [rect_arc(T, r[2], r[3], r[4], r[5]) for r in rects]
    # double sectors: real segment through zero (centre arg V) or imaginary segment through zero (centre arg V - pi/2)
    dbl = [segments && sarc[q][2] >= T(pi) && ((r[4] == 0 && r[5] == 0) || (r[2] == 0 && r[3] == 0)) for (q, r) in enumerate(rects)]
    dshift = [rects[q][4] == 0 && rects[q][5] == 0 ? zero(T) : T(pi) / 2 for q in eachindex(rects)]
    # corners = true: every device uses the corner enclosure (corner_support!), which needs no arc of arg S
    cor = [[(abs(complex(T(p), T(qq))), angle(complex(T(p), T(qq)))) for p in (r[2], r[3]) for qq in (r[4], r[5])] for r in rects]
    vmn = fill(T(lo), length(rects))
    psis(N) = [2 * T(pi) * (k - 1) / N for k in 1:N]
    ps, pr = psis(ndir), psis(ndir_rho)
    cw, cwr = cos(T(pi) / ndir), cos(T(pi) / ndir_rho)
    Jmax0 = Smax ./ T(lo)
    rho = [sum(abs(Zt[id[b], k_[q]]) * Jmax0[q] for q in 1:nd) for b in D]
    hist = T[maximum(rho)]
    phiJ = zeros(T, nd); wJ = fill(T(pi), nd); RJ = copy(Jmax0)
    acc = zeros(T, ndir_rho)
    for it in 1:iters
        for q in 1:nd
            v0 = abs(V0[k_[q]]); al = rho[q] < v0 ? asin(rho[q] / v0) : T(pi)
            vmn[q] = max(T(lo), v0 - rho[q]); RJ[q] = Smax[q] / vmn[q]
            if dbl[q]
                phiJ[q] = angle(V0[k_[q]]) - dshift[q]; wJ[q] = min(al, T(pi))
            else
                phiJ[q] = angle(V0[k_[q]]) - sarc[q][1]; wJ[q] = min(al + sarc[q][2], T(pi))
            end
        end
        new = similar(rho)
        for (p, b) in enumerate(D)
            fill!(acc, zero(T))
            for q in 1:nd
                zz = Zt[id[b], k_[q]]
                if corners
                    corner_support!(acc, pr, angle(V0[k_[q]]) + angle(zz), min(asin_or_pi(rho[q], abs(V0[k_[q]])), T(pi)), abs(zz) / vmn[q], cor[q])
                else
                    (dbl[q] ? double_support! : sector_support!)(acc, pr, phiJ[q] + angle(zz), wJ[q], abs(zz) * RJ[q])
                end
            end
            new[p] = maximum(acc) / cwr * safe
        end
        rho = min.(rho, new); push!(hist, maximum(rho))
        hist[end] > hist[end-1] * (1 - T(1e-4)) && break
    end
    out = Dict{String,Any}()
    for s in net.sources
        tt = id[s.bus]; y = ysd[s.name]
        a = CT(y * (Ed[s.name] - V0b[tt])); cs = CT.(-y .* RM.Z[tt, k_])
        term = (abs(a) + sum(abs.(cs) .* Jmax0)) * safe
        acc2 = abs(a) .* cos.(ps .- angle(a))
        for q in 1:nd
            if corners
                corner_support!(acc2, ps, angle(V0[k_[q]]) + angle(cs[q]), min(asin_or_pi(rho[q], abs(V0[k_[q]])), T(pi)), abs(cs[q]) / vmn[q], cor[q])
            else
                (dbl[q] ? double_support! : sector_support!)(acc2, ps, phiJ[q] + angle(cs[q]), wJ[q], abs(cs[q]) * RJ[q])
            end
        end
        sect = max(maximum(acc2), zero(T)) / cw * safe
        out[s.name] = (; sect, term, a0 = abs(a), a, cs)
    end
    full = count(q -> wJ[q] >= T(pi) && !dbl[q], 1:nd)
    # rho: final valid radii |V_b - V0_b| <= rho_b; vmn: valid lower bounds on |V_b|; V0D: no-load voltages; Smax: max |S_b|
    return (; out, rho_hist = hist, nd, full_sectors = full, double_sectors = count(dbl), wJ, RJ, D, rho, vmn, V0D = V0[k_], Smax)
end
