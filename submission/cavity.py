"""Integer-certified cavity lower bound for the four disjoint local saving events.

Uses:
  - A graded grid on [0, 2] with 1/(8n) refinement near zero and exact integer spacings.
  - Exact adaptive 0/1/2-bisection + weighted-harmonic geometric interpolation.
  - Precomputed integer kernel tables and O(1) cell lookup (zero floating-point operations).
  - Exact u-scaled Holder-compressed 4-corner geometric-mean tables for E (all disjoint
    (OUT x OUT), (V x OUT), (OUT x V) endpoint pairs and four guarded two-step hops).
  - Direct coordinatewise log-concave integration for Kind 1 in (x, b, t) and Kind 2 in (x, y, z).
"""
import operator
from bisect import bisect_right
from functools import lru_cache
from math import isqrt

SCALE = 100_000_000       # 1e8
UNIT = 100_000_000        # 1e8 represents 1/n
ROOT_LO = 141_421_356     # < sqrt(2) * 1e8
ROOT_HI = 141_421_357     # > sqrt(2) * 1e8
INV_R2_LO = 70_710_678    # < (1/sqrt(2)) * 1e8
RHO_HI = 94_280_905       # > (2*sqrt(2)/3) * 1e8
CG_NUM = 276_395_319      # < 2*sqrt(6/pi) * 1e8
C2_NUM = 47_746_482       # < 3/(2*pi) * 1e8
C3_NUM = 65_984_522       # < sqrt(108)/(2*pi)^(3/2) * 1e8
HALF_CDF_MESH = 2048
EXP_DEN = 1_000_000
EXP_CUT = 14_000_000
WEIGHT_CUT = 1_000
_mul = operator.mul


def ceil_div(a, b):
    return -(-a // b)


def geo_interp(a, b, num, den, _isqrt=isqrt):
    """Rigorous integer lower bound on a^(1-t) * b^t for t = num/den in [0, 1]."""
    if a == b or num <= 0:
        return a
    if num >= den:
        return b
    if not a or not b:
        return 0
    lo, hi = a, b
    mn = a if a < b else b
    diff = (b - a) if b > a else (a - b)
    if diff * 48 > mn:
        num2 = num << 1
        if num2 >= den:
            lo = _isqrt(lo * hi)
            num = num2 - den
        else:
            hi = _isqrt(lo * hi)
            num = num2
        if diff * 8 > mn:
            num2 = num << 1
            if num2 >= den:
                lo = _isqrt(lo * hi)
                num = num2 - den
            else:
                hi = _isqrt(lo * hi)
                num = num2
        if not num:
            return lo
        if not lo or not hi:
            return 0
    return (den * lo * hi) // ((den - num) * hi + num * lo)


@lru_cache(maxsize=None)
def exp_floor(i):
    """Rigorous integer lower bound on SCALE * exp(-i / EXP_DEN) for i >= 0."""
    if i <= 0:
        return SCALE
    if i > EXP_CUT:
        return 0
    h = 256 * EXP_DEN
    u = -(-(i * SCALE) // h)
    l = (i * SCALE) // h
    v = SCALE - u + (l * i) // (2 * h) - (-(-(-(-(u * i) // (2 * h)) * i) // (3 * h)))
    if v <= 0:
        return 0
    for _ in range(8):
        v = (v * v) // SCALE
    return v


def build_graded_grid(n):
    """Build exact graded grid on [0, 2] in units of UNIT (= 1/n), with n divisible by 8."""
    assert n % 8 == 0
    U = UNIT
    pts = [0, U // 8, U // 4, U // 2, U, (3 * U) // 2]
    k1 = (3 * n) // 8
    k2 = (3 * n) // 4
    k3 = (5 * n) // 4
    for k in range(2, k1 + 1):
        pts.append(k * U)
    for k in range(k1 + 2, k2 + 1, 2):
        pts.append(k * U)
    if pts[-1] < k2 * U:
        pts.append(k2 * U)
    for k in range(k2 + 4, k3 + 1, 4):
        pts.append(k * U)
    if pts[-1] < k3 * U:
        pts.append(k3 * U)
    for k in range(k3 + 8, 2 * n + 1, 8):
        pts.append(k * U)
    if pts[-1] < 2 * n * U:
        pts.append(2 * n * U)
    assert pts[-1] == 2 * n * U
    return pts


def cavity(n=64, nt1=16):
    """Return ((c_pos1, c_neg1), den1, (c_pos2, c_neg2), den2) as exact integer lower bounds."""
    if not __debug__:
        raise RuntimeError("Run the certificate with assertions enabled.")
    grid = build_graded_grid(n)
    side = len(grid)
    cside = side - 1
    size = side * side
    csize = cside * cside
    UNIT_N = UNIT * n
    bound = 2 * UNIT_N
    dg = [grid[k + 1] - grid[k] for k in range(cside)]
    mesh = HALF_CDF_MESH
    isqrt_ = isqrt
    geo_ = geo_interp
    U0125 = UNIT // 8

    V_QUANT = 40_000
    v_div = UNIT_N // V_QUANT
    assert UNIT_N % V_QUANT == 0
    assert all(g % 2 == 0 and g % U0125 == 0 for g in grid)

    # Exact O(1) lookup table for grid cell index, verified against bisect_right.
    n_lut = 16 * n
    idx_lut = [0] * (n_lut + 1)
    cur_j = 0
    for q_idx in range(n_lut + 1):
        y_q = q_idx * U0125
        while cur_j + 1 < cside and grid[cur_j + 1] <= y_q:
            cur_j += 1
        idx_lut[q_idx] = cur_j
    assert all(idx_lut[q] == min(cside - 1, bisect_right(grid, q * U0125) - 1) for q in range(n_lut + 1))

    @lru_cache(maxsize=None)
    def density(v):
        av = -v if v < 0 else v
        t = -(-av // v_div)
        q = -(-(3 * t * t) // 800)
        return (CG_NUM * exp_floor(q)) // UNIT

    mass = [0]
    m2 = mesh * mesh
    for i in range(1, 2 * mesh + 1):
        q = -(-(2 * EXP_DEN * (3 * i * i - 3 * i + 1)) // m2)
        mass.append(mass[-1] + (CG_NUM * exp_floor(q)) // (UNIT * mesh))
    max_mass_idx = 2 * mesh

    def half_mass(r):
        if r <= 0:
            return 0
        idx = (mesh * r) // UNIT_N
        return mass[max_mass_idx if idx > max_mass_idx else idx]

    def find_idx(y):
        if y < 0:
            return 0
        if y >= bound:
            return cside - 1
        return idx_lut[y // U0125]

    def locate(y_lo, y_hi):
        if y_lo < 0 or y_hi >= bound:
            return None
        j0 = idx_lut[y_lo // U0125]
        g0 = grid[j0]
        r0 = y_lo - g0
        d0 = dg[j0]
        if y_hi == y_lo:
            return (j0, r0, d0, -1, 0, 0)
        if y_hi < grid[j0 + 1]:
            return (j0, r0, d0, j0, y_hi - g0, d0)
        j1 = idx_lut[y_hi // U0125]
        return (j0, r0, d0, j1, y_hi - grid[j1], dg[j1])

    def eval_loc(mat, row_off, loc):
        j0, r0, d0, j1, r1, d1 = loc
        k0 = row_off + j0
        v0 = mat[k0]
        if not v0:
            return 0
        if j1 == j0:
            v0b = mat[k0 + 1]
            if not v0b:
                return 0
            if v0b >= v0:
                return geo_(v0, v0b, r0, d0) if r0 else v0
            return geo_(v0, v0b, r1, d0) if r1 else v0
        if r0:
            v0b = mat[k0 + 1]
            if not v0b:
                return 0
            a = geo_(v0, v0b, r0, d0)
            if not a:
                return 0
        else:
            a = v0
        if j1 < 0:
            return a
        k1 = row_off + j1
        v1 = mat[k1]
        if not v1:
            return 0
        if r1:
            v1b = mat[k1 + 1]
            if not v1b:
                return 0
            b = geo_(v1, v1b, r1, d1)
        else:
            b = v1
        return a if a < b else b

    def eval_batch(mats, row_off, loc):
        j0, r0, d0, j1, r1, d1 = loc
        k0 = row_off + j0
        if j1 == j0:
            k0p = k0 + 1
            out = []
            for mat in mats:
                v0 = mat[k0]
                if not v0:
                    out.append(0)
                    continue
                v0b = mat[k0p]
                if not v0b:
                    out.append(0)
                elif v0b >= v0:
                    out.append(geo_(v0, v0b, r0, d0) if r0 else v0)
                else:
                    out.append(geo_(v0, v0b, r1, d0) if r1 else v0)
            return tuple(out)
        if j1 < 0:
            if not r0:
                return tuple(mat[k0] for mat in mats)
            k0p = k0 + 1
            return tuple(geo_(mat[k0], mat[k0p], r0, d0) if (mat[k0] and mat[k0p]) else 0 for mat in mats)
        k0p = k0 + 1
        k1 = row_off + j1
        k1p = k1 + 1
        out = []
        for mat in mats:
            v0 = mat[k0]
            if not v0:
                out.append(0)
                continue
            if r0:
                v0b = mat[k0p]
                if not v0b:
                    out.append(0)
                    continue
                a = geo_(v0, v0b, r0, d0)
                if not a:
                    out.append(0)
                    continue
            else:
                a = v0
            v1 = mat[k1]
            if not v1:
                out.append(0)
                continue
            if r1:
                v1b = mat[k1p]
                if not v1b:
                    out.append(0)
                    continue
                b = geo_(v1, v1b, r1, d1)
            else:
                b = v1
            out.append(a if a < b else b)
        return tuple(out)

    ix_max1 = find_idx((7 * UNIT_N) // 10) + 1
    jb_max1 = find_idx((6 * UNIT_N) // 5) + 1
    ix_max2 = find_idx(UNIT_N // 2) + 1
    jy_max2 = find_idx(UNIT_N) + 1
    points = range(side)
    cpoints = range(cside)

    # 1. Precompute positive Gaussian transition kernels K and U-guarded kernels KU.
    kernels = []
    K_steps = []
    KU_steps = []
    for i in points:
        gi = grid[i]
        gi_half = gi // 2
        r_gi = (ROOT_LO * gi) // UNIT
        for j in points:
            gj = grid[j]
            base = -(-(ROOT_HI * gj) // UNIT) + gi_half
            base_u = (7 * gj) // 2 + r_gi
            krow = []
            ks = []
            kus = []
            ku0 = 0
            for k in points:
                gk = grid[k]
                val = density(base + gk)
                if not val:
                    break
                krow.append(val)
                if k > 0:
                    ks.append((k, val, dg[k - 1]))
                uval = (val * half_mass(base_u + (ROOT_LO * gk) // UNIT)) // SCALE
                if k == 0:
                    ku0 = uval
                else:
                    kus.append((k, uval, dg[k - 1]))
            kernels.append(krow)
            K_steps.append((krow[0], ks) if len(krow) >= 2 else None)
            KU_steps.append((ku0, kus) if len(krow) >= 2 else None)

    def apply_k_steps(mat, steps_tbl, i_lim=side):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
            iside = i * side
            for j in points:
                entry = steps_tbl[iside + j]
                if entry is None:
                    continue
                w0, ks = entry
                row = j * side
                v0 = mat[row]
                prev = (w0 * v0) // SCALE if v0 else 0
                tot = 0
                for k, wk, step_w in ks:
                    vk = mat[row + k]
                    if vk:
                        cur = (wk * vk) // SCALE
                        if prev and cur:
                            tot += isqrt_(prev * cur) * step_w
                        prev = cur
                    else:
                        prev = 0
                out[iside + j] = tot // UNIT_N
        return out

    def build_S_steps(strict, i_lim):
        mult = 2 if strict else 1
        tbl = [None] * (i_lim * side)
        for i in range(i_lim):
            gi_half = grid[i] // 2
            iside = i * side
            for j in points:
                gj = grid[j]
                ml = (ROOT_LO * gj) // UNIT - gi_half
                mh = -(-(ROOT_HI * gj) // UNIT) - gi_half
                c = -mult * ml if ml < 0 else 0
                if c >= bound:
                    continue
                k0 = idx_lut[c // U0125]
                rem = c - grid[k0]
                dg0 = dg[k0]
                w0 = density(c + mh if c + mh > 0 else 0)
                steps = []
                prev_c = c
                for kk in range(k0 + 1, side):
                    edge = grid[kk]
                    w = density(edge + mh)
                    if not w:
                        break
                    steps.append((kk, w, edge - prev_c))
                    prev_c = edge
                if steps:
                    tbl[iside + j] = (k0, rem, dg0, w0, steps)
        return tbl

    def apply_S_steps(mat, S_tbl, i_lim):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
            iside = i * side
            for j in points:
                entry = S_tbl[iside + j]
                if entry is None:
                    continue
                k0, rem, dg0, w0, steps = entry
                row = j * side
                if w0:
                    v0 = geo_(mat[row + k0], mat[row + k0 + 1], rem, dg0) if rem else mat[row + k0]
                    prev = (v0 * w0) // SCALE if v0 else 0
                else:
                    prev = 0
                tot = 0
                for kk, w, step_w in steps:
                    vk = mat[row + kk]
                    if vk:
                        cur = (w * vk) // SCALE
                        if prev and cur:
                            tot += isqrt_(prev * cur) * step_w
                        prev = cur
                    else:
                        prev = 0
                out[iside + j] = tot // UNIT_N
        return out

    # 2. Alternating crossing branches B_0 .. B_6, A, OUT, V_odd, and O_g.
    branches = [[half_mass((ROOT_LO * grid[j]) // UNIT + grid[i] // 2) for i in points for j in points]]
    for _ in range(6):
        branches.append(apply_k_steps(branches[-1], K_steps))

    S0_tbl = build_S_steps(False, ix_max1 + 1)
    S1_tbl_full = build_S_steps(True, side)

    def op_D_batch(f_list, g_list, pairs, i_lim=side):
        n_out = len(pairs)
        outs = [[0] * (i_lim * side) for _ in range(n_out)]
        d0_val = density(0)
        for i in range(i_lim):
            gi_half = grid[i] // 2
            iside = i * side
            for j in points:
                gj = grid[j]
                s0 = gi_half - (-(-(ROOT_HI * gj) // UNIT))
                if s0 <= 0:
                    continue
                s1 = gi_half - (ROOT_LO * gj) // UNIT
                row = j * side
                upper = bound if bound < s0 else s0
                kmax = find_idx(upper)
                prev_h = 0
                loc0 = locate(s0, s1)
                if loc0:
                    fv0 = eval_batch(f_list, row, loc0)
                    gv0 = eval_batch(g_list, row, loc0)
                    prev_vs = [((d0_val * fv0[fi]) // SCALE * gv0[gi_]) // SCALE for fi, gi_ in pairs]
                else:
                    prev_vs = [0] * n_out
                tots = [0] * n_out
                for kk in range(1, kmax + 2):
                    h = grid[kk] if kk <= kmax else upper
                    if h > upper:
                        h = upper
                    if h <= prev_h:
                        continue
                    dh = density(h)
                    if not dh:
                        break
                    lp = locate(s0 + h, s1 + h)
                    lm = locate(s0 - h, s1 - h)
                    if lp and lm:
                        fv = eval_batch(f_list, row, lp)
                        gv = eval_batch(g_list, row, lm)
                        step_w = h - prev_h
                        for p_idx, (fi, gi_) in enumerate(pairs):
                            cur = ((dh * fv[fi]) // SCALE * gv[gi_]) // SCALE
                            pv = prev_vs[p_idx]
                            if pv and cur:
                                tots[p_idx] += isqrt_(pv * cur) * step_w
                            prev_vs[p_idx] = cur
                    else:
                        prev_vs = [0] * n_out
                    prev_h = h
                    if h == upper:
                        break
                idx_ij = iside + j
                for p_idx in range(n_out):
                    outs[p_idx][idx_ij] = tots[p_idx] // UNIT_N
        return outs

    even = [branches[0], branches[2], branches[4], branches[6]]
    odd = [branches[1], branches[3], branches[5]]

    same_even = [apply_S_steps(mat, S0_tbl, ix_max1 + 1) for mat in even]
    A_mats = list(same_even) + op_D_batch(odd[:2], even[:3], [(a, b) for a in range(2) for b in range(3)], ix_max1 + 1)

    W1_1 = apply_S_steps(branches[1], S1_tbl_full, side)
    W1_3 = apply_S_steps(branches[3], S1_tbl_full, side)
    W1_5 = apply_S_steps(branches[5], S1_tbl_full, side)
    W2_11 = op_D_batch((branches[1],), (branches[1],), ((0, 0),), side)[0]

    # Precompute K_gm4_kerns and exact reachability masks for E_gm4 and T_mats.
    K_gm4_kerns = []
    for i in cpoints:
        i0 = i * side
        i1 = i0 + side
        for j in cpoints:
            k00 = kernels[i0 + j]
            k01 = kernels[i0 + j + 1]
            k10 = kernels[i1 + j]
            k11 = kernels[i1 + j + 1]
            nk = min(len(k00), len(k01), len(k10), len(k11))
            row_k = []
            for r in range(nk - 1):
                w = isqrt_(
                    isqrt_(isqrt_(k00[r] * k00[r + 1]) * isqrt_(k01[r] * k01[r + 1]))
                    * isqrt_(isqrt_(k10[r] * k10[r + 1]) * isqrt_(k11[r] * k11[r + 1]))
                )
                if w:
                    row_k.append((r, w * dg[r]))
            K_gm4_kerns.append(row_k)

    need_i = [0] * cside
    for i in range(ix_max1):
        need_i[i] = max(need_i[i], jb_max1)
    for i in range(ix_max2):
        need_i[i] = max(need_i[i], jy_max2)
    for _ in range(8):
        for i in cpoints:
            lim_j = need_i[i]
            ics = i * cside
            for j in range(lim_j):
                kl = K_gm4_kerns[ics + j]
                if kl:
                    rm = kl[-1][0] + 1
                    if rm > need_i[j]:
                        need_i[j] = rm

    t_max_j = [0] * side
    for i in range(ix_max1 + 1):
        iside = i * side
        for j in points:
            en = S1_tbl_full[iside + j]
            if en is not None:
                rm = en[4][-1][0] + 1
                if rm > t_max_j[j]:
                    t_max_j[j] = rm
    for _ in range(6):
        for i in points:
            iside = i * side
            for j in range(t_max_j[i]):
                en = K_steps[iside + j]
                if en is not None:
                    rm = en[1][-1][0] + 1
                    if rm > t_max_j[j]:
                        t_max_j[j] = rm

    v_i_lim = max(max(need_i) + 1, max(t_max_j))
    if v_i_lim > side:
        v_i_lim = side

    VU_steps = [None] * (v_i_lim * side)
    VW1_steps = [None] * (v_i_lim * side)
    VW13_steps = [None] * (v_i_lim * side)
    W_1_3 = (W1_1, W1_3)
    for i in range(v_i_lim):
        gi = grid[i]
        gi_half = gi // 2
        iside = i * side
        r_gi_hi = -(-(ROOT_HI * gi) // UNIT)
        for j in points:
            entry = S1_tbl_full[iside + j]
            if entry is None:
                continue
            k0, rem, dg0, w0, steps = entry
            gj = grid[j]
            ml = (ROOT_LO * gj) // UNIT - gi_half
            mh = -(-(ROOT_HI * gj) // UNIT) - gi_half
            c = -2 * ml if ml < 0 else 0
            row = j * side
            level = (7 * gj) // 2 - r_gi_hi
            u0 = (w0 * half_mass(level + (ROOT_LO * c) // UNIT)) // SCALE if w0 else 0
            us = [(kk, (w * half_mass(level + (ROOT_LO * grid[kk]) // UNIT)) // SCALE, step_w) for kk, w, step_w in steps]
            VU_steps[iside + j] = (k0, rem, dg0, u0, us)

            loc0 = locate(0 if 2 * ml + c < 0 else 2 * ml + c, 2 * mh + c)
            if w0 and loc0:
                g0_1, g0_3 = eval_batch(W_1_3, row, loc0)
                g0_1 = (w0 * g0_1) // SCALE
                g0_3 = (w0 * g0_3) // SCALE
            else:
                g0_1 = g0_3 = 0
            ws1 = []
            ws3 = []
            for kk, w, step_w in steps:
                edge = grid[kk]
                loc = locate(2 * ml + edge, 2 * mh + edge)
                if loc:
                    v1, v3 = eval_batch(W_1_3, row, loc)
                    ws1.append((kk, (w * v1) // SCALE if v1 else 0, step_w))
                    ws3.append((kk, (w * v3) // SCALE if v3 else 0, step_w))
                else:
                    ws1.append((kk, 0, step_w))
                    ws3.append((kk, 0, step_w))
            VW1_steps[iside + j] = (k0, rem, dg0, g0_1, ws1)
            VW13_steps[iside + j] = (k0, rem, dg0, g0_3, ws3)

    V_0U = apply_S_steps(branches[0], VU_steps, v_i_lim)
    V_2U = apply_S_steps(branches[2], VU_steps, v_i_lim)
    V_4U = apply_S_steps(branches[4], VU_steps, v_i_lim)
    V_0W1 = apply_S_steps(branches[0], VW1_steps, v_i_lim)
    V_2W1 = apply_S_steps(branches[2], VW1_steps, v_i_lim)
    V_0W13 = apply_S_steps(branches[0], VW13_steps, v_i_lim)
    del VU_steps, VW1_steps, VW13_steps

    # Guarded far-arm kernels for O_g (needed only up to ix_max2, i.e., x <= 17/32 in Kind 2).
    kw_i_lim = ix_max2 + 1
    KW1_steps = [None] * (kw_i_lim * side)
    KW13_steps = [None] * (kw_i_lim * side)
    for i in range(kw_i_lim):
        gi = grid[i]
        iside = i * side
        for j in points:
            kern = kernels[iside + j]
            nk = len(kern)
            if nk < 2:
                continue
            row = j * side
            gj = grid[j]
            c_lo = (2 * ROOT_LO * gj) // UNIT + gi
            c_hi = -(-(2 * ROOT_HI * gj) // UNIT) + gi
            loc0 = locate(c_lo, c_hi)
            if loc0:
                v1, v3 = eval_batch(W_1_3, row, loc0)
                w0_1 = (kern[0] * v1) // SCALE
                w0_3 = (kern[0] * v3) // SCALE
            else:
                w0_1 = w0_3 = 0
            s1 = []
            s3 = []
            for k in range(1, nk):
                gk = grid[k]
                loc = locate(c_lo + gk, c_hi + gk)
                if loc:
                    v1, v3 = eval_batch(W_1_3, row, loc)
                    s1.append((k, (kern[k] * v1) // SCALE if v1 else 0, dg[k - 1]))
                    s3.append((k, (kern[k] * v3) // SCALE if v3 else 0, dg[k - 1]))
                else:
                    s1.append((k, 0, dg[k - 1]))
                    s3.append((k, 0, dg[k - 1]))
            KW1_steps[iside + j] = (w0_1, s1)
            KW13_steps[iside + j] = (w0_3, s3)

    Og_mats = (
        [apply_k_steps(m_, KU_steps, kw_i_lim) for m_ in even]
        + [apply_k_steps(m_, KW1_steps, kw_i_lim) for m_ in even[:3]]
        + [apply_k_steps(branches[0], KW13_steps, kw_i_lim)]
    )
    del KW1_steps, KW13_steps, S0_tbl

    # 3. Exact Holder-compressed GM4 table E_gm4 for ALL disjoint endpoint pairs and 4 guarded hops.
    Kg_gm4_kerns = [[] for _ in range(csize)]
    for i in cpoints:
        lim_j = need_i[i]
        if not lim_j:
            continue
        a0, a1 = grid[i], grid[i + 1]
        ra0 = (ROOT_LO * a0) // UNIT
        ra1 = (ROOT_LO * a1) // UNIT
        ics = i * cside
        for j in range(lim_j):
            klist = K_gm4_kerns[ics + j]
            if not klist:
                continue
            b0, b1 = grid[j], grid[j + 1]
            u00 = (7 * b0) // 2 + ra0
            u01 = (7 * b1) // 2 + ra0
            u10 = (7 * b0) // 2 + ra1
            u11 = (7 * b1) // 2 + ra1
            c00_l = (2 * ROOT_LO * b0) // UNIT + a0
            c00_h = -(-(2 * ROOT_HI * b0) // UNIT) + a0
            c01_l = (2 * ROOT_LO * b1) // UNIT + a0
            c01_h = -(-(2 * ROOT_HI * b1) // UNIT) + a0
            c10_l = (2 * ROOT_LO * b0) // UNIT + a1
            c10_h = -(-(2 * ROOT_HI * b0) // UNIT) + a1
            c11_l = (2 * ROOT_LO * b1) // UNIT + a1
            c11_h = -(-(2 * ROOT_HI * b1) // UNIT) + a1
            r0 = j * side
            r1 = (j + 1) * side
            rmax = klist[-1][0] + 1
            u_arr = [
                (
                    half_mass(u00 + (ROOT_LO * grid[r]) // UNIT),
                    half_mass(u01 + (ROOT_LO * grid[r]) // UNIT),
                    half_mass(u10 + (ROOT_LO * grid[r]) // UNIT),
                    half_mass(u11 + (ROOT_LO * grid[r]) // UNIT),
                )
                for r in range(rmax + 1)
            ]
            w_arrs = []
            if c00_l < UNIT_N:
                locs_r = [
                    (
                        locate(c00_l + grid[r], c00_h + grid[r]),
                        locate(c01_l + grid[r], c01_h + grid[r]),
                        locate(c10_l + grid[r], c10_h + grid[r]),
                        locate(c11_l + grid[r], c11_h + grid[r]),
                    )
                    for r in range(rmax + 1)
                ]
                for wm in W_1_3:
                    wa = []
                    for l00, l01, l10, l11 in locs_r:
                        if l00 and l01 and l10 and l11:
                            v00 = eval_loc(wm, r0, l00)
                            if v00:
                                v01 = eval_loc(wm, r1, l01)
                                v10 = eval_loc(wm, r0, l10) if v01 else 0
                                v11 = eval_loc(wm, r1, l11) if v10 else 0
                                wa.append((v00, v01, v10, v11) if v11 else None)
                            else:
                                wa.append(None)
                        else:
                            wa.append(None)
                    w_arrs.append(wa)
            row_kg = []
            for r, w_dg in klist:
                ua, ub = u_arr[r], u_arr[r + 1]
                gtot = isqrt_(
                    isqrt_(isqrt_(ua[0] * ub[0]) * isqrt_(ua[1] * ub[1]))
                    * isqrt_(isqrt_(ua[2] * ub[2]) * isqrt_(ua[3] * ub[3]))
                )
                for wa in w_arrs:
                    wa0, wa1 = wa[r], wa[r + 1]
                    if wa0 and wa1:
                        gtot += isqrt_(
                            isqrt_(isqrt_(wa0[0] * wa1[0]) * isqrt_(wa0[1] * wa1[1]))
                            * isqrt_(isqrt_(wa0[2] * wa1[2]) * isqrt_(wa0[3] * wa1[3]))
                        )
                if gtot:
                    row_kg.append((r, (w_dg * gtot) // SCALE))
            Kg_gm4_kerns[ics + j] = row_kg

    def apply_gm4_kern(in_gm4, kerns_list):
        out = [0] * csize
        denom_k = SCALE * UNIT_N
        for i in cpoints:
            lim_j = need_i[i]
            if not lim_j:
                continue
            ics = i * cside
            for j in range(lim_j):
                klist = kerns_list[ics + j]
                if not klist:
                    continue
                jcs = j * cside
                tot = 0
                for r, w_dg in klist:
                    v = in_gm4[jcs + r]
                    if v:
                        tot += w_dg * v
                out[ics + j] = tot // denom_k
        return out

    # u-scaled T_gm4 covering ALL (OUT x OUT), (V x OUT), (OUT x V) pairs at once.
    h_cap_u = (13 * UNIT_N) // 10
    W_t = (W1_1, W1_3, W1_5, W2_11)
    V_t = (V_0U, V_2U, V_4U, V_0W1, V_2W1, V_0W13)
    WV_t = W_t + V_t
    nWt = len(W_t)
    nWV = len(WV_t)
    T_gm4 = [0] * csize
    for i in cpoints:
        lim_j = need_i[i]
        if not lim_j:
            continue
        a0, a1 = grid[i], grid[i + 1]
        a0h, a1h = a0 // 2, a1 // 2
        ra0 = (INV_R2_LO * a0) // UNIT
        ra1 = (INV_R2_LO * a1) // UNIT
        ics = i * cside
        for j in range(lim_j):
            b0, b1 = grid[j], grid[j + 1]
            M00_l = (ROOT_LO * b0) // UNIT + a0h
            if M00_l <= 0:
                continue
            M00_h = -(-(ROOT_HI * b0) // UNIT) + a0h
            M01_l = (ROOT_LO * b1) // UNIT + a0h
            M01_h = -(-(ROOT_HI * b1) // UNIT) + a0h
            M10_l = (ROOT_LO * b0) // UNIT + a1h
            M10_h = -(-(ROOT_HI * b0) // UNIT) + a1h
            M11_l = (ROOT_LO * b1) // UNIT + a1h
            M11_h = -(-(ROOT_HI * b1) // UNIT) + a1h
            U00 = h_cap_u if h_cap_u < M00_l else M00_l
            U01 = h_cap_u if h_cap_u < M01_l else M01_l
            U10 = h_cap_u if h_cap_u < M10_l else M10_l
            U11 = h_cap_u if h_cap_u < M11_l else M11_l
            U_gm = isqrt_(isqrt_(U00 * U01) * isqrt_(U10 * U11))
            if not U_gm:
                continue
            rb00 = (3 * b0) // 2 + ra0
            rb01 = (3 * b1) // 2 + ra0
            rb10 = (3 * b0) // 2 + ra1
            rb11 = (3 * b1) // 2 + ra1
            r0 = j * side
            r1 = (j + 1) * side
            nu = -(-U11 // (2 * UNIT))
            if nu < 6:
                nu = 6
            elif nu > 12:
                nu = 12
            do_plus = M00_h < bound
            prev_s = None
            tot = 0
            for kk in range(nu + 1):
                h00 = (kk * U00) // nu
                h00h = -(-(kk * U00) // nu)
                h01 = (kk * U01) // nu
                h01h = -(-(kk * U01) // nu)
                h10 = (kk * U10) // nu
                h10h = -(-(kk * U10) // nu)
                h11 = (kk * U11) // nu
                h11h = -(-(kk * U11) // nu)
                g3 = density(h11h)
                if not g3:
                    break
                cg = (density(h00h), density(h01h), density(h10h), g3)
                cum = (
                    half_mass(rb00 - (-(-(ROOT_HI * h00h) // UNIT))),
                    half_mass(rb01 - (-(-(ROOT_HI * h01h) // UNIT))),
                    half_mass(rb10 - (-(-(ROOT_HI * h10h) // UNIT))),
                    half_mass(rb11 - (-(-(ROOT_HI * h11h) // UNIT))),
                )
                cup = (
                    half_mass(rb00 + (ROOT_LO * h00) // UNIT),
                    half_mass(rb01 + (ROOT_LO * h01) // UNIT),
                    half_mass(rb10 + (ROOT_LO * h10) // UNIT),
                    half_mass(rb11 + (ROOT_LO * h11) // UNIT),
                )
                lm00 = locate(M00_l - h00h, M00_h - h00)
                lm01 = locate(M01_l - h01h, M01_h - h01) if lm00 else None
                lm10 = locate(M10_l - h10h, M10_h - h10) if lm01 else None
                lm11 = locate(M11_l - h11h, M11_h - h11) if lm10 else None
                if lm11:
                    b00 = eval_batch(WV_t, r0, lm00)
                    b01 = eval_batch(WV_t, r1, lm01)
                    b10 = eval_batch(WV_t, r0, lm10)
                    b11 = eval_batch(WV_t, r1, lm11)
                    cm = [
                        (b00[m_], b01[m_], b10[m_], b11[m_]) if (b00[m_] and b01[m_] and b10[m_] and b11[m_]) else None
                        for m_ in range(nWV)
                    ]
                else:
                    cm = None
                if kk == 0:
                    cp = cm
                elif do_plus and M11_h + h11h < bound:
                    lp00 = locate(M00_l + h00, M00_h + h00h)
                    lp01 = locate(M01_l + h01, M01_h + h01h) if lp00 else None
                    lp10 = locate(M10_l + h10, M10_h + h10h) if lp01 else None
                    lp11 = locate(M11_l + h11, M11_h + h11h) if lp10 else None
                    if lp11:
                        b00 = eval_batch(WV_t, r0, lp00)
                        b01 = eval_batch(WV_t, r1, lp01)
                        b10 = eval_batch(WV_t, r0, lp10)
                        b11 = eval_batch(WV_t, r1, lp11)
                        cp = [
                            (b00[m_], b01[m_], b10[m_], b11[m_]) if (b00[m_] and b01[m_] and b10[m_] and b11[m_]) else None
                            for m_ in range(nWV)
                        ]
                    else:
                        cp = None
                else:
                    cp = None
                cur_s = (cg, cum, cup, cm, cp)
                if prev_s is not None:
                    pg, pum, pup, pm, pp = prev_s
                    gw = isqrt_(
                        isqrt_(isqrt_(pg[0] * cg[0]) * isqrt_(pg[1] * cg[1]))
                        * isqrt_(isqrt_(pg[2] * cg[2]) * isqrt_(pg[3] * cg[3]))
                    )
                    if gw:
                        om = isqrt_(
                            isqrt_(isqrt_(pum[0] * cum[0]) * isqrt_(pum[1] * cum[1]))
                            * isqrt_(isqrt_(pum[2] * cum[2]) * isqrt_(pum[3] * cum[3]))
                        )
                        op = isqrt_(
                            isqrt_(isqrt_(pup[0] * cup[0]) * isqrt_(pup[1] * cup[1]))
                            * isqrt_(isqrt_(pup[2] * cup[2]) * isqrt_(pup[3] * cup[3]))
                        )
                        vm = 0
                        vp = 0
                        if pm and cm:
                            for m_ in range(nWt):
                                a_, b_ = pm[m_], cm[m_]
                                if a_ and b_:
                                    om += isqrt_(
                                        isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                        * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                    )
                            for m_ in range(nWt, nWV):
                                a_, b_ = pm[m_], cm[m_]
                                if a_ and b_:
                                    vm += isqrt_(
                                        isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                        * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                    )
                        if pp and cp:
                            for m_ in range(nWt):
                                a_, b_ = pp[m_], cp[m_]
                                if a_ and b_:
                                    op += isqrt_(
                                        isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                        * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                    )
                            for m_ in range(nWt, nWV):
                                a_, b_ = pp[m_], cp[m_]
                                if a_ and b_:
                                    vp += isqrt_(
                                        isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                        * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                    )
                        tot += (gw * ((om * (op + vp) + vm * op) // SCALE)) // SCALE
                prev_s = cur_s
            T_gm4[ics + j] = (U_gm * tot) // (nu * UNIT_N)

    E_gm4 = list(T_gm4)
    cur_gm4 = T_gm4
    for _ in range(4):
        mid_gm4 = apply_gm4_kern(cur_gm4, K_gm4_kerns)
        cur_gm4 = apply_gm4_kern(mid_gm4, Kg_gm4_kerns)
        E_gm4 = [a + b for a, b in zip(E_gm4, cur_gm4)]
    del K_gm4_kerns, Kg_gm4_kerns, T_gm4, cur_gm4, mid_gm4, kernels

    # 4. Point-level component tables for P_1 = S_1(E) (20 T pairs + 14 U-guarded hops = 34 components).
    m_mats = (W1_1, W1_3, W1_5, W2_11, V_0U, V_2U, V_4U, V_0W1, V_2W1, V_0W13)
    p_mats = (W1_1, W1_3, V_0U, V_2U)
    pair_specs = (
        (0, 0),            # 0: (U, U)
        (5, 0), (0, 3),    # 1,2: (V0U, U), (U, V0U)
        (1, 0), (0, 1),    # 3,4: (W1_1, U), (U, W1_1)
        (6, 0), (0, 4),    # 5,6: (V2U, U), (U, V2U)
        (2, 0), (0, 2),    # 7,8: (W1_3, U), (U, W1_3)
        (8, 0), (7, 0),    # 9,10: (V0W1, U), (V4U, U)
        (3, 0), (4, 0),    # 11,12: (W1_5, U), (W2_11, U)
        (9, 0), (10, 0),   # 13,14: (V2W1, U), (V0W13, U)
        (1, 3), (5, 1),    # 15,16: (W1_1, V0U), (V0U, W1_1)
        (1, 1),            # 17: (W1_1, W1_1)
        (2, 3), (5, 2),    # 18,19: (W1_3, V0U), (V0U, W1_3)
    )
    n_pairs = len(pair_specs)
    T_mats = [[0] * size for _ in range(n_pairs)]
    zero_m = (0,) * len(m_mats)
    zero_p = (0,) * len(p_mats)
    for i in points:
        lim_j = t_max_j[i]
        if not lim_j:
            continue
        gi = grid[i]
        gi_half = gi // 2
        iside = i * side
        r_gi = (INV_R2_LO * gi) // UNIT
        for j in range(lim_j):
            gj = grid[j]
            ml = (ROOT_LO * gj) // UNIT + gi_half
            if ml <= 0:
                continue
            mh = -(-(ROOT_HI * gj) // UNIT) + gi_half
            upper = h_cap_u if h_cap_u < ml else ml
            row = j * side
            r_base = (3 * gj) // 2 + r_gi
            kmax = find_idx(upper)
            prevs = [0] * n_pairs
            tots = [0] * n_pairs
            prev_h = 0
            hs = [grid[kk] for kk in range(kmax + 1)]
            if hs[-1] < upper:
                hs.append(upper)
            first = True
            for h in hs:
                dh = density(h)
                if not dh:
                    break
                lm = locate(ml - h, mh - h)
                u_m = half_mass(r_base - (-(-(ROOT_HI * h) // UNIT)))
                u_p = half_mass(r_base + (ROOT_LO * h) // UNIT)
                if lm is not None:
                    em = eval_batch(m_mats, row, lm)
                    dm = [(dh * u_m) // SCALE] + [(dh * v) // SCALE for v in em]
                else:
                    dm = [(dh * u_m) // SCALE] + list(zero_m)
                if h == 0:
                    vp = (u_p, em[0], em[1], em[4], em[5]) if lm is not None else ((u_p,) + zero_p)
                else:
                    lp = locate(ml + h, mh + h)
                    vp = ((u_p,) + eval_batch(p_mats, row, lp)) if lp is not None else ((u_p,) + zero_p)
                if first:
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_, b_ = dm[fi], vp[gi_]
                        prevs[t_idx] = (a_ * b_) // SCALE if (a_ and b_) else 0
                    first = False
                else:
                    step_w = h - prev_h
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_, b_ = dm[fi], vp[gi_]
                        cur = (a_ * b_) // SCALE if (a_ and b_) else 0
                        pv = prevs[t_idx]
                        if pv and cur:
                            tots[t_idx] += isqrt_(pv * cur) * step_w
                        prevs[t_idx] = cur
                prev_h = h
            idx_ij = iside + j
            for t_idx in range(n_pairs):
                T_mats[t_idx][idx_ij] = tots[t_idx] // UNIT_N

    K_T = [apply_k_steps(T_mats[k], K_steps) for k in (0, 1, 2, 3, 5, 7, 9, 15, 16)]
    h1_U = [apply_k_steps(m_, KU_steps) for m_ in K_T]
    h2_U = [apply_k_steps(apply_k_steps(h1_U[k], K_steps), KU_steps) for k in (0, 1, 3)]
    h3_U = [apply_k_steps(apply_k_steps(h2_U[0], K_steps), KU_steps), apply_k_steps(apply_k_steps(h2_U[1], K_steps), KU_steps)]
    E_pt = T_mats + h1_U + h2_U + h3_U
    P1_mats = [apply_S_steps(m_, S1_tbl_full, ix_max1 + 1) for m_ in E_pt]
    del T_mats, K_T, h1_U, h2_U, h3_U, E_pt, K_steps, KU_steps, S1_tbl_full

    def make_gm4_sum(mats, i_lim):
        out = [0] * (i_lim * cside)
        for i in range(i_lim):
            i0 = i * side
            i1 = i0 + side
            ics = i * cside
            for j in cpoints:
                tot = 0
                for mt in mats:
                    a, b_, c, d_ = mt[i0 + j], mt[i0 + j + 1], mt[i1 + j], mt[i1 + j + 1]
                    if a and b_ and c and d_:
                        tot += isqrt_(isqrt_(a * b_) * isqrt_(c * d_))
                out[ics + j] = tot
        return out

    BI_gm4 = make_gm4_sum(even, ix_max1)
    BO_gm4 = make_gm4_sum(odd, ix_max2)
    Og_gm4 = make_gm4_sum(Og_mats, ix_max2)

    # 5. Kind 1 (crossing degree 1: pos1, neg1) in coordinatewise affine (x, b, t).
    nA = len(A_mats)
    nP = len(P1_mats)
    d_ap = nA + nP
    AP_mats = tuple(A_mats + P1_mats)
    zero_ap = (0,) * d_ap

    def f2_opp(gx, gb):
        xq = -(-gx // v_div)
        bq = -(-gb // v_div)
        num = 9 * (xq * xq + bq * bq) + (-(-(12 * ROOT_HI * xq * bq) // UNIT))
        return (C2_NUM * exp_floor(-(-num // 3200))) // UNIT

    counts = [0, 0, 0, 0]
    old_sq = None
    for ix in range(ix_max1 + 1):
        gx = grid[ix]
        row_off = ix * side
        ap_cache = {}
        rows = []
        for jb in range(jb_max1 + 1):
            gb = grid[jb]
            w_f2 = f2_opp(gx, gb)
            if not w_f2:
                rows.append(None)
                continue
            s_lo = gb + (2 * ROOT_LO * gx) // UNIT
            jac = s_lo // (2 * n)
            if jac <= 0:
                rows.append(None)
                continue
            base_w = (w_f2 * jac) // SCALE
            if base_w * 3 < WEIGHT_CUT:
                rows.append(None)
                continue
            s_hi = gb + (-(-(2 * ROOT_HI * gx) // UNIT))
            row_kt = []
            for kt in range(nt1 + 1):
                h_hi = -(-(kt * s_hi) // (2 * nt1))
                w_all = (base_w * density(h_hi)) // SCALE
                if w_all < WEIGHT_CUT:
                    row_kt.append(None)
                    continue
                ylo = ((nt1 + kt) * s_lo) // (2 * nt1)
                yhi = -(-((nt1 + kt) * s_hi) // (2 * nt1))
                zlo = ((nt1 - kt) * s_lo) // (2 * nt1)
                zhi = -(-((nt1 - kt) * s_hi) // (2 * nt1))
                yk = (ylo, yhi)
                yv = ap_cache.get(yk)
                if yv is None:
                    loc = locate(ylo, yhi)
                    yv = zero_ap if loc is None else eval_batch(AP_mats, row_off, loc)
                    ap_cache[yk] = yv
                zk = (zlo, zhi)
                zv = ap_cache.get(zk)
                if zv is None:
                    loc = locate(zlo, zhi)
                    zv = zero_ap if loc is None else eval_batch(AP_mats, row_off, loc)
                    ap_cache[zk] = zv
                row_kt.append((w_all, yv, zv))
            rows.append(row_kt)

        lines = []
        for jb in range(jb_max1 + 1):
            rk = rows[jb]
            if rk is None:
                lines.append(None)
                continue
            lrow = []
            for kt in range(nt1):
                p0, p1 = rk[kt], rk[kt + 1]
                if p0 is None or p1 is None:
                    lrow.append(None)
                    continue
                w0, y0, z0 = p0
                w1, y1, z1 = p1
                wl = isqrt_(w0 * w1)
                if not wl:
                    lrow.append(None)
                    continue
                lrow.append((wl, list(map(isqrt_, map(_mul, y0, y1))), list(map(isqrt_, map(_mul, z0, z1)))))
            lines.append(lrow)

        squares = []
        for jb in range(jb_max1):
            l0, l1 = lines[jb], lines[jb + 1]
            if l0 is None or l1 is None:
                squares.append(None)
                continue
            srow = []
            for kt in range(nt1):
                p0, p1 = l0[kt], l1[kt]
                if p0 is None or p1 is None:
                    srow.append(None)
                    continue
                w0, y0, z0 = p0
                w1, y1, z1 = p1
                ws = isqrt_(w0 * w1)
                if not ws:
                    srow.append(None)
                    continue
                srow.append((ws, list(map(isqrt_, map(_mul, y0, y1))), list(map(isqrt_, map(_mul, z0, z1)))))
            squares.append(srow)

        if ix > 0:
            dx8 = dg[ix - 1] // U0125
            ics = (ix - 1) * cside
            for jb in range(jb_max1):
                bi_gm = BI_gm4[ics + jb]
                e_gm = E_gm4[ics + jb]
                if not bi_gm and not e_gm:
                    continue
                sr0, sr1 = old_sq[jb], squares[jb]
                if sr0 is None or sr1 is None:
                    continue
                vol = dx8 * (dg[jb] // U0125)
                s0 = 0
                s1 = 0
                for kt in range(nt1):
                    c0, c1 = sr0[kt], sr1[kt]
                    if c0 is None or c1 is None:
                        continue
                    w0, y0, z0 = c0
                    w1, y1, z1 = c1
                    w = isqrt_(w0 * w1)
                    if not w:
                        continue
                    if bi_gm:
                        ay = sum(map(isqrt_, map(_mul, y0[:nA], y1[:nA])))
                        az = sum(map(isqrt_, map(_mul, z0[:nA], z1[:nA]))) if ay else 0
                        if az:
                            s0 += (((((w * ay) // SCALE) * az) // SCALE) * bi_gm) // SCALE
                    if e_gm:
                        py = sum(map(isqrt_, map(_mul, y0[nA:], y1[nA:])))
                        pz = sum(map(isqrt_, map(_mul, z0[nA:], z1[nA:]))) if py else 0
                        if pz:
                            s1 += (((((w * py) // SCALE) * pz) // SCALE) * e_gm) // SCALE
                counts[0] += 3 * vol * s0
                counts[1] += 3 * vol * s1
        old_sq = squares

    # 6. Kind 2 (crossing degree 2: pos2, neg2) in (x, y, z).
    nA0 = len(same_even)
    nP1 = len(P1_mats)
    AP2 = tuple(same_even + P1_mats)
    z_ap2 = (0,) * (nA0 + nP1)

    def w_kind2(gx, gy, gz):
        xq = -(-gx // v_div)
        yq = -(-(gy + (-(-(RHO_HI * gx) // UNIT))) // v_div)
        zq = -(-(2 * gz + gy + (-(-(2 * ROOT_HI * gx) // UNIT))) // v_div)
        q = -(-(xq * xq + 9 * yq * yq + 3 * zq * zq) // 3200)
        return (C3_NUM * exp_floor(q)) // UNIT

    sums_set = {grid[jy] + grid[kz] for jy in range(jy_max2 + 1) for kz in range(jy, jy_max2 + 1)}

    def get_b_vals(ix):
        r = ix * side
        gx = grid[ix]
        r2_lo = (2 * ROOT_LO * gx) // UNIT
        r2_hi = -(-(2 * ROOT_HI * gx) // UNIT)
        out_c = {}
        for s in sums_set:
            loc = locate(r2_lo + s, r2_hi + s)
            out_c[s] = z_ap2 if loc is None else eval_batch(AP2, r, loc)
        return out_c

    w_prev_sq = None
    b_prev = get_b_vals(0)
    for ix in range(ix_max2 + 1):
        gx = grid[ix]
        wplane = [[w_kind2(gx, grid[jy], grid[kz]) for kz in range(jy_max2 + 1)] for jy in range(jy_max2 + 1)]
        wsq = [
            [
                isqrt_(isqrt_(wplane[jy][k] * wplane[jy][k + 1]) * isqrt_(wplane[jy + 1][k] * wplane[jy + 1][k + 1]))
                for k in range(jy_max2)
            ]
            for jy in range(jy_max2)
        ]
        if ix > 0:
            dx8 = dg[ix - 1] // U0125
            ics = (ix - 1) * cside
            b_cur = get_b_vals(ix)
            ap_gm_cache = {}

            def get_ap_gm(s00, s01, s10, s11):
                key = (s00, s01, s10, s11)
                res = ap_gm_cache.get(key)
                if res is not None:
                    return res
                p00, p01, p10, p11 = b_prev[s00], b_prev[s01], b_prev[s10], b_prev[s11]
                c00, c01, c10, c11 = b_cur[s00], b_cur[s01], b_cur[s10], b_cur[s11]
                sa = 0
                for idx_m in range(nA0):
                    a1, a2, a3, a4 = p00[idx_m], p01[idx_m], p10[idx_m], p11[idx_m]
                    if a1 and a2 and a3 and a4:
                        b1, b2, b3, b4 = c00[idx_m], c01[idx_m], c10[idx_m], c11[idx_m]
                        if b1 and b2 and b3 and b4:
                            sa += isqrt_(isqrt_(isqrt_(a1 * a2) * isqrt_(a3 * a4)) * isqrt_(isqrt_(b1 * b2) * isqrt_(b3 * b4)))
                sp = 0
                for idx_m in range(nA0, nA0 + nP1):
                    a1, a2, a3, a4 = p00[idx_m], p01[idx_m], p10[idx_m], p11[idx_m]
                    if a1 and a2 and a3 and a4:
                        b1, b2, b3, b4 = c00[idx_m], c01[idx_m], c10[idx_m], c11[idx_m]
                        if b1 and b2 and b3 and b4:
                            sp += isqrt_(isqrt_(isqrt_(a1 * a2) * isqrt_(a3 * a4)) * isqrt_(isqrt_(b1 * b2) * isqrt_(b3 * b4)))
                ap_gm_cache[key] = (sa, sp)
                return sa, sp

            for jy in range(jy_max2):
                dy8 = dg[jy] // U0125
                gy0, gy1 = grid[jy], grid[jy + 1]
                w0_row = w_prev_sq[jy]
                w1_row = wsq[jy]
                bi_y = BI_gm4[ics + jy]
                bo_y = BO_gm4[ics + jy]
                e_y = E_gm4[ics + jy]
                og_y = Og_gm4[ics + jy]
                for kz in range(jy_max2):
                    w = isqrt_(w0_row[kz] * w1_row[kz])
                    if w < WEIGHT_CUT:
                        continue
                    dz8 = dg[kz] // U0125
                    gz0, gz1 = grid[kz], grid[kz + 1]
                    a0, p1_ = get_ap_gm(gy0 + gz0, gy0 + gz1, gy1 + gz0, gy1 + gz1)
                    vol = dx8 * dy8 * dz8
                    bi_z = BI_gm4[ics + kz]
                    if a0:
                        counts[2] += 3 * vol * ((((w * a0) // SCALE) * (bi_y + (bo_y * bi_z) // SCALE)) // SCALE)
                    if p1_:
                        counts[3] += (
                            3
                            * vol
                            * ((((w * p1_) // SCALE) * ((bi_y * bi_z + e_y * Og_gm4[ics + kz] + og_y * E_gm4[ics + kz]) // SCALE)) // SCALE)
                        )
            b_prev = b_cur
        w_prev_sq = wsq

    den1 = 64 * SCALE * n * n * nt1
    den2 = 512 * SCALE * n * n * n
    return (counts[0], counts[1]), den1, (counts[2], counts[3]), den2


if __name__ == "__main__":
    import time

    t0 = time.time()
    c1, d1, c2, d2 = cavity(64, 16)
    dt = time.time() - t0
    p0, p1 = c1[0] / d1, c1[1] / d1
    p2, p3 = c2[0] / d2, c2[1] / d2
    print(
        f"[n=64] {dt:.2f}s | c1={c1}/{d1} c2={c2}/{d2} | "
        f"tot={(p0+p1+p2+p3)*100:.5f}% ({p0*100:.5f}%+{p1*100:.5f}%+{p2*100:.5f}%+{p3*100:.5f}%)"
    )
