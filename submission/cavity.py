"""Integer-certified cavity lower bound for the four disjoint local saving events.

Uses:
  - A graded grid on [0, 2] with 1/(8n) refinement near zero and exact integer spacings.
  - Exact two-bisection + weighted-harmonic geometric interpolation (zero t-rounding loss).
  - Exact u-scaled Holder-compressed 4-corner geometric-mean tables for E (all disjoint
    (OUT x OUT), (V x OUT), (OUT x V) endpoint pairs and four guarded two-step hops).
  - Direct coordinatewise log-concave integration for Kind 1 in (x, b, t) and Kind 2 in (x, y, z).
"""
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
    num2 = num << 1
    if num2 >= den:
        lo = _isqrt(lo * hi)
        num = num2 - den
    else:
        hi = _isqrt(lo * hi)
        num = num2
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
    if i > 20_000_000:
        return 0
    h = 256 * EXP_DEN
    u = ceil_div(i * SCALE, h)
    l = (i * SCALE) // h
    v = SCALE - u + (l * i) // (2 * h) - ceil_div(ceil_div(u * i, 2 * h) * i, 3 * h)
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
    k1 = (5 * n) // 8
    k2 = (9 * n) // 8
    k3 = (3 * n) // 2
    for k in range(2, k1 + 1):
        pts.append(k * U)
    for k in range(k1 + 2, k2 + 1, 2):
        pts.append(k * U)
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


def cavity(n=56, nt1=16):
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
    _bisect = bisect_right

    V_QUANT = 40_000
    v_div = UNIT_N // V_QUANT
    assert UNIT_N % V_QUANT == 0
    assert all(g % 2 == 0 and g % U0125 == 0 for g in grid)

    @lru_cache(maxsize=None)
    def density(v):
        t = ceil_div(abs(v), v_div)
        q = ceil_div(3 * t * t, 800)
        return (CG_NUM * exp_floor(q)) // UNIT

    mass = [0]
    m2 = mesh * mesh
    for i in range(1, 2 * mesh + 1):
        q = ceil_div(2 * EXP_DEN * (3 * i * i - 3 * i + 1), m2)
        mass.append(mass[-1] + (CG_NUM * exp_floor(q)) // (UNIT * mesh))

    def half_mass(r):
        if r <= 0:
            return 0
        idx = (mesh * r) // UNIT_N
        return mass[2 * mesh if idx > 2 * mesh else idx]

    def find_idx(y):
        j = _bisect(grid, y) - 1
        return cside - 1 if j >= cside else (0 if j < 0 else j)

    def locate(y_lo, y_hi):
        if y_lo < 0 or y_hi >= bound:
            return None
        j0 = _bisect(grid, y_lo) - 1
        if j0 >= cside:
            j0 = cside - 1
        r0 = y_lo - grid[j0]
        d0 = dg[j0]
        if y_hi == y_lo:
            return (j0, r0, d0, -1, 0, 0)
        j1 = _bisect(grid, y_hi) - 1
        if j1 >= cside:
            j1 = cside - 1
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

    def row_interp(mat, row_off, y_lo, y_hi):
        loc = locate(y_lo, y_hi)
        return eval_loc(mat, row_off, loc) if loc is not None else 0

    ix_max1 = find_idx((7 * UNIT_N) // 10) + 1
    jb_max1 = find_idx((6 * UNIT_N) // 5) + 1
    ix_max2 = find_idx(UNIT_N // 2) + 1
    jy_max2 = find_idx(UNIT_N) + 1
    points = range(side)
    cpoints = range(cside)

    # 1. Precompute positive Gaussian transition kernels on the graded grid.
    kernels = []
    for i in points:
        gi_half = grid[i] // 2
        for j in points:
            base = ceil_div(ROOT_HI * grid[j], UNIT) + gi_half
            krow = []
            for k in points:
                val = density(base + grid[k])
                if not val:
                    break
                krow.append(val)
            kernels.append(krow)

    # 2. Alternating crossing branches B_0 .. B_6.
    branches = [[half_mass((ROOT_LO * grid[j]) // UNIT + grid[i] // 2) for i in points for j in points]]
    for _ in range(6):
        old = branches[-1]
        mat = [0] * size
        for i in points:
            iside = i * side
            for j in points:
                kern = kernels[iside + j]
                nk = len(kern)
                if nk < 2:
                    continue
                row = j * side
                prev = (kern[0] * old[row]) // SCALE
                tot = 0
                for k in range(1, nk):
                    cur = (kern[k] * old[row + k]) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * dg[k - 1]
                    prev = cur
                mat[iside + j] = tot // UNIT_N
        branches.append(mat)

    def op_K(mat, i_lim=side):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
            iside = i * side
            for j in points:
                kern = kernels[iside + j]
                nk = len(kern)
                if nk < 2:
                    continue
                row = j * side
                prev = (kern[0] * mat[row]) // SCALE
                tot = 0
                for k in range(1, nk):
                    cur = (kern[k] * mat[row + k]) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * dg[k - 1]
                    prev = cur
                out[iside + j] = tot // UNIT_N
        return out

    def op_S(mat, strict, i_lim=side):
        out = [0] * (i_lim * side)
        mult = 2 if strict else 1
        for i in range(i_lim):
            gi_half = grid[i] // 2
            iside = i * side
            for j in points:
                gj = grid[j]
                ml = (ROOT_LO * gj) // UNIT - gi_half
                mh = ceil_div(ROOT_HI * gj, UNIT) - gi_half
                c = -mult * ml if ml < 0 else 0
                if c >= bound:
                    continue
                row = j * side
                k0 = find_idx(c)
                rem = c - grid[k0]
                v0 = geo_(mat[row + k0], mat[row + k0 + 1], rem, dg[k0]) if rem else mat[row + k0]
                prev = (v0 * density(c + mh if c + mh > 0 else 0)) // SCALE
                tot = 0
                for kk in range(k0 + 1, side):
                    edge = grid[kk]
                    w = density(edge + mh)
                    if not w:
                        break
                    cur = (w * mat[row + kk]) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * (edge - c)
                    prev, c = cur, edge
                out[iside + j] = tot // UNIT_N
        return out

    def op_D(f, g, i_lim=side):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
            gi_half = grid[i] // 2
            iside = i * side
            for j in points:
                gj = grid[j]
                s0 = gi_half - ceil_div(ROOT_HI * gj, UNIT)
                if s0 <= 0:
                    continue
                s1 = gi_half - (ROOT_LO * gj) // UNIT
                row = j * side
                upper = bound if bound < s0 else s0
                kmax = find_idx(upper)
                prev_h = 0
                loc0 = locate(s0, s1)
                prev_v = ((density(0) * eval_loc(f, row, loc0)) // SCALE * eval_loc(g, row, loc0)) // SCALE if loc0 else 0
                tot = 0
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
                    a = eval_loc(f, row, lp) if lp else 0
                    b = eval_loc(g, row, lm) if lm else 0
                    cur = ((dh * a) // SCALE * b) // SCALE
                    if prev_v and cur:
                        tot += isqrt_(prev_v * cur) * (h - prev_h)
                    prev_h, prev_v = h, cur
                    if h == upper:
                        break
                out[iside + j] = tot // UNIT_N
        return out

    def op_K_guard_U(mat, i_lim=side):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
            gi = grid[i]
            iside = i * side
            r_gi = (ROOT_LO * gi) // UNIT
            for j in points:
                kern = kernels[iside + j]
                nk = len(kern)
                if nk < 2:
                    continue
                row = j * side
                base_u = (7 * grid[j]) // 2 + r_gi
                prev = ((kern[0] * mat[row]) // SCALE * half_mass(base_u)) // SCALE
                tot = 0
                for k in range(1, nk):
                    cur = ((kern[k] * mat[row + k]) // SCALE * half_mass(base_u + (ROOT_LO * grid[k]) // UNIT)) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * dg[k - 1]
                    prev = cur
                out[iside + j] = tot // UNIT_N
        return out

    def op_K_guard_mat(mat, gd, i_lim=side):
        out = [0] * (i_lim * side)
        for i in range(i_lim):
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
                c_hi = ceil_div(2 * ROOT_HI * gj, UNIT) + gi
                prev = ((kern[0] * mat[row]) // SCALE * row_interp(gd, row, c_lo, c_hi)) // SCALE
                tot = 0
                for k in range(1, nk):
                    gk = grid[k]
                    gval = row_interp(gd, row, c_lo + gk, c_hi + gk)
                    cur = ((kern[k] * mat[row + k]) // SCALE * gval) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * dg[k - 1]
                    prev = cur
                out[iside + j] = tot // UNIT_N
        return out

    def make_V_U(bi_mat):
        out = [0] * size
        for i in points:
            gi = grid[i]
            gi_half = gi // 2
            iside = i * side
            r_gi_hi = ceil_div(ROOT_HI * gi, UNIT)
            for j in points:
                gj = grid[j]
                ml = (ROOT_LO * gj) // UNIT - gi_half
                mh = ceil_div(ROOT_HI * gj, UNIT) - gi_half
                c = -2 * ml if ml < 0 else 0
                if c >= bound:
                    continue
                row = j * side
                level = (7 * gj) // 2 - r_gi_hi
                k0 = find_idx(c)
                rem = c - grid[k0]
                b0 = geo_(bi_mat[row + k0], bi_mat[row + k0 + 1], rem, dg[k0]) if rem else bi_mat[row + k0]
                prev = ((density(c + mh if c + mh > 0 else 0) * b0) // SCALE * half_mass(level + (ROOT_LO * c) // UNIT)) // SCALE
                tot = 0
                for kk in range(k0 + 1, side):
                    edge = grid[kk]
                    w = density(edge + mh)
                    if not w:
                        break
                    cur = ((w * bi_mat[row + kk]) // SCALE * half_mass(level + (ROOT_LO * edge) // UNIT)) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * (edge - c)
                    prev, c = cur, edge
                out[iside + j] = tot // UNIT_N
        return out

    def make_V_mat(bi_mat, gd_mat):
        out = [0] * size
        for i in points:
            gi_half = grid[i] // 2
            iside = i * side
            for j in points:
                gj = grid[j]
                ml = (ROOT_LO * gj) // UNIT - gi_half
                mh = ceil_div(ROOT_HI * gj, UNIT) - gi_half
                c = -2 * ml if ml < 0 else 0
                if c >= bound:
                    continue
                row = j * side
                k0 = find_idx(c)
                rem = c - grid[k0]
                b0 = geo_(bi_mat[row + k0], bi_mat[row + k0 + 1], rem, dg[k0]) if rem else bi_mat[row + k0]
                g0 = row_interp(gd_mat, row, max(0, 2 * ml + c), 2 * mh + c)
                prev = ((density(c + mh if c + mh > 0 else 0) * b0) // SCALE * g0) // SCALE
                tot = 0
                for kk in range(k0 + 1, side):
                    edge = grid[kk]
                    w = density(edge + mh)
                    if not w:
                        break
                    gval = row_interp(gd_mat, row, 2 * ml + edge, 2 * mh + edge)
                    cur = ((w * bi_mat[row + kk]) // SCALE * gval) // SCALE
                    if prev and cur:
                        tot += isqrt_(prev * cur) * (edge - c)
                    prev, c = cur, edge
                out[iside + j] = tot // UNIT_N
        return out

    even = [branches[0], branches[2], branches[4], branches[6]]
    odd = [branches[1], branches[3], branches[5]]

    # Positive neighbor cover events A_0 and A.
    same_even = [op_S(mat, False, ix_max1 + 1) for mat in even]
    A_mats = list(same_even)
    for a in odd[:2]:
        for b in even[:3]:
            A_mats.append(op_D(a, b, ix_max1 + 1))

    # Disjoint OUT components (beyond U): W_1 = {S_1(B_1), S_1(B_3), S_1(B_5)}, W_2 = {D(B_1, B_1)}.
    W1_1 = op_S(branches[1], True)
    W1_3 = op_S(branches[3], True)
    W1_5 = op_S(branches[5], True)
    W2_11 = op_D(branches[1], branches[1])

    # Disjoint V_odd_in components (1 odd-arm neg child + 1 pos child in OUT).
    V_0U = make_V_U(branches[0])
    V_2U = make_V_U(branches[2])
    V_4U = make_V_U(branches[4])
    V_0W1 = make_V_mat(branches[0], W1_1)
    V_2W1 = make_V_mat(branches[2], W1_1)
    V_0W13 = make_V_mat(branches[0], W1_3)

    # Guarded far-arm events O_g (needed only up to ix_max2, i.e., x <= 29/56 in Kind 2).
    Og_mats = (
        [op_K_guard_U(m_, ix_max2 + 1) for m_ in even]
        + [op_K_guard_mat(m_, W1_1, ix_max2 + 1) for m_ in even[:3]]
        + [op_K_guard_mat(branches[0], W1_3, ix_max2 + 1)]
    )

    # 3. Exact Holder-compressed GM4 table E_gm4 for ALL disjoint endpoint pairs and 4 guarded hops.
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

    Kg_gm4_kerns = []
    W_guard = (W1_1, W1_3)
    for i in cpoints:
        a0, a1 = grid[i], grid[i + 1]
        ra0 = (ROOT_LO * a0) // UNIT
        ra1 = (ROOT_LO * a1) // UNIT
        ics = i * cside
        for j in cpoints:
            klist = K_gm4_kerns[ics + j]
            if not klist:
                Kg_gm4_kerns.append([])
                continue
            b0, b1 = grid[j], grid[j + 1]
            u00 = (7 * b0) // 2 + ra0
            u01 = (7 * b1) // 2 + ra0
            u10 = (7 * b0) // 2 + ra1
            u11 = (7 * b1) // 2 + ra1
            c00_l = (2 * ROOT_LO * b0) // UNIT + a0
            c00_h = ceil_div(2 * ROOT_HI * b0, UNIT) + a0
            c01_l = (2 * ROOT_LO * b1) // UNIT + a0
            c01_h = ceil_div(2 * ROOT_HI * b1, UNIT) + a0
            c10_l = (2 * ROOT_LO * b0) // UNIT + a1
            c10_h = ceil_div(2 * ROOT_HI * b0, UNIT) + a1
            c11_l = (2 * ROOT_LO * b1) // UNIT + a1
            c11_h = ceil_div(2 * ROOT_HI * b1, UNIT) + a1
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
                for wm in W_guard:
                    wa = []
                    for l00, l01, l10, l11 in locs_r:
                        v00 = eval_loc(wm, r0, l00) if l00 else 0
                        wa.append(
                            (
                                v00,
                                eval_loc(wm, r1, l01) if l01 else 0,
                                eval_loc(wm, r0, l10) if l10 else 0,
                                eval_loc(wm, r1, l11) if l11 else 0,
                            )
                            if v00
                            else None
                        )
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
            Kg_gm4_kerns.append(row_kg)

    def apply_gm4_kern(in_gm4, kerns_list):
        out = [0] * csize
        denom_k = SCALE * UNIT_N
        for i in cpoints:
            ics = i * cside
            for j in cpoints:
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
    nWt = len(W_t)
    nVt = len(V_t)
    T_gm4 = [0] * csize
    for i in cpoints:
        a0, a1 = grid[i], grid[i + 1]
        a0h, a1h = a0 // 2, a1 // 2
        ra0 = (INV_R2_LO * a0) // UNIT
        ra1 = (INV_R2_LO * a1) // UNIT
        ics = i * cside
        for j in cpoints:
            b0, b1 = grid[j], grid[j + 1]
            M00_l = (ROOT_LO * b0) // UNIT + a0h
            if M00_l <= 0:
                continue
            M00_h = ceil_div(ROOT_HI * b0, UNIT) + a0h
            M01_l = (ROOT_LO * b1) // UNIT + a0h
            M01_h = ceil_div(ROOT_HI * b1, UNIT) + a0h
            M10_l = (ROOT_LO * b0) // UNIT + a1h
            M10_h = ceil_div(ROOT_HI * b0, UNIT) + a1h
            M11_l = (ROOT_LO * b1) // UNIT + a1h
            M11_h = ceil_div(ROOT_HI * b1, UNIT) + a1h
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
            nu = ceil_div(U11, 2 * UNIT)
            if nu < 6:
                nu = 6
            elif nu > 16:
                nu = 16
            do_plus = M00_h < bound
            prev_s = None
            tot = 0
            for kk in range(nu + 1):
                h00 = (kk * U00) // nu
                h00h = ceil_div(kk * U00, nu)
                h01 = (kk * U01) // nu
                h01h = ceil_div(kk * U01, nu)
                h10 = (kk * U10) // nu
                h10h = ceil_div(kk * U10, nu)
                h11 = (kk * U11) // nu
                h11h = ceil_div(kk * U11, nu)
                g3 = density(h11h)
                if not g3:
                    break
                cg = (density(h00h), density(h01h), density(h10h), g3)
                cum = (
                    half_mass(rb00 - ceil_div(ROOT_HI * h00h, UNIT)),
                    half_mass(rb01 - ceil_div(ROOT_HI * h01h, UNIT)),
                    half_mass(rb10 - ceil_div(ROOT_HI * h10h, UNIT)),
                    half_mass(rb11 - ceil_div(ROOT_HI * h11h, UNIT)),
                )
                cup = (
                    half_mass(rb00 + (ROOT_LO * h00) // UNIT),
                    half_mass(rb01 + (ROOT_LO * h01) // UNIT),
                    half_mass(rb10 + (ROOT_LO * h10) // UNIT),
                    half_mass(rb11 + (ROOT_LO * h11) // UNIT),
                )
                lm00 = locate(M00_l - h00h, M00_h - h00)
                lm01 = locate(M01_l - h01h, M01_h - h01)
                lm10 = locate(M10_l - h10h, M10_h - h10)
                lm11 = locate(M11_l - h11h, M11_h - h11)
                m_ok = lm00 and lm01 and lm10 and lm11
                if do_plus and M11_h + h11h < bound:
                    lp00 = locate(M00_l + h00, M00_h + h00h)
                    lp01 = locate(M01_l + h01, M01_h + h01h)
                    lp10 = locate(M10_l + h10, M10_h + h10h)
                    lp11 = locate(M11_l + h11, M11_h + h11h)
                    p_ok = lp00 and lp01 and lp10 and lp11
                else:
                    p_ok = False
                cwm = []
                cwp = []
                cvm = []
                cvp = []
                for wm in W_t:
                    if m_ok:
                        v00 = eval_loc(wm, r0, lm00)
                        cwm.append(
                            (v00, eval_loc(wm, r1, lm01), eval_loc(wm, r0, lm10), eval_loc(wm, r1, lm11))
                            if v00
                            else None
                        )
                    else:
                        cwm.append(None)
                    if p_ok:
                        p00 = eval_loc(wm, r0, lp00)
                        cwp.append(
                            (p00, eval_loc(wm, r1, lp01), eval_loc(wm, r0, lp10), eval_loc(wm, r1, lp11))
                            if p00
                            else None
                        )
                    else:
                        cwp.append(None)
                for vm_ in V_t:
                    if m_ok:
                        v00 = eval_loc(vm_, r0, lm00)
                        cvm.append(
                            (v00, eval_loc(vm_, r1, lm01), eval_loc(vm_, r0, lm10), eval_loc(vm_, r1, lm11))
                            if v00
                            else None
                        )
                    else:
                        cvm.append(None)
                    if p_ok:
                        p00 = eval_loc(vm_, r0, lp00)
                        cvp.append(
                            (p00, eval_loc(vm_, r1, lp01), eval_loc(vm_, r0, lp10), eval_loc(vm_, r1, lp11))
                            if p00
                            else None
                        )
                    else:
                        cvp.append(None)
                cur_s = (cg, cum, cup, cwm, cwp, cvm, cvp)
                if prev_s is not None:
                    pg, pum, pup, pwm, pwp, pvm, pvp = prev_s
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
                        for m_ in range(nWt):
                            a_, b_ = pwm[m_], cwm[m_]
                            if a_ and b_:
                                om += isqrt_(
                                    isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                    * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                )
                            a_, b_ = pwp[m_], cwp[m_]
                            if a_ and b_:
                                op += isqrt_(
                                    isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                    * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                )
                        vm = 0
                        vp = 0
                        for m_ in range(nVt):
                            a_, b_ = pvm[m_], cvm[m_]
                            if a_ and b_:
                                vm += isqrt_(
                                    isqrt_(isqrt_(a_[0] * b_[0]) * isqrt_(a_[1] * b_[1]))
                                    * isqrt_(isqrt_(a_[2] * b_[2]) * isqrt_(a_[3] * b_[3]))
                                )
                            a_, b_ = pvp[m_], cvp[m_]
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

    # 4. Point-level component tables for P_1 = S_1(E) (20 T pairs + 16 guarded hops = 36 components).
    # Integer-indexed left/right tables for fast inner loop:
    # Left indices: 0=U, 1=W1_1, 2=W1_3, 3=W1_5, 4=W2_11, 5=V0U, 6=V2U, 7=V4U, 8=V0W1, 9=V2W1, 10=V0W13
    # Right indices: 0=U, 1=W1_1, 2=W1_3, 3=V0U, 4=V2U
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
    for i in points:
        gi = grid[i]
        gi_half = gi // 2
        iside = i * side
        r_gi = (INV_R2_LO * gi) // UNIT
        for j in points:
            gj = grid[j]
            ml = (ROOT_LO * gj) // UNIT + gi_half
            if ml <= 0:
                continue
            mh = ceil_div(ROOT_HI * gj, UNIT) + gi_half
            upper = bound if bound < ml else ml
            row = j * side
            r_base = (3 * gj) // 2 + r_gi
            kmax = find_idx(upper)
            prevs = [0] * n_pairs
            tots = [0] * n_pairs
            prev_h = 0
            for kk in range(kmax + 2):
                h = grid[kk] if kk <= kmax else upper
                if h > upper:
                    h = upper
                if kk > 0 and h <= prev_h:
                    continue
                dh = density(h)
                if not dh:
                    break
                lm = locate(ml - h, mh - h)
                lp = locate(ml + h, mh + h)
                u_m = half_mass(r_base - ceil_div(ROOT_HI * h, UNIT))
                u_p = half_mass(r_base + (ROOT_LO * h) // UNIT)
                if lm is not None:
                    dm = [(dh * u_m) // SCALE] + [(dh * eval_loc(mt, row, lm)) // SCALE for mt in m_mats]
                else:
                    dm = [(dh * u_m) // SCALE, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                if lp is not None:
                    vp = [u_p] + [eval_loc(mt, row, lp) for mt in p_mats]
                else:
                    vp = [u_p, 0, 0, 0, 0]
                step_w = h - prev_h
                if kk == 0:
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_ = dm[fi]
                        b_ = vp[gi_]
                        prevs[t_idx] = (a_ * b_) // SCALE if (a_ and b_) else 0
                else:
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_ = dm[fi]
                        b_ = vp[gi_]
                        cur = (a_ * b_) // SCALE if (a_ and b_) else 0
                        pv = prevs[t_idx]
                        if pv and cur:
                            tots[t_idx] += isqrt_(pv * cur) * step_w
                        prevs[t_idx] = cur
                prev_h = h
                if h == upper:
                    break
            idx_ij = iside + j
            for t_idx in range(n_pairs):
                T_mats[t_idx][idx_ij] = tots[t_idx] // UNIT_N

    K_T = [op_K(T_mats[k]) for k in (0, 1, 2, 3, 5, 7, 9, 15, 16)]
    h1_U = [op_K_guard_U(m_) for m_ in K_T]
    h1_W1 = [op_K_guard_mat(K_T[0], W1_1), op_K_guard_mat(K_T[1], W1_1)]
    h2_U = [op_K_guard_U(op_K(h1_U[k])) for k in (0, 1, 3)]
    h3_U = [op_K_guard_U(op_K(h2_U[0])), op_K_guard_U(op_K(h2_U[1]))]
    E_pt = T_mats + h1_U + h1_W1 + h2_U + h3_U
    P1_mats = [op_S(m_, True, ix_max1 + 1) for m_ in E_pt]

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
    d1 = 1 + 2 * d_ap
    w_row = (nt1 + 1) * d1
    stride1 = nt1 * d1
    zero_ap = (0,) * d_ap

    def f2_opp(gx, gb):
        xq = ceil_div(gx, v_div)
        bq = ceil_div(gb, v_div)
        num = 9 * (xq * xq + bq * bq) + ceil_div(12 * ROOT_HI * xq * bq, UNIT)
        return (C2_NUM * exp_floor(ceil_div(num, 3200))) // UNIT

    counts = [0, 0, 0, 0]
    old_sq = None
    for ix in range(ix_max1 + 1):
        gx = grid[ix]
        row_off = ix * side
        plane = []
        ap_cache = {}
        for jb in range(jb_max1 + 1):
            gb = grid[jb]
            w_f2 = f2_opp(gx, gb)
            if not w_f2:
                plane.extend((0,) * w_row)
                continue
            s_lo = gb + (2 * ROOT_LO * gx) // UNIT
            jac = s_lo // (2 * n)
            if jac <= 0:
                plane.extend((0,) * w_row)
                continue
            s_hi = gb + ceil_div(2 * ROOT_HI * gx, UNIT)
            base_w = (w_f2 * jac) // SCALE
            for kt in range(nt1 + 1):
                h_hi = ceil_div(kt * s_hi, 2 * nt1)
                w_all = (base_w * density(h_hi)) // SCALE
                if not w_all:
                    plane.extend((0,) * d1)
                    continue
                ylo = ((nt1 + kt) * s_lo) // (2 * nt1)
                yhi = ceil_div((nt1 + kt) * s_hi, 2 * nt1)
                zlo = ((nt1 - kt) * s_lo) // (2 * nt1)
                zhi = ceil_div((nt1 - kt) * s_hi, 2 * nt1)
                yk = (ylo, yhi)
                yv = ap_cache.get(yk)
                if yv is None:
                    loc = locate(ylo, yhi)
                    yv = zero_ap if loc is None else tuple(eval_loc(mt, row_off, loc) for mt in AP_mats)
                    ap_cache[yk] = yv
                zk = (zlo, zhi)
                zv = ap_cache.get(zk)
                if zv is None:
                    loc = locate(zlo, zhi)
                    zv = zero_ap if loc is None else tuple(eval_loc(mt, row_off, loc) for mt in AP_mats)
                    ap_cache[zk] = zv
                plane.append(w_all)
                plane.extend(yv)
                plane.extend(zv)

        lines = []
        for jb in range(jb_max1 + 1):
            st = jb * w_row
            lines.extend(isqrt_(a * b) for a, b in zip(plane[st : st + stride1], plane[st + d1 : st + w_row]))
        squares = [isqrt_(a * b) for a, b in zip(lines, lines[stride1:])]
        if ix > 0:
            dx8 = dg[ix - 1] // U0125
            ics = (ix - 1) * cside
            cube = [isqrt_(a * b) for a, b in zip(old_sq, squares)]
            for jb in range(jb_max1):
                bi_gm = BI_gm4[ics + jb]
                e_gm = E_gm4[ics + jb]
                if not bi_gm and not e_gm:
                    continue
                vol = dx8 * (dg[jb] // U0125)
                row_base = jb * stride1
                s0 = 0
                s1 = 0
                for kt in range(nt1):
                    off = row_base + kt * d1
                    w = cube[off]
                    if not w:
                        continue
                    ay = sum(cube[off + 1 : off + 1 + nA])
                    py = sum(cube[off + 1 + nA : off + 1 + d_ap])
                    az = sum(cube[off + 1 + d_ap : off + 1 + d_ap + nA])
                    pz = sum(cube[off + 1 + d_ap + nA : off + d1])
                    if bi_gm and ay and az:
                        s0 += (((((w * ay) // SCALE) * az) // SCALE) * bi_gm) // SCALE
                    if e_gm and py and pz:
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
        xq = ceil_div(gx, v_div)
        yq = ceil_div(gy + ceil_div(RHO_HI * gx, UNIT), v_div)
        zq = ceil_div(2 * gz + gy + ceil_div(2 * ROOT_HI * gx, UNIT), v_div)
        q = ceil_div(xq * xq + 9 * yq * yq + 3 * zq * zq, 3200)
        return (C3_NUM * exp_floor(q)) // UNIT

    sums_set = {grid[jy] + grid[kz] for jy in range(jy_max2 + 1) for kz in range(jy, jy_max2 + 1)}

    def get_b_vals(ix):
        r = ix * side
        gx = grid[ix]
        r2_lo = (2 * ROOT_LO * gx) // UNIT
        r2_hi = ceil_div(2 * ROOT_HI * gx, UNIT)
        out_c = {}
        for s in sums_set:
            loc = locate(r2_lo + s, r2_hi + s)
            out_c[s] = z_ap2 if loc is None else tuple(eval_loc(mt, r, loc) for mt in AP2)
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
                    b1, b2, b3, b4 = c00[idx_m], c01[idx_m], c10[idx_m], c11[idx_m]
                    if a1 and a2 and a3 and a4 and b1 and b2 and b3 and b4:
                        sa += isqrt_(isqrt_(isqrt_(a1 * a2) * isqrt_(a3 * a4)) * isqrt_(isqrt_(b1 * b2) * isqrt_(b3 * b4)))
                sp = 0
                for idx_m in range(nA0, nA0 + nP1):
                    a1, a2, a3, a4 = p00[idx_m], p01[idx_m], p10[idx_m], p11[idx_m]
                    b1, b2, b3, b4 = c00[idx_m], c01[idx_m], c10[idx_m], c11[idx_m]
                    if a1 and a2 and a3 and a4 and b1 and b2 and b3 and b4:
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
                    if not w:
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

    for n_val in (64, 56):
        t0 = time.time()
        c1, d1, c2, d2 = cavity(n_val, 16)
        dt = time.time() - t0
        p0, p1 = c1[0] / d1, c1[1] / d1
        p2, p3 = c2[0] / d2, c2[1] / d2
        print(
            f"[n={n_val}] {dt:.2f}s | c1={c1}/{d1} c2={c2}/{d2} | "
            f"tot={(p0+p1+p2+p3)*100:.5f}% ({p0*100:.5f}%+{p1*100:.5f}%+{p2*100:.5f}%+{p3*100:.5f}%)"
        )
