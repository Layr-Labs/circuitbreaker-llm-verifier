import operator
import time
from bisect import bisect_right
from fractions import Fraction
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
WEIGHT_CUT = 100
_mul = operator.mul


def lmlo(a, b, _isqrt=isqrt):
    if not a or not b:
        return 0
    if a == b:
        return a
    m = _isqrt(a * b)
    return (_isqrt(a * m) + _isqrt(b * m)) >> 1


def geo_interp(a, b, num, den, _isqrt=isqrt):
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


def cavity(n=72, nt1=18, nu_panels=20):
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
    lmlo_ = lmlo
    geo_ = geo_interp
    U0125 = UNIT // 8

    V_QUANT = 40_000
    v_div = UNIT_N // V_QUANT
    assert UNIT_N % V_QUANT == 0
    assert all(g % 2 == 0 and g % U0125 == 0 for g in grid)

    n_lut = 16 * n
    idx_lut = [0] * (n_lut + 1)
    cur_j = 0
    for q_idx in range(n_lut + 1):
        y_q = q_idx * U0125
        while cur_j + 1 < cside and grid[cur_j + 1] <= y_q:
            cur_j += 1
        idx_lut[q_idx] = cur_j

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
                            tot += lmlo_(prev, cur) * step_w
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
                            tot += lmlo_(prev, cur) * step_w
                        prev = cur
                    else:
                        prev = 0
                out[iside + j] = tot // UNIT_N
        return out

    # 2. Alternating crossing branches B_0 .. B_8
    branches = [[half_mass((ROOT_LO * grid[j]) // UNIT + grid[i] // 2) for i in points for j in points]]
    for _ in range(8):
        branches.append(apply_k_steps(branches[-1], K_steps))

    S0_tbl = build_S_steps(False, ix_max1 + 1)
    S1_tbl_full = build_S_steps(True, side)

    # Row-cached op_D_batch: since s0 + h = gi/2 - r_hi + h only depends on j and (gi/2 \pm h),
    # caching (s0+h, s1+h) per j avoids repeated locate + eval_batch across i!
    def op_D_batch(f_list, g_list, pairs, i_lim=side):
        n_out = len(pairs)
        outs = [[0] * (i_lim * side) for _ in range(n_out)]
        d0_val = density(0)
        same_fg = (f_list is g_list)
        for j in points:
            gj = grid[j]
            r_hi = -(-(ROOT_HI * gj) // UNIT)
            r_lo = (ROOT_LO * gj) // UNIT
            if grid[i_lim - 1] // 2 <= r_hi:
                continue
            row = j * side
            f_cache = {}
            g_cache = f_cache if same_fg else {}
            for i in range(i_lim):
                gi_half = grid[i] // 2
                s0 = gi_half - r_hi
                if s0 <= 0:
                    continue
                s1 = gi_half - r_lo
                upper = bound if bound < s0 else s0
                kmax = find_idx(upper)
                prev_h = 0
                fv0 = f_cache.get(gi_half)
                if fv0 is None:
                    loc0 = locate(s0, s1)
                    fv0 = eval_batch(f_list, row, loc0) if loc0 else None
                    f_cache[gi_half] = fv0
                if not same_fg:
                    gv0 = g_cache.get(gi_half)
                    if gv0 is None:
                        loc0 = locate(s0, s1)
                        gv0 = eval_batch(g_list, row, loc0) if loc0 else None
                        g_cache[gi_half] = gv0
                else:
                    gv0 = fv0
                if fv0 and gv0:
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
                    kp = gi_half + h
                    fv = f_cache.get(kp)
                    if kp not in f_cache:
                        lp = locate(s0 + h, s1 + h)
                        fv = eval_batch(f_list, row, lp) if lp else None
                        f_cache[kp] = fv
                    km = gi_half - h
                    gv = g_cache.get(km)
                    if km not in g_cache:
                        lm = locate(s0 - h, s1 - h)
                        gv = eval_batch(g_list, row, lm) if lm else None
                        g_cache[km] = gv
                    if fv and gv:
                        step_w = h - prev_h
                        for p_idx, (fi, gi_) in enumerate(pairs):
                            cur = ((dh * fv[fi]) // SCALE * gv[gi_]) // SCALE
                            pv = prev_vs[p_idx]
                            if pv and cur:
                                tots[p_idx] += lmlo_(pv, cur) * step_w
                            prev_vs[p_idx] = cur
                    else:
                        prev_vs = [0] * n_out
                    prev_h = h
                    if h == upper:
                        break
                idx_ij = i * side + j
                for p_idx in range(n_out):
                    outs[p_idx][idx_ij] = tots[p_idx] // UNIT_N
        return outs

    even = [branches[0], branches[2], branches[4], branches[6], branches[8]]
    odd = [branches[1], branches[3], branches[5], branches[7]]

    same_even = [apply_S_steps(mat, S0_tbl, ix_max1 + 1) for mat in even]
    pos_both = op_D_batch(odd[:3], even[:4], [(a, b) for a in range(3) for b in range(4)], ix_max1 + 1)

    W1_1 = apply_S_steps(branches[1], S1_tbl_full, side)
    W1_3 = apply_S_steps(branches[3], S1_tbl_full, side)
    W1_5 = apply_S_steps(branches[5], S1_tbl_full, side)
    W1_7 = apply_S_steps(branches[7], S1_tbl_full, side)
    W2_11, W2_13, W2_31, W2_33 = op_D_batch(odd[:2], odd[:2], [(0, 0), (0, 1), (1, 0), (1, 1)], side)

    # Precompute K_gm4_kerns and exact reachability masks
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
    for _ in range(10):
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
    VW2_steps = [None] * (v_i_lim * side)
    W_1_3_2 = (W1_1, W1_3, W2_11)
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
            c_p = -(-((-level + (UNIT_N // mesh)) * UNIT) // ROOT_LO) if level < (UNIT_N // mesh) else c
            if c_p <= c:
                u0 = (w0 * half_mass(level + (ROOT_LO * c) // UNIT)) // SCALE if w0 else 0
                us = [(kk, (w * half_mass(level + (ROOT_LO * grid[kk]) // UNIT)) // SCALE, step_w) for kk, w, step_w in steps]
                VU_steps[iside + j] = (k0, rem, dg0, u0, us)
            elif c_p < bound:
                kp = idx_lut[c_p // U0125]
                rem_p = c_p - grid[kp]
                dgp = dg[kp]
                wp = density(c_p + mh if c_p + mh > 0 else 0)
                up0 = (wp * half_mass(level + (ROOT_LO * c_p) // UNIT)) // SCALE if wp else 0
                us = []
                prev_cp = c_p
                for kk in range(kp + 1, side):
                    edge = grid[kk]
                    w = density(edge + mh)
                    if not w:
                        break
                    us.append((kk, (w * half_mass(level + (ROOT_LO * edge) // UNIT)) // SCALE, edge - prev_cp))
                    prev_cp = edge
                if us:
                    VU_steps[iside + j] = (kp, rem_p, dgp, up0, us)

            loc0 = locate(0 if 2 * ml + c < 0 else 2 * ml + c, 2 * mh + c)
            if w0 and loc0:
                g0_1, g0_3, g0_2 = eval_batch(W_1_3_2, row, loc0)
                g0_1 = (w0 * g0_1) // SCALE
                g0_3 = (w0 * g0_3) // SCALE
                g0_2 = (w0 * g0_2) // SCALE
            else:
                g0_1 = g0_3 = g0_2 = 0
            ws1, ws3, ws2 = [], [], []
            for kk, w, step_w in steps:
                edge = grid[kk]
                loc = locate(2 * ml + edge, 2 * mh + edge)
                if loc:
                    v1, v3, v2 = eval_batch(W_1_3_2, row, loc)
                    ws1.append((kk, (w * v1) // SCALE if v1 else 0, step_w))
                    ws3.append((kk, (w * v3) // SCALE if v3 else 0, step_w))
                    ws2.append((kk, (w * v2) // SCALE if v2 else 0, step_w))
                else:
                    ws1.append((kk, 0, step_w))
                    ws3.append((kk, 0, step_w))
                    ws2.append((kk, 0, step_w))
            VW1_steps[iside + j] = (k0, rem, dg0, g0_1, ws1)
            VW13_steps[iside + j] = (k0, rem, dg0, g0_3, ws3)
            VW2_steps[iside + j] = (k0, rem, dg0, g0_2, ws2)

    V_0U = apply_S_steps(branches[0], VU_steps, v_i_lim)
    V_2U = apply_S_steps(branches[2], VU_steps, v_i_lim)
    V_4U = apply_S_steps(branches[4], VU_steps, v_i_lim)
    V_6U = apply_S_steps(branches[6], VU_steps, v_i_lim)
    V_0W1 = apply_S_steps(branches[0], VW1_steps, v_i_lim)
    V_2W1 = apply_S_steps(branches[2], VW1_steps, v_i_lim)
    V_4W1 = apply_S_steps(branches[4], VW1_steps, v_i_lim)
    V_0W13 = apply_S_steps(branches[0], VW13_steps, v_i_lim)
    V_2W13 = apply_S_steps(branches[2], VW13_steps, v_i_lim)
    V_0W2 = apply_S_steps(branches[0], VW2_steps, v_i_lim)
    del VU_steps, VW1_steps, VW13_steps, VW2_steps

    # Guarded kernels for O_g and E_pt
    KW1_steps = [None] * (side * side)
    KW13_steps = [None] * (side * side)
    KW2_steps = [None] * ((ix_max2 + 1) * side)
    for j in points:
        row = j * side
        gj = grid[j]
        r2_lo = (2 * ROOT_LO * gj) // UNIT
        r2_hi = -(-(2 * ROOT_HI * gj) // UNIT)
        w_cache = {}
        for i in points:
            iside = i * side
            kern = kernels[iside + j]
            nk = len(kern)
            if nk < 2:
                continue
            gi = grid[i]
            do_w2 = (i <= ix_max2)
            v0 = w_cache.get(gi)
            if gi not in w_cache:
                loc0 = locate(r2_lo + gi, r2_hi + gi)
                v0 = eval_batch(W_1_3_2, row, loc0) if loc0 else None
                w_cache[gi] = v0
            if v0:
                w0_1 = (kern[0] * v0[0]) // SCALE
                w0_3 = (kern[0] * v0[1]) // SCALE
                w0_2 = (kern[0] * v0[2]) // SCALE if do_w2 else 0
            else:
                w0_1 = w0_3 = w0_2 = 0
            s1, s3 = [], []
            s2 = [] if do_w2 else None
            for k in range(1, nk):
                gik = gi + grid[k]
                vk = w_cache.get(gik)
                if gik not in w_cache:
                    loc = locate(r2_lo + gik, r2_hi + gik)
                    vk = eval_batch(W_1_3_2, row, loc) if loc else None
                    w_cache[gik] = vk
                if vk:
                    kk = kern[k]
                    dgk = dg[k - 1]
                    s1.append((k, (kk * vk[0]) // SCALE if vk[0] else 0, dgk))
                    s3.append((k, (kk * vk[1]) // SCALE if vk[1] else 0, dgk))
                    if do_w2:
                        s2.append((k, (kk * vk[2]) // SCALE if vk[2] else 0, dgk))
                else:
                    dgk = dg[k - 1]
                    s1.append((k, 0, dgk))
                    s3.append((k, 0, dgk))
                    if do_w2:
                        s2.append((k, 0, dgk))
            KW1_steps[iside + j] = (w0_1, s1)
            KW13_steps[iside + j] = (w0_3, s3)
            if do_w2:
                KW2_steps[iside + j] = (w0_2, s2)

    Og_mats = (
        [apply_k_steps(m_, KU_steps, ix_max2 + 1) for m_ in even]
        + [apply_k_steps(m_, KW1_steps, ix_max2 + 1) for m_ in even[:4]]
        + [apply_k_steps(m_, KW13_steps, ix_max2 + 1) for m_ in even[:2]]
        + [apply_k_steps(branches[0], KW2_steps, ix_max2 + 1)]
    )
    del KW2_steps, S0_tbl

    # 3. Holder-compressed GM4 table E_gm4 with corner caching and 5 guarded hops
    W_internal = (W1_1, W1_3, W1_5, W2_11)
    nW_int = len(W_internal)
    zero_wint = (0,) * nW_int

    @lru_cache(maxsize=None)
    def kg_corner_edges(i, j):
        a = grid[i]
        b = grid[j]
        ra = (ROOT_LO * a) // UNIT
        u_base = (7 * b) // 2 + ra
        c_l = (2 * ROOT_LO * b) // UNIT + a
        c_h = -(-(2 * ROOT_HI * b) // UNIT) + a
        row = j * side
        u_vals = [half_mass(u_base + (ROOT_LO * grid[r]) // UNIT) for r in points]
        u_edges = tuple(isqrt_(u_vals[r] * u_vals[r + 1]) for r in cpoints)
        if c_l >= UNIT_N:
            return u_edges, None
        w_vals = []
        for r in points:
            loc = locate(c_l + grid[r], c_h + grid[r])
            w_vals.append(eval_batch(W_internal, row, loc) if loc else zero_wint)
        w_edges = tuple(
            tuple(isqrt_(w_vals[r][m_] * w_vals[r + 1][m_]) for m_ in range(nW_int))
            for r in cpoints
        )
        return u_edges, w_edges

    Kg_gm4_kerns = [[] for _ in range(csize)]
    for i in cpoints:
        lim_j = need_i[i]
        if not lim_j:
            continue
        ics = i * cside
        for j in range(lim_j):
            klist = K_gm4_kerns[ics + j]
            if not klist:
                continue
            u00, w00 = kg_corner_edges(i, j)
            u01, w01 = kg_corner_edges(i, j + 1)
            u10, w10 = kg_corner_edges(i + 1, j)
            u11, w11 = kg_corner_edges(i + 1, j + 1)
            has_w = (w00 is not None and w01 is not None and w10 is not None and w11 is not None)
            row_kg = []
            for r, w_dg in klist:
                gtot = isqrt_(isqrt_(u00[r] * u01[r]) * isqrt_(u10[r] * u11[r]))
                if has_w:
                    wr00, wr01, wr10, wr11 = w00[r], w01[r], w10[r], w11[r]
                    for m_ in range(nW_int):
                        gtot += isqrt_(isqrt_(wr00[m_] * wr01[m_]) * isqrt_(wr10[m_] * wr11[m_]))
                if gtot:
                    row_kg.append((r, (w_dg * gtot) // SCALE))
            Kg_gm4_kerns[ics + j] = row_kg
    kg_corner_edges.cache_clear()

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

    h_cap_u = (13 * UNIT_N) // 10
    W_t = (W1_1, W1_3, W1_5, W1_7, W2_11, W2_13, W2_31, W2_33)
    V_t = (V_0U, V_2U, V_4U, V_6U, V_0W1, V_2W1, V_4W1, V_0W13, V_2W13, V_0W2)
    WV_t = W_t + V_t
    nWt = len(W_t)
    nWV = len(WV_t)
    zero_wv = (0,) * nWV
    nu = nu_panels

    @lru_cache(maxsize=None)
    def corner_time_edges(i, j):
        a = grid[i]
        b = grid[j]
        ah = a // 2
        M_l = (ROOT_LO * b) // UNIT + ah
        if M_l <= 0:
            return 0, 0, ()
        M_h = -(-(ROOT_HI * b) // UNIT) + ah
        U_c = h_cap_u if h_cap_u < M_l else M_l
        rb = (3 * b) // 2 + (INV_R2_LO * a) // UNIT
        h_u_max = ((rb - (UNIT_N // mesh)) * UNIT) // ROOT_HI
        U_u = U_c if U_c < h_u_max else (h_u_max if h_u_max > 0 else 0)
        row = j * side
        edges = []
        prev = None
        for kk in range(nu + 1):
            h_l = (kk * U_c) // nu
            h_h = -(-(kk * U_c) // nu)
            g_val = density(h_h)
            up_val = half_mass(rb + (ROOT_LO * h_l) // UNIT)
            lm = locate(M_l - h_h, M_h - h_l)
            bm = eval_batch(WV_t, row, lm) if lm else zero_wv
            if kk == 0:
                bp = bm
            elif M_h + h_h < bound:
                lp = locate(M_l + h_l, M_h + h_h)
                bp = eval_batch(WV_t, row, lp) if lp else zero_wv
            else:
                bp = zero_wv
            hu_l = (kk * U_u) // nu
            hu_h = -(-(kk * U_u) // nu)
            gu_val = density(hu_h) if U_u else 0
            um_val = half_mass(rb - (-(-(ROOT_HI * hu_h) // UNIT))) if U_u else 0
            upu_val = half_mass(rb + (ROOT_LO * hu_l) // UNIT) if U_u else 0
            if kk == 0:
                bpu = bm
            elif U_u and M_h + hu_h < bound:
                lpu = locate(M_l + hu_l, M_h + hu_h)
                bpu = eval_batch(WV_t, row, lpu) if lpu else zero_wv
            else:
                bpu = zero_wv
            cur = (g_val, up_val, bm, bp, gu_val, um_val, upu_val, bpu)
            if prev is not None:
                pg, pup, pbm, pbp, pgu, pum, pupu, pbpu = prev
                edges.append((
                    isqrt_(pg * g_val),
                    isqrt_(pup * up_val),
                    tuple(isqrt_(x * y) for x, y in zip(pbm, bm)),
                    tuple(isqrt_(x * y) for x, y in zip(pbp, bp)),
                    isqrt_(pgu * gu_val),
                    isqrt_(pum * um_val),
                    isqrt_(pupu * upu_val),
                    tuple(isqrt_(x * y) for x, y in zip(pbpu, bpu)),
                ))
            prev = cur
        return U_c, U_u, tuple(edges)

    T_gm4 = [0] * csize
    for i in cpoints:
        lim_j = need_i[i]
        if not lim_j:
            continue
        ics = i * cside
        for j in range(lim_j):
            U00, Uu00, e00 = corner_time_edges(i, j)
            U01, Uu01, e01 = corner_time_edges(i, j + 1)
            U10, Uu10, e10 = corner_time_edges(i + 1, j)
            U11, Uu11, e11 = corner_time_edges(i + 1, j + 1)
            if not (U00 and U01 and U10 and U11):
                continue
            U_gm = isqrt_(isqrt_(U00 * U01) * isqrt_(U10 * U11))
            Uu_gm = isqrt_(isqrt_(Uu00 * Uu01) * isqrt_(Uu10 * Uu11)) if (Uu00 and Uu01 and Uu10 and Uu11) else 0
            tot_c = 0
            tot_u = 0
            for (g00, up00, bm00, bp00, gu00, um00, upu00, bpu00), (g01, up01, bm01, bp01, gu01, um01, upu01, bpu01), (g10, up10, bm10, bp10, gu10, um10, upu10, bpu10), (g11, up11, bm11, bp11, gu11, um11, upu11, bpu11) in zip(e00, e01, e10, e11):
                gw = isqrt_(isqrt_(g00 * g01) * isqrt_(g10 * g11))
                if gw:
                    op = isqrt_(isqrt_(up00 * up01) * isqrt_(up10 * up11))
                    wm = 0
                    vm = 0
                    vp = 0
                    for m_ in range(nWt):
                        v = isqrt_(isqrt_(bm00[m_] * bm01[m_]) * isqrt_(bm10[m_] * bm11[m_]))
                        if v:
                            wm += v
                        v2 = isqrt_(isqrt_(bp00[m_] * bp01[m_]) * isqrt_(bp10[m_] * bp11[m_]))
                        if v2:
                            op += v2
                    for m_ in range(nWt, nWV):
                        v = isqrt_(isqrt_(bm00[m_] * bm01[m_]) * isqrt_(bm10[m_] * bm11[m_]))
                        if v:
                            vm += v
                        v2 = isqrt_(isqrt_(bp00[m_] * bp01[m_]) * isqrt_(bp10[m_] * bp11[m_]))
                        if v2:
                            vp += v2
                    tot_c += (gw * ((wm * (op + vp) + vm * op) // SCALE)) // SCALE
                if Uu_gm:
                    guw = isqrt_(isqrt_(gu00 * gu01) * isqrt_(gu10 * gu11))
                    if guw:
                        umu = isqrt_(isqrt_(um00 * um01) * isqrt_(um10 * um11))
                        if umu:
                            opu = isqrt_(isqrt_(upu00 * upu01) * isqrt_(upu10 * upu11))
                            for m_ in range(nWV):
                                v2 = isqrt_(isqrt_(bpu00[m_] * bpu01[m_]) * isqrt_(bpu10[m_] * bpu11[m_]))
                                if v2:
                                    opu += v2
                            tot_u += (guw * ((umu * opu) // SCALE)) // SCALE
            T_gm4[ics + j] = (U_gm * tot_c + Uu_gm * tot_u) // (nu * UNIT_N)
    corner_time_edges.cache_clear()

    E_gm4 = list(T_gm4)
    cur_gm4 = T_gm4
    for _ in range(5):
        mid_gm4 = apply_gm4_kern(cur_gm4, K_gm4_kerns)
        cur_gm4 = apply_gm4_kern(mid_gm4, Kg_gm4_kerns)
        E_gm4 = [a + b for a, b in zip(E_gm4, cur_gm4)]
    del K_gm4_kerns, Kg_gm4_kerns, T_gm4, cur_gm4, mid_gm4, kernels

    # 4. Point-level component tables for P_1 with row-cached T_mats
    m_mats = (W1_1, W1_3, W1_5, W2_11, V_0U, V_2U, V_4U, V_0W1, V_2W1, V_0W13, W2_13, W2_31, V_0W2)
    pair_specs = (
        (0, 0),            # 0: (U, U)
        (5, 0), (0, 5),    # 1,2: (V0U, U), (U, V0U)
        (1, 0), (0, 1),    # 3,4: (W1_1, U), (U, W1_1)
        (6, 0), (0, 6),    # 5,6: (V2U, U), (U, V2U)
        (2, 0), (0, 2),    # 7,8: (W1_3, U), (U, W1_3)
        (8, 0), (7, 0),    # 9,10: (V0W1, U), (V4U, U)
        (3, 0), (4, 0),    # 11,12: (W1_5, U), (W2_11, U)
        (9, 0), (10, 0),   # 13,14: (V2W1, U), (V0W13, U)
        (1, 5), (5, 1),    # 15,16: (W1_1, V0U), (V0U, W1_1)
        (1, 1),            # 17: (W1_1, W1_1)
        (2, 5), (5, 2),    # 18,19: (W1_3, V0U), (V0U, W1_3)
        (11, 0), (12, 0), (1, 6), (6, 1)
    )
    n_pairs = len(pair_specs)
    T_mats = [[0] * size for _ in range(n_pairs)]
    zero_m = (0,) * len(m_mats)
    for j in points:
        gj = grid[j]
        r_lo = (ROOT_LO * gj) // UNIT
        r_hi = -(-(ROOT_HI * gj) // UNIT)
        row = j * side
        r_base0 = (3 * gj) // 2
        m_cache = {}
        for i in points:
            if j >= t_max_j[i]:
                continue
            gi = grid[i]
            gi_half = gi // 2
            ml = r_lo + gi_half
            if ml <= 0:
                continue
            mh = r_hi + gi_half
            upper = h_cap_u if h_cap_u < ml else ml
            r_base = r_base0 + (INV_R2_LO * gi) // UNIT
            kmax = find_idx(upper)
            prevs = [0] * n_pairs
            tots = [0] * n_pairs
            prev_h = 0
            hs = [grid[kk] for kk in range(kmax + 1)]
            if hs[-1] < upper:
                hs.append(upper)
            h_u_max = ((r_base - (UNIT_N // mesh)) * UNIT) // ROOT_HI
            if 0 < h_u_max < upper:
                hs.append(h_u_max)
                hs.append(h_u_max >> 1)
                hs = sorted(set(hs))
            first = True
            for h in hs:
                dh = density(h)
                if not dh:
                    break
                km = gi_half - h
                em = m_cache.get(km)
                if km not in m_cache:
                    lm = locate(ml - h, mh - h)
                    em = eval_batch(m_mats, row, lm) if lm is not None else zero_m
                    m_cache[km] = em
                u_m = half_mass(r_base - (-(-(ROOT_HI * h) // UNIT)))
                u_p = half_mass(r_base + (ROOT_LO * h) // UNIT)
                dm = ((u_m,) + em)
                if h == 0:
                    dp = ((u_p,) + em)
                else:
                    kp = gi_half + h
                    ep = m_cache.get(kp)
                    if kp not in m_cache:
                        lp = locate(ml + h, mh + h)
                        ep = eval_batch(m_mats, row, lp) if lp is not None else zero_m
                        m_cache[kp] = ep
                    dp = ((u_p,) + ep)
                if first:
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_, b_ = dm[fi], dp[gi_]
                        prevs[t_idx] = (((dh * a_) // SCALE) * b_) // SCALE if (a_ and b_) else 0
                    first = False
                else:
                    step_w = h - prev_h
                    for t_idx, (fi, gi_) in enumerate(pair_specs):
                        a_, b_ = dm[fi], dp[gi_]
                        cur = (((dh * a_) // SCALE) * b_) // SCALE if (a_ and b_) else 0
                        pv = prevs[t_idx]
                        if pv and cur:
                            tots[t_idx] += lmlo_(pv, cur) * step_w
                        prevs[t_idx] = cur
                prev_h = h
            idx_ij = i * side + j
            for t_idx in range(n_pairs):
                T_mats[t_idx][idx_ij] = tots[t_idx] // UNIT_N

    K_T = [apply_k_steps(T_mats[k], K_steps) for k in (0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 15, 16, 17)]
    h1_U = [apply_k_steps(m_, KU_steps) for m_ in K_T]
    h1_W1 = [apply_k_steps(K_T[k], KW1_steps) for k in (0, 1, 3)]
    h1_W13 = [apply_k_steps(K_T[k], KW13_steps) for k in (0, 1)]
    K_h1 = [apply_k_steps(h1_U[k], K_steps) for k in (0, 1, 2, 3, 5)] + [apply_k_steps(h1_W1[0], K_steps)]
    h2_U = [apply_k_steps(m_, KU_steps) for m_ in K_h1]
    h2_W1 = [apply_k_steps(K_h1[0], KW1_steps)]
    K_h2 = [apply_k_steps(h2_U[k], K_steps) for k in (0, 1)]
    h3_U = [apply_k_steps(m_, KU_steps) for m_ in K_h2]
    E_pt = T_mats + h1_U + h1_W1 + h1_W13 + h2_U + h2_W1 + h3_U

    P1_S1 = [apply_S_steps(m_, S1_tbl_full, ix_max1 + 1) for m_ in E_pt]
    ee_both = op_D_batch(even[:3], even[:3], [(a, b) for a in range(3) for b in range(3)], ix_max1 + 1)
    e_sub1 = T_mats[:12] + h1_U[:4] + h2_U[:2]
    eo_list = [odd[0], odd[1]] + e_sub1
    eo_pairs = []
    for idx in range(len(e_sub1)):
        gi = 2 + idx
        eo_pairs.append((gi, 0))
        eo_pairs.append((0, gi))
    for idx in range(8):
        gi = 2 + (idx if idx < 6 else 12 + (idx - 6))
        eo_pairs.append((gi, 1))
        eo_pairs.append((1, gi))
    eo_both = op_D_batch(eo_list, eo_list, eo_pairs, ix_max1 + 1)
    del T_mats, K_T, h1_U, h1_W1, h1_W13, K_h1, h2_U, h2_W1, K_h2, h3_U, E_pt, K_steps, KU_steps, KW1_steps, KW13_steps, S1_tbl_full

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

    # 5. Kind 1 (crossing degree 1: pos1, neg1)
    # Split into S-type tables (active for both y and z) and D-type tables (active ONLY for z when ROOT_LO * z_lo < gx // 2)!
    nA_s = len(same_even)
    nP_s = len(P1_S1)
    d_ap_s = nA_s + nP_s
    AP_s = tuple(same_even + P1_S1)
    zero_ap_s = (0,) * d_ap_s

    nA_d = len(pos_both)
    nP_d = len(ee_both) + len(eo_both)
    d_ap_d = nA_d + nP_d
    AP_d = tuple(pos_both + ee_both + eo_both)
    zero_ap_d = (0,) * d_ap_d

    def f2_opp(gx, gb):
        xq = -(-gx // v_div)
        bq = -(-gb // v_div)
        num = 9 * (xq * xq + bq * bq) + (-(-(12 * ROOT_HI * xq * bq) // UNIT))
        return (C2_NUM * exp_floor(-(-num // 3200))) // UNIT

    counts = [0, 0, 0, 0]
    old_sq = None
    for ix in range(ix_max1 + 1):
        gx = grid[ix]
        gx_half = gx // 2
        row_off = ix * side
        aps_cache = {}
        apd_cache = {}
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
                yv = aps_cache.get(yk)
                if yv is None:
                    loc = locate(ylo, yhi)
                    yv = zero_ap_s if loc is None else eval_batch(AP_s, row_off, loc)
                    aps_cache[yk] = yv
                zk = (zlo, zhi)
                zv = aps_cache.get(zk)
                if zv is None:
                    loc = locate(zlo, zhi)
                    zv = zero_ap_s if loc is None else eval_batch(AP_s, row_off, loc)
                    aps_cache[zk] = zv
                if (ROOT_HI * zlo) // UNIT < gx_half:
                    zvd = apd_cache.get(zk)
                    if zvd is None:
                        loc = locate(zlo, zhi)
                        zvd = None if loc is None else eval_batch(AP_d, row_off, loc)
                        if zvd is not None and not any(zvd):
                            zvd = None
                        apd_cache[zk] = zvd
                else:
                    zvd = None
                row_kt.append((w_all, yv, zv, zvd))
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
                w0, y0, z0, zd0 = p0
                w1, y1, z1, zd1 = p1
                wl = isqrt_(w0 * w1)
                if not wl:
                    lrow.append(None)
                    continue
                zdl = list(map(isqrt_, map(_mul, zd0, zd1))) if (zd0 is not None and zd1 is not None) else None
                lrow.append((wl, list(map(isqrt_, map(_mul, y0, y1))), list(map(isqrt_, map(_mul, z0, z1))), zdl))
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
                w0, y0, z0, zd0 = p0
                w1, y1, z1, zd1 = p1
                ws = isqrt_(w0 * w1)
                if not ws:
                    srow.append(None)
                    continue
                zds = list(map(isqrt_, map(_mul, zd0, zd1))) if (zd0 is not None and zd1 is not None) else None
                srow.append((ws, list(map(isqrt_, map(_mul, y0, y1))), list(map(isqrt_, map(_mul, z0, z1))), zds))
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
                    w0, y0, z0, zd0 = c0
                    w1, y1, z1, zd1 = c1
                    w = isqrt_(w0 * w1)
                    if not w:
                        continue
                    has_zd = (zd0 is not None and zd1 is not None)
                    if bi_gm:
                        ay = sum(map(isqrt_, map(_mul, y0[:nA_s], y1[:nA_s])))
                        if ay:
                            az = sum(map(isqrt_, map(_mul, z0[:nA_s], z1[:nA_s])))
                            if has_zd:
                                az += sum(map(isqrt_, map(_mul, zd0[:nA_d], zd1[:nA_d])))
                            if az:
                                s0 += (((((w * ay) // SCALE) * az) // SCALE) * bi_gm) // SCALE
                    if e_gm:
                        py = sum(map(isqrt_, map(_mul, y0[nA_s:], y1[nA_s:])))
                        if py:
                            pz = sum(map(isqrt_, map(_mul, z0[nA_s:], z1[nA_s:])))
                            if has_zd:
                                pz += sum(map(isqrt_, map(_mul, zd0[nA_d:], zd1[nA_d:])))
                            if pz:
                                s1 += (((((w * py) // SCALE) * pz) // SCALE) * e_gm) // SCALE
                counts[0] += 3 * vol * s0
                counts[1] += 3 * vol * s1
        old_sq = squares

    # 6. Kind 2 (crossing degree 2: pos2, neg2) in (x, y, z)
    # In Kind 2, b = 2*sqrt(2)*x + y + z ALWAYS has m = sqrt(2)*b - x/2 >= 0, so D-type tables are identically 0!
    nA0 = len(same_even)
    nP1 = len(P1_S1)
    AP2 = tuple(same_even + P1_S1)
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
    U = Fraction(163309458029, 1250000000000)
    lam = Fraction(35, 71)
    R = Fraction(80089, 80080)
    t0 = time.time()
    c1, d1, c2, d2 = cavity(72, 18, 20)
    dt = time.time() - t0
    p0, p1 = Fraction(c1[0], d1), Fraction(c1[1], d1)
    p2, p3 = Fraction(c2[0], d2), Fraction(c2[1], d2)
    p = p0 + p1 + p2 + p3
    a0 = R * (U - (1 - lam) * p)
    coef = 1 + Fraction(1, 2) / a0
    print(f"n=72 nt1=18 nu=20 | {dt:5.2f}s | p={float(p)*100:.6f}% ({float(p0)*100:.5f}+{float(p1)*100:.5f}+{float(p2)*100:.5f}+{float(p3)*100:.5f}) | coef={float(coef):.7f}")
