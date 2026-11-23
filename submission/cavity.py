"""Lower bounds for the four local switching-gain integrals."""
from functools import lru_cache
from math import isqrt

SCALE = 1_000_000_000
UNIT = 1_000_000
ROOT_LO = 1_414_213
ROOT_HI = 1_414_214
RHO_LO = 942_809
RHO_HI = 942_810
INVROOT_LO = 707_106
QUARTER_ROOT_LO = 353_553
C_G_LO = 2_763_953
C_ROOT_LO = 659_845
HALF_CDF_MESH = 2048


def ceil_div(a, b):
    return -(-a // b)


def lmlo(a, b):
    if not a or not b:
        return 0
    if a == b:
        return a
    m = isqrt(a * b)
    return (isqrt(a * m) + isqrt(b * m)) // 2


def geomean(a, b, t):
    if not t or a == b:
        return a
    if t == UNIT:
        return b
    if not a or not b:
        return 0
    if b < a:
        a, b, t = b, a, UNIT - t
    if (b << 2) < 5 * a:
        w = (256 * t) // UNIT
        return (256 * a * b) // ((256 - w) * b + w * a)
    k = (64 * t) // UNIT
    if not k:
        return a
    v = a
    v = isqrt(v * (b if (k & 1) else a))
    v = isqrt(v * (b if (k & 2) else a))
    v = isqrt(v * (b if (k & 4) else a))
    v = isqrt(v * (b if (k & 8) else a))
    v = isqrt(v * (b if (k & 16) else a))
    v = isqrt(v * (b if (k & 32) else a))
    return v


@lru_cache(maxsize=None)
def exp_floor(i):
    if i > 4_000_000:
        return 0
    h = 6_400_000
    u = ceil_div(i * SCALE, h)
    l = (i * SCALE) // h
    v = SCALE - u + (l * i) // (2 * h) - ceil_div(ceil_div(u * i, 2 * h) * i, 3 * h)
    if v <= 0:
        return 0
    for _ in range(6):
        v = (v * v) // SCALE
    return v


def product(*values):
    acc = SCALE
    for v in values:
        acc = (acc * v) // SCALE
    return acc


def cavity(n=48):
    m = 2 * n
    side = m + 1
    size = side * side
    bound = m * UNIT
    mesh = HALF_CDF_MESH
    q_step = UNIT // 4
    h_step = UNIT // 2
    denom_un = UNIT * n

    @lru_cache(maxsize=None)
    def density(v):
        t = ceil_div(abs(v), n * 10)
        return (C_G_LO * exp_floor(ceil_div(6 * t * t, 100_000))) // 1_000_000

    mass = [0]
    for i in range(1, 2 * mesh + 1):
        e = exp_floor(ceil_div(200_000 * (3 * i * i - 3 * i + 1), mesh * mesh))
        mass.append(mass[-1] + (C_G_LO * e) // (1_000_000 * mesh))
    mass_end = mass[-1]

    def half_mass(r):
        if r <= 0:
            return 0
        idx = (mesh * r) // (n * UNIT)
        return mass[2 * mesh if idx > 2 * mesh else idx]

    def tail_mass(r):
        if r <= 0:
            return mass_end
        idx = ceil_div(mesh * r, n * UNIT)
        return 0 if idx >= 2 * mesh else mass_end - mass[idx]

    def lower_y(mat, row_off, yl, yh):
        if yl < 0 or yh >= bound:
            return 0
        j0, h0 = divmod(yl, UNIT)
        v0 = mat[row_off + j0]
        if not v0:
            return 0
        j1, h1 = divmod(yh, UNIT)
        if j0 == j1:
            if not h1:
                return v0
            v1 = mat[row_off + j0 + 1]
            if v1 >= v0:
                return geomean(v0, v1, h0) if h0 else v0
            return geomean(v0, v1, h1)
        a = geomean(v0, mat[row_off + j0 + 1], h0) if h0 else v0
        if not a:
            return 0
        v10 = mat[row_off + j1]
        if not v10:
            return 0
        b = geomean(v10, mat[row_off + j1 + 1], h1) if h1 else v10
        return a if a < b else b

    def lower(mat, xl, xh, yl, yh):
        if xl < 0 or yl < 0 or xh >= bound or yh >= bound:
            return 0
        i0, f0 = divmod(xl, UNIT)
        r0 = i0 * side
        a0 = lower_y(mat, r0, yl, yh)
        if not a0:
            return 0
        i1, f1 = divmod(xh, UNIT)
        if i0 == i1:
            if not f1:
                return a0
            b0 = lower_y(mat, r0 + side, yl, yh)
            return (geomean(a0, b0, f0) if f0 else a0) if b0 >= a0 else geomean(a0, b0, f1)
        b0 = lower_y(mat, r0 + side, yl, yh)
        v0 = geomean(a0, b0, f0) if f0 else a0
        if not v0:
            return 0
        r1 = i1 * side
        a1 = lower_y(mat, r1, yl, yh)
        if not a1:
            return 0
        v1 = geomean(a1, lower_y(mat, r1 + side, yl, yh), f1) if f1 else a1
        return v0 if v0 < v1 else v1

    points = range(side)
    d_q_grid = [density(k * q_step) for k in range(8 * n + 1)]

    kernels = []
    for i in points:
        half_i = i * h_step
        for j in points:
            base = ROOT_HI * j + half_i
            kr = []
            for k in points:
                val = density(base + k * UNIT)
                if not val:
                    break
                kr.append(val)
            kernels.append(kr)

    # 1. Branches B_0 .. B_6
    B0_mat = [half_mass(ROOT_LO * j + i * h_step) for i in points for j in points]
    B1_mat = []
    for i in points:
        half_i = i * h_step
        for j in points:
            kr = kernels[i * side + j]
            if not kr:
                B1_mat.append(0)
                continue
            row = j * side
            base_k = ROOT_HI * j + half_i
            half_j = j * h_step
            n_q = min(2, len(kr) - 1)
            tot4 = 0
            prev_q = (kr[0] * B0_mat[row]) // SCALE
            for q in range(1, 4 * n_q + 1):
                dq = q * q_step
                cur_q = (density(base_k + dq) * half_mass((ROOT_LO * dq) // UNIT + half_j)) // SCALE
                tot4 += lmlo(prev_q, cur_q)
                prev_q = cur_q
            tot = tot4 // 4
            prev = (kr[n_q] * B0_mat[row + n_q]) // SCALE
            for k in range(n_q + 1, len(kr)):
                cur = (kr[k] * B0_mat[row + k]) // SCALE
                tot += lmlo(prev, cur)
                prev = cur
            B1_mat.append(tot // n)

    branches = [B0_mat, B1_mat]
    for _ in range(5):
        old = branches[-1]
        mat = []
        for i in points:
            for j in points:
                kr = kernels[i * side + j]
                if not kr:
                    mat.append(0)
                    continue
                row = j * side
                prev = (kr[0] * old[row]) // SCALE
                tot = 0
                for k in range(1, len(kr)):
                    cur = (kr[k] * old[row + k]) // SCALE
                    tot += lmlo(prev, cur)
                    prev = cur
                mat.append(tot // n)
        branches.append(mat)

    even = branches[::2]  # B0, B2, B4, B6
    odd = branches[1::2]  # B1, B3, B5

    def op_same_b0():
        out = [0] * size
        for i in points:
            half_i = i * h_step
            for j in points:
                ml = ROOT_LO * j - half_i
                mh = ROOT_HI * j - half_i
                c = max(0, -ml)
                if c >= bound:
                    continue
                row = j * side
                half_j = j * h_step
                bval = half_mass((ROOT_LO * c) // UNIT + half_j)
                prev = (bval * density(max(c + mh, 0))) // SCALE
                tot = 0
                k_start = c // UNIT + 1
                if k_start < side:
                    span0 = k_start * UNIT - c
                    qs = span0 // 4
                    if qs > 0:
                        for q in range(1, 5):
                            dq = c + q * qs if q < 4 else k_start * UNIT
                            cur_q = (density(max(dq + mh, 0)) * half_mass((ROOT_LO * dq) // UNIT + half_j)) // SCALE
                            tot += lmlo(prev, cur_q) * (dq - c)
                            prev, c = cur_q, dq
                        k_start += 1
                for k in range(k_start, side):
                    edge = k * UNIT
                    w = density(edge + mh)
                    if not w:
                        break
                    cur = (w * B0_mat[row + k]) // SCALE
                    tot += lmlo(prev, cur) * (edge - c)
                    prev, c = cur, edge
                out[i * side + j] = tot // denom_un
        return out

    def op_same_batch(mats, strict):
        n_m = len(mats)
        outs = [[0] * size for _ in range(n_m)]
        for j in points:
            row = j * side
            r_lo = ROOT_LO * j
            r_hi = ROOT_HI * j
            act = [(mi, mm) for mi, mm in enumerate(mats) if any(mm[row : row + side])]
            if not act:
                continue
            for i in points:
                half_i = i * h_step
                ml = r_lo - half_i
                mh = r_hi - half_i
                c0 = max(0, -(2 if strict else 1) * ml)
                if c0 >= bound:
                    continue
                k0, rem = divmod(c0, UNIT)
                w0 = density(max(c0 + mh, 0))
                if not w0:
                    continue
                if rem:
                    prev = [(geomean(mm[row + k0], mm[row + k0 + 1], rem) * w0) // SCALE if mm[row + k0] else 0 for _, mm in act]
                else:
                    prev = [(mm[row + k0] * w0) // SCALE for _, mm in act]
                tots = [0] * len(act)
                c = c0
                for k in range(k0 + 1, side):
                    edge = k * UNIT
                    w = density(edge + mh)
                    if not w:
                        break
                    dt = edge - c
                    rk = row + k
                    for ai, (_, mm) in enumerate(act):
                        v = mm[rk]
                        if not v:
                            prev[ai] = 0
                            continue
                        cur = (w * v) // SCALE
                        p_ai = prev[ai]
                        if p_ai and cur:
                            tots[ai] += lmlo(p_ai, cur) * dt
                        prev[ai] = cur
                    c = edge
                cell = i * side + j
                for (mi, _), tot in zip(act, tots):
                    outs[mi][cell] = tot // denom_un
        return outs

    def op_both_batch(mats_list, pair_list):
        n_p = len(pair_list)
        n_m = len(mats_list)
        outs = [[0] * size for _ in range(n_p)]
        j_lim = (n * UNIT - 1) // ROOT_HI
        zeros_m = [0] * n_m
        for j in range(j_lim + 1):
            row = j * side
            half_j = j * h_step
            r_hi = ROOT_HI * j
            r_lo = ROOT_LO * j
            idx_min = ceil_div(r_hi, h_step)
            b0_tab = [0] * (4 * n + 1)
            v_tab = [zeros_m] * (4 * n + 1)
            for idx in range(idx_min, 4 * n + 1):
                yl = idx * h_step - r_hi
                yh = idx * h_step - r_lo
                if yh >= bound:
                    break
                b0_tab[idx] = half_mass((ROOT_LO * yl) // UNIT + half_j)
                v_tab[idx] = [lower_y(mm, row, yl, yh) for mm in mats_list]
            b0_zero = half_mass(half_j)
            v_zero = [lower_y(mm, row, 0, r_hi - r_lo) for mm in mats_list]

            i_start = (2 * r_hi) // UNIT + 1
            for i in range(i_start, side):
                s0 = i * h_step - r_hi
                s1 = i * h_step - r_lo
                prev = [0] * n_p
                tots = [0] * n_p
                left = 0
                lim = ceil_div(s0, h_step)
                for k in range(lim):
                    edge = k * h_step
                    dh = d_q_grid[2 * k]
                    if not dh:
                        break
                    ip = i + k
                    im = i - k
                    b0_p = b0_tab[ip] if ip <= 4 * n else half_mass((ROOT_LO * (s0 + edge)) // UNIT + half_j)
                    b0_m = b0_tab[im]
                    vp = v_tab[ip] if ip <= 4 * n else zeros_m
                    vm = v_tab[im]
                    dt = edge - left
                    for pi, (fi, gi) in enumerate(pair_list):
                        a = b0_p if fi == -1 else vp[fi]
                        if not a:
                            prev[pi] = 0
                            continue
                        b = b0_m if gi == -1 else vm[gi]
                        if not b:
                            prev[pi] = 0
                            continue
                        cur = ((dh * a) // SCALE * b) // SCALE
                        p_pi = prev[pi]
                        if p_pi and cur:
                            tots[pi] += lmlo(p_pi, cur) * dt
                        prev[pi] = cur
                    left = edge
                else:
                    dh = density(s0)
                    if dh:
                        b0_p = half_mass((ROOT_LO * (2 * s0)) // UNIT + half_j)
                        vp = [lower_y(mm, row, 2 * s0, s0 + s1) for mm in mats_list] if (s0 + s1) < bound else zeros_m
                        dt = s0 - left
                        for pi, (fi, gi) in enumerate(pair_list):
                            a = b0_p if fi == -1 else vp[fi]
                            if not a:
                                continue
                            b = b0_zero if gi == -1 else v_zero[gi]
                            if not b:
                                continue
                            cur = ((dh * a) // SCALE * b) // SCALE
                            p_pi = prev[pi]
                            if p_pi and cur:
                                tots[pi] += lmlo(p_pi, cur) * dt
                cell = i * side + j
                for mat, tot in zip(outs, tots):
                    mat[cell] = tot // denom_un
        return outs

    same_even = [op_same_b0()] + op_same_batch(even[1:], False)

    br_rest = branches[1:]
    br_pairs = [
        (0, -1), (0, 1), (0, 3), (2, -1), (2, 1), (2, 3), (4, -1), (4, 1),
        (0, 0), (0, 2), (2, 0),
        (-1, -1), (-1, 1), (1, -1), (1, 1), (-1, 3), (3, -1), (1, 3), (3, 1), (3, 3),
    ]
    br_both_outs = op_both_batch(br_rest, br_pairs)
    pos_both = br_both_outs[:8]
    w2_11, w2_13, w2_31 = br_both_outs[8:11]
    ee_both = br_both_outs[11:]
    positive = same_even + pos_both

    # Batched closed, private, rare1, rare3, wv_1, wv_3, vu_0, vu_2:
    closed = [0] * size
    private = [0] * size
    rare1 = [0] * size
    rare3 = [0] * size
    wv_1 = [0] * size
    wv_3 = [0] * size
    vu_0 = [0] * size
    vu_2 = [0] * size
    B2_mat = branches[2]
    B3_mat = branches[3]

    for i in points:
        inv_i = INVROOT_LO * i
        half_i = i * h_step
        for j in points:
            cell = i * side + j
            # 1) closed:
            a_val = (3 * UNIT * j) // 2 + inv_i
            h_max = min(bound, (a_val * UNIT) // ROOT_HI)
            if h_max > 0:
                dk_c = 1 if h_max <= 12 * UNIT else (2 if h_max <= 24 * UNIT else 4)
                step_c = dk_c * q_step
                k_lim_c = h_max // step_c
                ha0 = half_mass(a_val)
                prev_c = ((d_q_grid[0] * ha0) // SCALE * ha0) // SCALE
                tot_c = 0
                left_c = 0
                for ks in range(1, k_lim_c + 1):
                    k = ks * dk_c
                    edge = k * q_step
                    dh = d_q_grid[k]
                    if not dh:
                        prev_c = 0
                        break
                    cur_c = ((dh * half_mass(a_val - ceil_div(ROOT_HI * edge, UNIT))) // SCALE
                             * half_mass(a_val + (ROOT_LO * edge) // UNIT)) // SCALE
                    if prev_c and cur_c:
                        tot_c += lmlo(prev_c, cur_c)
                    prev_c, left_c = cur_c, edge
                tot_c *= step_c
                if left_c < h_max and prev_c:
                    dh = density(h_max)
                    if dh:
                        cur_c = ((dh * half_mass(a_val - ceil_div(ROOT_HI * h_max, UNIT))) // SCALE
                                 * half_mass(a_val + (ROOT_LO * h_max) // UNIT)) // SCALE
                        if cur_c:
                            tot_c += lmlo(prev_c, cur_c) * (h_max - left_c)
                closed[cell] = tot_c // denom_un

            # 2) Merged private, rare1, rare3, wv_1, wv_3, vu_0, vu_2:
            ml = ROOT_LO * j - half_i
            mh = ROOT_HI * j - half_i
            c_base = max(0, -2 * ml)
            if c_base >= bound:
                continue
            w0 = density(c_base + mh)
            if not w0:
                continue
            row = j * side
            half_j = j * h_step
            lev_base = (7 * UNIT * j) // 2 - ROOT_HI * i
            rad_base = (7 * UNIT * j) // 2 - ROOT_LO * i
            c_p = max(c_base, ceil_div(-lev_base * UNIT, ROOT_LO)) if lev_base < 0 else c_base
            k0, rem = divmod(c_base, UNIT)
            b0_0 = half_mass((ROOT_LO * c_base) // UNIT + half_j)
            b1_0 = geomean(B1_mat[row + k0], B1_mat[row + k0 + 1], rem) if rem else B1_mat[row + k0]
            b2_0 = geomean(B2_mat[row + k0], B2_mat[row + k0 + 1], rem) if rem else B2_mat[row + k0]
            b3_0 = geomean(B3_mat[row + k0], B3_mat[row + k0 + 1], rem) if rem else B3_mat[row + k0]
            hm0 = half_mass(lev_base + (ROOT_LO * c_base) // UNIT)
            tm0 = tail_mass(rad_base + ceil_div(ROOT_HI * c_base, UNIT))
            p_pr = (w0 * hm0) // SCALE
            p_r1 = ((w0 * b1_0) // SCALE * tm0) // SCALE
            p_r3 = ((w0 * b3_0) // SCALE * tm0) // SCALE
            p_wv1 = (p_pr * b1_0) // SCALE
            p_wv3 = (p_pr * b3_0) // SCALE
            p_vu0 = (p_pr * b0_0) // SCALE
            p_vu2 = (p_pr * b2_0) // SCALE
            t_pr = t_r1 = t_r3 = t_wv1 = t_wv3 = t_vu0 = t_vu2 = 0
            c_cur = c_base
            kp, rem_p = divmod(c_p, UNIT)
            for k in range(k0 + 1, side):
                edge = k * UNIT
                if k == kp + 1 and c_p < edge:
                    if c_p > c_cur:
                        wp = density(c_p + mh)
                        if wp:
                            dt = c_p - c_cur
                            b0_p = half_mass((ROOT_LO * c_p) // UNIT + half_j)
                            b1_p = geomean(B1_mat[row + kp], B1_mat[row + kp + 1], rem_p)
                            b2_p = geomean(B2_mat[row + kp], B2_mat[row + kp + 1], rem_p)
                            b3_p = geomean(B3_mat[row + kp], B3_mat[row + kp + 1], rem_p)
                            hmp = half_mass(lev_base + (ROOT_LO * c_p) // UNIT)
                            tmp = tail_mass(rad_base + ceil_div(ROOT_HI * c_p, UNIT))
                            c_pr = (wp * hmp) // SCALE
                            c_r1 = ((wp * b1_p) // SCALE * tmp) // SCALE
                            c_r3 = ((wp * b3_p) // SCALE * tmp) // SCALE
                            if p_r1 and c_r1: t_r1 += lmlo(p_r1, c_r1) * dt
                            if p_r3 and c_r3: t_r3 += lmlo(p_r3, c_r3) * dt
                            p_pr, p_r1, p_r3 = c_pr, c_r1, c_r3
                            p_wv1 = (c_pr * b1_p) // SCALE
                            p_wv3 = (c_pr * b3_p) // SCALE
                            p_vu0 = (c_pr * b0_p) // SCALE
                            p_vu2 = (c_pr * b2_p) // SCALE
                            c_cur = c_p
                    c_mid = (c_cur + edge) // 2
                    if c_cur < c_mid < edge:
                        wm = density(c_mid + mh)
                        if wm:
                            c_pr_m = (wm * half_mass(lev_base + (ROOT_LO * c_mid) // UNIT)) // SCALE
                            b0_m = half_mass((ROOT_LO * c_mid) // UNIT + half_j)
                            c_vu0_m = (c_pr_m * b0_m) // SCALE
                            if p_pr and c_pr_m: t_pr += lmlo(p_pr, c_pr_m) * (c_mid - c_cur)
                            if p_vu0 and c_vu0_m: t_vu0 += lmlo(p_vu0, c_vu0_m) * (c_mid - c_cur)
                            p_pr, p_vu0, c_p_last = c_pr_m, c_vu0_m, c_mid
                        else:
                            c_p_last = c_cur
                    else:
                        c_p_last = c_cur
                else:
                    c_p_last = c_cur
                w = density(edge + mh)
                if not w:
                    break
                dt = edge - c_cur
                dt_p = edge - c_p_last
                rk = row + k
                hm = half_mass(lev_base + ROOT_LO * k)
                tm = tail_mass(rad_base + ROOT_HI * k)
                b0_k, b1_k, b2_k, b3_k = B0_mat[rk], B1_mat[rk], B2_mat[rk], B3_mat[rk]
                c_pr = (w * hm) // SCALE
                c_r1 = ((w * b1_k) // SCALE * tm) // SCALE
                c_r3 = ((w * b3_k) // SCALE * tm) // SCALE
                c_wv1 = (c_pr * b1_k) // SCALE
                c_wv3 = (c_pr * b3_k) // SCALE
                c_vu0 = (c_pr * b0_k) // SCALE
                c_vu2 = (c_pr * b2_k) // SCALE
                if p_pr and c_pr: t_pr += lmlo(p_pr, c_pr) * dt_p
                if p_r1 and c_r1: t_r1 += lmlo(p_r1, c_r1) * dt
                if p_r3 and c_r3: t_r3 += lmlo(p_r3, c_r3) * dt
                if p_wv1 and c_wv1: t_wv1 += lmlo(p_wv1, c_wv1) * dt
                if p_wv3 and c_wv3: t_wv3 += lmlo(p_wv3, c_wv3) * dt
                if p_vu0 and c_vu0: t_vu0 += lmlo(p_vu0, c_vu0) * dt_p
                if p_vu2 and c_vu2: t_vu2 += lmlo(p_vu2, c_vu2) * dt
                p_pr, p_r1, p_r3, p_wv1, p_wv3, p_vu0, p_vu2 = c_pr, c_r1, c_r3, c_wv1, c_wv3, c_vu0, c_vu2
                c_cur = edge
            private[cell] = t_pr // denom_un
            rare1[cell] = t_r1 // denom_un
            rare3[cell] = t_r3 // denom_un
            wv_1[cell] = t_wv1 // denom_un
            wv_3[cell] = t_wv3 // denom_un
            vu_0[cell] = t_vu0 // denom_un
            vu_2[cell] = t_vu2 // denom_un

    weak1, weak3 = op_same_batch([B1_mat, B3_mat], True)

    def op_V_W_batch():
        outs = [[0] * size for _ in range(4)]
        for j in points:
            row = j * side
            half_j = j * h_step
            r2_lo = 2 * ROOT_LO * j
            r2_hi = 2 * ROOT_HI * j
            w1_d = [0] * (2 * m + 1)
            w3_d = [0] * (2 * m + 1)
            w2_d = [0] * (2 * m + 1)
            for d in range(-m, m + 1):
                sl = r2_lo + d * UNIT
                sh = r2_hi + d * UNIT
                if sh < 0 or sl >= bound or sh >= bound:
                    continue
                sl_pos = max(0, sl)
                idx_d = d + m
                w1_d[idx_d] = lower_y(weak1, row, sl_pos, sh)
                w3_d[idx_d] = lower_y(weak3, row, sl_pos, sh)
                w2_d[idx_d] = lower_y(w2_11, row, sl_pos, sh)
            sh_zero = r2_hi - r2_lo
            w1_z = lower_y(weak1, row, 0, sh_zero)
            w3_z = lower_y(weak3, row, 0, sh_zero)
            w2_z = lower_y(w2_11, row, 0, sh_zero)

            for i in points:
                half_i = i * h_step
                ml = ROOT_LO * j - half_i
                mh = ROOT_HI * j - half_i
                c = max(0, -2 * ml)
                if c >= bound:
                    continue
                if c > 0:
                    w1_0, w3_0, w2_0 = w1_z, w3_z, w2_z
                else:
                    idx_0 = m - i
                    w1_0, w3_0, w2_0 = w1_d[idx_0], w3_d[idx_0], w2_d[idx_0]
                b0_0 = half_mass((ROOT_LO * c) // UNIT + half_j)
                k0, rem = divmod(c, UNIT)
                b2_0 = geomean(B2_mat[row + k0], B2_mat[row + k0 + 1], rem) if rem else B2_mat[row + k0]
                d0 = density(c + mh)
                db0 = (d0 * b0_0) // SCALE
                db2 = (d0 * b2_0) // SCALE
                prev0 = (db0 * w1_0) // SCALE
                prev1 = (db2 * w1_0) // SCALE
                prev2 = (db0 * w3_0) // SCALE
                prev3 = (db0 * w2_0) // SCALE
                t0_v = t1_v = t2_v = t3_v = 0
                c_cur = c
                base_d = m - i
                for k in range(k0 + 1, side):
                    edge = k * UNIT
                    w = density(edge + mh)
                    if not w:
                        break
                    if r2_lo + (k - i) * UNIT >= bound:
                        break
                    idx_d = base_d + k
                    w1_k, w3_k, w2_k = w1_d[idx_d], w3_d[idx_d], w2_d[idx_d]
                    db0_k = (w * B0_mat[row + k]) // SCALE
                    db2_k = (w * B2_mat[row + k]) // SCALE
                    c0 = (db0_k * w1_k) // SCALE
                    c1 = (db2_k * w1_k) // SCALE
                    c2 = (db0_k * w3_k) // SCALE
                    c3 = (db0_k * w2_k) // SCALE
                    dt = edge - c_cur
                    if prev0 and c0: t0_v += lmlo(prev0, c0) * dt
                    if prev1 and c1: t1_v += lmlo(prev1, c1) * dt
                    if prev2 and c2: t2_v += lmlo(prev2, c2) * dt
                    if prev3 and c3: t3_v += lmlo(prev3, c3) * dt
                    prev0, prev1, prev2, prev3 = c0, c1, c2, c3
                    c_cur = edge
                cell = i * side + j
                outs[0][cell] = t0_v // denom_un
                outs[1][cell] = t1_v // denom_un
                outs[2][cell] = t2_v // denom_un
                outs[3][cell] = t3_v // denom_un
        return outs

    vw_01, vw_21, vw_03, vw_0_w2 = op_V_W_batch()

    # Row-precomputed tails (21 non-U-U endpoint pairs):
    all_t_mats = [private, weak1, rare1, rare3, weak3, w2_11, w2_13, w2_31,
                  vw_01, vw_21, vw_03, vw_0_w2, vu_0, vu_2, wv_1, wv_3]
    n_tm = len(all_t_mats)
    pair_idx = [
        (0, -1), (2, -1), (-1, 0), (8, -1), (1, 0), (3, -1), (12, 14),
        (9, -1), (10, -1), (5, -1), (4, 0), (12, 15), (13, 14),
        (6, -1), (7, -1), (0, 2), (8, 1), (-1, 2), (2, 2), (5, 0), (11, -1),
    ]
    n_pairs = len(pair_idx)
    tails = [[0] * size for _ in range(n_pairs)]
    zeros_tm = [0] * n_tm

    for j in points:
        row = j * side
        r_lo = ROOT_LO * j
        r_hi = ROOT_HI * j
        s_min = -(r_lo // q_step)
        s_max = (bound - 1 - r_hi) // q_step
        s_len = s_max - s_min + 1
        if s_len <= 0:
            continue
        s_odd_lim = (24 * UNIT - r_lo) // q_step
        s_tab = [
            [lower_y(mm, row, r_lo + s * q_step, r_hi + s * q_step) for mm in all_t_mats]
            if (not (s & 1) or s <= s_odd_lim) else zeros_tm
            for s in range(s_min, s_max + 1)
        ]
        vm_at_ml = [lower_y(mm, row, 0, r_hi - r_lo) for mm in all_t_mats]
        r_u_base = (3 * UNIT * j) // 2

        for i in points:
            half_i = 2 * i * q_step
            ml = r_lo + half_i
            mh = r_hi + half_i
            upper = min(bound, ml)
            if upper <= 0:
                continue
            r_u = r_u_base + INVROOT_LO * i
            dk = 1 if upper <= 12 * UNIT else (2 if upper <= 24 * UNIT else 4)
            step_t = dk * q_step
            k_upper = upper // step_t
            prev = [0] * n_pairs
            tots = [0] * n_pairs
            dh0 = d_q_grid[0]
            u0 = half_mass(r_u)
            s0_idx = 2 * i - s_min
            v0_m = s_tab[s0_idx] if 0 <= s0_idx < s_len else zeros_tm
            for ti, (fi, gi) in enumerate(pair_idx):
                a = u0 if fi == -1 else v0_m[fi]
                if not a:
                    continue
                b = u0 if gi == -1 else v0_m[gi]
                if not b:
                    continue
                prev[ti] = ((dh0 * a) // SCALE * b) // SCALE

            left = 0
            for ks in range(1, k_upper + 1):
                k = ks * dk
                edge = k * q_step
                dh = d_q_grid[k]
                if not dh:
                    prev = [0] * n_pairs
                    break
                up_u = half_mass(r_u + (ROOT_LO * edge) // UNIT)
                dn_u = half_mass(r_u - ceil_div(ROOT_HI * edge, UNIT))
                sm = 2 * i - k - s_min
                sp = 2 * i + k - s_min
                vm = s_tab[sm] if 0 <= sm < s_len else zeros_tm
                vp = s_tab[sp] if 0 <= sp < s_len else zeros_tm
                for ti, (fi, gi) in enumerate(pair_idx):
                    a = dn_u if fi == -1 else vm[fi]
                    if not a:
                        prev[ti] = 0
                        continue
                    b = up_u if gi == -1 else vp[gi]
                    if not b:
                        prev[ti] = 0
                        continue
                    cur = ((dh * a) // SCALE * b) // SCALE
                    p_ti = prev[ti]
                    if p_ti and cur:
                        tots[ti] += lmlo(p_ti, cur)
                    prev[ti] = cur
                left = edge
            tots = [tot * step_t for tot in tots]
            if left < upper:
                dh = density(upper)
                if dh:
                    up_u = half_mass(r_u + (ROOT_LO * upper) // UNIT)
                    dn_u = half_mass(r_u - ceil_div(ROOT_HI * upper, UNIT))
                    vm = vm_at_ml
                    vp = [lower_y(mm, row, ml + upper, mh + upper) for mm in all_t_mats] if (ml + upper) < bound else zeros_tm
                    dt = upper - left
                    for ti, (fi, gi) in enumerate(pair_idx):
                        a = dn_u if fi == -1 else vm[fi]
                        if not a:
                            continue
                        b = up_u if gi == -1 else vp[gi]
                        if not b:
                            continue
                        cur = ((dh * a) // SCALE * b) // SCALE
                        p_ti = prev[ti]
                        if p_ti and cur:
                            tots[ti] += lmlo(p_ti, cur) * dt
            idx_cell = i * side + j
            for mat, tot in zip(tails, tots):
                mat[idx_cell] = tot // denom_un

    def continue_arm_batch(mats, guard_mode=0):
        n_m = len(mats)
        outs = [[0] * size for _ in range(n_m)]
        for j in points:
            row = j * side
            act = [(mi, mm) for mi, mm in enumerate(mats) if any(mm[row : row + side])]
            if not act:
                continue
            jU = j * UNIT
            if guard_mode == 1:
                g_tab = [half_mass((7 * jU) // 2 + ROOT_LO * s) for s in range(2 * m + 1)]
            elif guard_mode == 2:
                g_tab = [lower_y(weak1, row, 2 * ROOT_LO * j + s * UNIT, 2 * ROOT_HI * j + s * UNIT) for s in range(2 * m + 1)]
            elif guard_mode == 3:
                g_tab = [lower_y(weak3, row, 2 * ROOT_LO * j + s * UNIT, 2 * ROOT_HI * j + s * UNIT) for s in range(2 * m + 1)]
            else:
                g_tab = None

            n_a = len(act)
            for i in points:
                kr = kernels[i * side + j]
                if not kr:
                    continue
                prev = [0] * n_a
                tots = [0] * n_a
                for k, w in enumerate(kr):
                    if guard_mode:
                        g_fac = g_tab[i + k]
                        if not g_fac:
                            prev = [0] * n_a
                            continue
                        wg = (w * g_fac) // SCALE
                        if not wg:
                            prev = [0] * n_a
                            continue
                    else:
                        wg = w
                    rk = row + k
                    for ai, (_, mm) in enumerate(act):
                        v = mm[rk]
                        if not v:
                            prev[ai] = 0
                            continue
                        val = (wg * v) // SCALE
                        p_ai = prev[ai]
                        if p_ai and val:
                            tots[ai] += lmlo(p_ai, val)
                        prev[ai] = val
                cell = i * side + j
                for (mi, _), tot in zip(act, tots):
                    outs[mi][cell] = tot // n
        return outs

    T_all = [closed] + tails

    protected_odd = (continue_arm_batch(even, 1)
                     + continue_arm_batch(even[:3], 2)
                     + continue_arm_batch(even[:2], 3))

    step1 = continue_arm_batch(T_all[:12], 0)
    two_U = continue_arm_batch(step1, 1)
    two_W1 = continue_arm_batch(step1[:6], 2)
    two_W3 = continue_arm_batch(step1[:2], 3)
    step3 = continue_arm_batch(two_U[:6], 0)
    four_U = continue_arm_batch(step3, 1)
    four_W1 = continue_arm_batch(step3[:2], 2)
    step5 = continue_arm_batch(four_U[:2], 0)
    six_U = continue_arm_batch(step5, 1)

    guarded_3d = T_all + two_U + two_W1 + two_W3 + four_U + four_W1 + six_U
    E_p = T_all[:16] + two_U[:8] + two_W1[:4] + two_W3[:2] + four_U[:4] + four_W1[:2] + six_U
    E_b1 = T_all[:12] + two_U[:4] + two_W1[:2] + four_U[:2]
    E_b3 = T_all[:6] + two_U[:2]

    protected = op_same_batch(E_p, True)
    eo_mats = [B1_mat, B3_mat] + E_b1
    eo_pairs = []
    for idx in range(len(E_b1)):
        gi = 2 + idx
        eo_pairs.append((gi, 0))
        eo_pairs.append((0, gi))
    for idx in range(len(E_b3)):
        gi = 2 + (idx if idx < 6 else 12 + (idx - 6))
        eo_pairs.append((gi, 1))
        eo_pairs.append((1, gi))
    eo_both = op_both_batch(eo_mats, eo_pairs)
    mixed_3d = protected + ee_both + eo_both
    og_3d = protected_odd

    # Kind 1 on [0, 1]^3
    b1_mats = even[1:] + guarded_3d
    b1_o = len(even) - 1

    y1_full = same_even + protected
    y1_fo = len(same_even)
    dim_yf = len(y1_full)
    w_yf = (n + 1) * dim_yf
    s_yf = n * dim_yf

    y1_both = pos_both + ee_both + eo_both
    y1_bo = len(pos_both)
    dim_yb = len(y1_both)
    w_yb = (n + 1) * dim_yb
    s_yb = n * dim_yb
    j_both = n // 6 + 2
    n10 = n * 10

    counts = [0, 0, 0, 0]
    old_w_sq = old_b1_sq = old_yf_sq = old_yb_sq = None

    for i in range(n + 1):
        i0_s = [0] * (2 * n + 1)
        f0_s = [0] * (2 * n + 1)
        i1_s = [0] * (2 * n + 1)
        f1_s = [0] * (2 * n + 1)
        xh_ok_s = [False] * (2 * n + 1)
        jac_s = [0] * (2 * n + 1)
        rh_s = [0] * (2 * n + 1)
        rl_s = [0] * (2 * n + 1)
        th_s = [0] * (2 * n + 1)
        tl_s = [0] * (2 * n + 1)
        xb2_s = [0] * (2 * n + 1)
        b1_pts = []
        for s in range(2 * n + 1):
            yz = s * UNIT
            num_x = QUARTER_ROOT_LO * i * s
            xl = num_x // n
            xh = ceil_div(num_x, n)
            i0_s[s], f0_s[s] = divmod(xl, UNIT)
            i1_s[s], f1_s[s] = divmod(xh, UNIT)
            xh_ok_s[s] = (xh < bound)
            jac_s[s] = (SCALE * QUARTER_ROOT_LO * s) // (UNIT * n)
            rh_s[s] = ceil_div(RHO_HI * xh, UNIT)
            rl_s[s] = (RHO_LO * xl) // UNIT
            th = ceil_div(2 * ROOT_HI * xh, UNIT)
            tl = (2 * ROOT_LO * xl) // UNIT
            th_s[s] = th
            tl_s[s] = tl
            xb = ceil_div(xh, n10)
            xb2_s[s] = xb * xb
            bl = yz - th
            bh = yz - tl
            b0_val = half_mass((ROOT_LO * bl) // UNIT + xl // 2)
            b1_pts.append([b0_val] + [lower(mm, xl, xh, bl, bh) for mm in b1_mats])
        b1_edge = [[isqrt(a * b) if a and b else 0 for a, b in zip(b1_pts[s], b1_pts[s + 1])] for s in range(2 * n)]
        b1_sq = [[isqrt(a * b) if a and b else 0 for a, b in zip(b1_edge[s], b1_edge[s + 1])] for s in range(2 * n - 1)]

        w_plane = [0] * ((n + 1) * (n + 1))
        yf_plane = [0] * ((n + 1) * w_yf)
        yb_plane = [0] * ((j_both + 1) * w_yb)
        for j in range(n + 1):
            y = j * UNIT
            row_w = j * (n + 1)
            row_yf = j * w_yf
            in_b = (j <= j_both)
            row_yb = j * w_yb
            for k in range(n + 1):
                s = j + k
                zz2_y = (2 * k + j) * UNIT
                yy = ceil_div(max(abs(y - rh_s[s]), abs(y - rl_s[s])), n10)
                zb = ceil_div(max(abs(zz2_y - th_s[s]), abs(zz2_y - tl_s[s])), n10)
                exp_arg = ceil_div(xb2_s[s] + 9 * yy * yy + 3 * zb * zb, 200_000)
                w = ((C_ROOT_LO * exp_floor(exp_arg)) // 1_000_000 * jac_s[s]) // SCALE
                w_plane[row_w + k] = w
                if w and xh_ok_s[s]:
                    pos_f = row_yf + k * dim_yf
                    i0, f0 = i0_s[s], f0_s[s]
                    i1, f1 = i1_s[s], f1_s[s]
                    c0 = i0 * side + j
                    c0_s = c0 + side
                    if i0 == i1:
                        if not f1:
                            yf_plane[pos_f : pos_f + dim_yf] = [mm[c0] for mm in y1_full]
                            if in_b:
                                pos_b = row_yb + k * dim_yb
                                yb_plane[pos_b : pos_b + dim_yb] = [mm[c0] for mm in y1_both]
                        else:
                            yf_plane[pos_f : pos_f + dim_yf] = [
                                (geomean(mm[c0], mm[c0_s], f0) if mm[c0_s] >= mm[c0] else geomean(mm[c0], mm[c0_s], f1))
                                if mm[c0] else 0
                                for mm in y1_full
                            ]
                            if in_b:
                                pos_b = row_yb + k * dim_yb
                                yb_plane[pos_b : pos_b + dim_yb] = [
                                    (geomean(mm[c0], mm[c0_s], f0) if mm[c0_s] >= mm[c0] else geomean(mm[c0], mm[c0_s], f1))
                                    if mm[c0] else 0
                                    for mm in y1_both
                                ]
                    else:
                        c1 = i1 * side + j
                        c1_s = c1 + side
                        yf_plane[pos_f : pos_f + dim_yf] = [
                            min(geomean(mm[c0], mm[c0_s], f0) if f0 else mm[c0],
                                geomean(mm[c1], mm[c1_s], f1) if f1 else mm[c1])
                            if (mm[c0] and mm[c1]) else 0
                            for mm in y1_full
                        ]
                        if in_b:
                            pos_b = row_yb + k * dim_yb
                            yb_plane[pos_b : pos_b + dim_yb] = [
                                min(geomean(mm[c0], mm[c0_s], f0) if f0 else mm[c0],
                                    geomean(mm[c1], mm[c1_s], f1) if f1 else mm[c1])
                                if (mm[c0] and mm[c1]) else 0
                                for mm in y1_both
                            ]

        w_lines = []
        for j in range(n + 1):
            st = j * (n + 1)
            w_lines.extend([isqrt(a * b) if a and b else 0 for a, b in zip(w_plane[st : st + n], w_plane[st + 1 : st + n + 1])])
        w_sq = [isqrt(a * b) if a and b else 0 for a, b in zip(w_lines, w_lines[n:])]

        yf_lines = []
        for j in range(n + 1):
            st = j * w_yf
            yf_lines.extend([isqrt(a * b) if a and b else 0 for a, b in zip(yf_plane[st : st + s_yf], yf_plane[st + dim_yf : st + w_yf])])
        yf_sq = [isqrt(a * b) if a and b else 0 for a, b in zip(yf_lines, yf_lines[s_yf:])]

        yb_lines = []
        for j in range(j_both + 1):
            st = j * w_yb
            yb_lines.extend([isqrt(a * b) if a and b else 0 for a, b in zip(yb_plane[st : st + s_yb], yb_plane[st + dim_yb : st + w_yb])])
        yb_sq = [isqrt(a * b) if a and b else 0 for a, b in zip(yb_lines, yb_lines[s_yb:])]

        if i:
            w_cube = [isqrt(a * b) if a and b else 0 for a, b in zip(old_w_sq, w_sq)]
            ev_b = [0] * (2 * n - 1)
            gr_b = [0] * (2 * n - 1)
            for s in range(2 * n - 1):
                bc = [isqrt(a * b) if a and b else 0 for a, b in zip(old_b1_sq[s], b1_sq[s])]
                ev_b[s] = sum(bc[: 1 + b1_o])
                gr_b[s] = sum(bc[1 + b1_o :])

            s_pos = [0] * (n * n)
            s_mix = [0] * (n * n)
            b_lim = j_both * n
            for idx in range(n * n):
                if not w_cube[idx]:
                    continue
                base_f = idx * dim_yf
                ycf = [isqrt(a * b) if a and b else 0 for a, b in zip(old_yf_sq[base_f : base_f + dim_yf], yf_sq[base_f : base_f + dim_yf])]
                sp = sum(ycf[:y1_fo])
                sm = sum(ycf[y1_fo:])
                if idx < b_lim:
                    base_b = idx * dim_yb
                    ycb = [isqrt(a * b) if a and b else 0 for a, b in zip(old_yb_sq[base_b : base_b + dim_yb], yb_sq[base_b : base_b + dim_yb])]
                    sp += sum(ycb[:y1_bo])
                    sm += sum(ycb[y1_bo:])
                s_pos[idx] = sp
                s_mix[idx] = sm

            for j in range(n):
                jn = j * n
                for k in range(n):
                    idx = jn + k
                    w = w_cube[idx]
                    if not w:
                        continue
                    t_idx = k * n + j
                    s = j + k
                    counts[0] += 3 * product(w, s_pos[idx], s_pos[t_idx], ev_b[s])
                    counts[1] += 3 * product(w, s_mix[idx], s_mix[t_idx], gr_b[s])
        old_w_sq, old_b1_sq, old_yf_sq, old_yb_sq = w_sq, b1_sq, yf_sq, yb_sq

    # Kind 2 on [0, 1/2] x [0, 1]^2
    nx2 = n // 2
    y2_mats = even + odd + guarded_3d + og_3d
    yo1 = len(even)
    yo2 = yo1 + len(odd)
    yo3 = yo2 + len(guarded_3d)
    yo4 = yo3 + len(og_3d)
    b2_mats = positive + mixed_3d
    bo1 = len(positive)
    bo2 = bo1 + len(mixed_3d)

    old_w_sq = old_y_edge = old_b_sq = None
    for i in range(nx2 + 1):
        x = i * UNIT
        row_i = i * side
        w_plane = [0] * ((n + 1) * (n + 1))
        xb = ceil_div(x, n10)
        xb2 = xb * xb
        rh = ceil_div(RHO_HI * x, UNIT)
        rl = (RHO_LO * x) // UNIT
        th = ceil_div(2 * ROOT_HI * x, UNIT)
        tl = (2 * ROOT_LO * x) // UNIT
        for j in range(n + 1):
            y = -j * UNIT
            yy = ceil_div(max(abs(y - rh), abs(y - rl)), n10)
            yy9 = xb2 + 9 * yy * yy
            row_off = j * (n + 1)
            for k in range(n + 1):
                zz2_y = -(2 * k + j) * UNIT
                zb = ceil_div(max(abs(zz2_y - th), abs(zz2_y - tl)), n10)
                exp_arg = ceil_div(yy9 + 3 * zb * zb, 200_000)
                w_plane[row_off + k] = (C_ROOT_LO * exp_floor(exp_arg)) // 1_000_000
        w_lines = []
        for j in range(n + 1):
            st = j * (n + 1)
            w_lines.extend([isqrt(a * b) if a and b else 0 for a, b in zip(w_plane[st : st + n], w_plane[st + 1 : st + n + 1])])
        w_sq = [isqrt(a * b) if a and b else 0 for a, b in zip(w_lines, w_lines[n:])]

        y_pts = [[mm[row_i + j] for mm in y2_mats] for j in range(n + 1)]
        y_edge = [[isqrt(a * b) if a and b else 0 for a, b in zip(y_pts[j], y_pts[j + 1])] for j in range(n)]

        bl0 = 2 * ROOT_LO * i
        bh0 = 2 * ROOT_HI * i
        b_pts = [[lower_y(mm, row_i, bl0 + s * UNIT, bh0 + s * UNIT) for mm in b2_mats] for s in range(2 * n + 1)]
        b_edge = [[isqrt(a * b) if a and b else 0 for a, b in zip(b_pts[s], b_pts[s + 1])] for s in range(2 * n)]
        b_sq = [[isqrt(a * b) if a and b else 0 for a, b in zip(b_edge[s], b_edge[s + 1])] for s in range(2 * n - 1)]

        if i:
            w_cube = [isqrt(a * b) if a and b else 0 for a, b in zip(old_w_sq, w_sq)]
            ev_y = [0] * n; od_y = [0] * n; gr_y = [0] * n; og_y = [0] * n
            for j in range(n):
                yc = [isqrt(a * b) if a and b else 0 for a, b in zip(old_y_edge[j], y_edge[j])]
                ev_y[j] = sum(yc[:yo1]); od_y[j] = sum(yc[yo1:yo2])
                gr_y[j] = sum(yc[yo2:yo3]); og_y[j] = sum(yc[yo3:yo4])
            pos_b = [0] * (2 * n - 1); mix_b = [0] * (2 * n - 1)
            for s in range(2 * n - 1):
                bc = [isqrt(a * b) if a and b else 0 for a, b in zip(old_b_sq[s], b_sq[s])]
                pos_b[s] = sum(bc[:bo1]); mix_b[s] = sum(bc[bo1:bo2])
            for j in range(n):
                jn = j * n
                ay, oy, sy, wy = ev_y[j], od_y[j], gr_y[j], og_y[j]
                for k in range(n):
                    w = w_cube[jn + k]
                    if not w:
                        continue
                    s = j + k
                    az, sz, wz = ev_y[k], gr_y[k], og_y[k]
                    counts[2] += 3 * product(w, pos_b[s], ay + product(oy, az))
                    counts[3] += 3 * product(w, mix_b[s], product(ay, az) + product(sy, wz) + product(wy, sz))
        old_w_sq, old_y_edge, old_b_sq = w_sq, y_edge, b_sq

    return counts, SCALE * (n ** 3)
