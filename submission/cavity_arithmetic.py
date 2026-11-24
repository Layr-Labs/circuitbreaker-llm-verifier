"""Directed four-root integration, adapted from xlib's PR 895 at cec7eb07dec8.

The point bank uses selected pairs. The compressed bank is complete.
"""
import operator
from array import array
from bisect import bisect_left, bisect_right
from functools import lru_cache
from math import isqrt, sumprod

from root_integration import arm_ranges, either_arm, companion_mean, mixed_arm_floors, row_lookup, square_cell_mean

from root_integration import (
    SCALE, UNIT, ROOT_LO, ROOT_HI, INV_R2_LO, CG_NUM, C2_NUM, C3_NUM,
    _edge_means, HALF_CDF_MESH, EXP_DEN, EXP_CUT, WEIGHT_CUT, build_graded_grid,
    curvature_floor, cubic_density_lookup, degree_two_weight_planes, exp_floor, packed_scaling, geo_interp, lmlo,
)

_mul = operator.mul


def cavity(n=96, nt1=18, nu_panels=18, progress=lambda name: None, *, positive_only=False, negative_only=False,
           centers_only=False, roots_only=False, endpoint_table=None, compute_endpoints=False):
    assert not (compute_endpoints and endpoint_table is not None)
    assert (not compute_endpoints and endpoint_table is None) or negative_only
    assert not centers_only or (negative_only and not roots_only)
    assert not (positive_only and negative_only)
    grid = build_graded_grid(n)
    side = len(grid)
    cside = side - 1
    size = side * side
    csize = cside * cside
    UNIT_N = UNIT * n
    bound = grid[-1]
    dg = [grid[k + 1] - grid[k] for k in range(cside)]
    mesh = HALF_CDF_MESH
    isqrt_ = isqrt
    lmlo_ = lmlo
    geo_ = geo_interp
    U0125 = UNIT // 8
    space_squared = UNIT_N**2
    h_cap_u = (13 * UNIT_N) // 10

    @lru_cache(maxsize=8192)
    def gaussian_span(width):
        return curvature_floor(width, width*width, space_squared)

    V_QUANT = 40_000
    v_div = UNIT_N // V_QUANT
    assert UNIT_N % V_QUANT == 0
    assert all(g % 2 == 0 and g % U0125 == 0 for g in grid)

    n_lut = bound // U0125
    idx_lut = [0] * (n_lut + 1)
    cur_j = 0
    for q_idx in range(n_lut + 1):
        y_q = q_idx * U0125
        while cur_j + 1 < cside and grid[cur_j + 1] <= y_q:
            cur_j += 1
        idx_lut[q_idx] = cur_j
        assert cur_j == min(cside - 1, bisect_right(grid, y_q) - 1)

    density = cubic_density_lookup(v_div)

    mass = [0]
    cumulative_numerator = 0
    m2 = mesh * mesh
    for i in range(1, 2 * mesh + 1):
        q = -(-(2 * EXP_DEN * (3 * i * i - 3 * i + 1)) // m2)
        cumulative_numerator += CG_NUM * exp_floor(q)
        mass.append(cumulative_numerator // (UNIT * mesh))
    max_mass_idx = 2 * mesh

    @lru_cache(maxsize=16384)
    def half_mass(r):
        if r <= 0:
            return 0
        if r >= 2 * UNIT_N:
            return mass[max_mass_idx]
        idx, remainder = divmod(mesh * r, UNIT_N)
        return mass[idx] + (mass[idx + 1] - mass[idx]) * remainder // UNIT_N

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

    ix_max1 = find_idx(((11 * UNIT_N) // 20) if positive_only else (UNIT_N // 2)) + 1
    center_axis = find_idx(68*UNIT_N//100)+1
    center_cells = find_idx(110*UNIT_N//100)+1
    single_limit = find_idx(95*UNIT_N//100)+2
    jb_max1 = find_idx((3 * UNIT_N) // 4) + 1
    ix_max2 = find_idx(UNIT_N // 2) + 1
    jy_max2 = find_idx(UNIT_N) + 1
    points = range(side)
    cpoints = range(cside)

    def quarter_nodes(left, right):
        if not left or not right:
            return (0, 0)
        midpoint = isqrt_(left * right)
        return isqrt_(left * midpoint), isqrt_(right * midpoint)

    density_nodes = lru_cache(maxsize=16384)(quarter_nodes)

    def prepare_panel_entry(entry, side_integral):
        if entry is None:
            return None
        if side_integral:
            k0, offset, width, previous, steps = entry
        else:
            previous, steps = entry
            offset = 0
        partial = None
        groups = []
        start = end = span = 0
        weights = []
        for k, value, length in steps:
            if previous and value:
                low, high = density_nodes(previous, value)
                panel_span = gaussian_span(length)
                if offset and k == k0 + 1:
                    partial = (k0, offset, width, panel_span*low, panel_span*high)
                else:
                    index = k - 1
                    if index != end or panel_span != span:
                        if weights:
                            groups.append((start,end,span,tuple(weights)))
                        start = index
                        span = panel_span
                        weights = []
                    weights.extend((low, high))
                    end = k
            previous = value
        if weights:
            groups.append((start,end,span,tuple(weights)))
        return (groups, partial) if groups or partial else None

    kernels = []
    K_steps = []
    K_panels = []
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
                if not positive_only:
                    uval = (val * half_mass(base_u + (ROOT_LO * gk) // UNIT)) // SCALE
                    if k == 0:
                        ku0 = uval
                    else:
                        kus.append((k, uval, dg[k - 1]))
            plain = (krow[0], ks) if len(krow) >= 2 else None
            if positive_only:
                K_panels.append(prepare_panel_entry(plain, side_integral=False))
            else:
                kernels.append(krow)
                K_steps.append(plain)
                KU_steps.append((ku0, kus) if len(krow) >= 2 else None)

    def quarter_panels(raw, *, side_integral=False):
        return [prepare_panel_entry(entry, side_integral) for entry in raw]

    def prepare_source(source):
        if not any(source):
            return None
        samples = [quarter_nodes(a, b) for a, b in zip(source, source[1:])]
        left, right = zip(*samples)
        first = next((i for i, (a, b) in enumerate(samples) if a or b), cside)
        stop = next((i + 1 for i in range(cside-1, first-1, -1) if left[i] or right[i]), first)
        return source, tuple(v for pair in samples for v in pair), first, stop

    panel_denominator = 2 * SCALE * UNIT_N

    def integrate_panel_entry(entry, source_samples):
        if entry is None or source_samples is None:
            return 0
        source, flat, first, stop = source_samples
        panels, partial = entry
        total = 0
        for start, end, span, weights in panels:
            if end > first and start < stop:
                total += span * sumprod(weights,flat[2*start:2*end])
        if partial:
            k, offset, width, wl, wr = partial
            start = geo_(source[k], source[k+1], offset, width)
            a, b = quarter_nodes(start, source[k+1])
            total += wl*a + wr*b
        return total // panel_denominator

    def apply_panels(mat, prepared, i_lim=side):
        out = [0] * (i_lim * side)
        for j in points:
            source_samples = prepare_source(mat[j*side:(j+1)*side])
            if source_samples is None:
                continue
            for ij in range(j, i_lim*side, side):
                out[ij] = integrate_panel_entry(prepared[ij], source_samples)
        return array("I", out)

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

    if not positive_only:
        K_panels = quarter_panels(K_steps)
        KU_panels = quarter_panels(KU_steps)
    del KU_steps
    branches = [[half_mass((ROOT_LO * grid[j]) // UNIT + grid[i] // 2) for i in points for j in points]]
    for _ in range(8):
        branches.append(apply_panels(branches[-1], K_panels))

    if positive_only:
        del K_panels
    S0_panels = None if negative_only else quarter_panels(
        build_S_steps(False, ix_max1 + 1), side_integral=True)
    if not positive_only:
        S1_tbl_full = build_S_steps(True, side)
        S1_panels = quarter_panels(S1_tbl_full, side_integral=True)

    def op_D_batch(f_list, g_list, pairs, i_lim=side):
        n_out = len(pairs)
        outs = [array("I",[0]) * (i_lim * side) for _ in range(n_out)]
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
                    weighted_f = tuple(d0_val*v//SCALE for v in fv0)
                    prev_vs = [weighted_f[fi]*gv0[gi_]//SCALE for fi, gi_ in pairs]
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
                        step_w = gaussian_span(h - prev_h)
                        weighted_f = tuple(dh*v//SCALE for v in fv)
                        for p_idx, (fi, gi_) in enumerate(pairs):
                            cur = weighted_f[fi]*gv[gi_]//SCALE
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

    if negative_only:
        same_even = []
        pos_both = []
    else:
        same_even = [apply_panels(mat, S0_panels, ix_max1 + 1) for mat in even]
        pos_both = op_D_batch(odd, even, [(a, b) for a in range(len(odd)) for b in range(len(even))], ix_max1 + 1)

    if positive_only:
        P1_S1 = []
        ee_both = []
        eo_both = []
        Og_mats = []
        E_gm4 = [0] * csize
        progress("positive-only branch families")
    else:
        W1_1 = apply_panels(branches[1], S1_panels, side)
        W1_3 = apply_panels(branches[3], S1_panels, side)
        W1_5 = apply_panels(branches[5], S1_panels, side)
        W1_7 = apply_panels(branches[7], S1_panels, side)
        W2_11, W2_13, W2_31, W2_33 = op_D_batch(odd[:2], odd[:2], [(0, 0), (0, 1), (1, 0), (1, 1)], side)

        def pack_kernel(row):
            if not row:
                return None
            index, weight = zip(*row)
            return array("I", index), array("Q", weight)

        # Precompute K_gm4_kerns and exact reachability masks
        K_gm4_kerns = []
        need_i = [0] * cside
        if endpoint_table is None:
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
                            row_k.append((r, w * gaussian_span(dg[r])))
                    K_gm4_kerns.append(pack_kernel(row_k))

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
                            rm = kl[0][-1] + 1
                            if rm > need_i[j]:
                                need_i[j] = rm

        t_max_j = [0] * side
        for i in range(max(center_axis,center_cells) + 1):
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

        v_i_lim = side if endpoint_table is not None else max(max(need_i) + 1, max(t_max_j))
        if v_i_lim > side:
            v_i_lim = side

        v_source_samples = [
            tuple(prepare_source(branches[k][j*side:(j+1)*side]) for k in (0, 2, 4, 6))
            for j in points
        ]
        v_tables = [
            [array("I", [0]) * (v_i_lim * side) for _ in range(count)]
            for count in (4, 3, 2, 1)
        ]
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
                u_entry = None
                if c_p <= c:
                    u0 = (w0 * half_mass(level + (ROOT_LO * c) // UNIT)) // SCALE if w0 else 0
                    us = [(kk, (w * half_mass(level + (ROOT_LO * grid[kk]) // UNIT)) // SCALE, step_w) for kk, w, step_w in steps]
                    u_entry = (k0, rem, dg0, u0, us)
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
                        u_entry = (kp, rem_p, dgp, up0, us)

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
                guard_entries = (
                    u_entry,
                    (k0, rem, dg0, g0_1, ws1),
                    (k0, rem, dg0, g0_3, ws3),
                    (k0, rem, dg0, g0_2, ws2),
                )
                for guard_entry, tables in zip(guard_entries, v_tables, strict=True):
                    prepared_guard = prepare_panel_entry(guard_entry, side_integral=True)
                    for table, source_samples in zip(tables, v_source_samples[j]):
                        table[iside + j] = integrate_panel_entry(prepared_guard, source_samples)

        del K_steps, S1_tbl_full, v_source_samples
        (V_0U, V_2U, V_4U, V_6U), (V_0W1, V_2W1, V_4W1), (V_0W13, V_2W13), (V_0W2,) = v_tables
        del v_tables

        # Guarded kernels for O_g and E_pt
        # 3. Holder-compressed GM4 table E_gm4 with corner caching and 5 guarded hops
        if endpoint_table is None and not centers_only:
            W_internal = (W1_1, W1_3, W1_5, W2_11)
            nW_int = len(W_internal)
            zero_wint = (0,) * nW_int

            @lru_cache(maxsize=2*side + 4)
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

            Kg_gm4_kerns = [None] * csize
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
                    for r, w_dg in zip(*klist, strict=True):
                        gtot = isqrt_(isqrt_(u00[r] * u01[r]) * isqrt_(u10[r] * u11[r]))
                        if has_w:
                            wr00, wr01, wr10, wr11 = w00[r], w01[r], w10[r], w11[r]
                            for m_ in range(nW_int):
                                gtot += isqrt_(isqrt_(wr00[m_] * wr01[m_]) * isqrt_(wr10[m_] * wr11[m_]))
                        if gtot:
                            row_kg.append((r, (w_dg * gtot) // SCALE))
                    Kg_gm4_kerns[ics + j] = pack_kernel(row_kg)
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
                        for r, w_dg in zip(*klist, strict=True):
                            v = in_gm4[jcs + r]
                            if v:
                                tot += w_dg * v
                        out[ics + j] = tot // denom_k
                return out

            W_t = (W1_1, W1_3, W1_5, W1_7, W2_11, W2_13, W2_31, W2_33)
            V_t = (V_0U, V_2U, V_4U, V_6U, V_0W1, V_2W1, V_4W1, V_0W13, V_2W13, V_0W2)
            WV_t = W_t + V_t
            nWt = len(W_t)
            nWV = len(WV_t)
            zero_wv = (0,) * nWV
            nu = nu_panels

            @lru_cache(maxsize=side)
            def terminal_lookup(j):
                return row_lookup(WV_t, j, grid, SCALE)

            @lru_cache(maxsize=2*side + 4)
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
                lookup = terminal_lookup(j)
                edges = []
                prev = None
                for kk in range(nu + 1):
                    h_l = (kk * U_c) // nu
                    h_h = -(-(kk * U_c) // nu)
                    g_val = density(h_h)
                    up_val = half_mass(rb + (ROOT_LO * h_l) // UNIT)
                    bm = lookup(M_l - h_h, M_h - h_l)
                    if kk == 0:
                        bp = bm
                    elif M_h + h_h < bound:
                        bp = lookup(M_l + h_l, M_h + h_h)
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
                        bpu = lookup(M_l + hu_l, M_h + hu_h)
                    else:
                        bpu = zero_wv
                    cur = (g_val, up_val, bm, bp, gu_val, um_val, upu_val, bpu)
                    if prev is not None:
                        pg, pup, pbm, pbp, pgu, pum, pupu, pbpu = prev
                        edges.append((
                            curvature_floor(isqrt_(pg * g_val), U_c * U_c, nu**2 * space_squared),
                            isqrt_(pup * up_val),
                            _edge_means(pbm, bm),
                            _edge_means(pbp, bp),
                            curvature_floor(isqrt_(pgu * gu_val), U_u * U_u, nu**2 * space_squared),
                            isqrt_(pum * um_val),
                            isqrt_(pupu * upu_val),
                            _edge_means(pbpu, bpu),
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
            terminal_lookup.cache_clear()

            E_gm4 = list(T_gm4)
            cur_gm4 = T_gm4
            for _ in range(5):
                mid_gm4 = apply_gm4_kern(cur_gm4, K_gm4_kerns)
                cur_gm4 = apply_gm4_kern(mid_gm4, Kg_gm4_kerns)
                E_gm4 = [a + b for a, b in zip(E_gm4, cur_gm4)]
            del K_gm4_kerns, Kg_gm4_kerns, T_gm4, cur_gm4, mid_gm4
            del e00, e01, e10, e11, W_t, V_t, WV_t, terminal_lookup, corner_time_edges
            progress("compressed guarded odd arms")
            if compute_endpoints:
                return E_gm4
        else:
            del K_gm4_kerns
            if endpoint_table is not None:
                assert len(endpoint_table) == csize
                E_gm4 = endpoint_table

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

        KW1_panels = quarter_panels(KW1_steps)
        KW13_panels = quarter_panels(KW13_steps)
        KW2_panels = quarter_panels(KW2_steps)
        del KW1_steps, KW13_steps, KW2_steps
        Og_mats = (
            [apply_panels(m_, KU_panels, ix_max2 + 1) for m_ in even]
            + [apply_panels(m_, KW1_panels, ix_max2 + 1) for m_ in even[:4]]
            + [apply_panels(m_, KW13_panels, ix_max2 + 1) for m_ in even[:2]]
            + [apply_panels(branches[0], KW2_panels, ix_max2 + 1)]
        )
        del KW2_panels, S0_panels
        progress("alternating branches and side guards")

        del kernels

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
            (11, 0), (12, 0), (1, 6), (6, 1),
            (0, 3), (0, 4), (0, 7), (0, 8), (0, 9), (0, 10),
            (1, 2), (2, 1), (0, 11), (0, 12), (13, 0), (0, 13),
            (1, 7), (7, 1), (1, 8), (8, 1),
        )
        # U and the W guards are outside. A pair cannot contain two V guards.
        outside_ids = {0, 1, 2, 3, 4, 11, 12}
        all_common_pairs = ((a, b) for a in range(11) for b in range(11) if a < 5 or b < 5)
        pair_specs += tuple(pair for pair in all_common_pairs if pair not in pair_specs)
        pair_specs = pair_specs[:44]
        assert len(set(pair_specs)) == len(pair_specs) == 44
        assert all(left in outside_ids or right in outside_ids for left, right in pair_specs)
        def point_endpoint_samples(columns, widths=None):
            result = []
            for values in columns:
                first = next((i for i, v in enumerate(values) if v), len(values))
                last = next((i for i in range(len(values)-1, first, -1) if values[i]), first)
                if first >= last:
                    result.append((first, first, [], []))
                    continue
                samples = tuple(map(quarter_nodes, values[first:last], values[first+1:last+1]))
                low, high = zip(*samples)
                if widths is not None:
                    spans = widths[first:last]
                    low = tuple(map(_mul, spans, low))
                    high = tuple(map(_mul, spans, high))
                result.append((first, last, low, high))
            return result

        n_pairs = len(pair_specs)
        T_mats = [array("I",[0]) * size for _ in range(n_pairs)]
        pack_row, lane_codec, lane_mask, fixed_reciprocal = packed_scaling(len(m_mats))

        for j in points:
            gj = grid[j]
            r_lo = (ROOT_LO * gj) // UNIT
            r_hi = -(-(ROOT_HI * gj) // UNIT)
            lookup = row_lookup(m_mats, j, grid, SCALE)
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
                left_events, right_events, used_nodes = [], [], []
                hs = [grid[kk] for kk in range(kmax + 1)]
                if hs[-1] < upper:
                    hs.append(upper)
                h_u_max = ((r_base - (UNIT_N // mesh)) * UNIT) // ROOT_HI
                if 0 < h_u_max < upper:
                    hs.append(h_u_max)
                    hs.append(h_u_max >> 1)
                    hs = sorted(set(hs))
                for h in hs:
                    dh = density(h)
                    if not dh:
                        break
                    km = gi_half - h
                    em = m_cache.get(km)
                    if em is None:
                        values = lookup(ml - h, mh - h)
                        em = values, pack_row(values)
                        m_cache[km] = em
                    u_m = half_mass(r_base - (-(-(ROOT_HI * h) // UNIT)))
                    u_p = half_mass(r_base + (ROOT_LO * h) // UNIT)
                    if h == 0:
                        dp = (u_p,) + em[0]
                    else:
                        kp = gi_half + h
                        ep = m_cache.get(kp)
                        if ep is None:
                            values = lookup(ml + h, mh + h)
                            ep = values, pack_row(values)
                            m_cache[kp] = ep
                        dp = (u_p,) + ep[0]
                    scaled = ((em[1] * (dh*fixed_reciprocal)) >> 96) & lane_mask
                    result = lane_codec.unpack(scaled.to_bytes(lane_codec.size,'little'))
                    left_events.append((dh*u_m//SCALE,) + result)
                    right_events.append(dp)
                    used_nodes.append(h)
                if len(used_nodes) < 2:
                    continue
                widths = [gaussian_span(right - left) for left, right in zip(used_nodes, used_nodes[1:])]
                # Split each product chord into its two factor chords.
                left_samples = point_endpoint_samples(zip(*left_events), widths)
                right_samples = point_endpoint_samples(zip(*right_events))
                cell = i*side + j
                for index, (left, right) in enumerate(pair_specs):
                    lf, ll, low_l, high_l = left_samples[left]
                    rf, rl, low_r, high_r = right_samples[right]
                    first = max(lf, rf)
                    stop = min(ll, rl)
                    if first >= stop:
                        continue
                    l0, l1 = first-lf, stop-lf
                    r0, r1 = first-rf, stop-rf
                    lower = sumprod(low_l[l0:l1],low_r[r0:r1])
                    upper = sumprod(high_l[l0:l1],high_r[r0:r1])
                    T_mats[index][cell] = (lower + upper) // panel_denominator

        k_indices = (0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 15, 16, 17, 8)
        K_T = [apply_panels(T_mats[k], K_panels) for k in k_indices]
        h1_U = [apply_panels(m_, KU_panels) for m_ in K_T]
        h1_W1 = [apply_panels(K_T[k], KW1_panels) for k in (0, 1, 3)]
        h1_W13 = [apply_panels(K_T[k], KW13_panels) for k in (0, 1)]
        K_h1 = [apply_panels(h1_U[k], K_panels) for k in (0, 1, 2, 3, 5)] + [apply_panels(h1_W1[0], K_panels)]
        h2_U = [apply_panels(m_, KU_panels) for m_ in K_h1]
        h2_W1 = [apply_panels(K_h1[0], KW1_panels)]
        K_h2 = [apply_panels(h2_U[k], K_panels) for k in (0, 1)]
        h3_U = [apply_panels(m_, KU_panels) for m_ in K_h2]
        E_pt = T_mats + h1_U + h1_W1 + h1_W13 + h2_U + h2_W1 + h3_U

        assert len(E_pt) == 74
        target_endpoints = E_pt
        center_i_lim = center_cells+1
        center_ids = (0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 15, 16, 44, 45, 46, 47)
        expanded = {j: apply_panels(E_pt[j], S1_panels, center_i_lim) for j in center_ids}
        P1_S1 = [expanded[j][:(center_axis + 1)*side] if j in expanded else
                 apply_panels(m_, S1_panels, center_axis + 1) for j, m_ in enumerate(E_pt)]
        center_guards = [expanded[j] if j in expanded else
                         table + array("I", [0]) * (center_i_lim*side-len(table))
                         for j, table in enumerate(P1_S1)]
        ee_both = op_D_batch(even, even, [(a, b) for a in range(len(even)) for b in range(len(even))], center_i_lim)
        if not centers_only:
            e_sub1 = T_mats[:12] + h1_U[:4] + h2_U[:2]
            eo_list = odd + e_sub1
            eo_pairs = []
            near_by_far = (range(18), (0, 1, 2, 3, 4, 5, 12, 13), range(7), range(4))
            for far, near_indices in enumerate(near_by_far):
                for near in near_indices:
                    index = len(odd) + near
                    eo_pairs.extend(((index, far), (far, index)))
            eo_both = op_D_batch(eo_list, eo_list, eo_pairs, ix_max1 + 1)

        del T_mats, K_T, h1_U, h1_W1, h1_W13, K_h1, h2_U, h2_W1, K_h2, h3_U, E_pt, K_panels, KU_panels, KW1_panels, KW13_panels, S1_panels
        del m_mats
        if not centers_only:
            del e_sub1, eo_list
        exp_floor.cache_clear()
        progress("point histories and mixed-arm neighbors")

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

    @lru_cache(maxsize=64)
    def selected_negative_lookup(index):
        return row_lookup(P1_S1, index, grid, SCALE)

    if not centers_only:
        BI_gm4 = make_gm4_sum(even, ix_max1)
        BI_ranges = arm_ranges(even, side, ix_max2, SCALE)
        BO_gm4 = make_gm4_sum(odd, ix_max2)
        Og_gm4 = make_gm4_sum(Og_mats, ix_max2)

        nA_s = len(same_even)
        AP_s = tuple(same_even + P1_S1)

        nA_d = len(pos_both)
        AP_d = tuple(pos_both + ee_both + eo_both)

        def f2_opp(gx, gb):
            xq = -(-gx // v_div)
            bq = -(-gb // v_div)
            num = 9 * (xq * xq + bq * bq) + (-(-(12 * ROOT_HI * xq * bq) // UNIT))
            return (C2_NUM * exp_floor(-(-num // 3200))) // UNIT

        counts = [0, 0, 0, 0]
        curvature_den1 = 48 * nt1 * nt1 * space_squared
        curvature_den2 = 12 * space_squared
        old_sq = None
        for ix in range(ix_max1 + 1):
            gx = grid[ix]
            gx_half = gx // 2
            lookup_s = selected_negative_lookup(ix) if negative_only else row_lookup(AP_s, ix, grid, SCALE)
            lookup_d = row_lookup(AP_d, ix, grid, SCALE)
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
                        yv = lookup_s(ylo, yhi)
                        aps_cache[yk] = yv
                    zk = (zlo, zhi)
                    zv = aps_cache.get(zk)
                    if zv is None:
                        zv = lookup_s(zlo, zhi)
                        aps_cache[zk] = zv
                    if (ROOT_HI * zlo) // UNIT < gx_half:
                        zvd = apd_cache.get(zk)
                        if zvd is None:
                            zvd = lookup_d(zlo, zhi)
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
                    zdl = _edge_means(zd0, zd1) if (zd0 is not None and zd1 is not None) else None
                    lrow.append((wl, _edge_means(y0, y1), _edge_means(z0, z1), zdl))
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
                    zds = _edge_means(zd0, zd1) if (zd0 is not None and zd1 is not None) else None
                    srow.append((ws, _edge_means(y0, y1), _edge_means(z0, z1), zds))
                squares.append(srow)

            if ix > 0:
                dx8 = dg[ix - 1] // U0125
                ics = (ix - 1) * cside
                for jb in range(jb_max1):
                    bi_gm = 0 if negative_only else BI_gm4[ics + jb]
                    e_gm = E_gm4[ics + jb]
                    if not bi_gm and not e_gm:
                        continue
                    sr0, sr1 = old_sq[jb], squares[jb]
                    if sr0 is None or sr1 is None:
                        continue
                    vol = dx8 * (dg[jb] // U0125)
                    dx2, db2 = dg[ix-1]**2, dg[jb]**2
                    variance_scale = 8*dx2 + db2
                    twice_mid_s = grid[jb] + grid[jb+1] + 2*ROOT_LO*(grid[ix-1] + gx)//UNIT
                    curvature_base = 36*nt1**2*(dx2 + db2) + 3*twice_mid_s**2 + 7*variance_scale
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
                        log_gain = curvature_base + 12*kt*(kt+1)*variance_scale
                        w = curvature_floor(w, log_gain, curvature_den1)
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

        progress("degree-one root cells")
        nA0 = len(same_even)
        AP2 = tuple(same_even + P1_S1)

        sums_set = {grid[jy] + grid[kz] for jy in range(jy_max2 + 1) for kz in range(jy, jy_max2 + 1)}

        def get_b_vals(ix):
            lookup = selected_negative_lookup(ix) if negative_only else row_lookup(AP2, ix, grid, SCALE)
            gx = grid[ix]
            r2_lo = (2 * ROOT_LO * gx) // UNIT
            r2_hi = -(-(2 * ROOT_HI * gx) // UNIT)
            out_c = {}
            for s in sums_set:
                out_c[s] = lookup(r2_lo + s, r2_hi + s)
            return out_c

        w_prev_sq = None
        b_prev = get_b_vals(0)
        planes = degree_two_weight_planes(grid, ix_max2, jy_max2, v_div, ROOT_HI, C3_NUM)
        for ix, wsq in enumerate(planes):
            if ix > 0:
                dx8 = dg[ix - 1] // U0125
                ics = (ix - 1) * cside
                b_cur = get_b_vals(ix)
                ap_gm_cache, ap_edges = {}, {}

                def ap_edge(l, r):
                    key = (l, r)
                    value = ap_edges.get(key)
                    if value is None:
                        value = (_edge_means(b_prev[l], b_prev[r]),
                                 _edge_means(b_cur[l], b_cur[r]))
                        ap_edges[key] = value
                    return value

                def get_ap_gm(s00, s01, s10, s11):
                    key = (s00, s01, s10, s11)
                    res = ap_gm_cache.get(key)
                    if res is not None:
                        return res
                    a,c = ap_edge(s00,s01)
                    b,d = ap_edge(s10,s11)
                    terms = [isqrt_(isqrt_(u*v)*isqrt_(x*y)) if u and v and x and y else 0
                             for u,v,x,y in zip(a,b,c,d)]
                    value = sum(terms[:nA0]), sum(terms[nA0:])
                    ap_gm_cache[key] = value
                    return value

                for jy in range(jy_max2):
                    dy8 = dg[jy] // U0125
                    gy0, gy1 = grid[jy], grid[jy + 1]
                    w0_row = w_prev_sq[jy]
                    w1_row = wsq[jy]
                    bi_y = BI_gm4[ics + jy]
                    bo_y = BO_gm4[ics + jy]
                    e_y = E_gm4[ics + jy]
                    og_y = Og_gm4[ics + jy]
                    for kz in range(jy, jy_max2):
                        w = isqrt_(w0_row[kz] * w1_row[kz])
                        if w < WEIGHT_CUT:
                            continue
                        dz8 = dg[kz] // U0125
                        log_gain = 33*dg[ix-1]**2 + 12*(dg[jy]**2 + dg[kz]**2)
                        w = curvature_floor(w, log_gain, curvature_den2)
                        gz0, gz1 = grid[kz], grid[kz + 1]
                        a0, p1_ = get_ap_gm(gy0 + gz0, gy0 + gz1, gy1 + gz0, gy1 + gz1)
                        vol = dx8 * dy8 * dz8 * (1 if kz == jy else 2)
                        bi_z = BI_gm4[ics + kz]
                        if a0:
                            cover = either_arm(SCALE, bi_y, bo_y, BI_ranges[ics + jy],
                                               bi_z, BO_gm4[ics + kz], BI_ranges[ics + kz])
                            counts[2] += 3 * vol * ((((w * a0) // SCALE) * cover) // SCALE)
                        if p1_:
                            counts[3] += (
                                3
                                * vol
                                * ((((w * p1_) // SCALE) * ((bi_y * bi_z + e_y * Og_gm4[ics + kz] + og_y * E_gm4[ics + kz]) // SCALE)) // SCALE)
                            )
                b_prev = b_cur
            w_prev_sq = wsq
    else:
        counts = [0, 0, 0, 0]

    if positive_only or roots_only:
        den1 = 64 * SCALE * n * n * nt1
        den2 = 512 * SCALE * n * n * n
        return (counts[0], counts[1]), den1, (counts[2], counts[3]), den2, [0, 0, 0, 0, 0], den2*SCALE

    root_limit = center_axis + 1
    double_guards = ee_both
    c1 = [0] * (single_limit * cside)
    for j, (v0, v1) in enumerate(zip(grid[:root_limit-1], grid[1:root_limit])):
        lower0 = (2 * ROOT_LO * v0) // UNIT
        upper0 = -(-(2 * ROOT_HI * v0)//UNIT)
        lower1 = (2 * ROOT_LO * v1) // UNIT
        upper1 = -(-(2 * ROOT_HI * v1)//UNIT)
        lo_lookup = selected_negative_lookup(j)
        hi_lookup = selected_negative_lookup(j+1)
        base = j*side
        next_base = (j+1)*side
        endpoint_means = []
        for k in range(cside):
            d0, d1 = grid[k], grid[k+1]
            left_total = right_total = whole_total = 0
            for table in target_endpoints:
                a, b = table[base+k], table[base+k+1]
                c, d = table[next_base+k], table[next_base+k+1]
                if a and b and c and d:
                    u, v = isqrt_(a*b), isqrt_(c*d)
                    if d1-d0 >= 4*UNIT:
                        left_total += isqrt_(isqrt_(a*u)*isqrt_(c*v))
                        right_total += isqrt_(isqrt_(u*b)*isqrt_(v*d))
                    else:
                        whole_total += isqrt_(u*v)
            if d1-d0 >= 4*UNIT:
                mid = (d0 + d1)//2
                endpoint_means.extend(((d0,mid,left_total),(mid,d1,right_total)))
            else:
                endpoint_means.append((d0,d1,whole_total))
        cached, edge_cache, plane_cache = {}, {}, {}
        r_hi0 = -(-(ROOT_HI*v0)//UNIT)
        r_hi1 = -(-(ROOT_HI*v1)//UNIT)

        def edge(left, right):
            key = left, right
            result = edge_cache.get(key)
            if result is None:
                a = lo_lookup(lower0+left,upper0+left)
                b = lo_lookup(lower0+right,upper0+right)
                c = hi_lookup(lower1+left,upper1+left)
                d = hi_lookup(lower1+right,upper1+right)
                result = (_edge_means(a,b),_edge_means(c,d))
                edge_cache[key] = result
            return result

        def plane(x, d0, d1):
            t = d0-x//2, d1-x//2
            value = plane_cache.get(t)
            if value is None:
                u,v = t
                a,b = density(u+r_hi0), density(v+r_hi0)
                c,d = density(u+r_hi1), density(v+r_hi1)
                value = isqrt_(isqrt_(a*b)*isqrt_(c*d))
                plane_cache[t] = value
            return value

        for i, (x0,x1) in enumerate(zip(grid[:single_limit-1],grid[1:single_limit])):
            total = 0
            for d0,d1,e in endpoint_means:
                if not e:
                    continue
                if lower0-x1+d0 < 0 or upper1-x0+d1 >= bound:
                    continue
                key = (d0-x1,d0-x0,d1-x1,d1-x0)
                sibling = cached.get(key)
                if sibling is None:
                    a,c = edge(key[0],key[1])
                    b,d = edge(key[2],key[3])
                    sibling = sum(isqrt_(isqrt_(u*v)*isqrt_(x*y))
                                  for u,v,x,y in zip(a,b,c,d)
                                  if u and v and x and y)
                    cached[key] = sibling
                if not sibling:
                    continue
                w = isqrt_(plane(x0,d0,d1)*plane(x1,d0,d1))
                total += w * e * sibling * gaussian_span(d1-d0)
            floor = total // (SCALE*SCALE*UNIT_N)
            c1[i*cside+j] = curvature_floor(floor, dg[i]**2 + 8*dg[j]**2, 4*space_squared)

    selected_negative_lookup.cache_clear()
    stable_c2 = make_gm4_sum(ee_both, center_cells)
    def opposite_half_width(x, v):
        r_lo = ROOT_LO * v // UNIT
        r_hi = -(-(ROOT_HI * v) // UNIT)
        return x // 2 - r_hi, x // 2 - r_lo

    mixed_c2, mixed_support = mixed_arm_floors(
        grid, root_limit, center_cells, target_endpoints, Og_mats,
        opposite_half_width, density, SCALE, UNIT_N)
    c2 = [a+b for a,b in zip(stable_c2,mixed_c2,strict=True)]
    center_guards = center_guards + ee_both
    cover_gm = make_gm4_sum(center_guards, center_cells)
    prev_companion = None
    shell_total = 0
    overlap_total = 0
    prev_plane = None
    pair_totals = [0,0,0]
    for i, x in enumerate(grid[:center_cells+1]):
        sum_lo = 2*ROOT_LO*x//UNIT
        sum_hi = -(-(2*ROOT_HI*x)//UNIT)
        supported = [
            (i,j) in mixed_support or any(guard[i*side+j] for guard in double_guards) or (i < single_limit and bool(c1[i*cside+j] or (j and c1[i*cside+j-1]) or (i and c1[(i-1)*cside+j]) or (i and j and c1[(i-1)*cside+j-1])))
            for j in range(root_limit)
        ]
        plane = [[0] * root_limit for _ in range(root_limit)]
        zero_cover = None
        companion = [[zero_cover] * root_limit for _ in range(root_limit)]
        cover_nonzero = [any(guard[i*side+j] for guard in center_guards) for j in range(root_limit)]
        cover_lookup = row_lookup(center_guards, i, grid, SCALE)
        cover_cache = {}
        x_term = 33*x*x
        for j, y in enumerate(grid[:root_limit]):
            if not supported[j] or not y:
                continue
            y_term = x_term + 12*y*y
            for k in range(j, root_limit):
                z = grid[k]
                magnitude = y + z
                if sum_lo <= magnitude:
                    break
                values = zero_cover
                if (cover_nonzero[k] or supported[k]) and 2*y+z < sum_lo:
                    values = cover_cache.get(magnitude)
                    if magnitude not in cover_cache:
                        values = cover_lookup(sum_lo-magnitude, sum_hi-magnitude)
                        if not any(values):
                            values = zero_cover
                        cover_cache[magnitude] = values
                    companion[j][k] = values
                if not supported[k] and values is zero_cover:
                    continue
                exponent = y_term + 12*z*(y+z) - (24*ROOT_LO*x*magnitude)//UNIT
                q = -(-(exponent*EXP_DEN)//(2*space_squared))
                weight = C3_NUM*exp_floor(q)//SCALE
                plane[j][k] = plane[k][j] = weight
        if prev_plane is not None:
            third_cache = {}
            for j in range(center_axis):
                u1 = c1[(i-1)*cside+j] if i < single_limit else 0
                u2 = c2[(i-1)*cside+j]
                if not (u1 or u2):
                    continue
                for k in range(j, center_axis):
                    v1 = c1[(i-1)*cside+k] if i < single_limit else 0
                    v2 = c2[(i-1)*cside+k]
                    cover_count = cover_gm[(i-1)*cside+k]
                    if not (v1 or v2 or cover_count):
                        continue
                    w = isqrt_(square_cell_mean(prev_plane,j,k)*square_cell_mean(plane,j,k))
                    gain = 33*dg[i-1]**2 + 12*(dg[j]**2 + dg[k]**2)
                    w = curvature_floor(w, gain, 12*space_squared)
                    volume = 3*(dg[i-1]//U0125)*(dg[j]//U0125)*(dg[k]//U0125)
                    pair_volume = volume if j == k else 2*volume
                    pair_totals[0] += pair_volume*(w*u1*v1//SCALE)
                    pair_totals[1] += pair_volume*(w*(u1*v2+u2*v1)//SCALE)
                    stable_j = stable_c2[(i-1)*cside+j]
                    stable_k = stable_c2[(i-1)*cside+k]
                    pair_totals[2] += pair_volume*(w*u2*v2//SCALE + w*stable_j*stable_k//SCALE)
                    if w and k > j and (cover_count or v1 or v2):
                        third = companion_mean(prev_companion, companion, j, k, third_cache)
                        if not third:
                            continue
                        prefix = w * (u1 + u2) * third
                        shell_total += volume * (prefix * cover_count // (SCALE * SCALE))
                        overlap_total += 2 * volume * (prefix * (v1 + v2) // (SCALE * SCALE))
        prev_plane = plane
        prev_companion = companion
    pair_totals.extend((shell_total, overlap_total))
    progress("isolated-center charge credits")
    den1 = 64 * SCALE * n * n * nt1
    den2 = 512 * SCALE * n * n * n
    return (counts[0], counts[1]), den1, (counts[2], counts[3]), den2, pair_totals, den2*SCALE
