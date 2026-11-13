#!/usr/bin/env python3
"""certify.py -- the certificate of Proposition C of cutwidth.md:  sup over all real t of G_d(t) <= 683/5000.

Declared experiment.  Python 3.12 standard library only; reads no file; input: one JSON value on stdin (ignored;
`{}`); output: exactly one JSON object on stdout.  "certified": true is printed only if every check below passed.

    echo '{}' | python3 -I -B certify.py

No math-library function enters any bound (math is not imported): every bound is computed with outward-rounded
interval arithmetic and with the enclosures of exp and of the normal cdf Phi defined in this file, or in exact
integer/rational arithmetic.  Every statement used is proved in cutwidth.md, sections 5 and 6 (section 6 names
the function of this file that implements each bound).

 (C1) the window energy f (cutwidth.md section 4) on {0,1}^10: range of 2f, flip sizes, mixed differences, the 48
      automorphisms, (k-1)^+ = sum_{i<j} b_i b_j - b_1 b_2 b_3.  Exhaustive, exact.
 (C2) the wave window (section 5.1): covariance w_d = (1 + d/3) 2^(-d/2) at tree distance d, exact in Q(sqrt2).
 (C3) sup_t G_d(t) <= 683/5000, G_d = V - (3/2) P2 + P3/2 (section 5.2): tails + a contiguous cover of [-5/4, 2]
      by leaves [a, b], on each of which   Vup(a,b) + P3up(b)/2 - (3/2) P2low(a,b) <= 683/5000   (Lemma 5.8),
      the comparison being exact (the binary64 bound, as a rational, against 683/5000).
 (C4) the numerical witnesses of section 7 (tables WITNESSES and GAMMA_WITNESSES below), exact rationals.
 (C5) Mills: Phi(-9) <= phi(9)/9 < 1.19e-19 < the binary64 number 1.2e-19 used by Phi_pt; exact rationals.
 (C6) global submodularity (Lemma 4.4), exact enumeration, with the sharpness witness lambda = 51/100.

Method for (C3) (section 6.3-6.5 of cutwidth.md; tau = tan psi):
  V(t), t >= 0:  (3/pi) int_0^sqrt3 g_t(tau) dtau,  g_t = exp(-al t^2) Phi(be t) / ((1+tau^2) sqrt(9+8tau^2)),
                 al = 9(1+tau^2)/(2(9+8tau^2)),  be = (1 - kap(1+tau^2))/sqrt(9+8tau^2),  kap = 6 sqrt2 - 8.
                 For t in [a,b] (0 <= a): g_t(tau) <= exp(-al a^2) Phi(max(be a, be b)) / (...), nonincreasing in
                 tau: its left Riemann sum on the cells [k/2048, (k+S)/2048], k = 0, S, 2S, ... < 3548
                 (3548/2048 > sqrt3) is an upper bound.  S = 16 first, S = 1 only where needed.
  V(t), t <= 0:  Phi(t) - Phi(c0 t) + (3/pi) int_0^sqrt3 exp(-al t^2) Phi(c0 t/sqrt(9+8tau^2)) / (...) dtau,
                 c0 = 3/(2 sqrt2); for t in [a,b], b <= 0: Phi(t) - Phi(c0 t) <= (c0-1)|a| phi(b), and the
                 integrand is <= [exp(-al b^2)/(...)](left end) * Phi(c0 b/sqrt(9+8tau^2))(right end) on a cell.
  P2(t) >= P2(a,b) = int_{u>=b} phi(u) K_a(u) du, K_a(u) = E[G(x1) G(x2)], G(x) = 1{x < M} H(kap_u - lam x),
                 H = (2Phi - 1)^+, M = 3(a - r1 u), kap_u = sqrt12(a - 5u/6), lam = 2 sqrt6/3, x1 = x,
                 x2 = -x/2 + (sqrt3/2) y, x, y independent standard normal.  Lower sums on the lattice
                 x_k = x0 + k sqrt3 h, y_l = sqrt3 x0 + l h, u_i = b + i (sqrt6/5) h: every H argument is
                 base + J sqrt2 h (J integer), every indicator is replaced by "J >= cut_i".
  P3(t) <= (1/r1 - 1) t phi(0) (1 - exp(-18 (1-r1)^2 t^2)) for t > 0, P3 = 0 for t <= 0 (Lemma 5.6).
  Tails: t <= -5/4: G_d <= Phi(-5/4);  t >= 2: G_d <= (7/2) Phi(-2)  (Lemma 5.7).
Deterministic; no randomness, no threads, no file access.
"""
import itertools, json, sys, time
from fractions import Fraction as Fr

# ====================================================================================================== (C4)
# Witnesses of cutwidth.md section 7.  Each entry: (label, function without arguments returning the truth value of
# an exact rational statement).  D("1.41421") is the exact rational 141421/100000; PI_Q = (lower, upper) are the
# exact rational bounds of pi from Machin's formula computed below.
D = Fr
GAMMA_D = Fr(683, 5000)                       # the target constant gamma_d of Proposition C

def exp_partial(z, terms):                    # partial sum of the exponential series: < e^z for rational z > 0
    s, t = Fr(0), Fr(1)
    for k in range(terms):
        s += t; t = t * z / (k + 1)
    return s

# ---- WITNESSES: only statements that do NOT depend on gamma_d
WITNESSES = [
    ("N1 1.41421^2 < 2 < 1.41422^2", lambda: D("1.41421") ** 2 < 2 < D("1.41422") ** 2),
    ("N1 3.3166^2 < 11, 9.7467^2 < 95", lambda: D("3.3166") ** 2 < 11 and D("9.7467") ** 2 < 95),
    ("N1 3.14159 < pi < 3.1416", lambda: D("3.14159") < PI_Q[0] and PI_Q[1] < D("3.1416")),
    ("N2 6*1.41422 + 60/3.3166 + 19 + 216/9.7467 < 67.738",
     lambda: 6 * D("1.41422") + 60 / D("3.3166") + 19 + 216 / D("9.7467") < D("67.738")),
    ("N2 3.1416*67.738/2 < 106.41", lambda: D("3.1416") * D("67.738") / 2 < D("106.41")),
    ("N2 eps^2 term: 45*2*1.41422/(2*3.14159) < 21", lambda: 45 * 2 * D("1.41422") / (2 * D("3.14159")) < 21),
    ("N2 w_d^2/(1 - w_d^2) = 8, 25/11, 1, 49/95 for d = 1..4",
     lambda: [x / (1 - x) for x in [((1 + Fr(d, 3)) ** 2) / 2 ** d for d in (1, 2, 3, 4)]]
     == [8, Fr(25, 11), 1, Fr(49, 95)]),
    ("N3 4*48 = 192 < 13.857^2", lambda: 4 * 48 < D("13.857") ** 2),
    ("N3 3 < 1.7321^2", lambda: 3 < D("1.7321") ** 2),
    ("N3 13.857 + (15/2)*1.7321 < 26.85", lambda: D("13.857") + Fr(15, 2) * D("1.7321") < D("26.85")),
    ("N3 (12 - 15/2)^2 * 3 = 243/4 < 192", lambda: (12 - Fr(15, 2)) ** 2 * 3 == Fr(243, 4) < 192),
    ("N3 26.85*1.41422/2 < 18.99", lambda: D("26.85") * D("1.41422") / 2 < D("18.99")),
    ("N4 45^2 < 17.96^2 * 2 * 3.14159", lambda: 45 ** 2 < D("17.96") ** 2 * 2 * D("3.14159")),
    ("N4 18.99 + 17.96 + 0.46 <= 38 and 30 < 0.46*2^35",
     lambda: D("18.99") + D("17.96") + D("0.46") <= 38 and 30 < D("0.46") * 2 ** 35),
    ("N5 106.41*9 = 957.69 < 958", lambda: D("106.41") * 9 < 958),
    ("N6 e^4.5 > 90 (partial sum of the exponential series)", lambda: exp_partial(Fr(9, 2), 30) > 90),
    ("N6 4/(90*2.5066) < 0.0178 and 2.5066^2 < 2*3.14159",
     lambda: Fr(4) / (90 * D("2.5066")) < D("0.0178") and D("2.5066") ** 2 < 2 * D("3.14159")),
    ("N7 2^40 < 22^9", lambda: Fr(2) ** 40 < Fr(22) ** 9),
]

# ---- GAMMA_WITNESSES: the statements of section 7 that DO depend on gamma_d (use GAMMA_D); same format.
GAMMA_WITNESSES = [
    ("N6 0.0178 < gamma_d = 683/5000", lambda: D("0.0178") < GAMMA_D == Fr(683, 5000)),
    ("N7 38/22 > 1.3634 = 6817/5000 = 3/2 - gamma_d",
     lambda: Fr(38, 22) > D("1.3634") == Fr(6817, 5000) == Fr(3, 2) - GAMMA_D),
    ("N8 1 + 1/(2 gamma_d) = 3183/683 and 3183/683 - 233/50 = 11/34150 > 0",
     lambda: 1 + 1 / (2 * GAMMA_D) == Fr(3183, 683) and Fr(3183, 683) - Fr(233, 50) == Fr(11, 34150) > 0),
]
# ============================================================================================================

# ------------------------------------------------------------------ fixed parameters of (C3)
H_LATTICE = 0.012            # h of the P2 lattice (a binary64 number; treated as the exact real it is)
X0, XMAX, ULEN = -3.5, 6.5, 1.0
KDEN = 64                    # base intervals [k/64, (k+1)/64], -80 <= k < 128
MAXDEPTH = 7                 # at most 7 bisections of a base interval
T_LO, T_HI = -1.25, 2.0
NT, TDEN = 3548, 2048        # tau cells: [0, 3548/2048], 3548/2048 > sqrt3
STRIDES = (16, 1)            # tau cell widths 16/2048 (first attempt) and 1/2048

# ------------------------------------------------------------------ outward rounding, intervals
EPS52, TINY = 2.0 ** -52, 5e-324
def up(x): return x + (abs(x) * EPS52 + TINY)
def dn(x): return x - (abs(x) * EPS52 + TINY)
def iadd(a, b): return (dn(a[0] + b[0]), up(a[1] + b[1]))
def isub(a, b): return (dn(a[0] - b[1]), up(a[1] - b[0]))
def imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (dn(min(p)), up(max(p)))
def idiv(a, b):
    assert b[0] > 0 or b[1] < 0
    p = (a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1])
    return (dn(min(p)), up(max(p)))
def pt(x): return (x, x)
def fr_iv(q):                       # enclosure of a rational (exact comparisons)
    f = float(q)
    lo, hi = f, f
    while Fr(lo) > q: lo = dn(lo)
    while Fr(hi) < q: hi = up(hi)
    return (lo, hi)
def isqrt(a):                       # enclosure of sqrt over a positive interval
    assert a[0] > 0
    def root(x, lower):
        r = x ** 0.5                # only a starting value; the loops below verify it
        assert r > 0
        if lower:
            while up(r * r) > x: r = dn(r)
        else:
            while dn(r * r) < x: r = up(r)
        return r
    return (root(a[0], True), root(a[1], False))

def machin_pi():                    # 16 atan(1/5) - 4 atan(1/239), alternating partial sums bracket atan
    def at(x, n): return sum(Fr((-1) ** k) * x ** (2 * k + 1) / (2 * k + 1) for k in range(n + 1))
    lo = 16 * at(Fr(1, 5), 21) - 4 * at(Fr(1, 239), 6)
    hi = 16 * at(Fr(1, 5), 22) - 4 * at(Fr(1, 239), 5)
    assert lo < hi
    return lo, hi
def inv_sqrt_2pi():
    f = 0.3989422804014327
    lo, hi = f, f
    while Fr(lo) ** 2 * 2 * PI_Q[1] > 1: lo = dn(lo)
    while Fr(hi) ** 2 * 2 * PI_Q[0] < 1: hi = up(hi)
    return (lo, hi)

PI_Q = PI = INVS2PI = S2 = S3 = S6 = S12 = R1 = LAM = KAP = C0 = THREE_OVER_PI = None
def init_constants():
    global PI_Q, PI, INVS2PI, S2, S3, S6, S12, R1, LAM, KAP, C0, THREE_OVER_PI
    PI_Q = machin_pi()
    PI = (fr_iv(PI_Q[0])[0], fr_iv(PI_Q[1])[1])
    INVS2PI = inv_sqrt_2pi()
    S2, S3, S6 = isqrt(pt(2.0)), isqrt(pt(3.0)), isqrt(pt(6.0))
    S12 = imul(pt(2.0), S3)
    R1 = idiv(imul(pt(2.0), S2), pt(3.0))
    LAM = idiv(imul(pt(2.0), S6), pt(3.0))
    KAP = isub(imul(pt(6.0), S2), pt(8.0))
    C0 = idiv(pt(3.0), imul(pt(2.0), S2))
    THREE_OVER_PI = idiv(pt(3.0), PI)
    assert Fr(NT, TDEN) ** 2 > 3 and KAP[0] > 0 and C0[0] > 1 and 0 < R1[0] and R1[1] < 1

# ------------------------------------------------------------------ exp and Phi
def exp_neg(z):
    """enclosure (lo, hi) of e^z for a float z <= 0"""
    assert z <= 0.0
    if z < -700.0: return (0.0, 1e-300)
    s, y = 0, z
    while y < -0.03125: y *= 0.5; s += 1
    p = 1.0
    for n in range(12, 0, -1): p = 1.0 + (y / n) * p
    lo, hi = dn(p * (1.0 - 1e-15)), up(p * (1.0 + 1e-15))
    for _ in range(s): lo, hi = dn(lo * lo), up(hi * hi)
    return (lo, hi)
def exp_iv(a):  # a = (lo, hi), hi <= 0
    return (exp_neg(a[0])[0], exp_neg(a[1])[1])
def halfphi_part(ax):
    """enclosure of Phi(ax) - 1/2 = phi(ax) S(ax) for a float 1e-150 <= ax <= 9"""
    x2 = ax * ax
    T, S, n = ax, ax, 0
    while True:
        T = T * x2 / (2 * n + 3); n += 1; S += T
        if x2 <= 0.5 * (2 * n + 3) and T <= 1e-17 * S: break
    d = 5.0 * (n + 2) * 1.1102230246251565e-16
    Slo, Shi = dn(S * (1.0 - d)), up(up(S + T * (1.0 + d)) * (1.0 + d))
    e = exp_iv((-up(x2) * 0.5, -dn(x2) * 0.5))
    lo = dn(dn(e[0] * INVS2PI[0]) * Slo)
    hi = up(up(e[1] * INVS2PI[1]) * Shi)
    return (lo, min(hi, 0.5))
def Phi_pt(x):
    """enclosure of Phi(x) for a float x"""
    ax = abs(x)
    if ax > 9.0: return (0.9999999999999998, 1.0) if x > 0 else (0.0, 1.2e-19)    # Phi(-9) < 1.19e-19: (C5)
    if ax < 1e-150: return (dn(0.5 - 0.4 * ax), up(0.5 + 0.4 * ax))
    lo, hi = halfphi_part(ax)
    if x >= 0: return (dn(0.5 + lo), min(1.0, up(0.5 + hi)))
    return (max(0.0, dn(0.5 - hi)), up(0.5 - lo))
def H_lo(z):
    """lower bound of H(z) = (2 Phi(z) - 1)^+ for a float z"""
    if z <= 1e-150: return 0.0
    if z > 9.0: return 0.9999999999999998
    return dn(2.0 * halfphi_part(z)[0])

def check_mills():
    """(C5): Phi(-9) <= phi(9)/9 <= (2/5) e^(-81/2) / 9 < 1.19e-19 < the float 1.2e-19; 1/sqrt(2 pi) < 2/5"""
    ok = Fr(2, 5) ** 2 * 2 * PI_Q[0] > 1
    ok = ok and Fr(2, 5) / (exp_partial(Fr(81, 2), 200) * 9) < Fr(119, 10 ** 21) < Fr(1.2e-19)
    ok = ok and 1 - Fr(0.9999999999999998) > Fr(24, 10 ** 20)
    return ok, "Phi(-9) <= phi(9)/9 < 1.19e-19 < float(1.2e-19); 1 - 0.9999999999999998 > 2.4e-19"

# ------------------------------------------------------------------ (C1), (C6), (C2)
def f2(s, lnum=1, lden=2):
    """lden * f on a window pattern s (0 = centre, 1..3 = neighbours, children of i: 2i+2, 2i+3); lambda = lnum/lden"""
    k = sum(s[i] * s[2 * i + 2] * s[2 * i + 3] for i in (1, 2, 3))
    return lden * s[0] - lden * s[0] * s[1] * s[2] * s[3] - lnum * (1 - s[0]) * max(k - 1, 0)

def window_autos():
    """the 48 automorphisms of the rooted window: permute the branches, swap the two children inside a branch"""
    res = []
    for perm in itertools.permutations((1, 2, 3)):
        for sw in itertools.product((0, 1), repeat=3):
            m = [0] * 10
            for i, j in zip((1, 2, 3), perm):
                m[i] = j
                m[2 * i + 2], m[2 * i + 3] = (2 * j + 2, 2 * j + 3) if not sw[i - 1] else (2 * j + 3, 2 * j + 2)
            res.append(m)
    return res

def check_f():
    """(C1): range, flips, mixed differences, symmetry of 2f; the identity for (k-1)^+"""
    pats = list(itertools.product((0, 1), repeat=10))
    vals = {s: f2(s) for s in pats}
    par = [None, 0, 0, 0, 1, 1, 2, 2, 3, 3]
    def path(x):
        r = [x]
        while par[r[-1]] is not None: r.append(par[r[-1]])
        return r
    def dist(x, y):
        px, py = path(x), path(y)
        c = [z for z in px if z in py][0]
        return px.index(c) + py.index(c)
    big, md_all, npairs = set(), set(), {}
    for p, q in itertools.combinations(range(10), 2):
        npairs[dist(p, q)] = npairs.get(dist(p, q), 0) + 1
        mds = set()
        for s in pats:
            if s[p] == 0 and s[q] == 0:
                s11 = list(s); s11[p] = s11[q] = 1
                s10 = list(s); s10[p] = 1
                s01 = list(s); s01[q] = 1
                mds.add(vals[tuple(s11)] - vals[tuple(s10)] - vals[tuple(s01)] + vals[s])
        md_all |= mds
        if max(abs(m) for m in mds) > 1: big.add((p, q))
    rng = (min(vals.values()), max(vals.values()))
    flip = [max(abs(vals[s] - vals[s[:x] + (1 - s[x],) + s[x + 1:]]) for s in pats) for x in range(10)]
    sym = all(vals[s] == vals[tuple(s[m[x]] for x in range(10))] for m in window_autos() for s in pats)
    ident = all(max(k - 1, 0) == b[0] * b[1] + b[0] * b[2] + b[1] * b[2] - b[0] * b[1] * b[2]
                for b in itertools.product((0, 1), repeat=3) for k in [sum(b)])
    want_big = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
    ok = (rng == (-2, 2) and flip == [3, 2, 2, 2, 1, 1, 1, 1, 1, 1] and min(md_all) == -2 and max(md_all) == 1
          and big == want_big and sym and ident and npairs == {1: 9, 2: 12, 3: 12, 4: 12})
    return ok, ("range of 2f %s, flip sizes %s, mixed differences of 2f in [%d, %d], |.| = 2 only at %s, "
                "48 automorphisms %s, (k-1)^+ identity %s"
                % (list(rng), flip, min(md_all), max(md_all), sorted(big), sym, ident))

def check_globsub():
    """(C6): exact maximum over all assignments of the summed mixed difference, distances 1..4 (and none beyond)"""
    adj, dep = {0: []}, {0: 0}
    n, frontier = 1, [0]
    while frontier:
        nxt = []
        for v in frontier:
            if dep[v] >= 6: continue
            for _ in range(3 if v == 0 else 2):
                w = n; n += 1
                adj[w] = [v]; adj[v].append(w); dep[w] = dep[v] + 1; nxt.append(w)
        frontier = nxt
    def bfs(s):
        d, fr = {s: 0}, [s]
        while fr:
            nf = []
            for v in fr:
                for w in adj[v]:
                    if w not in d: d[w] = d[v] + 1; nf.append(w)
            fr = nf
        return d
    D0 = bfs(0)
    def window(a):
        assert dep[a] <= 4 and len(adj[a]) == 3
        pos = [a] + list(adj[a])
        for y in adj[a]:
            assert len(adj[y]) == 3
            pos += [z for z in adj[y] if z != a]
        return pos
    ok, pos51 = True, 0
    report = []
    for d in range(1, 7):
        q = 1
        for _ in range(d - 1): q = adj[q][1]          # a child of q (adj[q][0] is the parent)
        assert D0[q] == d
        Dq = bfs(q)
        cent = [a for a in adj if D0.get(a, 99) <= 2 and Dq.get(a, 99) <= 2]
        if not cent:
            report.append("d=%d: no window" % d)
            ok = ok and d >= 5
            continue
        ok = ok and d <= 4
        wins = {a: window(a) for a in cent}
        owners = {}
        for a in cent:
            for x in wins[a]:
                if x not in (0, q): owners.setdefault(x, set()).add(a)
        shared = sorted(x for x in owners if len(owners[x]) >= 2)
        for lnum, lden in ((1, 2), (51, 100)):
            tabs = []
            for a in cent:
                w = wins[a]
                sh = [x for x in w if x in shared]
                pv = [x for x in w if x not in (0, q) and x not in shared]
                tmax, tmin = {}, {}
                for ssh in itertools.product((0, 1), repeat=len(sh)):
                    hi, lo = None, None
                    for spv in itertools.product((0, 1), repeat=len(pv)):
                        asg = dict(zip(sh, ssh)); asg.update(zip(pv, spv))
                        v = []
                        for e in range(4):
                            asg[0], asg[q] = e & 1, e >> 1
                            v.append(f2([asg[x] for x in w], lnum, lden))
                        m = v[3] - v[1] - v[2] + v[0]
                        hi = m if hi is None else max(hi, m)
                        lo = m if lo is None else min(lo, m)
                    tmax[ssh], tmin[ssh] = hi, lo
                tabs.append(([shared.index(x) for x in sh], tmax, tmin))
            best, worst = None, None
            for sig in itertools.product((0, 1), repeat=len(shared)):
                hi = sum(t[1][tuple(sig[i] for i in t[0])] for t in tabs)
                lo = sum(t[2][tuple(sig[i] for i in t[0])] for t in tabs)
                best = hi if best is None else max(best, hi)
                worst = lo if worst is None else min(worst, lo)
            if lden == 2: ok = ok and best <= 0
            elif best > 0: pos51 += 1
            if lden == 100: ok = ok and ((best == 2) if d == 2 else (best <= 0))
            report.append("d=%d lambda=%d/%d: %d windows, |U|=%d, %d shared, max %s, min %s"
                          % (d, lnum, lden, len(cent), len(owners) + 2, len(shared), Fr(best, lden), Fr(worst, lden)))
    ok = ok and pos51 == 1
    return ok, "; ".join(report)

def check_wave():
    """(C2): covariance of the construction of section 5.1, exact in Q(sqrt2)"""
    # numbers p + q sqrt2 as (p, q); variables as coefficient vectors over (U, x1, x2, x3, d1, d2, d3)
    def mul(a, b): return (a[0] * b[0] + 2 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def add(a, b): return (a[0] + b[0], a[1] + b[1])
    Z, ONE = (Fr(0), Fr(0)), (Fr(1), Fr(0))
    base = [[ONE if i == j else (Fr(-1, 2), Fr(0)) if 1 <= i <= 3 and 1 <= j <= 3 else Z for j in range(7)]
            for i in range(7)]
    for i in (4, 5, 6): base[i][i] = (Fr(1, 12), Fr(0))
    r1 = (Fr(0), Fr(2, 3))
    def var(**c):
        v = [Z] * 7
        for k, x in c.items(): v[int(k[1:])] = x
        return v
    pos = [var(c0=ONE)]
    for i in (1, 2, 3): pos.append(var(c0=r1, **{"c%d" % i: (Fr(1, 3), Fr(0))}))
    for i in (1, 2, 3):
        for sg in (1, -1):   # child = sqrt2 z_i - U/2 +- d_i = (5/6) U + (sqrt2/3) x_i +- d_i
            pos.append(var(c0=(Fr(5, 6), Fr(0)), **{"c%d" % i: (Fr(0), Fr(1, 3)), "c%d" % (i + 3): (Fr(sg), Fr(0))}))
    def cov(a, b):
        s = Z
        for i in range(7):
            for j in range(7):
                s = add(s, mul(mul(a[i], base[i][j]), b[j]))
        return s
    par = [None, 0, 0, 0, 1, 1, 2, 2, 3, 3]
    def path(x):
        r = [x]
        while par[r[-1]] is not None: r.append(par[r[-1]])
        return r
    def dist(x, y):
        px, py = path(x), path(y)
        common = [z for z in px if z in py][0]
        return px.index(common) + py.index(common)
    w = {0: (Fr(1), Fr(0)), 1: (Fr(0), Fr(2, 3)), 2: (Fr(5, 6), Fr(0)), 3: (Fr(0), Fr(1, 2)), 4: (Fr(7, 12), Fr(0))}
    ok = all(cov(pos[x], pos[y]) == w[dist(x, y)] for x in range(10) for y in range(10))
    ok = ok and all(w[d] == ((1 + Fr(d, 3)) / 2 ** (d // 2), Fr(0)) if d % 2 == 0 else
                    w[d] == (Fr(0), (1 + Fr(d, 3)) / 2 ** (d // 2 + 1)) for d in w)
    return ok, "covariance at tree distance d = (1+d/3) 2^(-d/2), d = 0..4, exact"

def check_witnesses():
    """(C4)"""
    bad = [n for n, c in WITNESSES + GAMMA_WITNESSES if not (c() if callable(c) else c)]
    return not bad, ("%d gamma-independent and %d gamma-dependent witnesses, %d failed%s"
                     % (len(WITNESSES), len(GAMMA_WITNESSES), len(bad), (": " + "; ".join(bad)) if bad else ""))

# ------------------------------------------------------------------ V
_vtab = None
def vtab():
    """per node tau_k = k/2048: enclosures of 1/((1+tau^2) sqrt(9+8tau^2)), al, be, c0/sqrt(9+8tau^2)"""
    global _vtab
    if _vtab is None:
        _vtab = []
        for k in range(NT + max(STRIDES) + 1):
            tau = k / TDEN                  # exact
            q = 9.0 + 8.0 * tau * tau       # exact
            one = 1.0 + tau * tau           # exact
            assert Fr(q) == 9 + 8 * Fr(k, TDEN) ** 2 and Fr(one) == 1 + Fr(k, TDEN) ** 2
            sq = isqrt(pt(q))
            A = idiv(pt(1.0), imul(pt(one), sq))                       # 1/((1+tau^2) sqrt(9+8tau^2))
            al = idiv(imul(pt(9.0), pt(one)), pt(2.0 * q))             # 9(1+tau^2)/(2(9+8tau^2))
            be = idiv(isub(pt(1.0), imul(KAP, pt(one))), sq)           # (1 - kap(1+tau^2))/sqrt(9+8tau^2)
            c0A = idiv(C0, sq)                                         # c0/sqrt(9+8tau^2)
            _vtab.append((A, al, be, c0A))
    return _vtab
def V_up(a, b, stride):
    """upper bound of sup of V on [a, b]; a, b floats with 0 <= a <= b or a <= b <= 0; tau cells of width stride/2048"""
    T = vtab()
    assert a <= b and (a >= 0.0 or b <= 0.0)
    w = stride / TDEN                       # exact
    s = 0.0
    if a >= 0.0:
        a2 = max(0.0, dn(a * a))
        for k in range(0, NT, stride):
            A, al, be, _ = T[k]
            e = exp_neg(-max(0.0, dn(al[0] * a2)))[1]
            arg = max(imul(be, pt(a))[1], imul(be, pt(b))[1])
            s += up(A[1] * e) * Phi_pt(arg)[1]
        return up(up(up(THREE_OVER_PI[1] * w) * s) * (1.0 + 1e-9)) + 1e-300
    b2 = max(0.0, dn(b * b))
    for k in range(0, NT, stride):
        A, al, _, _ = T[k]
        e = exp_neg(-max(0.0, dn(al[0] * b2)))[1]
        arg = imul(T[k + stride][3], pt(b))[1]
        s += up(A[1] * e) * Phi_pt(arg)[1]
    integral = up(up(up(THREE_OVER_PI[1] * w) * s) * (1.0 + 1e-9)) + 1e-300
    ph = up(exp_neg(-b2 * 0.5)[1] * INVS2PI[1])
    diff = up(up(isub(C0, pt(1.0))[1] * (-a)) * ph)
    return up(integral + diff)

def P3_up(b):
    """P3up(b) of Lemma 5.6 (0 for b <= 0)"""
    if b <= 0.0: return 0.0
    c1 = imul(isub(idiv(pt(1.0), R1), pt(1.0)), INVS2PI)[1]                # >= (1/r1 - 1) phi(0)
    c2 = imul(pt(18.0), imul(isub(pt(1.0), R1), isub(pt(1.0), R1)))[1]     # >= 18 (1 - r1)^2
    return up(up(c1 * b) * up(1.0 - exp_neg(-up(c2 * up(b * b)))[0]))

# ------------------------------------------------------------------ P2
class P2Grid:
    def __init__(self, h, x0=X0, xmax=XMAX, ulen=ULEN):
        self.h, self.x0, self.ulen = h, x0, ulen
        hh = pt(h)
        self.s2h = imul(S2, hh)                                    # sqrt2 h
        self.du = imul(idiv(S6, pt(5.0)), hh)                      # (sqrt6/5) h
        s3h = imul(S3, hh)
        nx = int((xmax - x0) / (1.7320508 * h)) + 1
        xs = [iadd(pt(x0), imul(pt(float(k)), s3h)) for k in range(nx + 1)]
        P = [Phi_pt(x[0])[0] for x in xs], [Phi_pt(x[1])[1] for x in xs]
        self.wx = [max(0.0, dn(P[0][k + 1] - P[1][k])) for k in range(nx)]      # <= Phi(x_{k+1}) - Phi(x_k)
        y0 = imul(S3, pt(x0))
        ny = int((2 * xmax + 2 - 1.7320508 * x0) / h) + 1
        ys = [iadd(y0, imul(pt(float(l)), hh)) for l in range(ny + 1)]
        Q = [Phi_pt(y[0])[0] for y in ys], [Phi_pt(y[1])[1] for y in ys]
        self.wy = [max(0.0, dn(Q[0][l + 1] - Q[1][l])) for l in range(ny)]      # <= Phi(y_{l+1}) - Phi(y_l)
        self.nx, self.ny, self.maxnu, self.maxpath = nx, ny, 0, 0

    def lower(self, a, b):
        """lower bound of P2(a, b) = int_{u >= b} phi(u) K_a(u) du (a <= b floats)"""
        g, h = self, self.h
        base = isub(imul(S12, isub(pt(a), idiv(imul(pt(5.0), pt(b)), pt(6.0)))), imul(LAM, pt(g.x0)))
        nu = 0
        while True:
            e = iadd(pt(b), imul(pt(float(nu + 1)), g.du))
            if e[0] - b >= g.ulen: break
            nu += 1
        es = [iadd(pt(b), imul(pt(float(i)), g.du)) for i in range(nu + 1)]     # enclosures of u_i
        cuts = []
        for i in range(1, nu + 1):
            M = imul(pt(3.0), isub(pt(a), imul(R1, es[i])))
            c = idiv(imul(idiv(pt(2.0), S3), isub(M, pt(g.x0))), pt(h))       # c_i
            hi = up(-c[0] - i)                                                  # >= -c_i - i
            cuts.append(int(hi // 1) + 1)                                       # J >= cut  =>  J > -c_i - i
        # the running sum R below keeps every J >= the current cut, so the cuts must not decrease in i
        assert all(cuts[i] >= cuts[i - 1] for i in range(1, len(cuts))), "cuts not monotone"
        Jmin = int(-base[1] / g.s2h[0] // 1) - 1          # indices J < Jmin are dropped (their terms are >= 0)
        H0 = {}
        def H0f(J):                                       # lower bound of H(base + J sqrt2 h)
            v = H0.get(J)
            if v is None:
                arg = iadd(base, imul(pt(float(J)), g.s2h))[0]
                v = H_lo(arg); H0[J] = v
            return v
        nmin, nmax = -1 - nu, g.nx
        R = [0.0] * (nmax - nmin + 1)
        wy, wx, ny = g.wy, g.wx, g.ny
        cur, total = nmax + 1, 0.0
        for i in range(nu, 0, -1):
            c = max(cuts[i - 1], Jmin)
            while cur > c:
                cur -= 1
                hv = H0f(cur)
                if hv == 0.0: continue
                lo, hi = max(cur, nmin), min(cur + ny - 1, nmax)
                if lo > hi: continue
                seg = wy[lo - cur:hi - cur + 1]
                o = lo - nmin
                R[o:o + len(seg)] = [r + hv * w for r, w in zip(R[o:o + len(seg)], seg)]
            J = 0.0
            for k in range(g.nx):
                JJ = -2 * (k + 1) - i
                if JJ < c: break
                g1 = H0f(JJ)
                if g1 == 0.0: break
                J += wx[k] * g1 * R[k - 1 - i - nmin]
            mass = dn(Phi_pt(es[i][0])[0] - Phi_pt(es[i - 1][1])[1])           # <= Phi(u_i) - Phi(u_{i-1})
            if mass > 0.0: total += mass * J
        g.maxnu = max(g.maxnu, nu)
        g.maxpath = max(g.maxpath, (nmax + 1 - cur) + g.nx + nu + 8)           # roundings on one path
        return max(0.0, total * (1.0 - 1e-9) - 1e-300)

# ------------------------------------------------------------------ cover
def frs(x):                                   # exact value of a float as a string p/q
    return str(Fr(x))
def dec_up(x, digits=9):                      # decimal string >= the float x, rounded up at the given digit
    q = Fr(x) * 10 ** digits
    n = -((-q.numerator) // q.denominator)
    return ("-" if n < 0 else "") + "%d.%0*d" % (abs(n) // 10 ** digits, digits, abs(n) % 10 ** digits)

def certify(G=GAMMA_D, h=H_LATTICE, strides=STRIDES, out=None):
    """(C3).  Returns (ok, record).  G: the target (a Fraction); every comparison with it is exact."""
    out = {} if out is None else out
    le = lambda x: Fr(x) <= G                # exact comparison of a binary64 number with the target
    grid = P2Grid(h)
    tail_lo = Phi_pt(T_LO)[1]                          # t <= -5/4: G_d <= Phi(t) <= Phi(-5/4)
    tail_hi = up(up(7.0 * Phi_pt(-T_HI)[1]) / 2.0)     # t >= 2:    G_d <= (7/2) Phi(-t) <= (7/2) Phi(-2)
    ok = le(tail_lo) and le(tail_hi)
    out["tails"] = {"t<=-5/4": {"bound": "Phi(-5/4)", "value_up": dec_up(tail_lo), "ok": le(tail_lo)},
                    "t>=2": {"bound": "(7/2)*Phi(-2)", "value_up": dec_up(tail_hi), "ok": le(tail_hi)}}
    k0, k1 = int(T_LO * KDEN), int(T_HI * KDEN)
    assert Fr(k0, KDEN) == Fr(T_LO) and Fr(k1, KDEN) == Fr(T_HI)
    leaves, cnt, failed = [], {"V_bounds_stride_%d" % s: 0 for s in strides}, None
    cnt["P2_bounds"] = 0
    for k in range(k0, k1):
        stack = [(k / KDEN, (k + 1) / KDEN, 0)]            # exact dyadic ends
        while stack and failed is None:
            a, b, d = stack.pop()
            p3h = up(P3_up(b) / 2.0)
            P, done = None, False
            for s in strides:
                U0 = up(V_up(a, b, s) + p3h); cnt["V_bounds_stride_%d" % s] += 1
                if P is None:
                    if le(U0):
                        leaves.append((a, b, U0, "V%d" % s)); done = True; break
                    P = dn(1.5 * grid.lower(a, b)); cnt["P2_bounds"] += 1
                U = up(U0 - P)
                if le(U):
                    leaves.append((a, b, U, "V%d+P2" % s)); done = True; break
            if done: continue
            if d >= MAXDEPTH:
                failed = (a, b, U); break
            m = 0.5 * (a + b)                              # exact (dyadic, depth <= 7)
            assert Fr(m) * 2 == Fr(a) + Fr(b)
            stack += [(m, b, d + 1), (a, m, d + 1)]
        if failed is not None: break
    out["counts"] = cnt
    out["counts"]["P2_lattice"] = {"h": repr(h), "x_cells": grid.nx, "y_cells": grid.ny, "max_u_cells": grid.maxnu,
                                   "max_roundings_on_a_path": grid.maxpath}
    assert grid.maxpath <= 100000 and NT + 8 <= 100000
    if failed is not None:
        out["failed_leaf"] = [frs(failed[0]), frs(failed[1]), dec_up(failed[2])]
        return False, out
    # the cover: recomputed from the list of accepted leaves, in exact rationals
    contiguous = (len(leaves) > 0 and Fr(leaves[0][0]) == Fr(-5, 4) and Fr(leaves[-1][1]) == 2
                  and all(Fr(leaves[i][1]) == Fr(leaves[i + 1][0]) for i in range(len(leaves) - 1))
                  and all(x[0] < x[1] and (x[0] >= 0.0 or x[1] <= 0.0) for x in leaves))
    all_le = all(Fr(x[2]) <= G for x in leaves)
    out["cover"] = {"from": "-5/4", "to": "2", "base_intervals": k1 - k0, "leaves": len(leaves),
                    "contiguous": contiguous, "every_leaf_bound_le_target": all_le}
    im = max(range(len(leaves)), key=lambda i: leaves[i][2])
    out["largest_certified_bound"] = dec_up(leaves[im][2])
    out["argmax_interval"] = [frs(leaves[im][0]), frs(leaves[im][1])]
    meth = {}
    for x in leaves: meth[x[3]] = meth.get(x[3], 0) + 1
    out["leaves_by_method"] = meth
    blocks = []
    for j in range(k0 // 16, k1 // 16):                    # summary of the whole cover in blocks of length 1/4
        L = [x for x in leaves if Fr(j, 4) <= Fr(x[0]) and Fr(x[1]) <= Fr(j + 1, 4)]
        blocks.append([str(Fr(j, 4)), str(Fr(j + 1, 4)), len(L), dec_up(max(x[2] for x in L), 6)])
    assert sum(r[2] for r in blocks) == len(leaves)
    out["block_table_columns"] = ["from", "to", "leaves", "largest_bound"]
    out["block_table"] = blocks
    lo = max(0, min(im - 30, len(leaves) - 60))
    out["leaf_table_columns"] = ["a", "b", "bound", "method"]
    out["leaf_table"] = [[frs(x[0]), frs(x[1]), dec_up(x[2], 8), x[3]] for x in leaves[lo:lo + 60]]
    return ok and contiguous and all_le, out

def main():
    t0 = time.time()
    res = {"experiment": "certify.py: sup_t G_d(t) <= gamma_d (Proposition C of cutwidth.md)",
           "target_gamma": str(GAMMA_D), "certified": False}
    try:
        raw = sys.stdin.read()
        try:
            res["input"] = json.loads(raw) if raw.strip() else None
        except ValueError:
            res["input"] = "unparsed"
        init_constants()
        res["constants"] = {"pi": [repr(PI[0]), repr(PI[1])], "1/sqrt(2pi)": [repr(INVS2PI[0]), repr(INVS2PI[1])],
                            "sqrt2": [repr(S2[0]), repr(S2[1])], "sqrt3": [repr(S3[0]), repr(S3[1])],
                            "sqrt6": [repr(S6[0]), repr(S6[1])]}
        checks, detail, allok = {}, {}, True
        for name, fn in (("C1", check_f), ("C2", check_wave), ("C4", check_witnesses), ("C5", check_mills),
                         ("C6", check_globsub)):
            ok, msg = fn()
            checks[name] = "ok" if ok else "FAIL"; detail[name] = msg; allok = allok and ok
        ok3, _ = certify(out=res)
        checks["C3"] = "ok" if ok3 else "FAIL"
        res["checks"], res["check_details"] = checks, detail
        res["certified"] = bool(allok and ok3)
    except BaseException as e:                       # any failure is reported in the JSON, never as a traceback
        res["certified"] = False
        res["error"] = "%s: %s" % (type(e).__name__, e)
    res["python"] = sys.version.split()[0]
    res["seconds"] = round(time.time() - t0, 2)
    sys.stdout.write(json.dumps(res, separators=(",", ":")) + "\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())
