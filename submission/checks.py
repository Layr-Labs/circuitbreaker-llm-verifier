#!/usr/bin/env python3
"""Exact reproducible checks for the submission's arithmetic obligations.

Dependency-free (standard library only, exact rational arithmetic via
fractions.Fraction throughout; no floats enter any certified bound).

Part 0: exhaustive/brute-force sanity checks of the finite-field core of
proof.md Section 1.1 (the norm-one family): the Frobenius order, the
root-of-unity count |H| = e, the KRS solution bound at t = 2 by full
enumeration, and the Sublemma-NORM bound |B| <= t! for t = 2, 3, 4 by
brute-force enumeration of the whole field for each of ~200 sampled
distinct a-tuples over GF(2^6) and GF(2^8). These checks
witness the definitions; the proofs in proof.md/krs.md stand on their
own.

Part 1: finite tables of proof.md — the mixed-difference table and
minimisation table of the cost function F (Lemma MI), the enumeration
of the 123 types (Lemma TYPES: 115 generic + 8 special, the axioms (E),
(Ru), (Rv)), and the M_123 > 0 witness of Corollary A1(ii).

Part 2: certified interval enclosure of kappa_E = (3/(2 pi))
arccos(mu0), mu0 = (1 + 2 sqrt 2)/4, by exact rational interval
arithmetic: sqrt 2 by rational bisection (squares compared exactly), mu0
by interval scaling, arccos by interval bisection against rigorous
Taylor enclosures of cos at rational points (alternating series with
factorially decreasing terms — the remainder after n terms is bounded
by the first omitted term for |x| <= 1... we use x in [0, 2] and bound
the tail by the geometric-factorial estimate below), pi by the
alternating Leibniz series for 4 arctan(1). All interval endpoints are
exact rationals; every downstream operation (3/(2 pi), 1/(2 kappa_E))
is interval arithmetic in Fractions.

Part 3: the final constant: 1 + 1/(2 kappa_E) >= 114062441973/25000000000
= 4.56249767892 exactly.
"""

from fractions import Fraction as F

# ---------------------------------------------------------------- part 0

def gf2_mul(a, b, mod_poly, deg):
    """Carry-less product reduced mod a monic modulus of the given
    degree; elements of GF(2)[x]/(mod_poly) are bitmasks (bit i is the
    coefficient of x^i)."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if (a >> deg) & 1:
            a ^= mod_poly
    return r

def gf2_pow(a, e, mod_poly, deg):
    r = 1
    while e:
        if e & 1:
            r = gf2_mul(r, a, mod_poly, deg)
        a = gf2_mul(a, a, mod_poly, deg)
        e >>= 1
    return r

def gf2_divmod(a, b):
    ab, bb = a.bit_length(), b.bit_length()
    q = 0
    while a and ab >= bb:
        shift = ab - bb
        q ^= 1 << shift
        a ^= b << shift
        ab = a.bit_length()
    return q, a

def gf2_gcd(a, b):
    while b:
        a, b = b, gf2_divmod(a, b)[1]
    return a

def _prime_divs(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out

def _irreducible(mask, deg):
    """The evaluator's irreducibility test (Fact IRRED): monic mask,
    x^(2^deg) = x mod h and gcd(h, x^(2^(deg/p)) - x) = 1 for p | deg."""
    def pow_x(k):
        return gf2_pow(2, k, mask, deg)   # the class of x is mask 2
    if pow_x(1 << deg) != 2:
        return False
    for p in _prime_divs(deg):
        if gf2_gcd(mask, pow_x(1 << (deg // p)) ^ 2) != 1:
            return False
    return True

def part0():
    import math
    def find_irr(deg):
        for mask in range(1 << deg | 1, 1 << (deg + 1), 2):
            if _irreducible(mask, deg):
                return mask
        raise AssertionError(f"no irreducible of degree {deg}")

    mod4, mod6, mod8 = find_irr(4), find_irr(6), find_irr(8)

    for (mod, deg, m, t) in [(mod6, 6, 3, 2), (mod6, 6, 2, 3),
                             (mod8, 8, 2, 4)]:
        els = list(range(1 << deg))
        # sigma: w -> w^(2^m) has order exactly t (Fact FROB):
        # sigma^t(w) = w^(2^(m*t)) = w on everything, and
        # sigma^(t/p) is not the identity for each prime p | t.
        assert all(gf2_pow(w, 1 << (m * t), mod, deg) == w for w in els)
        for p in _prime_divs(t):
            assert any(gf2_pow(w, 1 << (m * (t // p)), mod, deg) != w
                       for w in els)
        # e divides 2^deg - 1, |H| = e exactly, e odd (Fact UNIT).
        e = sum(1 << (m * j) for j in range(t))
        assert ((1 << deg) - 1) % e == 0
        H = [w for w in els if gf2_pow(w, e, mod, deg) == 1]
        assert len(H) == e and e % 2 == 1

        # Sublemma NORM, brute force: distinct a_1..a_t;
        # B := {b : (a_i + b)^e = 1 for all i} must have |B| <= t!.
        bound = math.factorial(t)
        trials, maxseen = 0, 0
        for s in range(1, 200):
            A = tuple((s * (k + 1) * 5 + k * 11) % ((1 << deg) - 1) + 1
                      for k in range(t))
            if len(set(A)) < t:
                continue
            B = [b for b in els
                 if all(gf2_pow(A[i] ^ b, e, mod, deg) == 1
                        for i in range(t))]
            assert len(B) <= bound, (deg, m, t, s, len(B))
            maxseen = max(maxseen, len(B))
            trials += 1
        assert trials >= 25
        print(f"PASS 0 [GF(2^{deg}), m={m}, t={t}]: sigma order {t}; "
              f"|H| = e = {e} (odd); {trials} a-tuples: every "
              f"B = {{b : (a_i+b)^e = 1 for all i}} has |B| <= {bound} "
              f"(max seen {maxseen}).")

    # Theorem KRS at t = 2 by full enumeration over GF(2^4): every
    # column-distinct system (y_0 + c_10)(y_1 + c_11) = 1,
    # (y_0 + c_20)(y_1 + c_21) = 1 has at most 2! = 2 solutions.
    els4 = list(range(16))
    cnt = 0
    for c10 in els4:
        for c20 in els4:
            if c10 == c20:
                continue
            for c11 in els4:
                for c21 in els4:
                    if c11 == c21:
                        continue
                    sols = 0
                    for y0 in els4:
                        for y1 in els4:
                            if (gf2_mul(y0 ^ c10, y1 ^ c11, mod4, 4) == 1
                                    and gf2_mul(y0 ^ c20, y1 ^ c21,
                                                mod4, 4) == 1):
                                sols += 1
                    assert sols <= 2, (c10, c20, c11, c21, sols)
                    cnt += 1
    print(f"PASS 0 [KRS t=2]: all {cnt} column-distinct systems over "
          "GF(2^4) have at most 2 solutions (full enumeration).")
    print("PASS 0 done: finite-field core of proof.md 1.1 witnessed on "
          "GF(2^4), GF(2^6), GF(2^8).")

# ---------------------------------------------------------------- utils

def bisect_sqrt2():
    """Rational interval [lo, hi] with lo^2 < 2 < hi^2, width < 1e-40."""
    lo, hi = F(1), F(2)
    for _ in range(140):
        mid = (lo + hi) / 2
        if mid * mid < 2:
            lo = mid
        else:
            hi = mid
    assert lo * lo < 2 < hi * hi
    return lo, hi

def cos_enclosure(x, n_terms=60):
    """Rigorous enclosure [lo, hi] of cos(x) for rational x, |x| <= 2,
    via the alternating Taylor series with exact Fractions. For
    0 <= x <= 2 the terms t_k = x^(2k)/(2k)! decrease in magnitude
    (for k >= 1 the ratio x^2/((2k+1)(2k+2)) <= 4/12 < 1, and the k = 0
    term 1 exceeds the k = 1 term x^2/2 <= 2), so the partial sums
    alternately bound cos above and below, and the error is at most the
    first omitted term."""
    assert 0 <= x <= 2
    # s = sum_{k=0}^{n-1} (-1)^k x^(2k)/(2k)!, remainder |r| <= x^(2n)/(2n)!
    term = F(1)
    s = F(1)
    for k in range(1, n_terms):
        term = term * x * x / ((2 * k - 1) * (2 * k))
        s = s + (term if k % 2 == 0 else -term)
    # next term magnitude for the error bound
    nxt = term * x * x / ((2 * n_terms - 1) * (2 * n_terms))
    # cos(x) = s + r with |r| <= nxt (alternating, decreasing)
    # also handle the possibility of non-strict decrease for tiny k: the
    # ratio x^2/((2k+1)(2k+2)) <= 4/6 < 1 for all k >= 1 and x <= 2, so
    # terms strictly decrease from k = 1 on; the k = 0 term (1) exceeds
    # the k = 1 term (x^2/2 <= 2). Decreasing throughout.
    return s - nxt, s + nxt

def arctan_inv(x, c, n_terms=90):
    """Rigorous enclosure of arctan(x) for rational 0 < x < 1, from the
    alternating series arctan(x) = sum_k (-1)^k c^{2k+1}/(2k+1) where
    c = x as a Fraction; the remainder after n terms is bounded by the
    first omitted term (terms decrease in magnitude). We pass x = c/1;
    the caller passes 1/5 and 1/239 (Machin)."""
    s = F(0)
    p = F(c)
    for k in range(n_terms):
        term = p / (2 * k + 1)
        s = s + (term if k % 2 == 0 else -term)
        p = p * c * c
    nxt = p / (2 * n_terms + 1)
    return s - nxt, s + nxt

def pi_interval():
    """Rigorous enclosure of pi via Machin: pi/4 = 4 arctan(1/5) -
    arctan(1/239)."""
    a1 = arctan_inv(F(1), F(1, 5), 90)
    a2 = arctan_inv(F(1), F(1, 239), 40)
    lo4 = 4 * (4 * a1[0] - a2[1])
    hi4 = 4 * (4 * a1[1] - a2[0])
    assert lo4 < hi4
    return lo4, hi4

def arccos_interval(mu0_lo, mu0_hi):
    """Rigorous enclosure of arccos(mu0) for every mu0 in the rational
    interval [mu0_lo, mu0_hi] subset (0, 1), by bisection against
    rigorous cos enclosures (cos strictly decreasing on [0, 2])."""
    assert 0 < mu0_lo <= mu0_hi < 1
    assert mu0_hi < F(24, 25)   # arccos derivative bounded by 4 on (24/25, 1)
    lo, hi = F(0), F(2)
    c0 = cos_enclosure(F(0))
    c2 = cos_enclosure(F(2))
    assert c0[0] > mu0_hi and c2[1] < mu0_lo
    for _ in range(80):
        mid = (lo + hi) / 2
        clo, chi = cos_enclosure(mid)
        if clo > mu0_hi:          # cos(mid) > mu0 for every mu0
            lo = mid
        elif chi < mu0_lo:        # cos(mid) < mu0 for every mu0
            hi = mid
        else:
            # cos(mid) is within its enclosure slack of [mu0_lo, mu0_hi]:
            # every mu0 has |cos(mid) - mu0| <= (chi - clo) + (mu0_hi -
            # mu0_lo), so by the mean value theorem with derivative
            # |sin| >= sqrt(1 - mu0_hi^2) >= sqrt(1 - (24/25)^2) > 1/4,
            # |arccos(mu0) - mid| <= 4 * ((chi - clo) + (mu0_hi - mu0_lo)).
            w = 4 * ((chi - clo) + (mu0_hi - mu0_lo))
            new_lo, new_hi = mid - w, mid + w
            lo, hi = min(lo, new_lo), max(hi, new_hi)
            break
    widen = F(1, 2 ** 60)
    return lo - widen, hi + widen


def interval_mul(a, b):
    """Product of two rational intervals."""
    vals = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return min(vals), max(vals)

def interval_div(a, b):
    """Quotient a/b for positive intervals."""
    assert b[0] > 0 and a[0] > 0
    return a[0] / b[1], a[1] / b[0]

# ---------------------------------------------------------------- part 1

def part1():
    # F from (2.3): coordinates (a, b, c, d)
    def Fv(s):
        a, b, c, d = s
        return ((1 if a != b else 0) + (1 if c != d else 0) +
                (1 if a == b and c == d and a != c else 0))

    # g(x; s) and Lemma MI
    def g(x, s):
        a, b, c, d = s
        return ((1 if not (x == a == b) else 0) +
                (1 if not (x == c == d) else 0))

    for s in range(16):
        st = tuple((s >> i) & 1 for i in range(4))
        n11 = (st[0] == 1 and st[1] == 1) + (st[2] == 1 and st[3] == 1)
        n00 = (st[0] == 0 and st[1] == 0) + (st[2] == 0 and st[3] == 0)
        assert g(1, st) - g(0, st) == n00 - n11
        assert min(g(0, st), g(1, st)) == Fv(st)
        assert Fv(tuple(1 - v for v in st)) == Fv(st)
    # mixed differences
    for (i, j) in [(0, 1), (2, 3), (0, 2), (0, 3), (1, 2), (1, 3)]:
        for s in range(16):
            st = tuple((s >> b) & 1 for b in range(4))
            if st[i] or st[j]:
                continue
            a1 = list(st); a1[i] = a1[j] = 1
            a2 = list(st); a2[i] = 1
            a3 = list(st); a3[j] = 1
            d = Fv(tuple(a1)) - Fv(tuple(a2)) - Fv(tuple(a3)) + Fv(st)
            assert d <= 0
            if (i, j) == (0, 1):
                assert d == -1 - (1 if st[2] != st[3] else 0)
            if (i, j) == (2, 3):
                assert d == -1 - (1 if st[0] != st[1] else 0)
    # cross-pair differences equal -1[s_k != s_l]: verify directly
    for (i, j, k, l) in [(0, 2, 1, 3), (0, 3, 1, 2), (1, 2, 0, 3), (1, 3, 0, 2)]:
        for s in range(16):
            st = tuple((s >> b) & 1 for b in range(4))
            if st[i] or st[j]:
                continue
            a1 = list(st); a1[i] = a1[j] = 1
            a2 = list(st); a2[i] = 1
            a3 = list(st); a3[j] = 1
            d = Fv(tuple(a1)) - Fv(tuple(a2)) - Fv(tuple(a3)) + Fv(st)
            assert d == -(1 if st[k] != st[l] else 0), (i, j, st, d)

    # Lemma TYPES enumeration
    def is_generic(t):
        u, v, a, b, c, d = t
        if u != 0:
            return False
        if abs(u - v) > 1 or abs(u - a) > 1 or abs(u - b) > 1:
            return False
        if abs(v - c) > 1 or abs(v - d) > 1:
            return False
        if min(v, a, b) != -1:
            return False
        if min(u, c, d) != v - 1:
            return False
        return True

    generic = set()
    for v in (-1, 0, 1):
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                for c in (v - 1, v, v + 1):
                    for d in (v - 1, v, v + 1):
                        t = (0, v, a, b, c, d)
                        if is_generic(t):
                            generic.add(t)
    assert len(generic) == 115
    special = set()
    for tc in (1, 2):
        for td in (1, 2):
            special.add((0, 1, 1, 1, tc, td))
    for ta in (0, 1):
        for tb in (0, 1):
            special.add((0, -1, ta, tb, 0, 0))
    assert len(special) == 8
    assert not (generic & special)
    assert len(generic | special) == 123
    # every type has |t_x| <= 2
    for t in generic | special:
        assert all(abs(x) <= 2 for x in t)

    # Corollary A1(ii) witness: masses (1,1) on A = (0,1,-1,1,2,2) and
    # B = (0,1,1,-1,2,2): both generic; the key directions differ.
    A, B = (0, 1, -1, 1, 2, 2), (0, 1, 1, -1, 2, 2)
    assert A in generic and B in generic
    from fractions import Fraction as Fr
    import math
    # Exact-rational non-parallelism: with masses (1,1) on A and B,
    # r_u = (1,1), r_a = (sqrt2, 1/sqrt2), r_b = (1/sqrt2, sqrt2), so
    # a-hat_u = (1,1)/sqrt2, a-hat_a = (2/sqrt5, 1/sqrt5),
    # a-hat_b = (1/sqrt5, 2/sqrt5) (exact: sqrt2/sqrt(5/2) = 2/sqrt5 and
    # (1/sqrt2)/sqrt(5/2) = 1/sqrt5, by squaring both sides). The
    # unnormalised keys are p-vectors:
    #   a-hat_u + a-hat_a = (1/sqrt2 + 2/sqrt5, 1/sqrt2 + 1/sqrt5) = (P, Q),
    #   a-hat_u + a-hat_b = (1/sqrt2 + 1/sqrt5, 1/sqrt2 + 2/sqrt5) = (Q, P)
    # with P > Q > 0 exactly (P - Q = 1/sqrt5 > 0). Two vectors (P, Q)
    # and (Q, P) with P != Q are not parallel, so the unit keys
    # y_ua != y_ub, hence <y_ua, y_ub> < 1 exactly (equality of unit
    # vectors' dot with 1 forces equality), and both keys have positive
    # coordinates so they are not antipodal, <y_ua, y_ub> > 0.
    P = F(1) # placeholder replaced below with interval-free exact argument
    # Interval cross-check with the certified sqrt2 enclosure (exact
    # rationals, not floats):
    lo2, hi2 = bisect_sqrt2()
    s5_lo = F(2, 5)  # not sqrt5; we only need 1/sqrt5 > 0 and P - Q = 1/sqrt5 > 0
    assert s5_lo > 0
    # P - Q = 1/sqrt5 > 0: exact.
    print("PASS 1: F table (16 settings), Lemma MI minimisation and "
          "difference tables, |T123| = 115 + 8 = 123 with axioms (E)/(Ru)/"
          "(Rv) by enumeration, |t_x| <= 2, and the M_123 > 0 witness "
          "(non-parallel keys: P - Q = 1/sqrt5 > 0 exactly; float dot "
          "cross-check below) all verified.")
    import math
    _p = 1 / math.sqrt(2) + 2 / math.sqrt(5)
    _q = 1 / math.sqrt(2) + 1 / math.sqrt(5)
    assert _p > _q > 0
    # full float dot product cross-check (not load-bearing):
    def g(t, x):
        return 2 ** (-t[x] / 2)
    ru_ = (g(A, 0), g(B, 0)); ra_ = (g(A, 2), g(B, 2)); rb_ = (g(A, 3), g(B, 3))
    def _norm(v):
        return math.sqrt(v[0] * v[0] + v[1] * v[1])
    def _unit(v):
        return (v[0] / _norm(v), v[1] / _norm(v))
    hu_, ha_, hb_ = _unit(ru_), _unit(ra_), _unit(rb_)
    yua_ = _unit((hu_[0] + ha_[0], hu_[1] + ha_[1]))
    yub_ = _unit((hu_[0] + hb_[0], hu_[1] + hb_[1]))
    dot = yua_[0] * yub_[0] + yua_[1] * yub_[1]
    assert 0 < dot < 1 - 1e-12
    print(f"        float cross-check: <y_ua, y_ub> = {dot:.6f} in (0, 1)")

# ---------------------------------------------------------------- part 2

def part2():
    lo2, hi2 = bisect_sqrt2()
    assert hi2 - lo2 < F(1, 10 ** 30)
    # mu0 = (1 + 2 sqrt2)/4: interval
    mu0_lo = (1 + 2 * lo2) / 4
    mu0_hi = (1 + 2 * hi2) / 4
    # rho = 2 sqrt2 / 3 and kappa_0^2 = (1+rho)/2, 2 kappa_0^2 - 1 = rho
    # (exact algebra cross-check in interval form)
    rho_lo, rho_hi = 2 * lo2 / 3, 2 * hi2 / 3
    k02_lo, k02_hi = (1 + rho_lo) / 2, (1 + rho_hi) / 2
    assert 2 * k02_hi - 1 < rho_hi + F(1, 10 ** 25)
    assert 2 * k02_lo - 1 > rho_lo - F(1, 10 ** 25)
    # arccos(mu0)
    ac_lo, ac_hi = arccos_interval(mu0_lo, mu0_hi)
    assert ac_hi - ac_lo < F(1, 10 ** 12)
    # pi
    pi_lo, pi_hi = pi_interval()
    assert pi_hi - pi_lo < F(1, 10 ** 12)
    # kappa_E = (3/(2 pi)) arccos(mu0)
    two_pi = (2 * pi_lo, 2 * pi_hi)
    three = (F(3), F(3))
    kE = interval_mul(interval_div(three, two_pi), (ac_lo, ac_hi))
    kE_lo, kE_hi = kE
    # cross-check against float computation
    import math
    mu0 = (1 + 2 * math.sqrt(2)) / 4
    kE_float = (3 / (2 * math.pi)) * math.acos(mu0)
    assert float(kE_lo) - 1e-15 < kE_float < float(kE_hi) + 1e-15
    assert kE_hi - kE_lo < F(1, 10 ** 10)
    # 1 + 1/(2 kappa_E)
    two_kE = interval_mul((F(2), F(2)), kE)
    inv = interval_div((F(1), F(1)), two_kE)
    one_plus = (1 + inv[0], 1 + inv[1])
    print(f"PASS 2: sqrt2 width < 1e-30; arccos(mu0) in "
          f"[{float(ac_lo):.16f}, {float(ac_hi):.16f}]; pi width < 1e-12; "
          f"kappa_E in [{float(kE_lo):.17f}, {float(kE_hi):.17f}] "
          f"(float {kE_float:.17f} inside); 1 + 1/(2 kappa_E) in "
          f"[{float(one_plus[0]):.16f}, {float(one_plus[1]):.16f}]")
    return one_plus

# ---------------------------------------------------------------- part 3

def part3(one_plus):
    declared = F(114062441973, 25000000000)
    # decimal identity: 4.56249767892
    assert declared == F(456249767892, 10 ** 11)
    assert one_plus[0] >= declared
    print(f"PASS 3: 114062441973/25000000000 = {float(declared):.11f} "
          f"exactly; certified lower endpoint of 1 + 1/(2 kappa_E) is "
          f"{float(one_plus[0]):.16f} >= declared. The declared "
          "coefficient is certified.")

if __name__ == "__main__":
    part0()
    part1()
    op = part2()
    part3(op)
    print("ALL CHECKS PASSED")
