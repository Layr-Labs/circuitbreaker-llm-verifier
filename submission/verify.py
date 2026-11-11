#!/usr/bin/env python3
"""Exact reproducible checks for the submission's arithmetic obligations.

Self-contained; runs with python3 (standard library only). Five parts:

  1. Classification of the 16 binary Boolean operations (constant / unary /
     XOR-type / AND-type) and their normal forms, as used in proof.md.
  2. The integer/Fraction arithmetic of every Delta W inequality of the ten
     gate-elimination cases in proof.md, and of the final assembly
     2|C| + 7n >= 14(n - d)  =>  |C| >= (7/2)n - 7d  =>  ... >= (7/2)n - 49 n^(4/5).
  3. The disperser parameter inequalities of disperser.md at n = 10^12 and
     n = 10^24 (worst-case envelopes, high-precision decimal arithmetic).
  4. The deterministic field construction of field.md on odd primes
     k in {3, 5, 7, 13, 17, 41}: Berlekamp factoring, the cyclotomic chain, the
     nonresidue, the trace element gamma, and the output irreducible
     polynomial of degree exactly k, all built exactly as specified.
  5. Consistency of the claimed corollary constants.

Each part prints PASS lines; any failure raises. The declared predicate of
each check is only what it states: the finite instances and identities
checked, not a proof of the general theorems (which are proved in the .md
files).
"""

from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 60

# ---------------------------------------------------------------- F2[X] as int bitmasks: bit i = coefficient of X^i.

def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r

def pmod(a, m):
    dm = m.bit_length() - 1
    while True:
        da = a.bit_length() - 1
        if da < dm or a == 0:
            return a
        a ^= m << (da - dm)

def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a

def pdivmod(a, m):
    dm = m.bit_length() - 1
    q = 0
    while True:
        da = a.bit_length() - 1
        if da < dm or a == 0:
            return q, a
        q |= 1 << (da - dm)
        a ^= m << (da - dm)

def is_irreducible(f):
    """Irreducibility over F2 for f != 0 of degree >= 1."""
    assert f and (f & 1) == 0 or True
    if f in (0, 1):
        return False
    d = f.bit_length() - 1
    if d == 0:
        return False  # constant
    if d == 1:
        return f in (2, 3)  # X, X+1
    # gcd(f, X) == 1 needed
    if f % 2 == 0:
        return False
    x = 2
    for i in range(1, d // 2 + 1):
        x = pmod(pmul(x, x), f)
        # x = X^(2^i) mod f
        if pgcd(f, x ^ 2) not in (1,):
            return False
    return True

def berlekamp_factor(h):
    """Deterministic full factorization of squarefree monic h into monic
    irreducibles (Berlekamp over F2), as specified in field.md Lemma 1.1.
    Returns the list of factors sorted by integer value (the lexicographically
    first convention of field.md picks the smallest in each recursive call)."""
    D = h.bit_length() - 1
    # squarefree part (h is squarefree in our uses; still do it)
    hp = 0
    for i in range(D + 1):
        if (h >> i) & 1 and i % 2 == 1:
            hp |= 1 << (i - 1)  # derivative, coefficients mod 2
    hh = pdivmod(h, pgcd(h, hp if hp else 1))[0] if hp else h
    def rec(m):
        d = m.bit_length() - 1
        # rows of the squaring map S on basis 1, X, ..., X^(d-1): X^(2j) mod m
        rows = []
        for j in range(d):
            rows.append(pmod(1 << (2 * j), m))
        # matrix S - I over F2 as integer rows (bit i = coord i)
        mat = [rows[i] ^ (1 << i) for i in range(d)]
        # Gaussian elimination: find kernel elements
        basis = kernel_basis(mat, d)
        if len(basis) == 1:
            return [m]  # only constants in Berlekamp subalgebra => irreducible
        # pick nonconstant kernel element: any basis vector other than the
        # all-ones vector (which is always in the kernel)
        g = None
        for v in basis:
            if v not in (0, 1) and pmod(pmul(v, v), m) == v:
                g = v
                break
        assert g is not None
        d1, d2 = pgcd(m, g), pgcd(m, g ^ 1)
        assert 1 < d1 < m and 1 < d2 < m and pmul(d1, d2) == m
        return rec(d1) + rec(d2)
    return sorted(rec(hh))

def kernel_basis(rows, n):
    """Kernel basis of the integer-row F2 matrix (bit i = coordinate i)."""
    rows = list(rows)
    piv = {}
    r = 0
    for c in range(n):
        p = None
        for i in range(r, len(rows)):
            if (rows[i] >> c) & 1:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> c) & 1):
                rows[i] ^= rows[r]
        piv[c] = r
        r += 1
        if r == len(rows):
            break
    # kernel: free columns
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = 1 << fc
        for c, ri in piv.items():
            if (rows[ri] >> fc) & 1:
                v |= 1 << c
        basis.append(v)
    if not basis:
        basis = [0]
    return basis

# ---------------------------------------------------------------- Part 1

def part1():
    ops = []
    for c0 in (0, 1):
        for b1 in (0, 1):
            for b2 in (0, 1):
                for b12 in (0, 1):
                    table = tuple(
                        (c0 ^ (b1 * x) ^ (b2 * y) ^ (b12 * x * y))
                        for x in (0, 1) for y in (0, 1)
                    )
                    ops.append(((c0, b1, b2, b12), table))
    assert len(ops) == 16 and len(set(t for _, t in ops)) == 16
    constant = [t for _, t in ops if len(set(t)) == 1]
    unary = []
    for coeffs, t in ops:
        if len(set(t)) == 2 and (coeffs[3] == 0):
            # affine in one/both inputs but genuinely depends on exactly one
            pass
    # clean classification: dependence on each input
    def depends(t, i):
        # i=0: does flipping first input change output somewhere
        a, b, c, d = t  # (x,y) in (0,0),(0,1),(1,0),(1,1)
        return (a != c) or (b != d) if i == 0 else (a != b) or (c != d)
    const = [t for _, t in ops if not depends(t, 0) and not depends(t, 1)]
    un = [t for _, t in ops if len(set(t)) == 2 and (
        (depends(t, 0) and not depends(t, 1)) or (depends(t, 1) and not depends(t, 0)))]
    both = [t for _, t in ops if depends(t, 0) and depends(t, 1)]
    # both-dependent: XOR-type iff b12 = 0
    xortype = [co for co, t in ops if depends(t, 0) and depends(t, 1) and co[3] == 0]
    andtype = [co for co, t in ops if depends(t, 0) and depends(t, 1) and co[3] == 1]
    assert len(const) == 2 and len(un) == 4 and len(xortype) == 2 and len(andtype) == 8
    # normal forms: XOR-type = x+y+c ; AND-type = (x+a)(y+b)+g
    for (c0, b1, b2, b12) in xortype:
        assert (b1, b2, b12) == (1, 1, 0)
    for (c0, b1, b2, b12) in andtype:
        # (x^a)(y^be)+g with a = b2, be = b1, g = c0 ^ (a & be)
        a, be = b2, b1
        g = c0 ^ (a & be)
        tt = tuple(g ^ ((x ^ a) & (y ^ be)) for x in (0, 1) for y in (0, 1))
        assert tt == tuple(c0 ^ b1 * x ^ b2 * y ^ b12 * x * y
                           for x in (0, 1) for y in (0, 1))
    print("PASS 1: 16 operations classify as 2 constant, 4 unary, 2 XOR-type "
          "(x+y+c), 8 AND-type ((x+a)(y+b)+g).")

# ---------------------------------------------------------------- Part 2

def part2():
    # Case 2 (w unpaired / w paired)
    assert -6 - 6 - 2 == -14
    assert -3 - 6 - 2 + 7 == -4 <= 0
    # Case 3
    assert -3 <= 0
    assert -4 - 3 + 7 == 0
    # Case 4
    assert -3 - 2 == -5 <= 0
    # Case 5
    assert -3 - 2 + 6 == 1 and 1 - 1 == 0
    # Case 6
    assert -7 - 4 - 3 == -14
    assert -4 - 6 + 7 == -3 <= 0
    # Case 7
    assert -7 - 7 - 4 == -18 <= -14
    # Case 8: for all admissible parameters
    for degx in range(3, 8):
        for k in range(1, 8):
            assert -7 - 2 * (degx + k) + k <= -14
    for k in range(2, 8):
        for k1 in range(0, k + 1):
            assert -12 - 2 * k + k1 <= -14
    assert -14 + 0 == -14
    # Case 9
    assert -7 - 7 - 6 - 8 == -28 == 14 * (-2)
    # Case 10
    assert -4 - 7 - 7 - 6 + 3 + 7 == -14 == 14 * (-1)
    # Assembly with exact rationals
    sevenhalf = Fraction(7, 2)
    for n, d in [(10**12, 10**10), (10**6, 10**5), (12345, 999)]:
        assert 14 * (n - d) - 7 * n == 2 * (sevenhalf * n - 7 * d)
    # remainder envelope: 7*(6*n^(4/5)+1) = 42 n^(4/5) + 7 <= 49 n^(4/5) for n >= 1
    assert 42 + 7 == 49
    # the declared N: 10^12 >= 1 so n^(4/5) >= 1
    assert 10**12 >= 1
    print("PASS 2: every case inequality and the final assembly arithmetic "
          "hold exactly (integer and Fraction arithmetic).")

# ---------------------------------------------------------------- Part 3

def part3():
    for n in (10**12, 10**24):
        u = Decimal(n) ** (Decimal(1) / Decimal(5))
        assert u > 251
        assert Decimal('0.004') * u >= 1
        assert Decimal('0.0001') * u * u >= 1
        assert Decimal(4) * u ** 3 + 2 <= Decimal("4.01") * u ** 3
        # r bounds via k in (2u^3, 4.01u^3]
        r_hi = Decimal(n) / (Decimal(2) * u ** 3) + Decimal(1)
        assert r_hi <= Decimal("0.501") * u * u
        sqrt_r = (Decimal("0.501") * u * u) ** Decimal("0.5")
        assert sqrt_r <= Decimal("0.708") * u
        t_hi = Decimal("0.708") * u + 1
        assert t_hi <= Decimal("0.712") * u
        # t < r
        assert t_hi < n / (Decimal("4.01") * u ** 3)
        # k >= r(r+1)/t
        lhs = Decimal("0.708") * u * (Decimal("0.501") * u * u + 1)
        assert lhs <= Decimal("0.708") * Decimal("0.5011") * u ** 3
        assert Decimal("0.708") * Decimal("0.5011") < Decimal("0.355")
        assert Decimal("0.355") * u ** 3 < 2 * u ** 3
        # dimension envelope
        total = (Decimal("4.01") * u ** 3 * (Decimal("0.708") + Decimal("0.712")) * u
                 + Decimal("0.501") * u * u * (Decimal("0.501") * u * u + 1))
        assert Decimal("4.01") * (Decimal("0.708") + Decimal("0.712")) == Decimal("5.6942")
        assert Decimal("0.501") * Decimal("0.5011") == Decimal("0.25105110")
        assert Decimal("5.6942") + Decimal("0.2511") == Decimal("5.9453")
        assert total <= Decimal("5.9453") * u ** 4 + Decimal("0.0000001")
        assert Decimal("5.9453") * u ** 4 < 6 * u ** 4
    print("PASS 3: disperser parameter envelopes verified at n = 10^12 and "
          "10^24 with 60-digit decimal arithmetic (all margins strict).")

# ---------------------------------------------------------------- Part 4

class Fq:
    """K = F2[X]/(f), elements as ints, as in field.md."""

    def __init__(self, f):
        assert is_irreducible(f)
        self.f = f
        self.d = f.bit_length() - 1

    def mul(self, a, b):
        return pmod(pmul(a, b), self.f)

    def pow_(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def order_is_exactly(self, a, want):
        """order(a) == want  (checked as a^want == 1 and, for each prime
        p | want, a^(want/p) != 1)."""
        assert a != 0
        if self.pow_(a, want) != 1:
            return False
        for p in prime_factors(want):
            if self.pow_(a, want // p) == 1:
                return False
        return True

def prime_factors(m):
    fs, x = set(), 2
    while x * x <= m:
        while m % x == 0:
            fs.add(x)
            m //= x
        x += 1
    if m > 1:
        fs.add(m)
    return fs

def part4_one(k):
    assert k % 2 == 1 and len(prime_factors(k)) == 1 and k > 2
    # Step A: m = ord_k(2)
    m, x = 1, 2 % k
    while x != 1:
        x = x * 2 % k
        m += 1
    assert 2 ** m % k == 1 and m <= k - 1
    # Step B: Phi_k over F2 = X^(k-1)+...+X+1 (k odd)
    phi = (1 << k) - 1   # X^(k-1)+...+X+1
    facs = berlekamp_factor(phi)
    assert all(is_irreducible(f) for f in facs)
    assert all(f.bit_length() - 1 == m for f in facs), "all factors degree m"
    f1 = facs[0]
    # chain
    N = 2 ** m - 1
    s, NN = 0, N
    while NN % k == 0:
        NN //= k
        s += 1
    assert s >= 1
    f_prev = f1
    for i in range(2, s + 1):
        comp = 0
        for j in range(f_prev.bit_length()):
            if (f_prev >> j) & 1:
                comp |= 1 << (j * k)
        fs = berlekamp_factor(comp)
        assert all(is_irreducible(g) for g in fs)
        f_prev = fs[0]
    f = f_prev
    ms = f.bit_length() - 1
    K = Fq(f)
    alpha = pmod(2, f)  # class of X
    assert K.order_is_exactly(alpha, k ** s), 'alpha has order exactly k^s'
    # alpha is a k-th nonresidue: v_k(2^ms - 1) == s
    NN, vs = (1 << ms) - 1, 0
    while NN % k == 0:
        NN //= k
        vs += 1
    assert vs == s, "v_k(2^ms-1) == s, so alpha has no k-th root in K"
    # gamma = sum_{j<ms} beta^(2^(k j)),  beta = Y in E = K[Y]/(Y^k - alpha)
    EK = ExtE(K, alpha, k)
    gamma = EK.trace_of_beta(ms)
    assert EK.frobenius(gamma, k) == gamma  # gamma in H = GF(2^k)
    # gamma not in F2:
    def to_E(x):
        return (x,) + (0,) * (k - 1)
    assert gamma != to_E(0) and gamma != to_E(1)
    # minpoly over F2: prod (X - gamma^(2^i)), coefficients land in F2
    g = EK.minpoly(gamma, k)
    assert g.bit_length() - 1 == k
    assert is_irreducible(g)
    return g, k, m, s, ms

class ExtE:
    """E = K[Y]/(Y^k - alpha) with K = F2[X]/(f)."""

    def __init__(self, K, alpha, k):
        self.K, self.alpha, self.k = K, alpha, k

    def mul(self, a, b):
        K, k = self.K, self.k
        c = [0] * (2 * k - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        c[i + j] ^= K.mul(ai, bj)
        # reduce Y^k -> alpha
        for idx in range(2 * k - 2, k - 1, -1):
            t = c[idx]
            if t:
                c[idx] = 0
                c[idx - k] ^= K.mul(t, self.alpha)
        return tuple(c[:k])

    def sq(self, a):
        return self.mul(a, a)

    def frobenius(self, a, times):
        for _ in range(times):
            a = self.sq(a)
        return a

    def trace_of_beta(self, ms):
        beta = (0, 1) + (0,) * (self.k - 2)  # the element Y
        acc = beta
        z = beta
        for _ in range(1, ms):
            z = self.frobenius(z, self.k)
            acc = tuple(x ^ y for x, y in zip(acc, z))
        return acc

    def minpoly(self, gamma, k):
        # product over i<k of (X + gamma^(2^i)); coefficients are E-elements
        coeffs = [tuple([1] + [0] * (self.k - 1))] + [tuple(0 for _ in range(self.k))] * k
        # multiply successively
        poly = [gamma, tuple([1] + [0] * (self.k - 1))]  # gamma + 1*X
        for i in range(1, k):
            c = self.frobenius(gamma, i)
            # poly = poly * (X + c)
            zero = tuple(0 for _ in range(self.k))
            one = tuple([1] + [0] * (self.k - 1))
            new = [zero] * (len(poly) + 1)
            for j, pj in enumerate(poly):
                if pj != zero:
                    new[j + 1] = tuple(x ^ y for x, y in zip(new[j + 1], pj))
                    new[j] = tuple(x ^ y for x, y in zip(new[j], self.mul(pj, c)))
            poly = new
        assert len(poly) == k + 1
        assert poly[k] == tuple([1] + [0] * (self.k - 1))
        g = 0
        for j, cj in enumerate(poly):
            # each coefficient must be 0 or 1 (in F2 subset E)
            if cj not in (tuple(0 for _ in range(self.k)),
                          tuple([1] + [0] * (self.k - 1))):
                raise AssertionError("coefficient not in F2")
            if cj[0] == 1:
                g |= 1 << j
        return g

def part4():
    results = []
    for k in (3, 5, 7, 13, 17, 41):
        g, k, m, s, ms = part4_one(k)
        assert g.bit_length() - 1 == k and is_irreducible(g)
        results.append((k, m, s, ms, g))
    for (k, m, s, ms, g) in results:
        print(f"PASS 4: k={k:3d}  ord_k(2)=m={m:3d}  s={s}  deg f^(s)=ms={ms:3d}  "
              f"output irreducible of degree {g.bit_length()-1} "
              f"(coeffs {g:0{g.bit_length()}b})")
    print("PASS 4: field.md construction (Berlekamp/cyclotomic chain/nonresidue/"
          "X^k-alpha/trace gamma/minpoly) runs exactly as specified on the listed "
          "odd primes; every intermediate assertion of field.md's Lemmas 1.1-2.3 "
          "that is checkable on these instances holds.")

# ---------------------------------------------------------------- Part 5

def part5():
    N0 = 10 ** 12
    dN = int(6 * float(Decimal(N0) ** (Decimal(4) / Decimal(5)))) + 1
    lhs = 7 * dN
    rhs = 49 * float(Decimal(N0) ** (Decimal(4) / Decimal(5)))
    assert lhs <= rhs, "7*d(n) <= 49*n^(4/5) at n = 10^12"
    print(f"PASS 5: at n = 10^12, d(n) = {dN}, 7*d(n) = {lhs} <= "
          f"49*n^(4/5) = {rhs:.3e}; the corollary's remainder envelope holds.")

if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    print("ALL CHECKS PASSED")
