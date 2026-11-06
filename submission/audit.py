"""Finite supplementary checks for this submission (proof.md Section 7).

Each check prints PASS/FAIL; the script exits nonzero on any failure.
These checks do not replace the general proofs in the .md files; they
guard against implementation and transcription slips on finite cases.
"""
import itertools
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from family import (
    evaluate, parameters, is_prime, field_mul, field_pow,
    irreducible_prime_degree,
)

random.seed(20261103)
failures = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), name, detail)
    if not ok:
        failures.append(name)


# --- 1. Prefix test exactness vs brute-force irreducibility ---------------
def brute_irreducible(poly, k):
    """poly encodes a monic degree-k polynomial; trial division over F_2."""
    for d in range(1, k // 2 + 1):
        for b in range(1 << d, 1 << (d + 1)):
            if b & 1 == 0:
                continue  # even b divisible by Z; Z covered separately
            a = poly
            while a and a.bit_length() >= b.bit_length():
                a ^= b << (a.bit_length() - b.bit_length())
            if a == 0:
                return False
    return True


ok = True
for k in (2, 3, 5, 7, 11, 13):
    for low in range(1 << k):
        poly = (1 << k) | low
        verdict = irreducible_prime_degree(poly, k)
        brute = (low % 2 == 1) and brute_irreducible(poly, k)
        if verdict != brute:
            ok = False
            print("  mismatch k=%d low=%d exact=%s brute=%s"
                  % (k, low, verdict, brute))
            break
    if not ok:
        break
check("prefix-test exactness vs brute factoring (k<=13, exhaustive)", ok)

# --- 2. Fast evaluator vs direct subset-sum for sigma_t -------------------
N = 331
while parameters(N) is None and N < 600:
    N += 1
guard_ok = parameters(N) is not None
check("guard example length found", guard_ok, f"N={N}")
if not guard_ok:
    print()
    print("FAILURES:", failures)
    sys.exit(1)
kp, mp, rp, tp = parameters(N)


def direct_eval(bits):
    k2, m2, r2, t2 = parameters(len(bits))
    modulus = (1 << k2) | sum(b << j for j, b in enumerate(bits[:k2]))
    if not irreducible_prime_degree(modulus, k2):
        return 0
    payload = list(bits[k2:]) + [0] * (r2 * k2 - m2)
    ys = []
    for index in range(r2):
        value = sum(b << j for j, b in
                    enumerate(payload[index * k2:(index + 1) * k2]))
        w = r2 + 2 - (index + 1)
        ys.append(field_pow(value, (1 << w) - 1, modulus))
    total = 0
    for subset in itertools.combinations(range(r2), t2):
        prod = 1
        for i in subset:
            prod = field_mul(prod, ys[i], modulus)
        total ^= prod
    return total & 1


agree = True
for trial in range(12):
    bits = [random.randint(0, 1) for _ in range(N)]
    if trial % 2 == 0:  # force a passing prefix half the time
        q = 0
        while not irreducible_prime_degree((1 << kp) | q, kp):
            q += 1
        for j in range(kp):
            bits[j] = (q >> j) & 1
    if evaluate(bits) != direct_eval(bits):
        agree = False
        break
check("fast evaluator = direct subset-sum sigma_t", agree, f"N={N}")

# --- 3. Guard threshold and loss growth ------------------------------------
first_pass = None
for n in range(1, 20001):
    if parameters(n) is not None:
        first_pass = n
        break
check("guards pass for some N <= 20000", first_pass is not None,
      f"first passing N={first_pass}")

# The theorem needs R(N)/N -> 0 with R(N) = O(N^(4/5)); the ratio is NOT
# monotone (t and r jump), so we verify an explicit growth bound instead.
worst = 0.0
for n in range(first_pass, 200001):
    p = parameters(n)
    if p is None:
        continue
    kk, mm, rr, tt = p
    d_n = (kk * rr) // tt + kk * tt + rr * (rr + 1) + 1
    loss = 3 * kk + 6 * d_n
    worst = max(worst, loss / n ** 0.8)
check("loss R(N) <= 60*N^(4/5) for N<=200000", worst <= 60.0,
      f"worst R(N)/N^(4/5) = {worst:.2f}")
p200k = parameters(200000)
k2, m2, r2, t2 = p200k
d200k = (k2 * r2) // t2 + k2 * t2 + r2 * (r2 + 1) + 1
check("loss ratio documented", True,
      f"R(200000)/200000 = {(3 * k2 + 6 * d200k) / 200000:.3f}")

# --- helpers for field checks ---------------------------------------------
def find_irreducible(k):
    for low in range(1, 1 << k):
        poly = (1 << k) | low
        if irreducible_prime_degree(poly, k):
            return poly
    raise RuntimeError(f"no irreducible of degree {k}")


def lin_eval(coeffs, z, k, modulus):
    """Evaluate sum_j coeffs[j] * z^(2^j); coeffs are field elements."""
    acc, zz = 0, z
    for j in range(k):
        if coeffs[j]:
            acc ^= field_mul(coeffs[j], zz, modulus)
        zz = field_mul(zz, zz, modulus)
    return acc


def mul_by_x_plus_v(coeffs, v, modulus):
    """Multiply a K-coefficient polynomial by (X + v), v in K."""
    out = [0] * (len(coeffs) + 1)
    for i, c in enumerate(coeffs):
        out[i + 1] ^= c
        out[i] ^= field_mul(c, v, modulus)
    return out


def poly_from_roots(roots, modulus):
    p = [1]
    for v in roots:
        p = mul_by_x_plus_v(p, v, modulus)
    return p


def max_cyclic_zero_run(bits):
    """Longest cyclic run of zeros in a list that is not all zero."""
    k = len(bits)
    best = 0
    for s in range(k):
        if bits[s] != 0:
            continue
        run = 0
        while bits[(s + run) % k] == 0:
            run += 1
        best = max(best, run)
    return best


def span(basis):
    space = {0}
    for b in basis:
        new = set()
        for a in space:
            for c in (0, b):
                new.add(a ^ c)
        if len(new) == 2 * len(space):
            space = new
        else:
            return None  # dependent
    return space


# --- 4. Claim 5: zero-run pattern of image polynomials --------------------
ok4 = True
detail = ""
for k in (5, 7, 13):
    modulus = find_irreducible(k)
    for d in (2, 3):
        for _ in range(8):
            basis = [random.randrange(1, 1 << k) for _ in range(d)]
            space = span(basis)
            if space is None:
                continue
            pv = poly_from_roots(sorted(space), modulus)
            if len(pv) - 1 != (1 << d):
                ok4 = False
                detail = f"deg P_V wrong (k={k},d={d})"
                continue
            # P_V must be linearized: support only at exponents 2^j
            if not all((i & (i - 1)) == 0 and i < (1 << k)
                       for i, c in enumerate(pv) if c):
                ok4 = False
                detail = f"P_V not linearized (k={k},d={d})"
                continue
            pv_lin = [0] * k
            for i, c in enumerate(pv):
                if c:
                    pv_lin[i.bit_length() - 1] = c
            img = {lin_eval(pv_lin, z, k, modulus)
                   for z in range(1 << k)}
            if len(img) != 1 << (k - d):
                ok4 = False
                detail = f"|image| wrong (k={k},d={d})"
                continue
            q = poly_from_roots(sorted(img), modulus)
            if len(q) - 1 != (1 << (k - d)):
                ok4 = False
                detail = f"deg Q wrong (k={k},d={d})"
                continue
            if not all((i & (i - 1)) == 0 and i < (1 << k)
                       for i, c in enumerate(q) if c):
                ok4 = False
                detail = f"Q not linearized (k={k},d={d})"
                continue
            q_lin = [0] * k
            for i, c in enumerate(q):
                if c:
                    q_lin[i.bit_length() - 1] = c
            coeffs = [1 if q_lin[j] else 0 for j in range(k)]
            run_max = max_cyclic_zero_run(coeffs)
            final_run = 0
            j = k - d + 1
            while coeffs[j % k] == 0:
                final_run += 1
                j += 1
            if run_max != d - 1 or final_run != d - 1:
                ok4 = False
                detail = (f"zero runs wrong (k={k},d={d}: max={run_max}, "
                          f"final={final_run})")
check("Claim 5 zero-run pattern (random subspaces, k=5,7,13)", ok4, detail)

# --- 5. Lemma D at the smallest non-vacuous parameters --------------------
# k=37, r=10, t=3: threshold = floor(370/3)+111+110 = 344; test D=345.
K, R, T = 37, 10, 3
modulus37 = find_irreducible(K)
DN = (K * R) // T + K * T + R * (R + 1) + 1
assert DN == 345, DN


def payload_f(a_blocks):
    e = [1] + [0] * T
    for index in range(R):
        w = R + 2 - (index + 1)
        y = field_pow(a_blocks[index], (1 << w) - 1, modulus37)
        for j in range(min(index + 1, T), 0, -1):
            e[j] ^= field_mul(y, e[j - 1], modulus37)
    return e[T] & 1


ok5 = True
for _ in range(25):
    gamma = [random.randrange(1 << K) for _ in range(R)]
    # build a full-rank DN-dimensional subspace by Gaussian echelon on
    # 370-bit rows (never enumerated: 2^345 points)
    pivots = {}      # pivot position -> row
    rows = []
    while len(rows) < DN:
        b = random.randrange(1, 1 << (K * R))
        r_ = b
        while r_:
            p = r_.bit_length() - 1
            if p in pivots:
                r_ ^= pivots[p]
            else:
                pivots[p] = r_
                rows.append(b)
                break
    base = payload_f(gamma)
    nonconst = False
    for b in rows[:80]:
        shifted = [gamma[i] ^ ((b >> (i * K)) & ((1 << K) - 1))
                   for i in range(R)]
        if payload_f(shifted) != base:
            nonconst = True
            break
    if not nonconst:
        ok5 = False
        break
check("Lemma D non-vacuous spot-check (k=37,r=10,t=3,D=345 of 370)", ok5)

print()
if failures:
    print("FAILURES:", failures)
    sys.exit(1)
print("ALL CHECKS PASSED")
