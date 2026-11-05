"""Exact evaluator for the BK/Shoup family of proof.md (small lengths).

Implements, with exact integer arithmetic, the family definition of
proof.md Section 1: for n >= 1, k is the least prime strictly above
2*n^(3/5) (decided exactly via k^5 > 32*n^3), r = ceil(n/k),
t = ceil(sqrt(r)); the input is zero-padded to r blocks of k bits, block
i read as an element phi_i of K = F_2[Z]/(q_k) (first bit = constant
coefficient), y_i = phi_i^(2^(i+1)-1), and f_n(x) is the constant
coefficient of the elementary symmetric polynomial Sym_r^t(y_1..y_r),
computed by the descending recurrence of proof.md Section 2.

The proof fixes the irreducible polynomial q_k as the output of the
deterministic construction of proof.md Sections 1 and 1A (V. Shoup,
Math. Comp. 54 (1990), 435-447, specialized to p = 2 and prime degree
and reproduced with proof there; deterministic, time polynomial in k).
This file does NOT implement that construction; for these finite
evaluations it fixes the concrete deterministic stand-in rule
"numerically smallest monic irreducible of degree k over F_2". Every
quantity audit.py checks (parameters, threshold arithmetic, LP and
ledger constants, evaluator well-definedness) is independent of which
fixed irreducible is chosen, so the stand-in is sound for the audit's
scope; the polynomial-time premise of the proof rests on proof.md
Section 1A, not on this file.

Scope: exact evaluation is exponential in n if all inputs are enumerated
and is intended for small n only; nothing here proves any asymptotic or
all-length statement.
"""
from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def parameters(n):
    """(k, r, t) of proof.md Section 1, decided with integers only."""
    if n < 1:
        raise ValueError("parameters are defined for n >= 1")
    k = 2
    while not (is_prime(k) and k ** 5 > 32 * n ** 3):
        k += 1
    r = -(-n // k)
    t = isqrt(r)
    if t * t < r:
        t += 1
    return k, r, t


# --- F_2[Z] arithmetic: polynomials are bit-packed ints (bit j = coeff of Z^j).

def _pmul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
    return out


def _pmod(a, p):
    dp = p.bit_length() - 1
    while a and a.bit_length() - 1 >= dp:
        a ^= p << (a.bit_length() - 1 - dp)
    return a


def _is_irreducible(p):
    d = p.bit_length() - 1
    if d < 1:
        return False
    # trial division by every polynomial of degree 1..d/2
    for deg in range(1, d // 2 + 1):
        for cand in range(1 << deg, 1 << (deg + 1)):
            if cand & 1 == 0:
                continue
            if _pmod(p, cand) == 0:
                return False
    return True


def irreducible(k):
    """Numerically smallest monic irreducible of degree k over F_2."""
    for cand in range((1 << k) | 1, 1 << (k + 1), 2):
        if _is_irreducible(cand):
            return cand
    raise AssertionError("no irreducible found")


# --- K = F_2[Z]/(q_k) arithmetic.

def _fmul(a, b, q):
    return _pmod(_pmul(a, b), q)


def _fpow(a, e, q):
    out = 1
    while e:
        if e & 1:
            out = _fmul(out, a, q)
        a = _fmul(a, a, q)
        e >>= 1
    return out


def evaluate(n, x):
    """f_n(x) for x an integer with 0 <= x < 2^n. Returns 0 or 1."""
    if n == 0:
        return 0
    k, r, t = parameters(n)
    q = irreducible(k)
    mask = (1 << k) - 1
    ys = []
    for i in range(1, r + 1):
        block = (x >> ((i - 1) * k)) & mask
        ys.append(_fpow(block, (1 << (i + 1)) - 1, q))
    # descending recurrence E_j <- E_j + y_i * E_{j-1} (proof.md Section 2)
    E = [0] * (t + 1)
    E[0] = 1
    for i, y in enumerate(ys, start=1):
        for j in range(min(i, t), 0, -1):
            E[j] ^= _fmul(y, E[j - 1], q)
    return E[t] & 1


def evaluate_direct(n, x):
    """Sym_r^t by its definition (sum over t-subsets); cross-check only."""
    from itertools import combinations
    if n == 0:
        return 0
    k, r, t = parameters(n)
    q = irreducible(k)
    mask = (1 << k) - 1
    ys = []
    for i in range(1, r + 1):
        block = (x >> ((i - 1) * k)) & mask
        ys.append(_fpow(block, (1 << (i + 1)) - 1, q))
    total = 0
    for subset in combinations(range(r), t):
        prod = 1
        for idx in subset:
            prod = _fmul(prod, ys[idx], q)
        total ^= prod
    return total & 1
