#!/usr/bin/env python3
"""Exact-arithmetic checks for the 4n bound submission.

Everything here is bookkeeping, not mathematics.  The checks confirm:

  1. the classification of the sixteen binary operations used by the counting
     argument, and which of them are affine;
  2. every row of the accounting table in gate-elimination.md Section 9,
     recomputed from its stated changes in s, q, m and in parameter weights;
  3. the arithmetic of Theorem 10.1:  W_0 < 10(n-d)  is equivalent to
     |C_0| < 4n - 10d, and 0 <= W_t <= W_0 + 10(a_t - n) makes a_t <= d
     impossible;
  4. the parameter arithmetic of family.md Lemma 5.1 (k prime, r, t, d(n) < n,
     t < r) and the Ben-Sasson--Kopparty parameter hypotheses for Theorem 2.2;
  5. the exact rationals of bound.json, the identity of the two structured
     claims, and the implication from the declared guarantee to the declared
     projection, via the repository's expressions.py;
  6. the eight-file, 64 KiB, 256 KiB submission limits.

Run:  python3 submission/checks.py            (from the repository root)
      python3 checks.py                       (from submission/)

It does not verify the case analysis, the invariant, or the disperser property.
"""

from fractions import Fraction
from math import isqrt
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

FAILURES = []


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    if not condition:
        FAILURES.append(name)
    print(f"{status}  {name}" + (f"   [{detail}]" if detail else ""))


# ---------------------------------------------------------------------------
# integer roots and exact thresholds
# ---------------------------------------------------------------------------

def iroot(x, k):
    """floor(x ** (1/k)) for nonnegative integers, exactly."""
    if x < 0:
        raise ValueError("negative radicand")
    if x == 0:
        return 0
    lo, hi = 0, 1
    while hi ** k <= x:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** k <= x:
            lo = mid
        else:
            hi = mid
    return lo


def is_prime(m):
    if m < 2:
        return False
    if m % 2 == 0:
        return m == 2
    f = 3
    while f * f <= m:
        if m % f == 0:
            return False
        f += 2
    return True


def least_prime_above_two_n_to_3_5(n):
    """Least prime k with k**5 > 32 n**3, i.e. k > 2 n^(3/5)."""
    bound = 32 * n ** 3
    k = iroot(bound, 5)
    while not (k ** 5 > bound and is_prime(k)):
        k += 1
    return k


def floor_six_n_to_4_5(n):
    """floor(6 n^(4/5)) = floor((6**5 n**4)^(1/5))."""
    return iroot(6 ** 5 * n ** 4, 5)


# ---------------------------------------------------------------------------
# 1. the sixteen operations
# ---------------------------------------------------------------------------

def classify(b1, b2, b12):
    if b12 == 1:
        return "AND-type"
    if b1 == 0 and b2 == 0:
        return "constant"
    if b1 == 1 and b2 == 1:
        return "XOR-type"
    return "unary"


def check_operations():
    counts = {"constant": 0, "unary": 0, "XOR-type": 0, "AND-type": 0}
    affine = 0
    seen = set()
    for c in (0, 1):
        for b1 in (0, 1):
            for b2 in (0, 1):
                for b12 in (0, 1):
                    table = tuple(
                        (c + b1 * t1 + b2 * t2 + b12 * t1 * t2) % 2
                        for t1 in (0, 1) for t2 in (0, 1)
                    )
                    seen.add(table)
                    counts[classify(b1, b2, b12)] += 1
                    if b12 == 0:
                        affine += 1
    check("all 16 truth tables are distinct", len(seen) == 16, f"{len(seen)} tables")
    check("2 constant, 4 unary, 2 XOR-type, 8 AND-type",
          counts == {"constant": 2, "unary": 4, "XOR-type": 2, "AND-type": 8},
          str(counts))
    check("8 affine operations (constant + unary + XOR-type)",
          affine == 8 and (counts["constant"] + counts["unary"] + counts["XOR-type"]) == 8,
          f"affine={affine}")


# ---------------------------------------------------------------------------
# 2. accounting table of gate-elimination.md Section 9
# ---------------------------------------------------------------------------
# Each row: name, ds, dq, dm, list of unpaired-weight changes (each in {0,6},
# or -6 for a parameter that stops being an unpaired used parameter).
# dW = ds + 4*dq + sum(weights), and da = dm - dq.

ROWS = [
    ("D1/D2/D3 deletion",              -1,  0,  0, []),
    ("D4 unit parameter deletion",     -1,  0,  0, []),
    ("AR reparameterisation",           0,  0,  0, []),
    ("move gate into affine system",    0,  0,  0, []),
    ("Case 1, deg(p) = 0",              0, -1, -1, [0]),    # p becomes unused
    ("Case 1, deg(p) >= 2",            -2, -1, -1, [6]),    # r becomes used
    ("Case 2, shared child",           -1, -1, -1, [0]),    # survivor unused
    ("Case 3, paired param in system", -2, -1, -1, [6]),    # p becomes used/isolated
    ("Case 4, both inputs unpaired",   -2,  1,  0, [-6, -6]),
    ("Case 5, k = 1",                  -2, -1, -2, [-6]),
    ("Case 5, k = 2",                  -2, -2, -3, []),
]

EXPECTED_DW = {
    "D1/D2/D3 deletion": -1,
    "D4 unit parameter deletion": -1,
    "AR reparameterisation": 0,
    "move gate into affine system": 0,
    "Case 1, deg(p) = 0": -4,
    "Case 1, deg(p) >= 2": 0,
    "Case 2, shared child": -5,
    "Case 3, paired param in system": 0,
    "Case 4, both inputs unpaired": -10,
    "Case 5, k = 1": -12,
    "Case 5, k = 2": -10,
}


def check_accounting():
    for name, ds, dq, dm, weights in ROWS:
        da = dm - dq
        dW = ds + 4 * dq + sum(weights)
        check(f"{name}: da = dm - dq in {{0,-1}}", da in (0, -1), f"da={da}")
        check(f"{name}: dW <= 10*da", dW <= 10 * da, f"dW={dW}, 10*da={10 * da}")
        check(f"{name}: matches the derived bound", dW == EXPECTED_DW[name],
              f"dW={dW}, expected={EXPECTED_DW[name]}")
    neg = {name for name, ds, dq, dm, w in ROWS if dm - dq == -1}
    need = {name for name in EXPECTED_DW if EXPECTED_DW[name] <= -10}
    check("rows with da = -1 are exactly the rows needing dW <= -10", neg == need,
          f"{sorted(neg)}")


# ---------------------------------------------------------------------------
# 3. Theorem 10.1 arithmetic
# ---------------------------------------------------------------------------

def check_theorem_arithmetic():
    bad_equiv = []
    bad_sum = []
    for n in (16807, 10 ** 5, 10 ** 6, 10 ** 7):
        d = floor_six_n_to_4_5(n) + 1
        if not d < n:
            bad_equiv.append((n, "d<n", d))
            continue
        bound = 4 * n - 10 * d
        # W_0 <= |C_0| + 6n  <  10(n-d)   <==>   |C_0| < 4n - 10d
        for c in range(max(0, bound - 4), bound + 4):
            if (c + 6 * n < 10 * (n - d)) != (c < bound):
                bad_equiv.append((n, c))
        # premise |C_0| < 4n - 10d, i.e. W_0 < 10(n-d); then for a_t <= d,
        # 0 <= W_t <= W_0 + 10(a_t - n) < 10(a_t - d) <= 0 is impossible.
        # Only meaningful when the hypothesis is satisfiable, i.e. bound > 0.
        if bound <= 0:
            continue
        for c0 in (0, max(0, bound - 1)):
            w0 = c0 + 6 * n
            if not w0 < 10 * (n - d):
                bad_equiv.append((n, c0, "w0"))
                continue
            for a_t in (d - 5, d - 1, d):
                upper = w0 + 10 * (a_t - n)
                if upper >= 0:
                    bad_sum.append((n, a_t, w0, upper))
    check("W_0 < 10(n-d) is equivalent to |C_0| < 4n - 10d",
          not bad_equiv, f"{len(bad_equiv)} mismatches: {bad_equiv[:2]}")
    check("a_t <= d is refuted by 0 <= W_t <= W_0 + 10(a_t - n)", not bad_sum,
          f"{len(bad_sum)} violations: {bad_sum[:2]}")


# ---------------------------------------------------------------------------
# 4. family.md Lemma 5.1 and the BK parameter hypotheses
# ---------------------------------------------------------------------------

def check_family(n_max=100000):
    bad = []
    first_d_lt_n = None
    for n in range(64, n_max + 1):
        k = least_prime_above_two_n_to_3_5(n)
        if not (k ** 5 > 32 * n ** 3):
            bad.append(("k > 2n^3/5", n, k))
        if not (k ** 5 < 4 ** 5 * n ** 3):
            bad.append(("k < 4n^3/5", n, k))
        r = (n + k - 1) // k
        t = isqrt(r - 1) + 1
        if r < 2:
            bad.append(("r >= 2", n, r))
        if not (1 <= t <= r):
            bad.append(("1 <= t <= r", n, t, r))
        if r >= 3 and not (t < r):
            bad.append(("t < r", n, t, r))
        d = floor_six_n_to_4_5(n) + 1
        if d < n and first_d_lt_n is None:
            first_d_lt_n = n
        if n >= 16807 and not (d < n):
            bad.append(("d(n) < n", n, d))
        # BK Theorem 2.2 hypothesis: k >= r(r+1)/t, i.e. k*t >= r(r+1)
        if not (k * t >= r * (r + 1)):
            bad.append(("k*t >= r(r+1)", n, k, t, r))
        # and kr/t + kt + r(r+1) <= 6 n^(4/5), compared exactly via 5th powers
        s = Fraction(k * r, t) + k * t + r * (r + 1)
        lhs = s.numerator ** 5
        rhs = (6 ** 5) * (n ** 4) * s.denominator ** 5
        if not (lhs <= rhs):
            bad.append(("kr/t+kt+r(r+1) <= 6n^(4/5)", n, float(s), 6 * n ** 0.8))
    check(f"family parameters valid for 64 <= n <= {n_max}", not bad,
          f"{len(bad)} failures; first={bad[:3]}")
    check("d(n) < n for every n >= 16807 (and at n = 16807)", not bad,
          f"{len(bad)} failures")
    check("the exact first index with d(n) < n is 7781", first_d_lt_n == 7781,
          str(first_d_lt_n))


# ---------------------------------------------------------------------------
# 5. bound.json, claims, and the projection implication
# ---------------------------------------------------------------------------

def check_declaration():
    bound = json.load(open(os.path.join(HERE, "bound.json")))
    sub = json.load(open(os.path.join(HERE, "submission.json")))

    g, p = bound["guarantee"], bound["projection"]
    check("guarantee kind is polynomial_log", g["kind"] == "polynomial_log", g["kind"])
    check("declared coefficient is 4/1",
          (g["coefficient"]["numerator"], g["coefficient"]["denominator"]) == ("4", "1"),
          str(g["coefficient"]))
    check("declared exponent is 1/1",
          (g["exponent"]["numerator"], g["exponent"]["denominator"]) == ("1", "1"),
          str(g["exponent"]))
    check("log and loglog exponents are 0/1",
          (g["log_exponent"]["numerator"], g["log_exponent"]["denominator"]) == ("0", "1")
          and (g["loglog_exponent"]["numerator"], g["loglog_exponent"]["denominator"]) == ("0", "1"))
    check("projection is the finite pair (4/1, 1/1)",
          p["kind"] == "finite"
          and (p["coefficient"]["numerator"], p["coefficient"]["denominator"]) == ("4", "1")
          and (p["exponent"]["numerator"], p["exponent"]["denominator"]) == ("1", "1"),
          json.dumps(p, sort_keys=True))

    try:
        import expressions
    except Exception as exc:                                   # pragma: no cover
        check("import expressions.py", False, repr(exc))
        return
    try:
        scores = expressions.validate_bound(bound)
        check("expressions.validate_bound accepts bound.json", True, str(scores))
    except Exception as exc:
        check("expressions.validate_bound accepts bound.json", False, repr(exc))
        return
    expected_scores = [{"numerator": "4", "denominator": "1"},
                       {"numerator": "1", "denominator": "1"}]
    check("transport scores are the reduced rationals 4/1 and 1/1",
          scores == expected_scores, str(scores))

    statement = expressions.joint_statement(bound)
    claims = sub["claims"]
    check("at least two structured claims", len(claims) >= 2, str(len(claims)))
    check("first two statements are identical",
          claims[0]["statement"] == claims[1]["statement"])
    check("first statement is the canonical joint statement",
          claims[0]["statement"] == statement)
    ev = claims[0]["evidence"]
    check("identical sorted evidence lists containing bound.json and proof.md",
          claims[0]["evidence"] == claims[1]["evidence"] == sorted(ev)
          and "bound.json" in ev and "proof.md" in ev,
          str(ev))
    check("every evidence file exists in submission/",
          all(os.path.exists(os.path.join(HERE, f)) for f in ev), str(ev))
    check("claim scores are the two transport coordinates",
          claims[0]["score"] == {"numerator": "4", "denominator": "1"}
          and claims[1]["score"] == {"numerator": "1", "denominator": "1"},
          f"{claims[0]['score']} / {claims[1]['score']}")
    check("canonical claims_for output matches the file",
          claims[:2] == expressions.claims_for(bound, ev))
    check("remaining claims are non-empty strings",
          all(isinstance(c, str) and c.strip() for c in claims[2:]),
          f"{len(claims) - 2} string claims")
    for i, c in enumerate(claims[2:]):
        check(f"string claim {i + 2} within the 8000-char claim limit", len(c) <= 8000,
              str(len(c)))


# ---------------------------------------------------------------------------
# 6. submission limits
# ---------------------------------------------------------------------------

def check_limits():
    files = sorted(f for f in os.listdir(HERE)
                   if not f.startswith(".") and f != "__pycache__")
    sizes = {f: os.path.getsize(os.path.join(HERE, f)) for f in files}
    total = sum(sizes.values())
    check("at most eight files in submission/", len(files) <= 8, f"{len(files)}: {files}")
    check("each file at most 64 KiB", all(v <= 65536 for v in sizes.values()),
          str({k: v for k, v in sizes.items() if v > 65536}))
    check("aggregate at most 256 KiB", total <= 262144, f"{total} bytes")


def main():
    print("== 1. sixteen operations ==")
    check_operations()
    print("\n== 2. accounting table ==")
    check_accounting()
    print("\n== 3. Theorem 10.1 arithmetic ==")
    check_theorem_arithmetic()
    print("\n== 4. family parameters ==")
    check_family()
    print("\n== 5. declaration ==")
    check_declaration()
    print("\n== 6. limits ==")
    check_limits()
    print()
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s): {FAILURES}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
