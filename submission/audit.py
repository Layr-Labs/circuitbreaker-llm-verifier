"""Supplementary exact finite checks for the (31/10, 1) submission.

Every check below is exact (integer or Fraction arithmetic) and finite.
These computations verify: the parameter recipe of proof.md Section 1
and its Section 2 consequences that are exactly decidable (minimality
of k, the Bertrand consequence k < 4*n^(3/5), the chain link
r*(r+1)/t < k, the disperser threshold d(n) and its vacuity range); the
exact rational arithmetic of Section 7 (the six quantities of
Theorem 6.1, the LP optimality certificate with weights 7/10, 1/10,
1/5, the final assembly giving 31/10), of Lemma 5.4 (delta_p = 69/5),
and of the Case 5.2 ledger in elimination.md (Delta_mu >= 173/5 over
t = 2 substitutions); and the well-definedness of the family evaluator
of proof.md Sections 1-2 at small lengths, including agreement of the
descending recurrence with the subset definition of Sym_r^t.

No finite run proves the asymptotic theorem, the BK disperser estimate
of disperser.md, or any all-length statement; those rest on the prose
proof alone. This file only certifies that the finite, exactly
decidable arithmetic the prose relies on is correct as stated.
"""
import importlib.util
import random
from fractions import Fraction
from math import isqrt
from pathlib import Path


def load_family():
    path = Path('/snapshot/files/family.py')
    if not path.exists():
        path = Path(__file__).with_name('family.py')
    spec = importlib.util.spec_from_file_location('candidate_family', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    if not ok:
        raise AssertionError(f"FAILED: {name} {detail}")


def prev_prime(k):
    p = k - 1
    while p >= 2:
        if all(p % d for d in range(2, isqrt(p) + 1)):
            return p
        p -= 1
    return None


def run():
    f = load_family()

    # ---- Section 1 parameters and Section 2 exactly-decidable consequences
    for n in range(1, 3001):
        k, r, t = f.parameters(n)
        assert f.is_prime(k)
        assert k ** 5 > 32 * n ** 3                      # k > 2 n^{3/5}
        pp = prev_prime(k)
        assert pp is None or pp ** 5 <= 32 * n ** 3      # k is the LEAST such prime
        assert (r - 1) * k < n <= r * k                  # r = ceil(n/k)
        assert (t - 1) ** 2 < r <= t ** 2                # t = ceil(sqrt(r))
        assert k ** 5 < 1024 * n ** 3                    # Bertrand: k < 4 n^{3/5}
        assert Fraction(r * (r + 1), t) < k              # Section 2 chain link
        assert r * r <= t * t * (r + 1)                  # r(r+1)/t <= (r+1)^{3/2}
    check("parameters + Section 2 chain (n=1..3000)", True)

    # ---- disperser threshold d(n) = floor(6 n^{4/5}) + 1
    def d_of(n):
        m = int(round((7776 * n ** 4) ** 0.2))
        while (m + 1) ** 5 <= 7776 * n ** 4:
            m += 1
        while m ** 5 > 7776 * n ** 4:
            m -= 1
        return m + 1

    for n in range(1, 9001):
        assert (7776 * n ** 4 >= n ** 5) == (n <= 7776)        # exact vacuity
    for n in range(7777, 20001):
        assert 1 <= d_of(n) <= n                               # threshold nontrivial
    assert d_of(7776) >= 7776                                  # boundary consistency
    check("threshold d(n): vacuity n<=7776 exact; d(n)<=n for 7777..20000", True)

    # ---- Section 7 rational arithmetic (exact)
    aI, aQ, ap = Fraction(48, 5), Fraction(9, 5), Fraction(1, 5)
    six = [aI / 3, 2 - 2 * ap + aQ, 4 - 4 * ap, 3 + ap, 5 - aQ,
           (5 - 2 * ap + aQ) / 2]
    assert six == [Fraction(16, 5), Fraction(17, 5), Fraction(16, 5),
                   Fraction(16, 5), Fraction(16, 5), Fraction(16, 5)]
    u = min(six)
    delta = aI + u
    assert u == Fraction(16, 5) and delta == Fraction(64, 5)
    assert [x - u for x in six] == [0, Fraction(1, 5), 0, 0, 0, 0]  # slacks
    cert = (Fraction(7, 10) * (3 + ap) + Fraction(1, 10) * (5 - aQ)
            + Fraction(1, 5) * ((5 - 2 * ap + aQ) / 2))
    assert cert == Fraction(31, 10) + ap / 2                        # LP certificate
    assert u - ap / 2 == Fraction(31, 10)                           # objective
    assert delta - aI - ap / 2 == Fraction(31, 10)                  # final assembly
    check("Section 7: six quantities, LP certificate (7/10,1/10,1/5), 31/10", True)

    # ---- Lemma 5.4 pair-elimination decrement
    delta_p = min(3 * aI / 2, 3 - 3 * ap + aI + aQ)
    assert delta_p == Fraction(69, 5) and delta_p > delta
    check("Lemma 5.4: delta_p = 69/5 > delta = 64/5", True)

    # ---- Case 5.2 ledger (elimination.md): Delta_mu >= 6 + 3 aI - ap over t = 2
    ledger = 6 + 3 * aI - ap
    assert ledger == Fraction(173, 5)
    assert ledger >= 2 * delta and ledger / 2 >= delta
    check("Case 5.2 ledger: Delta_mu = 173/5 >= 2*delta; per-substitution >= delta", True)

    # ---- evaluator well-definedness at small lengths
    assert f.evaluate(0, 0) == 0
    for n in range(1, 13):
        vals = {f.evaluate(n, x) for x in range(1 << n)}
        assert vals == {0, 1}, (n, vals)                 # total and nonconstant
    rng = random.Random(7)
    for n in (5, 8, 10, 12):
        for _ in range(60):
            x = rng.randrange(1 << n)
            assert f.evaluate(n, x) == f.evaluate_direct(n, x)
    check("evaluator: total, nonconstant (n=1..12); recurrence == subset definition", True)

    return {name: ok for name, ok, _ in RESULTS}


if __name__ == "__main__":
    run()
    for name, ok, _ in RESULTS:
        print(f"{'PASS' if ok else 'FAIL'}  {name}")
    print(f"AUDIT: PASS ({len(RESULTS)} checks, all exact)")
