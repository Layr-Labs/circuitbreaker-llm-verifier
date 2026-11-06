"""Exact bounded supplementary predicates; universal claims rest on written proofs."""
import importlib.util
import itertools
import json
from fractions import Fraction
from math import comb, isqrt
from pathlib import Path
from random import Random


def load_family():
    path = Path('/snapshot/files/family.py')
    if not path.exists():
        path = Path(__file__).with_name('family.py')
    spec = importlib.util.spec_from_file_location('candidate_family', path)
    f = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f)
    return f


def literal_dimensions(k, r, t):
    H = sum(range(2, r+2))
    T = sum(max(k-j, 0) for j in range(2, t+1))
    U = sum(max(k-j, 0) for j in range(1, r+1))
    old = H+max(T, U//(t+1))+1
    if k < 2*r+1 or 2*k < (t+1)**2:
        return old, old
    h = next(h for h in range(2, r+1) if h*(t+1) >= 2*r)
    F = Fraction(r*(2*k-2*r-1+h), 2*(t+1))
    new = H+max(T, F.numerator//F.denominator)+1
    assert Fraction(U, t+1)-F == Fraction(r*(r-h), 2*(t+1))
    assert new <= old
    return old, new


def oracle(n, refined):
    best = None
    data = None
    for k in range(2, n):
        if any(k%d == 0 for d in range(2, isqrt(k)+1)):
            continue
        m = n-k
        r = (m+k-1)//k
        if r < 3 or k <= r+1:
            continue
        for t in range(2, r):
            old, new = literal_dimensions(k, r, t)
            d = new if refined else old
            B = sum(range(2, t+2))
            if not (B <= d < m):
                continue
            key = (4*k+8*d-4*B-3, k, t)
            if best is None or key < best:
                best = key
                data = (k, m, r, t)
    return best, data


def run():
    f = load_family()
    prefixes = 0
    for k in (2, 3, 5, 7):
        for low in range(1 << k):
            q = (1 << k) | low
            reducible = any(f.poly_mod(q, (1 << d) | a) == 0
                            for d in range(1, k//2+1)
                            for a in range(1 << d))
            assert f.irreducible_prime_degree(q, k) == (not reducible)
            prefixes += 1

    distributions = 0
    for m in range(3, 9):
        for t in range(2, m):
            h = next(h for h in range(2, m+1) if h*(t+1) >= 2*m)
            high = tuple(range(m-h, m))
            low = tuple(range(m-h))
            lam = Fraction(h*(t+1), m)-2
            assert 0 <= lam <= 1 and len(low) >= t-1
            coefficients = [Fraction(0) for _ in range(m)]
            marked_weight = Fraction(0)
            mass = Fraction(0)
            for extra_high, case_mass in ((0, 1-lam), (1, lam)):
                if not case_mass:
                    continue
                for j in high:
                    for hs in itertools.combinations([i for i in high if i != j], extra_high):
                        for ls in itertools.combinations(low, t-1-extra_high):
                            weight = case_mass/h/comb(h-1, extra_high)/comb(len(low), t-1-extra_high)
                            chosen = (j,)+hs+ls
                            assert len(chosen) == len(set(chosen)) == t
                            mass += weight
                            marked_weight += weight*(j+2)
                            for i in chosen:
                                coefficients[i] += weight
                            coefficients[j] += weight
            assert mass == 1
            assert coefficients == [Fraction(t+1, m)]*m
            assert marked_weight == Fraction(sum(i+2 for i in high), h)
            distributions += 1

    dimension_cases = 0
    for r in range(3, 41):
        for k in (r+2, 2*r, 2*r+1, 2*r+3):
            for t in range(2, r):
                old, new = literal_dimensions(k, r, t)
                d, B = f.certificate_data(k, r, t)
                assert d == new and B == sum(range(2, t+2))
                assert new <= old
                dimension_cases += 1

    vectors = 0
    obligations = 0
    for r in range(3, 6):
        k = 2*r+1
        for t in range(2, r):
            if 2*k < (t+1)**2:
                continue
            _, dimension = literal_dimensions(k, r, t)
            for ds in itertools.product(range(k+1), repeat=r):
                vectors += 1
                if sum(ds) < dimension:
                    continue
                obligations += 1
                positive = [i for i, d in enumerate(ds) if d > i+2]
                assert any(sum(ds[i]-(i+2) for i in chosen)+max(ds[i] for i in chosen) > k+1
                           for chosen in itertools.combinations(positive, t))

    comparisons = []
    for n, expected_old, expected_new in ((200, 669, 661), (1000, 2533, 2493), (10000, 16809, 16601)):
        old, _ = oracle(n, False)
        new, params = oracle(n, True)
        assert (old[0], new[0]) == (expected_old, expected_new)
        assert f.parameters(n) == params
        k, _, r, t = params
        d, B = f.certificate_data(k, r, t)
        assert new == (4*k+8*d-4*B-3, k, t)
        comparisons.append({'N': n, 'old_loss': old[0], 'new_loss': new[0]})

    recurrence_cases = 0
    for values in itertools.product(range(8), repeat=3):
        q = 0b1011
        e = [1, 0, 0]
        for value in values:
            for j in (2, 1):
                e[j] ^= f.field_mul(value, e[j-1], q)
        literal = 0
        for a, b in itertools.combinations(values, 2):
            literal ^= f.field_mul(a, b, q)
        assert e[2] == literal
        recurrence_cases += 1

    assert f.parameters(32) == (7, 25, 4, 2)
    k, m, r, t = f.parameters(32)
    q = 131
    assert f.irreducible_prime_degree(q, k)
    prefix = [(q >> j) & 1 for j in range(k)]
    rng = Random(20261106)
    for _ in range(16):
        payload = [rng.randrange(2) for _ in range(m)]
        padded = payload+[0]*(r*k-m)
        ys = []
        for i in range(r):
            value = sum(bit << j for j, bit in enumerate(padded[i*k:(i+1)*k]))
            ys.append(f.field_pow(value, (1 << (r+1-i))-1, q))
        expected = 0
        for subset in itertools.combinations(ys, t):
            product = 1
            for value in subset:
                product = f.field_mul(product, value, q)
            expected ^= product
        assert f.evaluate(prefix+payload) == (expected & 1)
    for n in range(25):
        assert f.evaluate([0]*n) == 0
    return {'all_passed': True, 'prime_degree_prefixes': prefixes,
            'marked_distributions': distributions, 'dimension_cases': dimension_cases,
            'exhaustive_dimension_vectors': vectors, 'isolation_obligations': obligations,
            'global_comparisons': comparisons, 'symmetric_recurrences': recurrence_cases,
            'valid_prefix_evaluator_cases': 16, 'zero_inputs': 25}


if __name__ == '__main__':
    json.load(__import__('sys').stdin)
    print(json.dumps(run(), sort_keys=True))
