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


def literal_integer_bound(k, r, t, witness=False):
    best = sum(k-j for j in range(2, t+1))
    best_dims = [k if i < t-1 else i+2 for i in range(r)]
    for count in range(t, r+1):
        for A in range(t-1, k-count-1):
            upper = min(A//(t-1), (k-t-A)//2,
                        (2*(t-1)*k-t*(t-1)-2*t*A)//(2*(t-1)))
            if upper < 1:
                continue
            # Small witness cases enumerate EVERY pivot, not only the optimizer.
            pivots = range(1, upper+1) if witness else (upper,)
            for x in pivots:
                value = A+sum(min(x, (k-j-A)//2) for j in range(t, count+1))
                if value <= best:
                    continue
                best = value
                if witness:
                    top = [x]*(t-1)
                    left = A-sum(top)
                    for j in range(1, t):
                        add = min(left, k-A-x-j-top[j-1])
                        assert add >= 0
                        top[j-1] += add
                        left -= add
                    assert left == 0
                    deficits = top+[min(x, (k-j-A)//2) for j in range(t, count+1)]
                    best_dims = [d+i+2 for i, d in enumerate(deficits)]+[i+2 for i in range(count, r)]
    return best, best_dims


def oracle(n, mode):
    best = None
    data = None
    for k in range(2, n):
        if any(k%d == 0 for d in range(2, isqrt(k)+1)):
            continue
        m = n-k
        r = (m+k-1)//k
        if r < 3 or k <= r+1:
            continue
        H = sum(range(2, r+2))
        for t in range(2, r):
            old, previous = literal_dimensions(k, r, t)
            if mode == 'old':
                d = old
            elif mode == 'weighted':
                d = previous
            else:
                capacity = sum(k-j for j in range(2, t+1))
                if H+capacity+1 >= m:
                    continue
                d = previous if previous == H+capacity+1 else H+literal_integer_bound(k, r, t)[0]+1
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
    attaining_vectors = 0
    for r in range(3, 21):
        for k in (r+2, 2*r, 2*r+1, 2*r+3):
            for t in range(2, r):
                old, previous = literal_dimensions(k, r, t)
                d, B = f.certificate_data(k, r, t)
                assert B == sum(range(2, t+2)) and d <= previous <= old
                if r <= 10:
                    E, ds = literal_integer_bound(k, r, t, witness=True)
                    assert d == sum(range(2, r+2))+E+1
                    assert sum(ds) == d-1 and all(0 <= x <= k for x in ds)
                    positive = [i for i, x in enumerate(ds) if x > i+2]
                    assert not any(sum(ds[i]-(i+2) for i in chosen)+max(ds[i] for i in chosen) > k+1
                                   for chosen in itertools.combinations(positive, t))
                    attaining_vectors += 1
                dimension_cases += 1

    vectors = 0
    maxima_candidates = 0
    exact_maxima = 0
    for r in range(3, 6):
        k = 2*r+1
        for t in range(2, r):
            actual = -1
            for ds in itertools.product(range(k+1), repeat=r):
                vectors += 1
                total = sum(ds)
                if total <= actual:
                    continue
                maxima_candidates += 1
                positive = [i for i, d in enumerate(ds) if d > i+2]
                good = any(sum(ds[i]-(i+2) for i in chosen)+max(ds[i] for i in chosen) > k+1
                           for chosen in itertools.combinations(positive, t))
                if not good:
                    actual = total
            assert f.certificate_data(k, r, t)[0] == actual+1
            exact_maxima += 1

    comparisons = []
    for n, expected in ((500, (1461, 1437, 1429)), (1000, (2533, 2493, 2485)),
                        (2000, (4513, 4449, 4433)), (10000, (16809, 16601, 16577))):
        old, _ = oracle(n, 'old')
        previous, _ = oracle(n, 'weighted')
        new, params = oracle(n, 'exact')
        assert (old[0], previous[0], new[0]) == expected
        assert f.parameters(n) == params
        k, _, r, t = params
        d, B = f.certificate_data(k, r, t)
        assert new == (4*k+8*d-4*B-3, k, t)
        comparisons.append({'N': n, 'old_loss': old[0], 'weighted_loss': previous[0], 'new_loss': new[0]})

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
            'exhaustive_dimension_vectors': vectors, 'maxima_candidates': maxima_candidates,
            'exact_maximum_cases': exact_maxima, 'attaining_vectors': attaining_vectors,
            'global_comparisons': comparisons, 'symmetric_recurrences': recurrence_cases,
            'valid_prefix_evaluator_cases': 16, 'zero_inputs': 25}


if __name__ == '__main__':
    json.load(__import__('sys').stdin)
    print(json.dumps(run(), sort_keys=True))
