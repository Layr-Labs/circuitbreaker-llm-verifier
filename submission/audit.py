"""Bounded supplementary checks, not a proof of the circuit theorem."""
import importlib.util
import itertools
import json
from pathlib import Path


def run():
    path = Path('/snapshot/files/family.py')
    if not path.exists():
        path = Path(__file__).with_name('family.py')
    spec = importlib.util.spec_from_file_location('family', path)
    f = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f)
    cases = 0
    for k in (2, 3, 5, 7):
        for low in range(1 << k):
            q = (1 << k) | low
            # Literal divisibility by every possible monic proper factor.
            reducible = any(f.poly_mod(q, (1 << d) | a) == 0
                            for d in range(1, k//2+1)
                            for a in range(1 << d))
            assert f.irreducible_prime_degree(q, k) == (not reducible)
            cases += 1
    capacities = 0
    for k in range(5, 40):
        for r in range(3, k-1):
            for t in range(2, r):
                d, _ = f.certificate_data(k, r, t)
                literal = r*(r+3)//2 + max(
                    sum(max(k-j, 0) for j in range(2, t+1)),
                    sum(max(k-j, 0) for j in range(1, r+1))//(t+1)) + 1
                assert d == literal
                capacities += 1
    # Division-free symmetric recurrence versus literal subset enumeration.
    symmetric = 0
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
        symmetric += 1
    for n in range(25):
        assert f.evaluate([0]*n) == 0
    # Exercise a nonexceptional complete evaluator path, not just zero guards.
    k, m, r, t = f.parameters(94)
    assert (k, m, r, t) == (31, 63, 3, 2)
    q = (1 << 31) | (1 << 3) | 1
    assert f.irreducible_prime_degree(q, k)
    prefix = [(q >> j) & 1 for j in range(k)]
    complete = 0
    for mask in range(8):
        payload = [((j*j+3*j+mask) >> (mask % 3)) & 1 for j in range(m)]
        bits = prefix+payload
        padded = payload+[0]*(r*k-m)
        ys = []
        for i in range(r):
            a = sum(bit << j for j, bit in enumerate(padded[i*k:(i+1)*k]))
            ys.append(f.field_pow(a, (1 << (r+1-i))-1, q))
        total = 0
        for subset in itertools.combinations(ys, t):
            product = 1
            for a in subset:
                product = f.field_mul(product, a, q)
            total ^= product
        assert f.evaluate(bits) == (total & 1)
        complete += 1
    return {'all_passed': True, 'prime_degree_prefixes': cases,
            'literal_capacity_cases': capacities,
            'symmetric_recurrences': symmetric, 'zero_inputs': 25,
            'valid_prefix_evaluator_cases': complete}


if __name__ == '__main__':
    json.load(__import__('sys').stdin)
    print(json.dumps(run(), sort_keys=True))
