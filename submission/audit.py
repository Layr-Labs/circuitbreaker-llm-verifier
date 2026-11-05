"""Supplementary exact finite checks; no finite run proves the asymptotic theorem."""
import importlib.util
import itertools
import json
from fractions import Fraction
from pathlib import Path
import sys
import random


def load_family():
    path = Path('/snapshot/files/family.py')
    if not path.exists():
        path = Path(__file__).with_name('family.py')
    spec = importlib.util.spec_from_file_location('candidate_family', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def table_value(table, x, y):
    return (table >> (2 * x + y)) & 1


def matrix_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit in pivots:
                row ^= pivots[bit]
            else:
                pivots[bit] = row
                break
    return len(pivots)


def mixed_audit(f, rng):
    obstruction_cases = 0
    zero_cases = 0
    for degree in [4, 8]:
        field = f.tower_field(degree)
        half, nu = degree // 2, 1 << (degree // 2)
        for z in range(1, 1 << degree):
            z3, z7 = field.power(z, 3), field.power(z, 7)
            for mu in range(1, 1 << half):
                a = field.mul(mu, field.mul(nu, z7)) >> half
                b = field.mul(mu, z3) >> half
                if a == 0:
                    assert b != 0
                    zero_cases += 1
                obstruction_cases += 1
    trace_cases = 0
    for degree in [4, 8, 16]:
        field = f.tower_field(degree)
        half, nu = degree // 2, 1 << (degree // 2)
        for _ in range(100):
            x, mu = rng.randrange(1 << degree), rng.randrange(1, 1 << half)
            value = field.power(x, 3) ^ field.mul(nu, field.power(x, 7))
            relative = value ^ field.power(value, 1 << half)
            assert relative == f.mixed_cubic(field, x)
            assert field.trace(field.mul(mu, value)) == field.base.trace(field.base.mul(mu, relative))
            trace_cases += 1
    field = f.tower_field(16)
    base, half, nu = field.base, 8, 1 << 8
    inverse_nu = field.power(nu, (1 << 16) - 2)
    inverse_seven = pow(7, -1, (1 << 16) - 1)
    rank_cases = 0
    formula_cases = 0
    ranks = []
    for index in range(4):
        target, mu = rng.randrange(1, 1 << half), rng.randrange(1, 1 << half)
        z = field.power(field.mul(inverse_nu, target), inverse_seven)
        z3, z7 = field.power(z, 3), field.power(z, 7)
        lam = field.mul(mu, field.mul(nu, z7))
        B = field.mul(mu, z3) >> half
        assert 0 < lam < 1 << half and B
        xi = 0 if index % 2 == 0 else rng.randrange(1, 1 << half)
        h = rng.randrange(1 << half) | (xi << half)
        c = field.mul(z, h)
        actual = [base.trace(base.mul(mu, f.mixed_cubic(field, c ^ field.mul(z, y))))
                  for y in range(1 << half)]
        C3 = base.mul(lam, base.power(xi, 4)) ^ B
        C5, C6 = base.mul(lam, base.power(xi, 2)), base.mul(lam, xi)
        quadratic = [base.trace(base.mul(C3, base.power(y, 3))
                                ^ base.mul(C5, base.power(y, 5))
                                ^ base.mul(C6, base.power(y, 6)))
                     for y in range(1 << half)]
        constant = actual[0] ^ quadratic[0]
        linear = [actual[1 << j] ^ quadratic[1 << j] ^ constant for j in range(half)]
        for y in range(1 << half):
            assert actual[y] ^ quadratic[y] == constant ^ (sum(linear[j] for j in range(half) if y >> j & 1) & 1)
            formula_cases += 1
        rows = [sum((actual[(1 << i) ^ (1 << j)] ^ actual[1 << i]
                     ^ actual[1 << j] ^ actual[0]) << j for j in range(half))
                for i in range(half)]
        value = matrix_rank(rows)
        assert value >= half - (4 if xi else 2)
        ranks.append(value)
        rank_cases += 1
    return {'mixed_subfield_cases': obstruction_cases,
            'mixed_zero_cubic_cases': zero_cases,
            'mixed_trace_cases': trace_cases, 'mixed_rank_cases': rank_cases,
            'mixed_rank_values': ranks, 'mixed_quadratic_formula_cases': formula_cases}


def run():
    f = load_family()
    results = {}
    # All AND-type operations and every possible second B2 operation.
    and_tables = [v for v in range(16) if v.bit_count() in (1, 3)]
    points = list(itertools.product((0, 1), repeat=2))
    common_line_cases = 0
    for g, h in itertools.product(and_tables, range(16)):
        assert any(
            table_value(g, *p) == table_value(g, *q)
            and table_value(h, *p) == table_value(h, *q)
            for p, q in itertools.combinations(points, 2)
        )
        common_line_cases += 1
    results['common_child_affine_line_cases'] = common_line_cases

    # Composition h(x,g(x,y)) always remains one B2 table.
    composition_cases = 0
    for g, h in itertools.product(and_tables, range(16)):
        composed = sum(table_value(h, x, table_value(g, x, y)) << (2*x+y)
                       for x, y in points)
        assert 0 <= composed < 16
        for x, y in points:
            assert table_value(composed, x, y) == table_value(h, x, table_value(g, x, y))
        composition_cases += 1
    results['composed_B2_tables'] = composition_cases

    # Quadratic reparameterization is a bijection on all 3-bit inputs.
    qr_cases = 0
    for g in and_tables:
        for c in (0, 1):
            image = set()
            for x, y, w in itertools.product((0, 1), repeat=3):
                z = w ^ table_value(g, x, y) ^ c
                assert table_value(g, x, y) ^ z ^ c == w
                image.add((x, y, z))
                qr_cases += 1
            assert len(image) == 8
    results['quadratic_pivot_assignments'] = qr_cases

    # Arithmetic only for the ternary-pair potential W=s+4I+2q.
    # Universal applicability and preservation of domains are proof obligations.
    rows = [(-1,0,0,0), (0,0,-1,0), (-2,1,-1,0),
            (-1,0,-1,0), (-3,-1,-1,-1), (-4,0,-2,-1)]
    for ds,di,dq,da in rows:
        assert ds+4*di+2*dq <= 8*da
    for children in range(1,33):
        assert -(1+children)-8+2 <= -8
    results['ternary_potential_rows_checked'] = len(rows)+32

    rng = random.Random(2026110502)
    root_cases = 0
    for n in range(1000):
        for degree in range(1, 12):
            r = f.floor_root(n, degree)
            assert r ** degree <= n < (r + 1) ** degree
            c = f.ceil_root(n, degree)
            assert c ** degree >= n and (c == 0 or (c - 1) ** degree < n)
            root_cases += 1
    
    field_cases = 0
    for degree in [1, 2, 4, 8, 16]:
        tower = f.tower_field(degree)
        assert tower.trace(tower.trace_one) == 1
        for a in range(1, min(1 << degree, 300)):
            assert tower.power(a, (1 << degree) - 1) == 1
            assert tower.power(a, 1 << degree) == a
            assert tower.trace(a) in [0, 1]
            field_cases += 1
        for _ in range(300):
            a, b, c = [rng.randrange(1 << degree) for _ in range(3)]
            assert tower.mul(a, b) == tower.mul(b, a)
            assert tower.mul(a, b ^ c) == tower.mul(a, b) ^ tower.mul(a, c)
            assert tower.mul(tower.mul(a, b), c) == tower.mul(a, tower.mul(b, c))
            field_cases += 1
    
    irreducible_cases = 0
    for degree in range(1, 10):
        for low in range(1 << degree):
            modulus = (1 << degree) | low
            trial = all(f.poly_mod(modulus, divisor) != 0
                        for d in range(1, degree // 2 + 1)
                        for divisor in range(1 << d, 1 << (d + 1)))
            assert f.irreducible(modulus) == trial, modulus
            irreducible_cases += 1
    
    design_cases = 0
    for ell in range(2, 13):
        for intersection in range(1, 4):
            count = min(20, (1 << (ell - 1).bit_length()) ** (intersection + 1))
            dimension, sets = f.polynomial_design(ell, count, intersection)
            assert all(len(set(s)) == ell and max(s) < dimension for s in sets)
            assert all(len(set(a) & set(b)) <= intersection for a, b in itertools.combinations(sets, 2))
            design_cases += 1
    for ell, count, intersection, alphabet_bits in [(4, 8, 2, 1), (5, 8, 3, 1), (4, 10, 2, 2)]:
        dimension, sets = f.greedy_design(ell, count, intersection, alphabet_bits)
        assert all(len(set(a) & set(b)) <= intersection for a, b in itertools.combinations(sets, 2))
        design_cases += 1
    
    code_cases = 0
    for n in [1, 2, 7, 16, 47, 100]:
        code = f.LinearCode(n, 3, 2)
        for _ in range(50):
            a, b = rng.randrange(1 << n), rng.randrange(1 << n)
            index = rng.randrange(1 << code.ell)
            assert code.bit(a ^ b, index) == code.bit(a, index) ^ code.bit(b, index)
            code_cases += 1
    
    results.update(root_cases=root_cases, field_cases=field_cases,
                   irreducible_cases=irreducible_cases, design_cases=design_cases,
                   code_linearity_cases=code_cases)
    delayed_root_cases = 0
    for n in range(500):
        for exponent in range(12):
            degree = 1 << exponent
            root = f.floor_root_power_two(n, exponent)
            assert root ** degree <= n < (root + 1) ** degree
            delayed_root_cases += 1
    results['delayed_degree_root_cases'] = delayed_root_cases
    subset_cases = 0
    for count in range(11):
        for size in range(count + 1):
            assert set(f.index_subsets(count, size)) == set(itertools.combinations(range(count), size))
            subset_cases += 1
    results['subset_enumeration_cases'] = subset_cases
    results.update(mixed_audit(f, rng))
    evaluations = 0
    for n in [0, 1, 2, 3, 17, 100, 1024, 10000]:
        assert f.FamilyShape.build(n) is None
        assert f.evaluate([0] * n) == 0
        assert f.evaluate([j & 1 for j in range(n)]) == 0
        evaluations += 2
    results['actual_guard_evaluations'] = evaluations
    return {'all_passed': True, 'checks': results,
            'scope': 'Finite exact evidence only. All full-family evaluations use the documented zero guard. Universal extractor properties, uniformity, actual-length charge and asymptotics are proved in the Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
