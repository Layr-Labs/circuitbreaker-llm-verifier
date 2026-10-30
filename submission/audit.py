"""Supplementary exact finite checks; no finite run proves the asymptotic theorem."""
import importlib.util
import itertools
import json
from fractions import Fraction
from pathlib import Path
import sys


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

    # The universal inequality is proved in elimination.md; this bounded
    # grid checks the arithmetic and the exact exceptional configuration.
    potential_cases = 0
    for u in range(2, 13):
        for h in range(1, 13):
            for h1 in range(h + 1):
                bound = -7 - 2*(u+h) + h1 - (1 if u == 2 else 0)
                if (u, h, h1) == (2, 1, 1):
                    assert bound == -13
                else:
                    assert bound <= -14
                potential_cases += 1
    assert -4+3+(7-6)-7-7 == -14
    assert -7-7-6-8 == -28
    for numerator in range(33):
        sigma = Fraction(numerator, 8)
        assert min(4-sigma, 3+sigma) <= Fraction(7, 2)
    results['potential_grid_cases'] = potential_cases

    # Cross-check prime-degree irreducibility against exhaustive trial
    # polynomial division, for every monic polynomial at these degrees.
    irreducibility_cases = 0
    irreducible_counts = {}
    for k in (2, 3, 5, 7):
        count = 0
        for low in range(1 << k):
            q = (1 << k) | low
            trial = not any(f.poly_mod(q, divisor) == 0
                            for degree in range(1, k//2 + 1)
                            for divisor in range(1 << degree, 1 << (degree+1)))
            assert f.irreducible_prime_degree(q, k) == trial
            count += trial
            irreducibility_cases += 1
        assert count == ((1 << k) - 2) // k
        irreducible_counts[str(k)] = count
    results['prefix_polynomials_checked'] = irreducibility_cases
    results['irreducible_counts'] = irreducible_counts

    # Weighted block-isolation check. For nested allowed destinations i<=j,
    # all suffix-capacity inequalities are necessary for an allocation.
    # No unequal equal-cardinality pair should satisfy them.
    block_cases = 0
    for r in range(2, 9):
        for t in range(1, r):
            sets = list(itertools.combinations(range(1, r+1), t))
            for target in sets:
                for source in sets:
                    possible = all(
                        sum(r+2-i for i in target if i >= threshold)
                        <= sum(r+2-j for j in source if j >= threshold)
                        for threshold in range(1, r+1)
                    )
                    assert possible == (target == source)
                    block_cases += 1
    results['block_capacity_pairs'] = block_cases

    # Exercise the exact family code at actual lengths, and compare the
    # symmetric recurrence against direct enumeration on small valid cases.
    evaluation_cases = 0
    tested_parameters = []
    for n in (0, 1, 2, 8, 16, 32, 64, 128, 192, 256):
        params = f.parameters(n)
        if params is None:
            assert f.evaluate([0]*n) == 0
            evaluation_cases += 1
            continue
        k, m, r, t = params
        assert (k >= 2 and f.is_prime(k) and k**5 > 32*n**3
                and k+m == n and r*k >= m and t*t >= r and t < r
                and k*t >= r*(r+1))
        assert f.evaluate([0]*n) == 0  # Z^k is an invalid prefix.
        # Bounded test fixture selection only. The language itself never
        # performs this search, and the theorem does not depend on it.
        modulus = next((1 << k) | low for low in range(1, 8192, 2)
                       if low < (1 << k)
                       and f.irreducible_prime_degree((1 << k) | low, k))
        prefix = [(modulus >> j) & 1 for j in range(k)]
        for seed in (0, 1, 2, 3):
            payload = [((j*j + 3*j + seed) >> (j % 5)) & 1 for j in range(m)]
            padded = payload + [0]*(r*k-m)
            ys = []
            for i in range(r):
                value = sum(padded[i*k+j] << j for j in range(k))
                ys.append(f.field_pow(value, (1 << (r-i+1))-1, modulus))
            total = 0
            for selected in itertools.combinations(range(r), t):
                term = 1
                for i in selected:
                    term = f.field_mul(term, ys[i], modulus)
                total ^= term
            assert f.evaluate(prefix+payload) == (total & 1)
            evaluation_cases += 1
        tested_parameters.append({'N': n, 'k': k, 'm': m, 'r': r, 't': t})
    results['family_evaluation_cases'] = evaluation_cases
    results['family_parameters_tested'] = tested_parameters
    # Exact constants in the new universal analytic argument. Finite
    # samples exercise ceilings and the sufficient bound; they are not
    # extrapolated to establish the all-length theorem.
    a = Fraction(401, 100)
    sqrt_upper = Fraction(801, 400)
    assert sqrt_upper**2 - a == Fraction(1, 160000)
    envelope = (14*sqrt_upper + Fraction(7, 4)
                + Fraction(21, 2)*a/Fraction(10)
                + Fraction(7, 2)/Fraction(100) + Fraction(7, 10000))
    assert envelope == Fraction(42539, 1250) and envelope < 35
    assert Fraction(7, 2)*10 == 35
    results['quantitative_envelope'] = str(envelope)
    quantitative_cases = 0
    sample_lengths = sorted(set(range(100000-5, 100000+6))
                            | {10**power for power in range(5, 13)})
    for n in sample_lengths:
        k, m, r, t = f.parameters(n)
        assert r*k < n <= (r+1)*k
        assert t*t >= r and t < r and k*t >= r*(r+1)
        d = (k*r)//t + k*t + r*(r+1) + 1
        twice_loss = 7*k + 14*d
        # Both sides are nonnegative; fifth powers preserve order.
        assert twice_loss**5 <= 70**5 * n**4
        if n >= 100000:
            assert (100*k)**5 < 401**5 * n**3
            assert r >= 24
        quantitative_cases += 1
    results['quantitative_parameter_cases'] = quantitative_cases
    results['quantitative_parameter_range'] = [min(sample_lengths), max(sample_lengths)]
    for epsilon in (Fraction(1, 2), Fraction(1, 10), Fraction(1, 100)):
        threshold = (Fraction(35)/epsilon)**5
        assert threshold.denominator == 1
        assert epsilon**5 * int(threshold) == 35**5
    results['explicit_threshold_cases'] = 3
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
