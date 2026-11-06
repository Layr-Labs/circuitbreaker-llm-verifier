"""Bounded exact checks for the supplied 4n-o(n) proof, not a proof by computation."""
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
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def table_value(table, x, y):
    return (table >> (2 * x + y)) & 1


def run():
    f = load_family()
    results = {}

    # On pr=0, every binary table is constant on at least one of the three
    # lines p=0, r=0, p+r=1. This supports Case 2 only.
    pair_lines = (
        ((0, 0), (0, 1)),
        ((0, 0), (1, 0)),
        ((0, 1), (1, 0)),
    )
    pair_line_cases = 0
    for table in range(16):
        assert any(
            len({table_value(table, *point) for point in line}) == 1
            for line in pair_lines
        )
        pair_line_cases += 1
    results['pair_line_truth_tables'] = pair_line_cases

    # Exact arithmetic for every potential estimate in elimination.md.
    # The literals correspond, in order, to s, 4q, and 6I changes.
    potential = {
        'case1_unused': -4,
        'case1_degree_at_least_2': -4 - 2 + 6,
        'case2_common_child': -4 - 1,
        'case3_affine_system': -4 - 2 + 6,
        'case4_pair_minimal_and': -2 + 4 - 6 - 6,
        'case5_one_paired_input': -2 - 4 - 6,
        'case5_two_paired_inputs': -2 - 8,
    }
    expected = {
        'case1_unused': -4,
        'case1_degree_at_least_2': 0,
        'case2_common_child': -5,
        'case3_affine_system': 0,
        'case4_pair_minimal_and': -10,
        'case5_one_paired_input': -12,
        'case5_two_paired_inputs': -10,
    }
    assert potential == expected
    for name in ('case1_unused', 'case1_degree_at_least_2',
                 'case2_common_child', 'case3_affine_system'):
        assert Fraction(potential[name]) <= 0
    for name in ('case4_pair_minimal_and', 'case5_one_paired_input',
                 'case5_two_paired_inputs'):
        assert Fraction(potential[name]) <= -10
    results['potential_estimates'] = potential

    # The integer dimension boundary and final coefficient arithmetic.
    dimension_cases = 0
    for a in range(65):
        for d in range(65):
            if a > d:
                assert a - 1 >= d
                dimension_cases += 1
    results['integer_dimension_cases'] = dimension_cases
    accounting_cases = 0
    for n in range(1, 65):
        for d in range(n + 2):
            bound = 4 * n - 10 * d
            for size in range(max(0, bound)):
                # At a=d, W <= W0+10(d-n) <= size+6n+10(d-n) < 0.
                assert size + 6 * n + 10 * (d - n) < 0
                accounting_cases += 1
    assert Fraction(4) + Fraction(6) == Fraction(10)
    results['final_accounting_cases'] = accounting_cases

    # Compare the prime-degree Frobenius/gcd prefix test with exhaustive
    # polynomial trial division at these bounded degrees.
    irreducibility_cases = 0
    irreducible_counts = {}
    for k in (2, 3, 5, 7):
        count = 0
        for low in range(1 << k):
            q = (1 << k) | low
            trial = not any(
                f.poly_mod(q, divisor) == 0
                for degree in range(1, k // 2 + 1)
                for divisor in range(1 << degree, 1 << (degree + 1))
            )
            assert f.irreducible_prime_degree(q, k) == trial
            count += trial
            irreducibility_cases += 1
        assert count == ((1 << k) - 2) // k
        irreducible_counts[str(k)] = count
    results['prefix_polynomials_checked'] = irreducibility_cases
    results['irreducible_counts'] = irreducible_counts

    # Exhaustively test the suffix-capacity form of block isolation. This is
    # finite corroboration of, not a replacement for, disperser.md section 7.
    block_cases = 0
    for r in range(2, 9):
        for t in range(1, r):
            sets = list(itertools.combinations(range(1, r + 1), t))
            for target in sets:
                for source in sets:
                    possible = all(
                        sum(r + 2 - i for i in target if i >= threshold)
                        <= sum(r + 2 - j for j in source if j >= threshold)
                        for threshold in range(1, r + 1)
                    )
                    assert possible == (target == source)
                    block_cases += 1
    results['block_capacity_pairs'] = block_cases

    # Exercise every guard encountered at selected actual lengths. Whenever
    # the guards pass, compare the dynamic symmetric recurrence with direct
    # subset enumeration on several payloads and an exact valid prefix.
    evaluation_cases = 0
    tested_parameters = []
    for n in (0, 1, 2, 8, 16, 32, 64, 128, 192, 256):
        params = f.parameters(n)
        if params is None:
            assert f.evaluate([0] * n) == 0
            evaluation_cases += 1
            continue
        k, m, r, t = params
        assert (k >= 2 and f.is_prime(k) and k ** 5 > 32 * n ** 3
                and k + m == n and r * k >= m and t * t >= r
                and t < r and k * t >= r * (r + 1))
        assert f.evaluate([0] * n) == 0  # Z^k is not squarefree.
        # This bounded search belongs to the test only, not to the language.
        modulus = next(
            (1 << k) | low for low in range(1, 8192, 2)
            if low < (1 << k)
            and f.irreducible_prime_degree((1 << k) | low, k)
        )
        prefix = [(modulus >> j) & 1 for j in range(k)]
        for seed in (0, 1, 2, 3):
            payload = [((j * j + 3 * j + seed) >> (j % 5)) & 1
                       for j in range(m)]
            padded = payload + [0] * (r * k - m)
            ys = []
            for i in range(r):
                value = sum(padded[i * k + j] << j for j in range(k))
                ys.append(f.field_pow(value, (1 << (r - i + 1)) - 1,
                                      modulus))
            total = 0
            for selected in itertools.combinations(range(r), t):
                term = 1
                for i in selected:
                    term = f.field_mul(term, ys[i], modulus)
                total ^= term
            assert f.evaluate(prefix + payload) == (total & 1)
            evaluation_cases += 1
        tested_parameters.append({'N': n, 'k': k, 'm': m, 'r': r, 't': t})
    results['family_evaluation_cases'] = evaluation_cases
    results['family_parameters_tested'] = tested_parameters

    return {
        'all_passed': True,
        'checks': results,
        'scope': ('Bounded exact audit only. Universal circuit coverage, '
                  'affine algebra, uniformity, and asymptotics are proved '
                  'in the supplied Markdown files.'),
    }


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
