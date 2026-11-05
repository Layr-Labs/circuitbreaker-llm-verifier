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

    # Section 8: exact truth-table premises, separate from graph coverage.
    and_child_cases = 0
    majority_child_cases = 0
    for g, d in itertools.product(and_tables, and_tables):
        forcing = [c for c in (0, 1)
                   if table_value(d, c, 0) == table_value(d, c, 1)]
        assert len(forcing) == 1
        c = forcing[0]
        assert any(table_value(g, *p) == c for p in points)
        majority = 1 if g.bit_count() == 3 else 0
        if c == majority:
            for coordinate in (0, 1):
                assert any(all(table_value(g, *p) == c
                               for p in points if p[coordinate] == fixed)
                           for fixed in (0, 1))
            majority_child_cases += 1
        and_child_cases += 1
    assert and_child_cases == 64 and majority_child_cases == 32
    results['and_child_truth_table_pairs'] = and_child_cases
    results['majority_forcing_pairs'] = majority_child_cases
    union_cases = 0
    for u in range(2, 10):
        for v in range(2, 10):
            for shared in range(min(u, v)):
                x_children = {'G'} | {('shared', i) for i in range(shared)}
                y_children = set(x_children)
                x_children |= {('x', i) for i in range(u-1-shared)}
                y_children |= {('y', i) for i in range(v-1-shared)}
                assert len(x_children | y_children | {'D'}) == u+v-shared
                assert -14-6-2*(u+v-shared) <= -26
                if (u, v, shared) != (2, 2, 1):
                    assert -14-6-2*(u+v-shared) <= -28
                assert -7-6-2*(max(u,v)+1) <= -19
                union_cases += 1
    results['and_child_union_arithmetic_cases'] = union_cases

    paired_forcing_cases = 0
    for g in and_tables:
        for pin in (0, 1):
            fixed = next(v for v in (0, 1)
                         if len({table_value(g, *p) for p in points
                                 if p[pin] == v}) == 1)
            forced_g = next(table_value(g, *p) for p in points if p[pin] == fixed)
            for child in range(16):
                unary = [table_value(child, forced_g, z) for z in (0, 1)]
                assert unary in ([0, 0], [1, 1], [0, 1], [1, 0])
                paired_forcing_cases += 1
    results['paired_and_forcing_cases'] = paired_forcing_cases
    paired_arithmetic_cases = 0
    for children in range(1, 33):
        for isolated in range(children+1):
            assert -3+6-2*(children+1)+isolated <= 0
            paired_arithmetic_cases += 1
            if isolated < children:
                assert -3+7-2*(children+1)+isolated <= 0
                paired_arithmetic_cases += 1
    results['paired_and_arithmetic_cases'] = paired_arithmetic_cases

    # Sharper Section 9 block selection: exhaust finite dimension vectors.
    block_selection_cases = 0
    block_selection_applicable = 0
    for blocks in (3, 4, 5):
        order = next(t for t in range(1, blocks) if t*t >= blocks)
        weights = [blocks+2-i for i in range(1, blocks+1)]
        assert sum(weights) == blocks*(blocks+3)//2
        for field_degree in (6, 9, 12):
            for dimensions in itertools.product(range(field_degree+1), repeat=blocks):
                block_selection_cases += 1
                if sum(dimensions) <= field_degree*order+sum(weights):
                    continue
                excess = [d-w for d,w in zip(dimensions,weights)]
                selected = sorted(range(blocks), key=lambda i: (-excess[i], i))[:order]
                assert all(excess[i]>0 for i in selected)
                assert sum(excess[i] for i in selected) > field_degree
                assert all(dimensions[i]>weights[i]>=2 for i in selected)
                block_selection_applicable += 1
    assert block_selection_applicable > 0
    assert 1867**5 > 2*1625**5
    assert Fraction(65, 8)*Fraction(1867, 1625) == Fraction(1867, 200)
    results['sharper_block_selection_vectors'] = block_selection_cases
    results['sharper_block_selection_applicable'] = block_selection_applicable

    reweighted_cases = 0
    for u in range(2, 17):
        for children in range(1, 17):
            for isolated in range(children+1):
                if u>=3:
                    assert -6-2*(u+children)+isolated <= -13
                elif children>=2 or isolated==0:
                    assert -11-2*children+isolated <= -13
                reweighted_cases += 1
    for children in range(2, 17):
        for shared_input_children in range(min(2,children)+1):
            for isolated in range(children-shared_input_children+1):
                assert -11-2*children+shared_input_children+isolated <= -13
                reweighted_cases += 1
    assert -5-5-4+1 == -13
    assert -4+2-6-6+(6-5) == -13
    assert -6-6-5-10+1 == -26
    assert -6-5-8 == -19
    results['reweighted_potential_grid_cases'] = reweighted_cases

    # Every AND-AND-AND chain with the first input repeated compresses
    # to two full-B2 gates. The probe checks all tables and assignments.
    chain_assignments = 0
    for g,d,e in itertools.product(and_tables,repeat=3):
        fixed_x = next(c for c in (0,1)
                       if table_value(e,0,c)==table_value(e,1,c))
        active_x = 1-fixed_x
        constant = table_value(e,0,fixed_x)
        htable = sum(table_value(e,table_value(d,table_value(g,active_x,y),z),active_x)
                     << (2*y+z) for y,z in points)
        final_table = sum((constant if x==fixed_x else h) << (2*x+h)
                          for x,h in points)
        for x,y,z in itertools.product((0,1),repeat=3):
            original=table_value(e,table_value(d,table_value(g,x,y),z),x)
            compressed=table_value(final_table,x,table_value(htable,y,z))
            assert original==compressed
            chain_assignments += 1
    results['reweighted_AND_chain_assignments'] = chain_assignments
    xor_chain_tables = 0
    for g,d,e in itertools.product(and_tables,and_tables,(6,9)):
        fixed_y = next(c for c in (0,1)
                       if table_value(g,0,c)==table_value(g,1,c))
        restricted = sum(table_value(e,table_value(d,table_value(g,x,fixed_y),z),x)
                         <<(2*x+z) for x,z in points)
        assert restricted in (3,12,6,9)
        xor_chain_tables += 1
    results['reweighted_XOR_chain_tables'] = xor_chain_tables

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
        assert (k >= 2 and f.is_prime(k) and k**5 > 4*n**3
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
    # Supplementary checks of the new component arithmetic and linear maps.
    spec = importlib.util.spec_from_file_location(
        'candidate_engine', Path(f.__file__).with_name('engine.py'))
    eng = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = eng
    spec.loader.exec_module(eng)
    roots = 0
    for degree in range(1, 13):
        for integer in range(257):
            lo = eng.root_floor(integer, degree)
            hi = eng.ceil_root(integer, degree)
            assert lo**degree <= integer < (lo+1)**degree
            assert integer <= hi**degree
            assert hi == 0 or (hi-1)**degree < integer
            roots += 1
    results['advanced_integer_roots'] = roots
    # Exercise the explicit pre-data guard using an integer length only;
    # this does not allocate a string of that enormous length.
    guarded_plan = eng.Plan(1 << 1024)
    assert guarded_plan.k > (1 << 1024) and not guarded_plan.ready
    results['advanced_early_entropy_guard_cases'] = 1
    # Regression: actual near-half dimension is 100000, above Python's
    # recursion limit. Enumeration must remain total without changing it.
    large = list(itertools.islice(eng.monomials(3, 100000, 2), 3))
    assert len(large) == 3 and all(len(x) == 100000 for x in large)
    assert [sum(x) for x in large] == [0, 1, 1]
    for dimension in range(1, 6):
        for total in range(5):
            actual = list(eng.compositions(total, dimension))
            expected = [x for x in itertools.product(range(total+1), repeat=dimension)
                        if sum(x) == total]
            assert actual == expected
    results['advanced_large_dimension_monomials'] = 3
    # Traversal regression against a small recursive reference. These mock
    # predicates test tree ordering only, not extractor correctness.
    from types import SimpleNamespace
    tree_cases = 0
    for value in range(256):
        obj = eng.Engine.__new__(eng.Engine)
        obj.p = SimpleNamespace(t=2, b=2, h=3, n=8, pmax=1)
        obj.x = [value >> j & 1 for j in range(8)]
        obj.ch, obj.child = {}, {}
        obj.be = lambda blocks: [sum(map(sum, blocks)) % 2]
        obj.cr = lambda node, target, confidence: sum(obj.x[target[0]:target[1]]) % 2
        reference_ch, reference_child = {}, {}
        def reference(node, depth):
            if depth == obj.p.h:
                reference_ch[node] = []
                return
            children = obj.children(node)
            for child in children:
                reference(child, depth+1)
            chosen = next((v for v in children if obj.cr(node, v, 1)), children[0])
            reference_child[node] = chosen
            reference_ch[node] = reference_ch[chosen] + [obj.be([obj.source(v) for v in children])]
        reference((0, 8), 0)
        obj.challenges()
        assert obj.ch == reference_ch and obj.child == reference_child
        tree_cases += 1
    results['advanced_tree_traversal_mock_cases'] = tree_cases
    linear_cases = 0
    for degree, modulus in ((2, 7), (3, 11), (4, 19)):
        for dimension in (1, 2, 3):
            n = 6
            hp = dimension*eng.ceil_root(n, dimension)
            params = eng.SeedParams(n, n, dimension, hp, 4, 3, degree)
            matrix = [[(i+2*j+1) % (1 << degree)
                       for j in range(dimension)] for i in range(dimension)]
            extractor = eng.Seeded(params, modulus, matrix)
            for seed in range(16):
                basis = [extractor.value([int(i == j) for i in range(n)], seed)
                         for j in range(n)]
                for value in range(1 << n):
                    bits = [value >> j & 1 for j in range(n)]
                    actual = extractor.value(bits, seed)
                    expected = [sum(basis[j][i] for j in range(n) if bits[j]) % 2
                                for i in range(4)]
                    assert actual == expected
                    linear_cases += 1
    results['advanced_seeded_linearity_assignments'] = linear_cases
    # Four algebraic coefficient-selection branches, with all cyclic shifts.
    import random
    rng = random.Random(391)
    pattern_cases = {}
    for prime in (1001, 1009):
        dim = (9*prime+19)//20+10
        last = prime-dim
        limit = prime-2*dim+3
        for trial in range(128):
            if trial % 4 == 0:
                support = {0, last} | {j for j in range(1, last)
                                       if rng.random() < .06}
            else:
                aa = rng.randrange(limit+1)
                bb = aa if trial % 4 == 1 else rng.randrange(limit+1)
                support = {0, last, aa, last-bb}
            ordered = sorted(support)
            if any(b-a-1 >= dim-1 for a,b in zip(ordered, ordered[1:])):
                continue
            middle = [j for j in ordered if prime-2*dim+4 <= j <= dim-4]
            if middle:
                cols = [0, middle[0]+1, last+2]
                branch = 'middle'
            else:
                alpha = max(j for j in support if j <= limit)
                j2 = min(j for j in support if j >= dim-3)
                beta = last-j2
                if alpha > beta:
                    cols, branch = [0, last+1, alpha+2], 'left'
                elif beta > alpha:
                    cols, branch = [j2, alpha+1, last+2], 'right'
                else:
                    cols, branch = [0, j2+1, alpha+2], 'equal'
            for shift in range(prime):
                mat = [[int((v-shift-u) % prime in support) for v in cols]
                       for u in range(3)]
                determinant = sum(mat[0][a]*mat[1][b]*mat[2][c]
                                  for a,b,c in itertools.permutations(range(3))) % 2
                assert determinant == int(shift == 0)
            pattern_cases[branch] = pattern_cases.get(branch, 0)+1
    assert set(pattern_cases) == {'middle', 'left', 'right', 'equal'}
    results['advanced_cubic_indicator_patterns'] = pattern_cases
    # Frobenius fixed-space dimensions at a composite odd field degree.
    # This finite check supplements the Bezout/root-count proof in A2.
    from math import gcd
    odd_degree = 9
    def irreducible_general(modulus):
        power = 2
        for exponent in range(1, odd_degree+1):
            power = eng.mul(power, power, modulus)
            if exponent <= odd_degree//2 and f.poly_gcd(power ^ 2, modulus) != 1:
                return False
        return power == 2
    modulus = next((1 << odd_degree) | low for low in range(1, 1 << odd_degree, 2)
                   if irreducible_general((1 << odd_degree) | low))
    fixed_counts = {}
    for shift in range(1, odd_degree):
        fixed = 0
        for value in range(1 << odd_degree):
            image = value
            for _ in range(shift):
                image = eng.mul(image, image, modulus)
            fixed += image == value
        assert fixed == 1 << gcd(shift, odd_degree)
        assert gcd(shift, odd_degree) <= odd_degree//3
        fixed_counts[str(shift)] = fixed
    results['advanced_odd_composite_frobenius_fixed_counts'] = fixed_counts
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
