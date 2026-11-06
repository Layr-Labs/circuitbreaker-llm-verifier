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

    # Verify legal promised-pair lines for every binary truth table.
    pair_lines = (((0,0),(0,1)), ((0,0),(1,0)), ((0,1),(1,0)))
    for table in range(16):
        assert any(table_value(table,*p)==table_value(table,*q)
                   for p,q in pair_lines)
    results['promised_pair_tables'] = 16
    # Exact whole-step arithmetic for V=s+2q+4I. Universal
    # coverage and distinct-object deletions are proved in elimination.md.
    steps = [('unused',0,0,-1,0), ('two_children',0,-2,-1,1),
             ('common',0,-1,-1,0), ('row',0,-2,-1,1),
             ('D',-1,-2,1,-2), ('E1',-1,-2,-1,-1),
             ('E2',-1,-4,-2,0)]
    for name,da,ds,dq,di in steps:
        assert ds+2*dq+4*di <= 8*da
    results['promised_potential_cases'] = len(steps)

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

    # Selection estimate: for small dimension vectors satisfying the
    # new strict threshold, check the selected high blocks meet the
    # coefficient-isolation hypothesis. No extrapolation is asserted.
    selection_cases = 0
    for r in range(2,6):
        for k in range(2,9):
            for t in range(1,r):
                threshold = max(k*t, Fraction(k*r,t))+Fraction(r*(r+3),2)
                for dims in itertools.product(range(k+1), repeat=r):
                    if sum(dims) <= threshold:
                        continue
                    high = sorted((v-(r+1-i) for i,v in enumerate(dims) if v>r+1-i),reverse=True)
                    assert len(high)>=t
                    assert sum(high[:t])>k
                    selection_cases += 1
    # Nonvacuous arbitrary-sized examples for the same scalar implication.
    for r in range(2,15):
        k=10*r*r
        for t in range(2,r):
            dims=[k]*r
            threshold=max(k*t,Fraction(k*r,t))+Fraction(r*(r+3),2)
            if sum(dims)>threshold:
                assert sum(sorted((v-(r+1-i) for i,v in enumerate(dims)),reverse=True)[:t])>k
                selection_cases += 1
    assert selection_cases>0
    results['refined_selection_cases'] = selection_cases

    full_selection_cases=0
    for r in range(2,6):
        H=Fraction(r*(r+3),2)
        for k in range(2,9):
            for t in range(1,r):
                capacity=sum(max(k-j,0) for j in range(2,t+1))
                U=sum(max(k-j,0) for j in range(1,r+1))
                threshold=H+max(capacity,Fraction(U,t+1))
                for dims in itertools.product(range(k+1),repeat=r):
                    if sum(dims)<=threshold:
                        continue
                    candidates=[(v-(r+1-i),v) for i,v in enumerate(dims) if v>r+1-i]
                    assert len(candidates)>=t
                    assert any(sum(delta for delta,v in chosen)+max(v for delta,v in chosen)>k+1
                               for chosen in itertools.combinations(candidates,t))
                    full_selection_cases+=1
    assert full_selection_cases>0
    results['full_hypothesis_selection_cases']=full_selection_cases
    sharp_dimension_cases=0
    for r in range(2,41):
        t=1
        while t*t<r:
            t+=1
        H=r*(r+3)//2
        for k in range(2,81):
            u=min(t,k-1)
            capacity=(u-1)*k-u*(u+1)//2+1
            assert capacity==sum(max(k-j,0) for j in range(2,t+1))
            U=sum(max(k-j,0) for j in range(1,r+1))
            sharp=H+max(capacity,U//(t+1))+1
            base=k*t+H+1
            assert (k*t+t)//(t+1)<=base-sharp<=k+t*(t+1)//2
            sharp_dimension_cases+=1
    results['exact_sharp_dimension_cases']=sharp_dimension_cases

    # Exact certificates for the new irrational leading-loss comparison.
    assert 23**5>2*20**5
    assert 25*979**2<24*1000**2
    assert Fraction(23,16)<Fraction(19,20)*(Fraction(979,1000)+Fraction(625,1152))
    results['direct_degree_constant_certificates']=3

    # Exact terminal ledger arithmetic; inertness persistence and the
    # excluded one-pair case are universal arguments in elimination.md.
    ledger_cases=0
    for d in range(1,65):
        for q in range(d+1):
            for inert in range(min(d-1,d-q)+1):
                if q==1 and inert==d-1:
                    continue
                n=4*d+17
                lower=4*n-4*d-2*q-4*inert
                assert lower>=4*n-8*d+4
                ledger_cases+=1
    results['inert_terminal_ledger_cases']=ledger_cases

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
        old_degree = 2
        while old_degree**5<=4*n**3 or not f.is_prime(old_degree):
            old_degree+=1
        assert k==old_degree
        assert k**5>4*n**3
        assert (k >= 2 and f.is_prime(k)
                and k+m == n and r*k >= m and t*t >= r and t < r
                )
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
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
