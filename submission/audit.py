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

    def anf_degree(table,width):
        coefficients=[(table>>a)&1 for a in range(1<<width)]
        for bit in range(width):
            for a in range(1<<width):
                if(a>>bit)&1:
                    coefficients[a]^=coefficients[a^(1<<bit)]
        return max((a.bit_count() for a,c in enumerate(coefficients) if c),default=0)

    # All Boolean functions on F_8: compare reduced field-polynomial
    # binary weight with exact binary algebraic normal-form degree.
    modulus=0b1011
    indicators=[]
    for a in range(8):
        poly=[f.field_pow(a,7-j,modulus) for j in range(8)]
        poly[0]^=1
        indicators.append(poly)
    for table in range(256):
        poly=[0]*8
        for a in range(8):
            if(table>>a)&1:
                poly=[x^y for x,y in zip(poly,indicators[a])]
        weight=max((j.bit_count() for j,c in enumerate(poly) if c),default=0)
        assert weight==anf_degree(table,3)
    results['field_weight_binary_degree_cases']=256

    # Deterministic small straight-line probes for the multiplication
    # degree bound. Universal induction is supplied in elimination.md.
    state=18317
    def next_int(limit):
        nonlocal state
        state=(1103515245*state+12345)&((1<<31)-1)
        return state%limit
    for trial in range(1000):
        values=[sum(((a>>i)&1)<<a for a in range(16)) for i in range(4)]
        multiplications=0
        for step in range(5):
            op=next_int(16);x=values[next_int(len(values))];y=values[next_int(len(values))]
            out=sum(table_value(op,(x>>a)&1,(y>>a)&1)<<a for a in range(16))
            values.append(out)
            multiplications+=op.bit_count() in(1,3)
            assert anf_degree(out,4)<=multiplications+1
    results['multiplication_degree_probe_circuits']=1000

    degree_ledger_cases=0
    for d in range(1,33):
        for B in range(1,d+1):
            for q in range(d+1):
                for inert in range(min(d-B,d-q)+1):
                    n=4*d+17
                    lower=4*n-4*d+(B-1)-2*q-4*inert
                    assert lower>=4*n-8*d+3*B-1
                    degree_ledger_cases+=1
    results['degree_sensitive_terminal_cases']=degree_ledger_cases

    # Three allowed lines: their highest-degree parities cannot all
    # equal one. Exhaust all small tables on a pair times free bits.
    parity_cases=0
    for free_bits in range(3):
        width=1<<free_bits
        for table in range(1<<(3*width)):
            parities=[]
            for a,b in ((0,1),(0,2),(1,2)):
                parities.append(sum((table>>(p*width+j))&1 for p in (a,b) for j in range(width))%2)
            assert parities[0]^parities[1]^parities[2]==0
            assert not all(parities)
            parity_cases+=1
    results['three_line_degree_parity_tables']=parity_cases
    strong_terminal_cases=0
    for d in range(1,33):
        for B in range(1,d+1):
            for q in range(d+1):
                cap=min(d-B-(1 if q else 0),d-q)
                for inert in range(cap+1):
                    r=d-inert;a0=r-q
                    charged=max(B-1,2*B-a0-1)
                    n=4*d+17
                    lower=4*n-4*d+charged-2*q-4*inert
                    assert lower>=4*n-8*d+min(5*B-1,4*B+1)
                    strong_terminal_cases+=1
    results['strong_degree_terminal_cases']=strong_terminal_cases

    # Independently compute adjacent-order certificates from literal sums.
    # These are bounded checks of the implemented selector, not its limit proof.
    from math import isqrt
    order_cases = 0
    for k in range(2, 36):
        for r in range(3, 56):
            lower = isqrt(r)
            upper = lower + (lower*lower != r)
            scores = {}
            for order in {lower, upper}:
                H = r*(r+3)//2
                T = sum(max(k-j, 0) for j in range(2, order+1))
                U = sum(max(k-j, 0) for j in range(1, r+1))
                d = H + max(T, U//(order+1)) + 1
                B = order*(order+3)//2
                scores[order] = 8*d - 4*B - 1
                assert f.certificate_score(k, r, order) == scores[order]
                order_cases += 1
            chosen = min(scores, key=lambda order: (scores[order], order))
            assert scores[chosen] <= scores[upper]
    results['adjacent_order_certificate_cases'] = order_cases
    square_cases = 0
    for u in range(9, 80):
        k, r = 2*u**3, u*u+1
        H = r*(r+3)//2
        U = r*k-r*(r+1)//2
        dimensions = []
        for order in (u, u+1):
            T = (order-1)*k-order*(order+1)//2+1
            assert (order+1)*T >= U
            dimensions.append(H+T+1)
        assert dimensions[1]-dimensions[0] == k-u-1
        assert f.certificate_score(k,r,u+1)-f.certificate_score(k,r,u) == 8*k-12*u-16 > 0
        square_cases += 1
    results['square_adjacent_saving_cases'] = square_cases

    # Independent literal-sum branch oracle, including every prefix bit.
    def literal_branch(n, k):
        if k >= n:
            return None
        m = n-k
        r = (m+k-1)//k
        lo = isqrt(r)
        hi = lo+(lo*lo != r)
        if hi >= r:
            return None
        choices = []
        for t in {lo,hi}:
            H = r*(r+3)//2
            T = sum(max(k-j,0) for j in range(2,t+1))
            U = sum(max(k-j,0) for j in range(1,r+1))
            d = H+max(T,U//(t+1))+1
            B = t*(t+3)//2
            choices.append((4*k+8*d-4*B-1,k,t))
        _,_,t = min(choices)
        return k,m,r,t

    def literal_parameters(n):
        if not n:
            return None
        B = 1
        while B**5 <= 4*n**3:
            B += 1
        above = B
        while not f.is_prime(above):
            above += 1
        plus = literal_branch(n,above)
        if plus is None:
            return None
        options = [plus]
        below = B-1
        while below>=2 and not f.is_prime(below):
            below -= 1
        if below>=2:
            minus=literal_branch(n,below)
            if minus is not None:
                options.append(minus)
        def key(params):
            k,_,r,t=params
            H=r*(r+3)//2
            T=sum(max(k-j,0) for j in range(2,t+1))
            U=sum(max(k-j,0) for j in range(1,r+1))
            d=H+max(T,U//(t+1))+1
            return 4*k+8*d-2*t*(t+3)-1,k,t
        return min(options,key=key)

    for n in range(401):
        assert f.parameters(n)==literal_parameters(n)
    results['dual_prime_parameter_cases']=401
    prime_length_cases=0
    common_block_cases=0
    strict_dual_cases=0
    for q in range(101,1000):
        if not f.is_prime(q):
            continue
        lo,hi=0,q*q
        while lo+1<hi:
            mid=(lo+hi)//2
            if 4*mid**3>=q**5:
                hi=mid
            else:
                lo=mid
        n=hi
        assert q**5<4*n**3<(q+1)**5
        p=q+1
        while not f.is_prime(p):
            p+=1
        plus=literal_branch(n,p)
        minus=literal_branch(n,q)
        actual=f.parameters(n)
        assert actual==literal_parameters(n)
        prime_length_cases+=1
        if plus[2]==minus[2]:
            r,t=plus[2],plus[3]
            delta=p-q
            old=4*p+f.certificate_score(p,r,t)
            same_order=4*q+f.certificate_score(q,r,t)
            assert old-same_order >= (8*t-20)*delta
            current=4*actual[0]+f.certificate_score(actual[0],actual[2],actual[3])
            assert current<=same_order
            common_block_cases+=1
            if t>=3:
                assert current<old
                strict_dual_cases+=1
    results['prime_derived_actual_lengths']=prime_length_cases
    results['common_block_dual_cases']=common_block_cases
    results['strict_dual_certificate_cases']=strict_dual_cases

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
        assert f.parameters(n)==literal_parameters(n)
        assert old_degree**5>4*n**3
        assert (k >= 2 and f.is_prime(k)
                and k+m == n and r*k >= m and t in {isqrt(r), isqrt(r)+1} and t < r
                )
        lower = isqrt(r)
        upper = lower + (lower*lower != r)
        assert t == min({lower,upper}, key=lambda order: (f.certificate_score(k,r,order),order))
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
