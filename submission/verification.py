"""Finite diagnostics for the proof; these do not certify the universal theorem."""
import itertools
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a

def multiply(a, b, q):
    c = 0
    while b:
        if b & 1:
            c ^= a
        a <<= 1
        b >>= 1
    return remainder(c, q)

def field_test(q, k):
    if not q & 1 or q.bit_count() % 2 == 0:
        return False
    z = 2
    for _ in range(k):
        z = multiply(z, z, q)
    return z == 2

def trial_irreducible(q, k):
    for degree in range(1, k//2+1):
        for tail in range(1 << degree):
            if remainder(q, (1 << degree) | tail) == 0:
                return False
    return True

def span_basis(values):
    pivots = {}
    for x in values:
        while x:
            h = x.bit_length()-1
            if h not in pivots:
                pivots[h] = x
                break
            x ^= pivots[h]
    return list(pivots.values())

def evaluate(coeff, x, q):
    value = 0
    for a in coeff:
        value ^= multiply(a, x, q)
        x = multiply(x, x, q)
    return value

def kernel_polynomial(basis, q):
    coeff = [1]
    for v in basis:
        c = evaluate(coeff, v, q)
        assert c
        new = [0]*(len(coeff)+1)
        for j, a in enumerate(coeff):
            new[j] ^= multiply(c, a, q)
            new[j+1] ^= multiply(a, a, q)
        coeff = new
    return coeff

def image_polynomial(basis, k, q):
    p = kernel_polynomial(basis, q)
    image_basis = span_basis(evaluate(p, 1 << j, q) for j in range(k))
    assert len(image_basis) == k-len(basis)
    coeff = kernel_polynomial(image_basis, q)
    return coeff+[0]*(k-len(coeff))

def frobenius(a, s, q):
    for _ in range(s):
        a = multiply(a, a, q)
    return a

def permanent(coeff, rows, columns, k, q):
    value = 0
    for columns_in_order in itertools.permutations(columns):
        term = 1
        for u, v in zip(rows, columns_in_order):
            term = multiply(term, frobenius(coeff[(v-u) % k], u, q), q)
        value ^= term
    return value

def coefficient_test(profile, k, q, rng):
    profile = sorted(profile, key=lambda de: -de[0])
    polynomials = []
    for d, e in profile:
        basis = []
        while len(basis) < d:
            basis = span_basis(basis+[rng.randrange(1, 1 << k)])
        a = image_polynomial(basis, k, q)
        assert a[0] and a[k-d] == 1
        assert not any(a[k-d+1:])
        for start in range(k):
            assert any(a[(start+j) % k] for j in range(d))
        for start in range(k):
            if not any(a[(start+j) % k] for j in range(d-1)):
                assert start == k-d+1
        polynomials.append(a)
    threshold = k-profile[0][0]+1
    assert sum(d-e-1 for d,e in profile) >= threshold
    L = 0
    active = True
    columns = []
    zero_branch = 0
    for (d,e), a in zip(profile, polynomials):
        if not active:
            B = [0]+list(range(k-d+1,k-d+e+1))
        elif a[L]:
            B = [L]+list(range(k-d+1,k-d+e+1))
        else:
            zero_branch += 1
            b = max(j for j in range(L) if a[j])
            c = min(j for j in range(L+1,k-d+1) if a[j])
            gap = c-b-1
            assert gap <= d-2
            s = min(gap,e)
            B = [c]+list(range(b+1,b+s+1))+list(range(k-d+s+1,k-d+e+1))
        assert len(B) == e+1 and len(set(B)) == e+1
        columns.append(B)
        if active:
            L += d-e-1
            active = L < threshold
    terms = []
    for shift in range(k):
        term = 1
        for (d,e), a, B in zip(profile,polynomials,columns):
            minor = permanent(a,range(e+1),[(v-shift) % k for v in B],k,q)
            term = multiply(term,frobenius(minor,shift,q),q)
        terms.append(term)
    assert terms[0] and not any(terms[1:])
    return zero_branch


def inverse(rows, size):
    a = [row | (1 << (size+i)) for i,row in enumerate(rows)]
    for col in range(size):
        pivot = next((j for j in range(col,size) if (a[j] >> col) & 1),None)
        if pivot is None:
            return None
        a[col],a[pivot] = a[pivot],a[col]
        for j in range(size):
            if j != col and (a[j] >> col) & 1:
                a[j] ^= a[col]
    return [row >> size for row in a]

def solve(inv, rhs):
    return sum(((row & rhs).bit_count() & 1) << j for j,row in enumerate(inv))

def matvec(rows, vector):
    return sum(((row & vector).bit_count() & 1) << j for j,row in enumerate(rows))

def dropbit(value, index):
    return (value & ((1 << index)-1)) | ((value >> (index+1)) << index)

def sparse_system_checks():
    systems = exchanges = short_deletions = unit_deletions = assignments = 0
    for r in (1,2,3):
        for m in (1,2,3):
            for rows in itertools.product(range(1 << r),repeat=r):
                inv = inverse(rows,r)
                if inv is None:
                    continue
                choices = [
                    [v for v in range(1 << m) if v.bit_count()+row.bit_count() <= 3]
                    for row in rows
                ]
                for pars in itertools.product(*choices):
                    systems += 1
                    constants = range(1 << r) if r <= 2 else (0,(1 << r)-1)
                    for x in range(m):
                        nx = sum(((v >> x) & 1) << i for i,v in enumerate(pars))
                        dependence = solve(inv,nx)
                        for j in range(r):
                            if not (dependence >> j) & 1:
                                continue
                            swapped_a = [
                                (v & ~(1 << j)) | (((nx >> i) & 1) << j)
                                for i,v in enumerate(rows)
                            ]
                            swapped_n = [
                                (v & ~(1 << x)) | (((rows[i] >> j) & 1) << x)
                                for i,v in enumerate(pars)
                            ]
                            assert inverse(swapped_a,r) is not None
                            assert all(a.bit_count()+n.bit_count() <= 3 for a,n in zip(swapped_a,swapped_n))
                            exchanges += 1
                            for c in constants:
                                seen = set()
                                for u in range(1 << m):
                                    v = solve(inv,matvec(pars,u)^c)
                                    new_u = (u & ~(1 << x)) | (((v >> j) & 1) << x)
                                    new_v = (v & ~(1 << j)) | (((u >> x) & 1) << j)
                                    assert matvec(swapped_a,new_v) == matvec(swapped_n,new_u)^c
                                    seen.add(new_u); assignments += 1
                                assert len(seen) == 1 << m
                            if nx.bit_count() == 1:
                                e = nx.bit_length()-1
                                reduced = [dropbit(v,j) for i,v in enumerate(swapped_a) if i != e]
                                assert inverse(reduced,r-1) is not None
                                unit_deletions += 1
                    for e,(ar,nr) in enumerate(zip(rows,pars)):
                        if ar.bit_count()+nr.bit_count() > 2:
                            continue
                        j = (ar & -ar).bit_length()-1
                        new_a=[]; new_n=[]
                        for i,(a,n) in enumerate(zip(rows,pars)):
                            if i == e:
                                continue
                            if (a >> j) & 1:
                                a ^= ar; n ^= nr
                            new_a.append(dropbit(a,j)); new_n.append(n)
                        new_inv = inverse(new_a,r-1)
                        assert new_inv is not None
                        assert all(a.bit_count()+n.bit_count() <= 3 for a,n in zip(new_a,new_n))
                        for c in constants:
                            new_c = sum((((c >> i) & 1)^(((rows[i] >> j)&1)*((c >> e)&1))) << k
                                        for k,i in enumerate(i for i in range(r) if i != e))
                            for u in range(1 << m):
                                v=solve(inv,matvec(pars,u)^c)
                                assert matvec(new_a,dropbit(v,j)) == matvec(new_n,u)^new_c
                                assignments += 1
                        short_deletions += 1
    return {
        'sparse_fair_systems_r_m_1_through_3':systems,
        'parameter_auxiliary_exchanges':exchanges,
        'short_row_deletions':short_deletions,
        'unit_parameter_cofactor_deletions':unit_deletions,
        'explicit_assignment_equivalences':assignments,
        'constants_scope':'all constants for r<=2, zero and all-ones constants for r=3'
    }

def assert_source(table, m, pairs):
    assert len(table) == 1 << m and len(set(table)) == len(table)
    used=[x for pair in pairs for x in pair]
    assert len(used) == len(set(used))
    coeff=list(table)
    for bit in range(m):
        for u in range(1 << m):
            if (u >> bit) & 1:
                coeff[u] ^= coeff[u ^ (1 << bit)]
    allowed={sum(1 << x for x in pair) for pair in pairs}
    assert all(v == 0 or mask.bit_count() <= 1 or mask in allowed for mask,v in enumerate(coeff))

def source_restriction(table,m,pairs,x,coeff,const,remove_pair=None):
    rest=[i for i in range(m) if i != x]
    new_index={old:i for i,old in enumerate(rest)}
    new_table=[]
    for u in range(1 << (m-1)):
        old=sum(((u >> j)&1) << i for j,i in enumerate(rest))
        bit=(coeff & old).bit_count() % 2 ^ const
        old |= bit << x
        new_table.append(table[old])
    new_pairs=[tuple(new_index[i] for i in pair) for pair in pairs if pair != remove_pair]
    assert_source(new_table,m-1,new_pairs)

def coordinate_checks():
    qr=ar=uf=pf=0
    m=4
    for x,y,z in itertools.permutations(range(m),3):
        for alpha,beta,gamma,c in itertools.product((0,1),repeat=4):
            def gate(u):
                return ((((u >> x)&1)^alpha)*(((u >> y)&1)^beta))^gamma
            table=[]
            for u in range(1 << m):
                old=u ^ ((gate(u)^c) << z)
                assert gate(old)^((old >> z)&1)^c == ((u >> z)&1)
                table.append(old)
            pairs=[tuple(sorted((x,y)))]
            assert_source(table,m,pairs); qr += 1
            for v in range(m):
                if v in (x,y):
                    continue
                for coeff in range(1 << m):
                    if (coeff >> v)&1:
                        continue
                    for const in (0,1):
                        changed=[table[u ^ ((((coeff & u).bit_count()&1)^const) << v)] for u in range(1 << m)]
                        assert_source(changed,m,pairs); ar += 1
                        source_restriction(table,m,pairs,v,coeff,const); uf += 1
            for p,r in ((x,y),(y,x)):
                for const in (0,1):
                    source_restriction(table,m,pairs,p,0,const,remove_pair=pairs[0]); pf += 1
                    source_restriction(table,m,pairs,p,1 << r,const,remove_pair=pairs[0]); pf += 1
    return {'quadratic_reparameterizations':qr,'affine_reparameterizations':ar,
            'unpaired_affine_restrictions':uf,'paired_restrictions':pf,
            'scope':'all Boolean assignments to four parameters, all signed bottleneck tables, and stated subsequent coordinate operations'}

def resource_assignment_checks():
    # A generic integral max flow, independent of the moment argument in D7.
    def feasible(r,I,J):
        if not I:
            return True
        n=2*r+2; source=2*r; sink=source+1
        cap=[[0]*n for _ in range(n)]
        demand=sum(r+2-i for i in I)
        for j in J:
            cap[source][j-1]=r+2-j
            for i in I:
                if i <= j:
                    cap[j-1][r+i-1]=demand
        for i in I:
            cap[r+i-1][sink]=r+2-i
        flow=0
        while True:
            parent=[-1]*n; parent[source]=source; queue=[source]
            for u in queue:
                for v,c in enumerate(cap[u]):
                    if c and parent[v] < 0:
                        parent[v]=u; queue.append(v)
                if parent[sink] >= 0:
                    break
            if parent[sink] < 0:
                return flow == demand
            v=sink; amount=demand
            while v != source:
                amount=min(amount,cap[parent[v]][v]); v=parent[v]
            v=sink
            while v != source:
                u=parent[v]; cap[u][v]-=amount; cap[v][u]+=amount; v=u
            flow += amount
    count=0
    for r in range(1,9):
        for t in range(r+1):
            subsets=list(itertools.combinations(range(1,r+1),t))
            for I in subsets:
                for J in subsets:
                    assert feasible(r,I,J) == (I == J)
                    count += 1
    return count

def run():
    count=0
    for k in (3,5,7):
        for tail in range(1 << k):
            q=(1 << k)|tail
            assert field_test(q,k) == trial_irreducible(q,k)
            count += 1
    lines=[[(0,0),(0,1)],[(1,0),(1,1)],[(0,0),(1,0)],
           [(0,1),(1,1)],[(0,0),(1,1)],[(0,1),(1,0)]]
    for a in range(16):
        if a.bit_count() not in (1,3):
            continue
        for b in range(16):
            assert any(len({(a >> (2*x+y))&1 for x,y in line}) == 1 and
                       len({(b >> (2*x+y))&1 for x,y in line}) == 1 for line in lines)
    rng=random.Random(20261017); zero_branches=0; tests=0
    for k,q,profile,repeats in [
        (7,0b10000011,[(4,1),(4,1)],500),
        (7,0b10000011,[(5,2),(4,1)],20),
        (11,0b100000000101,[(6,3)]*3,20),
        (11,0b100000000101,[(8,3)],20)]:
        assert field_test(q,k)
        for _ in range(repeats):
            zero_branches += coefficient_test(profile,k,q,rng); tests += 1
    assert zero_branches > 0
    result={
        'field_validity_comparisons':count,'and_other_truth_table_pairs':128,
        'independent_subspace_coefficient_instances':tests,
        'internal_zero_column_branches_exercised':zero_branches,
        'sparse_affine_systems':sparse_system_checks(),
        'paired_source_coordinate_operations':coordinate_checks(),
        'resource_assignment_pairs_r_1_through_8':resource_assignment_checks()
    }
    claim=json.loads((ROOT/'submission.json').read_text())['claims'][0]
    assert claim['score'] == {'numerator':'7','denominator':'2'}
    assert math.gcd(7,2) == 1
    for evidence in claim['evidence']:
        p=ROOT/evidence
        assert p.is_file() and p.resolve().is_relative_to(ROOT)
    files=[p for p in ROOT.rglob('*') if p.is_file()]
    assert len(files) <= 32 and all(p.stat().st_size <= 65536 for p in files)
    assert sum(p.stat().st_size for p in files) <= 262144
    result['submission_file_bytes']={str(p.relative_to(ROOT)):p.stat().st_size for p in files}
    result['scope']='Finite diagnostics only. They do not certify the universal circuit theorem; its proof is in proof.md.'
    print(json.dumps(result,sort_keys=True))

if __name__ == '__main__':
    run()
