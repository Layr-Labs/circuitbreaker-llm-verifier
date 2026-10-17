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

def run():
    count = 0
    for k in (3,5,7):
        for tail in range(1 << k):
            q = (1 << k) | tail
            assert field_test(q,k) == trial_irreducible(q,k)
            count += 1
    binary_lines = [
        [(0,0),(0,1)],[(1,0),(1,1)],
        [(0,0),(1,0)],[(0,1),(1,1)],
        [(0,0),(1,1)],[(0,1),(1,0)]
    ]
    and_tables = [a for a in range(16) if a.bit_count() in (1,3)]
    for a in and_tables:
        for b in range(16):
            assert any(len({(a >> (2*x+y)) & 1 for x,y in line}) == 1 and
                       len({(b >> (2*x+y)) & 1 for x,y in line}) == 1
                       for line in binary_lines)
    rng = random.Random(20261017)
    zero_branches = 0
    tests = 0
    for k, q, profile, repeats in [
        (7,0b10000011,[(4,1),(4,1)],500),
        (7,0b10000011,[(5,2),(4,1)],20),
        (11,0b100000000101,[(6,3)]*3,20),
        (11,0b100000000101,[(8,3)],20),
    ]:
        assert field_test(q,k)
        for _ in range(repeats):
            zero_branches += coefficient_test(profile,k,q,rng)
            tests += 1
    assert zero_branches > 0, "diagnostic did not exercise the internal-zero branch"
    j = json.loads((ROOT/'submission.json').read_text())
    claim = j['claims'][0]
    score = claim['score']
    assert score == {'numerator':'31','denominator':'10'}
    assert math.gcd(int(score['numerator']),int(score['denominator'])) == 1
    for evidence in claim['evidence']:
        p = ROOT/evidence
        assert p.is_file() and p.resolve().is_relative_to(ROOT)
    files = [p for p in ROOT.rglob('*') if p.is_file()]
    assert len(files) <= 32
    assert all(p.stat().st_size <= 65536 for p in files)
    assert sum(p.stat().st_size for p in files) <= 262144
    print(json.dumps({
        'field_validity_comparisons':count,
        'and_other_table_pairs':128,
        'independent_subspace_coefficient_instances':tests,
        'internal_zero_column_branches_exercised':zero_branches,
        'submission_file_bytes':{str(p.relative_to(ROOT)):p.stat().st_size for p in files},
        'scope':'Finite diagnostics only; universal coverage is the written proof.'
    },sort_keys=True))

if __name__ == '__main__':
    run()
