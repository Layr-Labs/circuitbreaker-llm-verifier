"""Supplementary exact finite checks; no finite run proves the asymptotic theorem."""
import base64
import struct
import time
import zlib
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


def apply_matrix(rows, value):
    return sum(((row & value).bit_count() & 1) << j
               for j, row in enumerate(rows))


def matrix_pivot_audit():
    """Exhaust all invertible binary systems through three auxiliaries.

    Check a role exchange by solving both systems on every RHS, independently
    of the determinant proof. These are finite examples, not a size cutoff
    in the universal row lemma.
    """
    systems = pivots = assignments = 0
    for n in range(1, 4):
        for rows in itertools.product(range(1 << n), repeat=n):
            solutions = {apply_matrix(rows, v): v for v in range(1 << n)}
            if len(solutions) != 1 << n:
                continue
            systems += 1
            for nx in range(1, 1 << n):
                coefficient = solutions[nx]
                for column in range(n):
                    if not (coefficient >> column) & 1:
                        continue
                    # Replace auxiliary Q's column by the old x column Nx.
                    swapped = tuple((row & ~(1 << column)) |
                                    (((nx >> j) & 1) << column)
                                    for j, row in enumerate(rows))
                    new_solutions = {apply_matrix(swapped, v): v
                                     for v in range(1 << n)}
                    assert len(new_solutions) == 1 << n
                    pivots += 1
                    for rhs in range(1 << n):
                        for w in (0, 1):
                            new_rhs = rhs ^ apply_matrix(rows, w << column)
                            new_values = new_solutions[new_rhs]
                            x = (new_values >> column) & 1
                            old_values = solutions[rhs ^ (nx if x else 0)]
                            assert ((old_values >> column) & 1) == w
                            assert ((old_values ^ new_values) & ~(1 << column)) == 0
                            assignments += 1
    return dict(invertible_systems=systems, legal_pivots=pivots,
                solution_assignments=assignments)


def promised_pair_audit():
    points = ((0, 0), (0, 1), (1, 0))
    lines = tuple(itertools.combinations(points, 2))
    constant_lines = 0
    for table in range(16):
        witnesses = [(p, q) for p, q in lines
                     if table_value(table, *p) == table_value(table, *q)]
        assert witnesses
        for p, q in witnesses:
            assert p[0] == q[0] == 0 or p[1] == q[1] == 0 or (
                p[0] ^ p[1] == q[0] ^ q[1] == 1)
            constant_lines += 1
    normalization_assignments = 0
    for alpha, beta, gamma in itertools.product((0, 1), repeat=3):
        for x, y in points:
            old_x, old_y = x ^ alpha, y ^ beta
            assert ((old_x ^ alpha) & (old_y ^ beta)) ^ gamma == gamma
            normalization_assignments += 1
    # E1 and E2: fixing a partner to zero respects its old promise.
    # Every listed point lies in the old promised domain and forces G.
    forcing_assignments = 0
    for alpha, beta, gamma in itertools.product((0, 1), repeat=3):
        for y in (0, 1):
            x, hx, hy = alpha, 0, 0
            assert x * hx == y * hy == 0
            assert ((x ^ alpha) & (y ^ beta)) ^ gamma == gamma
            forcing_assignments += 1
    steps = {'A_unused': (0, -2), 'A_two_children': (0, -2-2+4),
             'B_shared': (0, -2-1), 'C_row': (0, -2-2+4),
             'D': (-1, -2+2-4-4), 'E1': (-1, -3-2-4),
             'E2': (-1, -4-4)}
    for delta_a, delta_w in steps.values():
        assert delta_w <= 8 * delta_a
    # Check only the numeric terminal inequalities of elimination.md section 3.
    # The universal proof establishes why real terminal states satisfy the guards.
    terminal_cases = 0
    for d in range(1, 25):
        for q in range(d+1):
            cap = d-1 if q == 0 else d-2 if q == 1 else d-q
            for inert in range(max(-1, cap)+1):
                assert -4*d-2*q-4*inert >= -8*d+4
                terminal_cases += 1
    refined_terminal_cases = 0
    for d in range(2, 25):
        for q in range(d+1):
            cap = d-1 if q == 0 else d-2 if q in (1, 2) else d-q
            for inert in range(cap+1):
                # Section 4 proves these two strengthened terminal obligations.
                slack = int(q == 0 and inert == d-1)
                objects = int(q == 2 and inert == d-2)
                assert -4*d-2*q-4*inert+slack+objects >= -8*d+5
                refined_terminal_cases += 1
    assert 9*8//2 > (1 << 5)-1
    return dict(refined_terminal_ledger_cases=refined_terminal_cases,
                terminal_ledger_cases=terminal_cases,
                tables=16, constant_line_witnesses=constant_lines,
                normalization_assignments=normalization_assignments,
                forcing_assignments=forcing_assignments, potential_steps=steps)


def terminal_cover_audit(path):
    started=time.monotonic()
    n,d=8,3
    legal=[a for a in range(256) if not(a&1 and a&4)
           and not(a&2 and a&8) and not(a&16 and a&32)]
    legal_set=set(legal);index={v:i for i,v in enumerate(legal)};planes=[]
    def extend(offset,points,pivots,low):
        if len(pivots)==d:
            planes.append(tuple(index[v] for v in points));return
        k=len(pivots)
        for pivot in range(low,n-(d-k)+1):
            if offset>>pivot&1:continue
            free=[j for j in range(pivot) if j not in pivots]
            for bits in range(1<<len(free)):
                v=(1<<pivot)|sum(((bits>>j)&1)<<b for j,b in enumerate(free))
                more=[x^v for x in points]
                if all(x in legal_set for x in more):
                    extend(offset,points+more,pivots+[pivot],pivot+1)
    for a in legal:extend(a,[a],[],0)
    assert len(legal)==108 and len(planes)==14391
    def trees(labels):
        if len(labels)==1:return [labels[0]]
        out=[];first=labels[0];rest=labels[1:]
        for k in range(len(rest)):
            for more in itertools.combinations(rest,k):
                left=(first,)+more;right=tuple(x for x in labels if x not in left)
                for a in trees(left):
                    for b in trees(right):out.append((a,b))
        return out
    def serial(tree,perm):
        if isinstance(tree,int):return str(perm[tree])
        a,b=serial(tree[0],perm),serial(tree[1],perm)
        return '('+','.join(sorted((a,b)))+')'
    permutations=[]
    for a,b,e in itertools.product((0,1),repeat=3):
        p=list(range(7))
        if a:p[1],p[2]=p[2],p[1]
        if b:p[3],p[4]=p[4],p[3]
        if e:p[5],p[6]=p[6],p[5]
        permutations.append(p)
    all_trees=trees(tuple(range(7)));seen=set();representatives=[]
    for tree in all_trees:
        key=min(serial(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(all_trees)==10395 and len(representatives)==1710
    certificate=json.loads(Path(path).read_text())
    assert set(certificate)=={'version','encoding','data'}
    assert certificate['version']==1 and certificate['encoding']=='zlib-base64-terminal-covers-v1'
    compressed=base64.b64decode(certificate['data'],validate=True)
    decoder=zlib.decompressobj();raw=decoder.decompress(compressed,200000)
    assert decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail
    assert raw[:4]==b'TC61' and struct.unpack_from('<H',raw,4)[0]==1710
    covers=[];pos=6
    for _ in representatives:
        count=raw[pos];pos+=1
        assert 1<=count<=255
        ids=list(struct.unpack_from('<'+'H'*count,raw,pos));pos+=2*count
        assert all(0<=i<len(planes) for i in ids)
        covers.append(ids)
    assert pos==len(raw)
    assignments=64*3**6;full=(1<<assignments)-1
    repeat=full//((1<<64)-1)
    signs=[sum(((j>>i)&1)<<j for j in range(64))*repeat for i in range(6)]
    op_masks=[]
    for i in range(6):
        stride=64*3**i
        a=(full//((1<<(3*stride))-1))*((1<<stride)-1)
        op_masks.append((a,a<<stride,a<<(2*stride)))
    total_planes=0
    for tree,ids in zip(representatives,covers):
        operations=[]
        def compile_tree(t):
            if isinstance(t,int):return t
            a,b=compile_tree(t[0]),compile_tree(t[1])
            operations.append((a,b));return 6+len(operations)
        root=compile_tree(tree);assert root==12 and len(operations)==6
        points={j for pi in ids for j in planes[pi]};values={}
        for j in points:
            a=legal[j]
            leaves=[signs[i] if not(a>>i&1) else full^signs[i] for i in range(6)]
            leaves.extend(full if a>>i&1 else 0 for i in (6,7))
            v=[leaves[0]&leaves[1]]+leaves[2:]
            for (left,right),(aa,oo,xx) in zip(operations,op_masks):
                x,y=v[left],v[right]
                v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
            values[j]=v[root]
        covered=0
        for pi in ids:
            zero=one=full
            for j in planes[pi]:
                one &= values[j];zero &= full^values[j]
            covered |= zero|one
        assert covered==full
        total_planes +=len(ids)
    return {'passed':True,'labeled_component_trees':10395,'shape_orbits':1710,
            'domain_points':108,'affine_three_flats':14391,
            'normalized_assignments_per_shape':assignments,
            'normalized_assignments_verified':assignments*len(representatives),
            'cover_plane_witnesses':total_planes,'certificate_raw_bytes':len(raw),
            'elapsed_seconds':time.monotonic()-started}



def run():
    f = load_family()
    certificate_path = Path('/snapshot/files/terminal_cover.txt')
    if not certificate_path.exists():
        certificate_path = Path(__file__).with_name('terminal_cover.txt')
    results = {"terminal_read_once_cover": terminal_cover_audit(certificate_path),
               "affine_role_exchanges": matrix_pivot_audit(),
               "promised_pair_reductions": promised_pair_audit()}
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

    # Retained historical arithmetic from the older 7/2 argument.
    # This grid does not certify the replacement promised-domain gate proof.
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
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
