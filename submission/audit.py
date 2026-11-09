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


def affine_flat_indices(legal,n,d):
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
    return planes


def decode_terminal_blob(envelope,version,encoding,maximum):
    assert set(envelope)=={'version','encoding','data'}
    assert envelope['version']==version and envelope['encoding']==encoding
    compressed=base64.b64decode(envelope['data'],validate=True)
    decoder=zlib.decompressobj();raw=decoder.decompress(compressed,maximum)
    assert decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail
    return raw


def transform_coloring(value,permutation,flip):
    return sum((((value>>i)&1)^flip)<<j for i,j in enumerate(permutation))


def canonical_colorings(group):
    remaining=set(range(4096));representatives=[]
    while remaining:
        value=min(remaining)
        orbit={transform_coloring(value,p,f) for p,f in group}
        assert value==min(orbit)
        representatives.append(value);remaining.difference_update(orbit)
    return representatives


def ramsey84_normalized_roots():
    legal=[a for a in range(256) if not(a&1 and a&4)
           and not(a&2 and a&8) and not(a&16 and a&32)]
    planes=affine_flat_indices(legal,8,3)
    assert len(legal)==108 and len(planes)==14391
    ids={};values=[]
    for a in legal:
        key=(int(not(a&1) and not(a&2)),a>>2)
        if key not in ids:ids[key]=len(ids)+1
        values.append(ids[key])
    assert len(ids)==84
    edges=sorted({tuple(sorted({values[j] for j in plane})) for plane in planes})
    assert len(edges)==11349
    base=[]
    for edge in edges:base.extend([list(edge),[-v for v in edge]])
    small=[a for a in range(16) if not(a&1 and a&2)]
    small_planes=[sum(1<<small.index(x) for x in points)
        for points in itertools.combinations(small,4)
        if points[0]^points[1]^points[2]^points[3]==0]
    bad={v for v in range(4096) if not any(v&m in (0,m) for m in small_planes)}
    assert len(small_planes)==39 and len(bad)==48
    core=238;orbit=set();stabilizer=[]
    for pair in itertools.permutations((0,1,2)):
        for r1,r2 in itertools.permutations((1,2,3),2):
            for shifts in itertools.product(range(4),repeat=3):
                permutation=[]
                for a in small:
                    b=a>>2
                    free=((b&r1).bit_count()&1)|(((b&r2).bit_count()&1)<<1)
                    free ^=shifts[a&3]
                    permutation.append(small.index(pair[a&3]|(free<<2)))
                for flip in (0,1):
                    result=transform_coloring(core,permutation,flip);orbit.add(result)
                    if result==core:stabilizer.append((permutation,flip))
    assert orbit==bad and len(stabilizer)==96
    zeros=[i for i in range(12) if not(core>>i&1)]
    ones=[i for i in range(12) if core>>i&1]
    allowed=set()
    for iz,io in itertools.product([-1]+zeros,[-1]+ones):
        code=4095^core
        for j in (iz,io):
            if j>=0:code ^=1<<j
        allowed.add(code)
    assert len(allowed)==49
    cases=[];remaining=set(allowed)
    while remaining:
        code=min(remaining)
        current={transform_coloring(code,p,f) for p,f in stabilizer}
        assert current<=allowed
        cases.append(code);remaining.difference_update(current)
    assert cases==[1809,1811,1813,3857]
    core0=[ids[(0,t<<2)] for t in small]
    core1=[ids[(1,t<<2)] for t in small]
    A=[ids[(0,1|(t<<2))] for t in small]
    B=[ids[(0,2|(t<<2))] for t in small]
    translations=[[small.index(a^(shift<<2)) for a in small] for shift in range(4)]
    canonical_B=canonical_colorings([(p,0) for p in translations])
    assert len(canonical_B)==1072
    outputs=[]
    for code in cases:
        current=[(p,f) for p,f in stabilizer if transform_coloring(code,p,f)==code]
        group=[]
        for p,f in current:
            for shift in range(4):
                group.append(([small.index(small[j]^(shift<<2)) for j in p],f))
        canonical_A=canonical_colorings(group)
        clauses=[list(q) for q in base]
        for j,v in enumerate(core0):clauses.append([v if core>>j&1 else -v])
        for j,v in enumerate(core1):clauses.append([v if code>>j&1 else -v])
        def prefix_exclusions(codes,variables):
            def visit(xs,bit,prefix):
                if not xs:
                    clauses.append([-variables[j] if b else variables[j] for j,b in prefix]);return
                if len(xs)==1<<(bit+1):return
                assert bit>=0
                for b in (0,1):
                    visit([x for x in xs if (x>>bit&1)==b],bit-1,prefix+[(bit,b)])
            visit(codes,11,[])
        prefix_exclusions(canonical_A,A);prefix_exclusions(canonical_B,B)
        previous=None;nv=84
        for j in range(11,-1,-1):
            a,b=A[j],B[j];prefix=[] if previous is None else [-previous]
            clauses.append(prefix+[-a,b])
            if j:
                nv+=1
                clauses.extend([prefix+[a,b,nv],prefix+[-a,-b,nv]])
                previous=nv
        assert nv==95
        units={};work=[tuple(q) for q in clauses]
        while True:
            new=[];found={}
            for clause in work:
                if any(abs(x) in units and units[abs(x)]==(x>0) for x in clause):continue
                clause=tuple(x for x in clause if abs(x) not in units)
                assert clause
                if len(clause)==1:
                    variable=abs(clause[0]);value=clause[0]>0
                    assert variable not in found or found[variable]==value
                    found[variable]=value
                else:new.append(clause)
            work=new
            if not found:break
            for variable,value in found.items():
                assert variable not in units or units[variable]==value
            units.update(found)
        roots=sorted(set(tuple(sorted(q,key=lambda x:(abs(x),x))) for q in work),
                     key=lambda q:(len(q),q))
        outputs.append((code,roots,len(canonical_A),len(units)))
    assert [len(x[1]) for x in outputs]==[5130,5049,4773,4408]
    assert [x[2] for x in outputs]==[270,364,262,49]
    return outputs


class ReverseUnitPropagation:
    def __init__(self,roots,nv=95):
        self.nv=nv;self.clauses=[]
        self.watches={v:[] for v in range(-nv,nv+1) if v}
        self.units=set();self.visits=0
        for clause in roots:self.add(clause)

    def add(self,clause):
        assert clause and all(1<=abs(v)<=self.nv for v in clause)
        assert len({abs(v) for v in clause})==len(clause)
        clause=list(clause);index=len(self.clauses);self.clauses.append(clause)
        if len(clause)==1:self.units.add(clause[0])
        else:
            self.watches[clause[0]].append(index)
            self.watches[clause[1]].append(index)

    def proves(self,claim):
        assert all(1<=abs(v)<=self.nv for v in claim)
        assert len({abs(v) for v in claim})==len(claim)
        assignment=[0]*(self.nv+1);queue=[]
        def assign(literal):
            variable=abs(literal);sign=1 if literal>0 else -1
            if assignment[variable]:return assignment[variable]==sign
            assignment[variable]=sign;queue.append(literal);return True
        for literal in self.units:
            if not assign(literal):return True
        for literal in claim:
            if not assign(-literal):return True
        q=0
        while q<len(queue):
            false=-queue[q];q+=1;watch=self.watches[false];i=0
            while i<len(watch):
                self.visits+=1;index=watch[i];clause=self.clauses[index]
                if clause[0]==false:clause[0],clause[1]=clause[1],clause[0]
                assert clause[1]==false
                other=clause[0]
                value=assignment[abs(other)]*(1 if other>0 else -1)
                if value==1:i+=1;continue
                for k in range(2,len(clause)):
                    literal=clause[k]
                    if assignment[abs(literal)]*(1 if literal>0 else -1)!=-1:break
                else:
                    if value==-1:return True
                    assert assign(other);i+=1;continue
                clause[1],clause[k]=clause[k],clause[1]
                watch[i]=watch[-1];watch.pop()
                self.watches[clause[1]].append(index)
        return False


def ramsey84_rup_audit(envelope):
    started=time.monotonic()
    raw=decode_terminal_blob(envelope,2,'zlib-base64-rup-ramsey84-nohints',200000)
    assert raw[:5]==b'RU84\x04';position=5
    def integer():
        nonlocal position
        value=0
        for shift in range(0,35,7):
            assert position<len(raw)
            b=raw[position];position+=1;value |= (b&127)<<shift
            if b<128:return value
        raise AssertionError('oversized integer')
    cases=[]
    for code,roots,canonical_A,unit_count in ramsey84_normalized_roots():
        assert position+2<=len(raw) and struct.unpack_from('<H',raw,position)[0]==code
        position+=2;count=integer();assert 1<=count<=len(roots)
        selected=[];index=0
        for _ in range(count):
            delta=integer();assert delta>0;index+=delta
            assert index<=len(roots);selected.append(roots[index-1])
        checker=ReverseUnitPropagation(selected)
        lemmas=integer();assert 1<=lemmas<=50000
        for j in range(lemmas):
            size=integer();assert 0<=size<=95 and position+size<=len(raw)
            encoded=raw[position:position+size];position+=size
            assert all(v<190 for v in encoded)
            claim=[((v>>1)+1)*(-1 if v&1 else 1) for v in encoded]
            assert checker.proves(claim)
            if claim:
                assert j<lemmas-1;checker.add(claim)
            else:assert j==lemmas-1
        assert not claim
        cases.append({'code':code,'root_clauses':len(roots),'selected_roots':count,
                      'derived_clauses':lemmas,'watch_visits':checker.visits,
                      'canonical_A_values':canonical_A,'forced_units':unit_count})
    assert position==len(raw)
    return {'passed':True,'cases':cases,'raw_bytes':len(raw),
            'derived_clauses':sum(x['derived_clauses'] for x in cases),
            'elapsed_seconds':time.monotonic()-started}


def readonce_seven_audit(envelope):
    started=time.monotonic()
    legal=[a for a in range(128) if not(a&1 and a&2) and not(a&4 and a&8)]
    planes=affine_flat_indices(legal,7,3)
    assert len(legal)==72 and len(planes)==5049
    def trees(labels):
        if len(labels)==1:return [labels[0]]
        out=[];first=labels[0];rest=labels[1:]
        for k in range(len(rest)):
            for more in itertools.combinations(rest,k):
                left=(first,)+more;right=tuple(x for x in labels if x not in left)
                for a in trees(left):
                    for b in trees(right):out.append((a,b))
        return out
    def serial(tree,permutation):
        if isinstance(tree,int):return str(permutation[tree])
        a,b=serial(tree[0],permutation),serial(tree[1],permutation)
        return '('+','.join(sorted((a,b)))+')'
    permutations=[]
    for swap,a,b in itertools.product((0,1),repeat=3):
        for free in itertools.permutations((4,5,6)):
            p=[0,1,2,3]+list(free)
            if a:p[0],p[1]=p[1],p[0]
            if b:p[2],p[3]=p[3],p[2]
            if swap:p[:4]=p[2:4]+p[:2]
            permutations.append(p)
    all_trees=trees(tuple(range(7)));representatives=[];seen=set()
    for tree in all_trees:
        key=min(serial(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(permutations)==48 and len(all_trees)==10395 and len(representatives)==350
    raw=decode_terminal_blob(envelope,1,'zlib-base64-readonce-seven-two-pairs',64000)
    assert raw[:4]==b'RO73' and struct.unpack_from('<H',raw,4)[0]==350
    position=6;covers=[]
    for _ in representatives:
        assert position<len(raw);count=raw[position];position+=1
        assert 1<=count<=255 and position+2*count<=len(raw)
        ids=list(struct.unpack_from('<'+'H'*count,raw,position));position+=2*count
        assert all(0<=i<len(planes) for i in ids);covers.append(ids)
    assert position==len(raw)
    assignments=16*3**6;full=(1<<assignments)-1;repeat=full//65535
    signs=[sum(((v>>i)&1)<<v for v in range(16))*repeat for i in range(4)]
    op_masks=[]
    for i in range(6):
        stride=16*3**i
        a=full//((1<<(3*stride))-1)*((1<<stride)-1)
        op_masks.append((a,a<<stride,a<<(2*stride)))
    for tree,ids in zip(representatives,covers):
        operations=[]
        def compile_tree(t):
            if isinstance(t,int):return t
            a,b=compile_tree(t[0]),compile_tree(t[1])
            operations.append((a,b));return 6+len(operations)
        root=compile_tree(tree);assert root==12 and len(operations)==6
        values={}
        for j in {j for i in ids for j in planes[i]}:
            a=legal[j]
            v=[signs[i] if not(a>>i&1) else full^signs[i] for i in range(4)]
            v.extend(full if a>>i&1 else 0 for i in (4,5,6))
            for (left,right),(aa,oo,xx) in zip(operations,op_masks):
                x,y=v[left],v[right]
                v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
            values[j]=v[root]
        covered=0
        for i in ids:
            zero=one=full
            for j in planes[i]:one &=values[j];zero &=full^values[j]
            covered |=zero|one
        assert covered==full
    return {'passed':True,'shape_orbits':350,'assignments_each':assignments,
            'assignments_verified':350*assignments,
            'plane_witnesses':sum(map(len,covers)),'raw_bytes':len(raw),
            'elapsed_seconds':time.monotonic()-started}


def terminal_certificate_audit(path):
    certificate=json.loads(Path(path).read_text())
    assert set(certificate)=={'version','ramsey84','readonce7'} and certificate['version']==1
    return {'ramsey84':ramsey84_rup_audit(certificate['ramsey84']),
            'readonce7':readonce_seven_audit(certificate['readonce7'])}




def ramsey6_roots():
    """Regenerate every normalized coloring case; no stored CNF is trusted."""
    flats6=affine_flat_indices(list(range(64)),6,3)
    assert len(flats6)==11160
    roots=[]
    for flat in flats6:
        for omit in range(8):
            clause=tuple(v+1 for j,v in enumerate(flat) if j!=omit)
            roots.extend((clause,tuple(-v for v in clause)))
    roots.extend((-(a+1),) if a<6 else (a+1,) for a in range(8))
    roots.extend(((-9,),(-17,),(-33,)))
    assert len(roots)==178571
    yield 'maximum-six',64,roots
    del roots

    flats5=affine_flat_indices(list(range(32)),5,3)
    assert len(flats5)==620
    assignment_count=1<<16;full=(1<<assignment_count)-1
    masks=[]
    for j in range(16):
        stride=1<<j
        masks.append((((1<<stride)-1)<<stride)
                     *(full//((1<<(2*stride))-1)))
    case_zero_sets=({0},{0,1},{0,1,2},{0,1,2,4})
    expected_patterns=(3136,5760,6232,5296)
    expected_groups=(168,48,24,24)
    expected_orbits=(6,18,28,30)
    for case,zeros in enumerate(case_zero_sets,1):
        core=[0]*7+[1]+[int(j not in zeros) for j in range(8)]
        allowed=full
        for flat in flats5:
            ones=zeroes=full
            for j in flat:
                if j<16:
                    ones &= full if core[j] else 0
                    zeroes &= full if not core[j] else 0
                else:
                    ones &= masks[j-16]
                    zeroes &= full^masks[j-16]
            allowed &=full^(ones|zeroes)
        patterns=set()
        while allowed:
            bit=(allowed&-allowed).bit_length()-1
            allowed &=allowed-1;patterns.add(bit)
        assert len(patterns)==expected_patterns[case-1]
        group=[]
        for columns in itertools.permutations(range(1,8),3):
            linear=[0]*8
            for v in range(8):
                for j in range(3):
                    if v>>j&1:linear[v]^=columns[j]
            if len(set(linear))!=8:continue
            shift=7^linear[7]
            for shear in range(8):
                permutation=[(linear[v&7]^shift^(shear if v&8 else 0))
                             |(v&8) for v in range(16)]
                if all(core[permutation[v]]==core[v] for v in range(16)):
                    group.append(permutation)
        assert len(group)==expected_groups[case-1]
        seen=set();representatives=[];orbit_sizes=[]
        for pattern in sorted(patterns):
            if pattern in seen:continue
            orbit=set()
            for permutation in group:
                image=transform_coloring(pattern,permutation,0)
                for shift in range(16):
                    orbit.add(sum(((image>>j)&1)<<(j^shift)
                                  for j in range(16)))
            assert orbit<=patterns
            seen |=orbit;orbit_sizes.append(len(orbit))
            representatives.append(min(v for v in orbit if not(v&1)))
        assert seen==patterns and len(representatives)==expected_orbits[case-1]
        assert sum(orbit_sizes)==len(patterns)
        for index,pattern in enumerate(representatives):
            fixed=core+[(pattern>>j)&1 for j in range(16)]
            roots={(-1,)}
            for flat in flats6:
                for color in (0,1):
                    if any(j<32 and fixed[j]!=color for j in flat):continue
                    roots.add(tuple((1 if color==0 else -1)*(j-31)
                                    for j in flat if j>=32))
            roots=sorted(roots,key=lambda clause:(len(clause),clause))
            assert roots and all(clause for clause in roots)
            yield 'core-seven-%d-%d'%(case,index),32,roots


def ramsey6_audit(path):
    started=time.monotonic()
    raw=decode_terminal_blob(json.loads(path.read_text()),1,
                             'zlib-base64-rup-ramsey6',100000)
    assert raw[:5]==b'R6UP\x53';position=5
    def integer():
        nonlocal position
        value=0
        for shift in range(0,35,7):
            assert position<len(raw)
            byte=raw[position];position+=1;value|=(byte&127)<<shift
            if byte<128:return value
        raise AssertionError('oversized Ramsey integer')
    cases=[]
    for name,nv,roots in ramsey6_roots():
        count=integer();assert 1<=count<=len(roots)
        selected=[];index=0
        for _ in range(count):
            delta=integer();assert delta>0;index+=delta
            assert index<=len(roots);selected.append(roots[index-1])
        checker=ReverseUnitPropagation(selected,nv)
        lemmas=integer();assert 1<=lemmas<=50000
        for j in range(lemmas):
            size=integer();assert size<=nv and position+size<=len(raw)
            data=raw[position:position+size];position+=size
            assert all(byte<2*nv for byte in data)
            clause=[((byte>>1)+1)*(-1 if byte&1 else 1) for byte in data]
            assert len({abs(v) for v in clause})==len(clause)
            assert checker.proves(clause)
            if clause:
                assert j<lemmas-1;checker.add(clause)
            else:assert j==lemmas-1
        assert not clause
        cases.append({'case':name,'root_clauses':len(roots),
                      'selected_roots':count,'derived_clauses':lemmas})
    assert len(cases)==83 and position==len(raw)
    return {'passed':True,'cases':len(cases),'raw_bytes':len(raw),
            'selected_roots':sum(r['selected_roots'] for r in cases),
            'derived_clauses':sum(r['derived_clauses'] for r in cases),
            'elapsed_seconds':time.monotonic()-started}


def run():
    f = load_family()
    certificate_path = Path('/snapshot/files/finite_certificates.txt')
    if not certificate_path.exists():
        certificate_path = Path(__file__).with_name('finite_certificates.txt')
    terminal_results = terminal_certificate_audit(certificate_path)
    ramsey_path=certificate_path.with_name("ramsey6_certificate.txt")
    six_variable_results=ramsey6_audit(ramsey_path)
    results = {"terminal_certificates": terminal_results,
               "six_variable_ramsey": six_variable_results,
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
