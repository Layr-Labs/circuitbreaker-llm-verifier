"""Exact finite witnesses and the unchanged total uniform family evaluator."""
import base64
import struct
import time
import lzma
import itertools
import json
from pathlib import Path
import sys

from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def parameters(n):
    if n == 0:
        return None
    target = 32 * n**3
    lo, hi = 0, 1
    while hi**5 <= target:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**5 > target:
            hi = mid
        else:
            lo = mid
    k = hi
    while not is_prime(k):
        k += 1
    if k >= n:
        return None
    m = n - k
    r = (m + k - 1) // k
    t = isqrt(r)
    if t * t != r:
        t += 1
    if t >= r or k * t < r * (r + 1):
        return None
    return k, m, r, t


def poly_mod(a, b):
    if b == 0:
        raise ZeroDivisionError("zero polynomial divisor")
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_mod(a, b)
    return a


def field_mul(a, b, modulus):
    degree = modulus.bit_length() - 1
    high_bit = 1 << degree
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & high_bit:
            a ^= modulus
    return out


def field_pow(a, exponent, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = field_mul(result, a, modulus)
        exponent >>= 1
        if exponent:
            a = field_mul(a, a, modulus)
    return result


def irreducible_prime_degree(modulus, k):
    """Exact for prime k, as established in proof.md section 2."""
    if modulus.bit_length() != k + 1 or not is_prime(k):
        return False
    x = 2
    for _ in range(k):
        x = field_mul(x, x, modulus)
    return x == 2 and poly_gcd(modulus, 0b110) == 1


def evaluate(bits):
    """Return F_N on the sequence of N bits, in written input order."""
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("the input alphabet is {0,1}")
    params = parameters(len(bits))
    if params is None:
        return 0
    k, m, r, t = params
    modulus = (1 << k) | sum(bit << j for j, bit in enumerate(bits[:k]))
    if not irreducible_prime_degree(modulus, k):
        return 0
    payload = list(bits[k:]) + [0] * (r * k - m)
    symmetric = [1] + [0] * t
    for index in range(r):
        block = payload[index * k:(index + 1) * k]
        value = sum(bit << j for j, bit in enumerate(block))
        # i=index+1, e_i=r+1-i=r-index, w_i=e_i+1.
        power = (1 << (r - index + 1)) - 1
        value = field_pow(value, power, modulus)
        for j in range(min(index + 1, t), 0, -1):
            symmetric[j] ^= field_mul(value, symmetric[j - 1], modulus)
    return symmetric[t] & 1


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
    compressed=base64.b85decode(envelope['data'])
    assert base64.b85encode(compressed).decode()==envelope['data']
    decoder=lzma.LZMADecompressor(memlimit=16*1024*1024)
    raw=decoder.decompress(compressed,max_length=maximum)
    assert decoder.eof and not decoder.unused_data
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
    raw=decode_terminal_blob(envelope,2,'xz-base85-rup-ramsey84-nohints',200000)
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
    trees=binary_trees;serial=tree_string
    permutations=[]
    for swap,a,b in itertools.product((0,1),repeat=3):
        for free in itertools.permutations((4,5,6)):
            p=[0,1,2,3]+list(free)
            if a:p[0],p[1]=p[1],p[0]
            if b:p[2],p[3]=p[3],p[2]
            if swap:p[:4]=p[2:4]+p[:2]
            permutations.append(tuple(p))
    all_trees=trees(tuple(range(7)));representatives=[];seen=set()
    for tree in all_trees:
        key=min(serial(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(permutations)==48 and len(all_trees)==10395 and len(representatives)==350
    clear_tree_caches();del all_trees
    raw=decode_terminal_blob(envelope,1,'xz-base85-readonce-seven-two-pairs',64000)
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
        op_masks.append((a|(a<<stride),(a<<stride)|(a<<(2*stride))))
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
            for (left,right),(ma,mb) in zip(operations,op_masks):
                x,y=v[left],v[right]
                v.append(((x&y)&ma)^((x^y)&mb))
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


def ramsey6_audit(envelope):
    started=time.monotonic()
    raw=decode_terminal_blob(envelope,1,
                             'xz-base85-rup-ramsey6',100000)
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

from functools import lru_cache

@lru_cache(maxsize=None)
def binary_trees(labels):
    if len(labels)==1:return (labels[0],)
    output=[];first=labels[0];rest=labels[1:]
    for size in range(len(rest)):
        for selected in itertools.combinations(rest,size):
            left=(first,)+selected
            right=tuple(v for v in labels if v not in left)
            for a in binary_trees(left):
                for b in binary_trees(right):output.append((a,b))
    return tuple(output)

@lru_cache(maxsize=100000)
def tree_string(tree,permutation):
    if isinstance(tree,int):return str(permutation[tree])
    return '('+','.join(sorted((tree_string(tree[0],permutation),
                                tree_string(tree[1],permutation))))+')'

def clear_tree_caches():
    binary_trees.cache_clear();tree_string.cache_clear()

class FiniteReader:
    def __init__(self,envelope,encoding,magic,maximum=500000):
        self.raw=decode_terminal_blob(envelope,1,encoding,maximum)
        assert self.raw[:len(magic)]==magic
        self.position=len(magic)
    def number(self,size):
        assert size in (1,2,4) and self.position+size<=len(self.raw)
        value=int.from_bytes(self.raw[self.position:self.position+size],'little')
        self.position+=size;return value
    def common(self,plane_count):
        count=self.number(2);assert count==64
        values=[self.number(2) for _ in range(count)]
        assert len(set(values))==count and all(v<plane_count for v in values)
        return values
    def hints(self,plane_count,size=2):
        count=self.number(size);assert count<=100
        values=[self.number(2) for _ in range(count)]
        assert len(set(values))==count and all(v<plane_count for v in values)
        return values
    def finish(self):assert self.position==len(self.raw)


def check_common_cover(values,flats,common,hints,full):
    """Exhaust all parameter assignments; the fixed small scan trusts no result."""
    left=full;used=set();cache={};position=0;steps=0
    while left:
        bit=(left&-left).bit_length()-1;choice=None
        for candidate in common:
            if candidate in used:continue
            flat=flats[candidate];color=(values[flat[0]]>>bit)&1
            if all(((values[j]>>bit)&1)==color for j in flat[1:]):
                choice=candidate;break
        if choice is None:
            assert position<len(hints)
            choice=hints[position];position+=1
        assert choice not in used
        if choice not in cache:
            zero=one=full
            for j in flats[choice]:
                one &=values[j];zero &=full^values[j]
            cache[choice]=zero|one
        mask=cache[choice]&left
        assert mask>>bit&1
        left &=full^mask;used.add(choice);steps+=1
        assert steps<=len(common)+len(hints)
    assert position==len(hints)
    return steps


def check_union_cover(values,flats,common,hints,full):
    left=full;steps=0
    for choice in common:
        zero=one=full
        for j in flats[choice]:one &=values[j];zero &=full^values[j]
        left &=full^(zero|one);steps+=1
        if not left:break
    for choice in hints:
        zero=one=full
        for j in flats[choice]:one &=values[j];zero &=full^values[j]
        mask=(zero|one)&left;assert mask
        left &=full^mask;steps+=1
    assert not left
    return steps


def repeated_component_shapes(kind):
    permutations=[]
    if kind=='cross':
        for swap in (0,1):
            for free in itertools.permutations((3,4,5)):
                p=[0,1,2,*free,6,6]
                if swap:p[1],p[2]=p[2],p[1]
                permutations.append(tuple(p))
        forbidden=(1,2)
    else:
        assert kind=='mixed';forbidden=(1,2,3)
        for a,b in itertools.product((0,1),repeat=2):
            p=list(range(7))+[6]
            if a:p[2],p[3]=p[3],p[2]
            if b:p[4],p[5]=p[5],p[4]
            permutations.append(tuple(p))
    seen=set();shapes=[]
    for size in range(1,6):
        for subset in itertools.combinations(range(6),size):
            if size==1 and subset[0] in forbidden:continue
            images=[tuple(sorted(p[v] for v in subset)) for p in permutations]
            if subset!=min(images):continue
            stabilizer=[p for p,image in zip(permutations,images) if image==subset]
            outside=tuple(sorted(tuple(i for i in range(6) if i not in subset)+(6,7)))
            for inside in binary_trees(subset):
                for outer in binary_trees(outside):
                    key=min(tree_string(inside,p)+'|'+tree_string(outer,p)
                            for p in stabilizer)
                    if key in seen:continue
                    seen.add(key);shapes.append((inside,outer,size))
    if kind=='mixed':
        for outer in binary_trees(tuple(range(7))):
            key='extra|'+min(tree_string(outer,p) for p in permutations)
            if key not in seen:seen.add(key);shapes.append((-1,outer,0))
    assert len(shapes)==(4208 if kind=='cross' else 12468)
    clear_tree_caches()
    return shapes


def near_readonce_audit(envelope,kind):
    started=time.monotonic()
    legal=[a for a in range(128) if not(a&1 and a&4) and not(a&2 and a&8)]
    flats=affine_flat_indices(legal,7,3);assert len(flats)==5049
    if kind=='mixed':
        permutation=(1,3,2,4,0,5,6)
        legal=[sum(((a>>j)&1)<<permutation[j] for j in range(7)) for a in legal]
    reader=FiniteReader(envelope,'xz-base85-near-'+kind,
                        b'NU6C' if kind=='cross' else b'NU6M')
    common=reader.common(len(flats));shapes=repeated_component_shapes(kind)
    assert reader.number(4)==len(shapes)
    prepared={};base_arrays={};bottom_cache={};previous=None;cached=None
    total_assignments=total_hints=total_steps=0
    for inside,outside,size in shapes:
        signbits=3 if kind=='cross' else (5 if size==0 else 4)
        if signbits not in prepared:
            signs=1<<signbits;count=signs*3**6;full=(1<<count)-1
            sign_masks=[sum(((v>>i)&1)<<v for v in range(signs))
                        *(full//((1<<signs)-1)) for i in range(signbits)]
            operation_masks=[]
            for i in range(6):
                stride=signs*3**i
                mask=full//((1<<(3*stride))-1)*((1<<stride)-1)
                operation_masks.append((mask|(mask<<stride),
                                        (mask<<stride)|(mask<<(2*stride))))
            prepared[signbits]=(count,full,sign_masks,operation_masks)
        count,full,sm,om=prepared[signbits];total_assignments+=count
        operations=[]
        def compile_tree(tree):
            if isinstance(tree,int):return tree
            a=compile_tree(tree[0]);b=compile_tree(tree[1])
            operations.append((a,b));return 7+len(operations)
        inner=compile_tree(inside) if size else -1
        cut=len(operations);root=compile_tree(outside)
        assert len(operations)==6
        if signbits not in base_arrays:
            rows=[]
            for point in legal:
                if kind=='cross':
                    values=[full if not(point&1) and not(point&2) else 0]
                    values +=[full if point>>i&1 else 0 for i in (2,3,4,5,6)]
                    values[1]^=sm[0];values[2]^=sm[1]
                    values +=[0,0]
                else:
                    g=full if point&1 and not(point&2) else 0
                    if not size:g=(sm[4]^(full if point&2 else 0)) if point&1 else 0
                    values=[g]+[full if point>>i&1 else 0 for i in (2,3,4,5,6)]
                    for j in range(3):values[j+1]^=sm[j]
                    values +=[0,0]
                    if not size:values[6]=(full if point&1 else 0)^sm[3]
                rows.append(values)
            base_arrays[signbits]=[list(column) for column in zip(*rows)]
        def node_value(values,a,b,index):
            key=(signbits,index,min(a,b),max(a,b))
            limit=7 if not size else 6
            simple=a<limit and b<limit
            if simple and key in bottom_cache:return bottom_cache[key]
            ma,mb=om[index]
            output=[((x&y)&ma)^((x^y)&mb) for x,y in zip(values[a],values[b])]
            if simple:bottom_cache[key]=output
            return output
        if previous!=(inside,signbits):
            cached=base_arrays[signbits][:]
            for index,(a,b) in enumerate(operations[:cut]):
                cached.append(node_value(cached,a,b,index))
            previous=(inside,signbits)
        values=cached[:]
        if size:
            values[6]=values[inner]
            values[7]=[v^sm[2 if kind=='cross' else 3] for v in values[inner]]
        for index,(a,b) in enumerate(operations[cut:],cut):
            values.append(node_value(values,a,b,index))
        hints=reader.hints(len(flats));total_hints+=len(hints)
        total_steps+=check_union_cover(values[root],flats,common,hints,full)
    reader.finish()
    return {'passed':True,'shapes':len(shapes),'assignments':total_assignments,
            'hinted_flats':total_hints,'cover_steps':total_steps,
            'elapsed_seconds':time.monotonic()-started}


def readonce_eight_audit(envelope):
    started=time.monotonic()
    legal=[a for a in range(256) if not(a&1 and a&4) and not(a&2 and a&8)
           and not(a&16 and a&32) and not(a&64 and a&128)]
    flats=affine_flat_indices(legal,8,3);assert len(flats)==2106
    reader=FiniteReader(envelope,'xz-base85-readonce-eight-auto',b'RA84')
    common=reader.common(len(flats));permutations=[]
    for a,b,e,swap in itertools.product((0,1),repeat=4):
        p=list(range(7))
        if a:p[1],p[2]=p[2],p[1]
        if b:p[3],p[4]=p[4],p[3]
        if e:p[5],p[6]=p[6],p[5]
        if swap:p[3:7]=p[5:7]+p[3:5]
        permutations.append(tuple(p))
    seen=set();representatives=[]
    for tree in binary_trees(tuple(range(7))):
        key=min(tree_string(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(representatives)==880 and reader.number(2)==880
    clear_tree_caches();count=64*3**6;full=(1<<count)-1
    signs=[sum(((v>>i)&1)<<v for v in range(64))*(full//((1<<64)-1))
           for i in range(6)]
    masks=[]
    for i in range(6):
        stride=64*3**i;mask=full//((1<<(3*stride))-1)*((1<<stride)-1)
        masks.append((mask|(mask<<stride),(mask<<stride)|(mask<<(2*stride))))
    rows=[]
    for point in legal:
        row=[full if not(point&1) and not(point&2) else 0]
        row +=[(full if point>>(i+2)&1 else 0)^signs[i] for i in range(6)]
        rows.append(row)
    base=[list(column) for column in zip(*rows)];bottom={}
    total_hints=total_steps=0
    for tree in representatives:
        operations=[]
        def compile_tree(tree):
            if isinstance(tree,int):return tree
            a=compile_tree(tree[0]);b=compile_tree(tree[1])
            operations.append((a,b));return 6+len(operations)
        root=compile_tree(tree);assert len(operations)==6
        values=base[:]
        for index,(a,b) in enumerate(operations):
            key=(index,min(a,b),max(a,b));simple=a<7 and b<7
            if simple and key in bottom:result=bottom[key]
            else:
                ma,mb=masks[index]
                result=[((x&y)&ma)^((x^y)&mb) for x,y in zip(values[a],values[b])]
                if simple:bottom[key]=result
            values.append(result)
        hints=reader.hints(len(flats),1);total_hints+=len(hints)
        total_steps+=check_common_cover(values[root],flats,common,hints,full)
    reader.finish()
    return {'passed':True,'shapes':880,'assignments':880*count,
            'hinted_flats':total_hints,'cover_steps':total_steps,
            'elapsed_seconds':time.monotonic()-started}


def graft_new_leaf(tree,label):
    yield ('n',tree,label)
    if isinstance(tree,int):return
    if tree[0]=='o':
        _,table,a,b=tree
        for value in graft_new_leaf(a,label):yield ('o',table,value,b)
        for value in graft_new_leaf(b,label):yield ('o',table,a,value)
    else:
        _,a,b=tree
        for value in graft_new_leaf(a,label):yield ('n',value,b)
        for value in graft_new_leaf(b,label):yield ('n',a,value)


def tagged_tree_string(tree,permutation):
    if isinstance(tree,int):return str(permutation[tree])
    label='o'+str(tree[1]) if tree[0]=='o' else 'n'
    return label+'('+','.join(sorted((tagged_tree_string(tree[-2],permutation),
                                      tagged_tree_string(tree[-1],permutation))))+')'


def extension_shapes():
    a=('o',14,0,2);b=('o',14,1,3);c=('o',14,4,5)
    originals=(('o',6,a,('o',6,b,c)),('o',6,('o',6,a,b),c))
    group=[]
    for swap,free,partners in itertools.product((0,1),repeat=3):
        p=list(range(9))
        if swap:p[0],p[1]=p[1],p[0];p[2],p[3]=p[3],p[2]
        if free:p[4],p[5]=p[5],p[4]
        if partners:p[7],p[8]=p[8],p[7]
        group.append(tuple(p))
    seen=set();output=[]
    for tree in originals:
        for first in graft_new_leaf(tree,6):
            for second in graft_new_leaf(first,7):
                for third in graft_new_leaf(second,8):
                    key=min(tagged_tree_string(third,p) for p in group)
                    if key not in seen:seen.add(key);output.append(third)
    assert len(output)==1450
    return output


def extension_values(legal):
    prepared={};bases={};bottom={}
    for tree in extension_shapes():
        operations=[]
        def compile_tree(tree):
            if isinstance(tree,int):return tree,tree<6
            a,old_a=compile_tree(tree[-2]);b,old_b=compile_tree(tree[-1])
            if tree[0]=='o':table=tree[1];options=None
            else:
                assert not old_b
                table=None;options=(4,6,14) if old_a else (2,4,6,8,14)
            operations.append((a,b,table,options))
            return 8+len(operations),old_a or old_b
        root,_=compile_tree(tree);assert len(operations)==8
        choices=tuple(options for _,_,_,options in operations if options)
        assert len(choices)==3
        if choices not in prepared:
            count=4
            for options in choices:count*=len(options)
            full=(1<<count)-1;stride=4;coefficients=[]
            for options in choices:
                k=len(options);mask=full//((1<<(k*stride))-1)*((1<<stride)-1)
                blocks=[mask<<(j*stride) for j in range(k)]
                m1=m2=m3=0
                for table,block in zip(options,blocks):
                    if table&4:m1|=block
                    if table&2:m2|=block
                    if ((table>>1)^(table>>2)^(table>>3))&1:m3|=block
                coefficients.append((m1,m2,m3));stride*=k
            prepared[choices]=(count,full,coefficients)
        count,full,coefficients=prepared[choices]
        if count not in bases:
            repeat=full//15;rows=[]
            for point in legal:
                row=[full if point>>j&1 else 0 for j in range(6)]
                alpha=0 if point&64 else 1;beta=0 if point&256 else 1
                row.append((1<<(2*alpha+beta))*repeat)
                row.extend((full if point&128 else 0,full if point&512 else 0))
                rows.append(row)
            bases[count]=[list(column) for column in zip(*rows)]
        values=bases[count][:];new_index=0
        for a,b,table,options in operations:
            if table is not None:
                assert table in (6,14)
                key=(count,table,min(a,b),max(a,b));simple=a<6 and b<6
                if simple and key in bottom:result=bottom[key]
                elif table==6:result=[x^y for x,y in zip(values[a],values[b])]
                else:result=[x|y for x,y in zip(values[a],values[b])]
                if simple:bottom[key]=result
            else:
                m1,m2,m3=coefficients[new_index];new_index+=1
                if m1==full:
                    result=[x^(y&m2)^((x&y)&m3) for x,y in zip(values[a],values[b])]
                else:
                    result=[(x&m1)^(y&m2)^((x&y)&m3) for x,y in zip(values[a],values[b])]
            values.append(result)
        yield count,full,values[root]


def exception_extension_audit(envelope):
    started=time.monotonic()
    legal=[a for a in range(1024) if not(a&1 and a&2)
           and not(a&64 and a&128) and not(a&256 and a&512)]
    index={v:j for j,v in enumerate(legal)};assert len(legal)==432
    reader=FiniteReader(envelope,'xz-base85-exception-extension',b'EA10')
    number=reader.number(2);assert 1<=number<=5000
    flats=[];seen=set()
    for _ in range(number):
        offset=reader.number(2);basis=[reader.number(2) for _ in range(4)]
        assert offset<1024 and all(0<v<1024 for v in basis)
        points=[offset]
        for v in basis:points +=[x^v for x in points]
        assert len(set(points))==16 and all(v in index for v in points)
        key=tuple(sorted(points));assert key not in seen;seen.add(key)
        flats.append(tuple(index[v] for v in points))
    common=reader.common(len(flats));assert reader.number(2)==1450
    total_assignments=total_hints=steps=cases=0
    for count,full,values in extension_values(legal):
        hints=reader.hints(len(flats));total_assignments+=count;total_hints+=len(hints)
        steps+=check_common_cover(values,flats,common,hints,full);cases+=1
    assert cases==1450;reader.finish()
    return {'passed':True,'shapes':cases,'assignments':total_assignments,
            'verified_flats':len(flats),'hinted_flats':total_hints,'cover_steps':steps,
            'elapsed_seconds':time.monotonic()-started}


def one_pair_counterexamples():
    """Classify every six-input read-once exception, including full truth tables."""
    legal=[a for a in range(64) if not(a&1 and a&2)]
    indexed=affine_flat_indices(legal,6,3)
    flats=[tuple(legal[j] for j in flat) for flat in indexed]
    assert len(flats)==1770
    group=[tuple(([1,0] if swap else [0,1])+list(p))
           for swap in (0,1) for p in itertools.permutations(range(2,6))]
    seen=set();trees=[]
    for tree in binary_trees(tuple(range(6))):
        key=min(tree_string(tree,p) for p in group)
        if key not in seen:seen.add(key);trees.append(tree)
    assert len(trees)==41;clear_tree_caches()
    count=4*3**5;full=(1<<count)-1
    signs=[sum(((v>>i)&1)<<v for v in range(4))*(full//15) for i in range(2)]
    masks=[]
    for i in range(5):
        stride=4*3**i;m=full//((1<<(3*stride))-1)*((1<<stride)-1)
        masks.append((m|(m<<stride),(m<<stride)|(m<<(2*stride))))
    canonical=lambda a:int(bool(a&1 or a&4)^bool(a&2 or a&8)^bool(a&16 or a&32))
    orbit=set();truth_full=(1<<64)-1
    for shift in range(16):
        truth=sum(canonical(a^(shift<<2))<<a for a in range(64))
        orbit.update((truth,truth^truth_full))
    originals=(((0,2),((1,3),(4,5))),(((0,2),(1,3)),(4,5)))
    exceptions=[]
    for tree in trees:
        operations=[]
        def compile_tree(tree):
            if isinstance(tree,int):return tree
            a=compile_tree(tree[0]);b=compile_tree(tree[1])
            operations.append((a,b));return 5+len(operations)
        root=compile_tree(tree);outputs=[]
        for point in range(64):
            values=[full if point>>j&1 else 0 for j in range(6)]
            values[0]^=signs[0];values[1]^=signs[1]
            for (a,b),(ma,mb) in zip(operations,masks):
                x,y=values[a],values[b];values.append(((x&y)&ma)^((x^y)&mb))
            outputs.append(values[root])
        left=full
        for flat in flats:
            zero=one=full
            for j in flat:one &=outputs[j];zero &=full^outputs[j]
            left &=full^(one|zero)
        while left:
            bit=(left&-left).bit_length()-1;left &=left-1
            truth=sum(((v>>bit)&1)<<j for j,v in enumerate(outputs))
            assert truth in orbit and tree in originals
            value=bit//4;tables=[]
            for _ in range(5):tables.append((8,14,6)[value%3]);value//=3
            assert value==0
            # The three disjoint cherries are AND/OR; both joins are XOR.
            for (a,b),table in zip(operations,tables):
                assert table in ((8,14) if a<6 and b<6 else (6,))
            exceptions.append({'tree':tree,'signs':bit%4,'tables':tuple(tables)})
    assert len(exceptions)==16
    assert all(sum(row['tree']==tree for row in exceptions)==8 for tree in originals)
    return exceptions


def cleanup_legal(partner):
    base=[a for a in range(128) if not(a&1 and a&2) and not(a&4 and a&8)]
    remaining=[j for j in range(2,6) if j!=partner]
    permutation=(0,1,partner,6,*remaining)
    return [sum(((a>>j)&1)<<permutation[j] for j in range(7)) for a in base]


def counter_tag(tree,operations):
    if isinstance(tree,int):return tree
    a=counter_tag(tree[0],operations);b=counter_tag(tree[1],operations)
    index=len(operations);operations.append((a,b));return ('o',index,a,b)


def counter_masks(cases,full):
    repeat=full//65535
    sign0=sum((row['signs']&1)<<i for i,row in enumerate(cases))*repeat
    sign1=sum(((row['signs']>>1)&1)<<i for i,row in enumerate(cases))*repeat
    # Low three bits select one of eight exceptions, bit3 is the extra sign.
    sign0 |=sign0<<8;sign1 |=sign1<<8
    flip=65280*repeat
    old=[]
    for j in range(5):
        tables=[row['tables'][j] for row in cases]
        if set(tables)=={6}:old.append(None)
        else:
            assert set(tables)<={8,14}
            value=sum(int(t==14)<<i for i,t in enumerate(tables))
            old.append((value|(value<<8))*repeat)
    return sign0,sign1,flip,old


def zero_table_masks(options,stride,full):
    k=len(options);mask=full//((1<<(k*stride))-1)*((1<<stride)-1)
    m1=m2=m3=0
    for i,table in enumerate(options):
        block=mask<<(i*stride)
        if table&4:m1|=block
        if table&2:m2|=block
        if ((table>>1)^(table>>2)^(table>>3))&1:m3|=block
    return m1,m2,m3


def reverse_A_values(exceptions):
    originals=list(dict.fromkeys(row['tree'] for row in exceptions))
    for partner in range(2,6):
        legal=cleanup_legal(partner)
        for original in originals:
            cases=[row for row in exceptions if row['tree']==original]
            tree=counter_tag(original,[]);seen=set();shapes=[]
            for first in graft_new_leaf(tree,6):
                for second in graft_new_leaf(first,6):
                    key=tagged_tree_string(second,tuple(range(7)))
                    if key not in seen:seen.add(key);shapes.append(second)
            assert len(shapes)==77;prepared={}
            for tree in shapes:
                operations=[]
                def compile_tree(tree):
                    if isinstance(tree,int):return tree,tree<6
                    a,old_a=compile_tree(tree[-2]);b,old_b=compile_tree(tree[-1])
                    if tree[0]=='o':slot=tree[1];options=None
                    else:
                        assert not old_b;slot=-1
                        options=(4,6,12,14) if old_a else (2,4,6,8,14)
                    operations.append((a,b,slot,options))
                    return 6+len(operations),old_a or old_b
                root,_=compile_tree(tree);assert len(operations)==7
                choices=tuple(options for _,_,_,options in operations if options)
                assert len(choices)==2
                if choices not in prepared:
                    count=16
                    for options in choices:count*=len(options)
                    full=(1<<count)-1;stride=16;new=[]
                    for options in choices:
                        new.append(zero_table_masks(options,stride,full));stride*=len(options)
                    prepared[choices]=(count,full,new,counter_masks(cases,full))
                count,full,new,(s0,s1,flip,old)=prepared[choices]
                rows=[]
                for point in legal:
                    row=[full if point>>j&1 else 0 for j in range(7)]
                    row[0]^=s0;row[1]^=s1;row[partner]^=flip;rows.append(row)
                values=[list(column) for column in zip(*rows)];new_index=0
                for a,b,slot,options in operations:
                    if slot>=0:
                        mask=old[slot]
                        if mask is None:result=[x^y for x,y in zip(values[a],values[b])]
                        else:result=[(x&y)^((x^y)&mask) for x,y in zip(values[a],values[b])]
                    else:
                        m1,m2,m3=new[new_index];new_index+=1
                        if m1==full:result=[x^(y&m2)^((x&y)&m3) for x,y in zip(values[a],values[b])]
                        else:result=[(x&m1)^(y&m2)^((x&y)&m3) for x,y in zip(values[a],values[b])]
                    values.append(result)
                yield partner,count,full,values[root]


def reverse_C_values(exceptions):
    originals=list(dict.fromkeys(row['tree'] for row in exceptions))
    for partner in range(2,6):
        legal=cleanup_legal(partner);others=[j for j in range(6) if j!=partner]
        for original in originals:
            cases=[row for row in exceptions if row['tree']==original]
            tree=counter_tag(original,[])
            for coupled in (False,True):
                choices=others if coupled else [-1]+others
                count=16*len(choices)*(1 if coupled else 4)
                full=(1<<count)-1;s0,s1,flip,old=counter_masks(cases,full)
                choice_masks=[]
                stride=16;k=len(choices)
                mask=full//((1<<(k*stride))-1)*((1<<stride)-1)
                for i in range(k):choice_masks.append(mask<<(i*stride))
                rows=[]
                for point in legal:
                    row=[full if point>>j&1 else 0 for j in range(7)]
                    selected=0
                    for choice,mask in zip(choices,choice_masks):
                        if choice>=0 and point>>choice&1:selected |=mask
                        if coupled and point&64:row[choice]^=mask
                    row[partner]=(full if point>>partner&1 else 0)^selected^flip
                    if coupled and point&64:row[partner]^=full
                    row[0]^=s0;row[1]^=s1;rows.append(row)
                base=[list(column) for column in zip(*rows)]
                new=zero_table_masks((4,6,12,14),16*k,full) if not coupled else None
                shapes=(tree,) if coupled else tuple(graft_new_leaf(tree,6))
                assert len(shapes)==(1 if coupled else 11)
                for shape in shapes:
                    operations=[]
                    def compile_tree(tree):
                        if isinstance(tree,int):return tree
                        a=compile_tree(tree[-2]);b=compile_tree(tree[-1])
                        operations.append((a,b,tree[1] if tree[0]=='o' else -1))
                        return 6+len(operations)
                    root=compile_tree(shape);values=base[:]
                    for a,b,slot in operations:
                        if slot>=0:
                            mask=old[slot]
                            if mask is None:result=[x^y for x,y in zip(values[a],values[b])]
                            else:result=[(x&y)^((x^y)&mask) for x,y in zip(values[a],values[b])]
                        else:
                            m1,m2,m3=new;assert m1==full
                            result=[x^(y&m2)^((x&y)&m3) for x,y in zip(values[a],values[b])]
                        values.append(result)
                    yield partner,count,full,values[root]


def reverse_cleanup_audit(envelope,mode,exceptions):
    started=time.monotonic();assert mode in ('A','C')
    legal=[a for a in range(128) if not(a&1 and a&2) and not(a&4 and a&8)]
    flats=affine_flat_indices(legal,7,3);assert len(flats)==5049
    reader=FiniteReader(envelope,'xz-base85-reverse-'+mode,b'RC7'+mode.encode())
    common=reader.common(len(flats));expected=616 if mode=='A' else 96
    assert reader.number(2)==expected
    generator=reverse_A_values if mode=='A' else reverse_C_values
    count_cases=total_assignments=total_hints=steps=0
    for partner,count,full,values in generator(exceptions):
        hints=reader.hints(len(flats));total_hints+=len(hints)
        total_assignments+=count;steps+=check_common_cover(values,flats,common,hints,full)
        count_cases+=1
    assert count_cases==expected;reader.finish()
    return {'passed':True,'cases':count_cases,'assignments':total_assignments,
            'hinted_flats':total_hints,'cover_steps':steps,
            'elapsed_seconds':time.monotonic()-started}


def packet(name):
    path=Path('/snapshot/files')/name
    if not path.exists():path=Path(__file__).with_name(name)
    return json.loads(path.read_text())


def certificate_parameters(n):
    selected=parameters(n)
    if selected is None:return None
    k,m,r,t=selected
    h=r*(r+3)//2;capacity=sum(max(k-j,0) for j in range(2,t+1))
    total=sum(max(k-j,0) for j in range(1,r+1))
    dimension=h+max(capacity,total//(t+1))+1
    degree=t*(t+3)//2;coarse=(k*r)//t+k*t+r*(r+1)+1
    return k,m,r,t,dimension,degree,coarse


def run():
    first=packet('finite_certificates.txt');second=packet('refinement_certificates.txt')
    assert set(first)=={'version','ramsey84','readonce7','readonce8'} and first['version']==1
    assert set(second)=={'version','ramsey6','near_cross','near_mixed','extension','reverse_A','reverse_C'} and second['version']==1
    checks={}
    checks['ramsey84']=ramsey84_rup_audit(first['ramsey84'])
    checks['readonce7']=readonce_seven_audit(first['readonce7'])
    checks['ramsey6']=ramsey6_audit(second['ramsey6'])
    checks['near_cross']=near_readonce_audit(second['near_cross'],'cross')
    checks['near_mixed']=near_readonce_audit(second['near_mixed'],'mixed')
    checks['readonce8']=readonce_eight_audit(first['readonce8'])
    checks['extension']=exception_extension_audit(second['extension'])
    exceptions=one_pair_counterexamples();checks['classified_exceptions']=len(exceptions)
    checks['reverse_A']=reverse_cleanup_audit(second['reverse_A'],'A',exceptions)
    checks['reverse_C']=reverse_cleanup_audit(second['reverse_C'],'C',exceptions)
    for n in (0,1,2,8,32,128,256):assert evaluate([0]*n)==0
    checks['total_evaluator_zero_inputs']=7
    examples=[]
    for n in (128,192,256,1024,16384,262144,4194304,67108864):
        data=certificate_parameters(n)
        if data is None:continue
        k,m,r,t,d,b,coarse=data
        assert k>r+1 and 4<=b<=d
        assert d<=coarse and 3*k+8*d-4*b-3<4*k+8*coarse-8
        examples.append({'N':n,'k':k,'r':r,'t':t,'D':d,'B':b,
                         'loss':3*k+8*d-4*b-3,'old_loss':4*k+8*coarse-8})
    checks['certificate_parameter_examples']=examples
    isolated=0
    for k in (2,3,5,7,11):
        for low in range(1<<k):
            q=(1<<k)|low
            if not irreducible_prime_degree(q,k):continue
            assert q.bit_count()%2==1
            for j in range(k):assert not irreducible_prime_degree(q^(1<<j),k)
            isolated+=1
    checks['isolated_valid_prefixes']=isolated
    return {'all_passed':True,'checks':checks,
            'scope':'Exact finite obligations only. Universal reductions, eligibility and all-length asymptotics are supplied in the Markdown proofs.'}


if __name__=='__main__':
    request=json.load(sys.stdin)
    if request!={}:raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(),sort_keys=True))
