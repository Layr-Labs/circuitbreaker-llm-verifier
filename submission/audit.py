"""Exact scoped component checks. Universal claims rest on supplied proofs."""
import json,random,itertools,time
from fractions import Fraction
import importlib.util
from pathlib import Path
import sys
sys.dont_write_bytecode=True
family_path=Path('/snapshot/files/family.py')
if not family_path.exists():family_path=Path(__file__).with_name('family.py')
spec=importlib.util.spec_from_file_location('candidate_family',family_path)
f=importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)

def rank(rows):
    rows=list(rows);out=0
    while rows:
        a=rows.pop()
        if not a:continue
        pivot=a.bit_length()-1;out+=1
        rows=[b^a if b>>pivot&1 else b for b in rows]
    return out

def xor_basis(basis,mask):
    out=0
    for i,b in enumerate(basis):
        if mask>>i&1:out^=b
    return out

def restrict_form(rows,basis):
    return tuple(sum(((sum(((row&v).bit_count()%2)<<i for i,row in enumerate(rows))&u).bit_count()%2)<<j
                     for j,v in enumerate(basis)) for u in basis)

def run():
    rng=random.Random(20261106);result={};root_cases=0
    for n in range(257):
        for degree in (1,2,3,5,8,10,100,1000):
            a=f.root(n,degree);c=f.ceil_root(n,degree)
            assert a**degree<=n<(a+1)**degree
            assert c==a+(a**degree<n);root_cases+=1
    result['exact_integer_roots']=root_cases
    fields=[]
    for k in (3,5,7,11,13,17,29,37):
        q=f.irreducible_prime(k);assert q.bit_length()==k+1
        value=2
        for _ in range(k):value=f.mul(value,value,q)
        assert value==2 and f.gcd(q,0b110)==1
        fields.append([k,q])
    result['prime_field_constructors']=fields
    linear_cases=0
    for k in (3,5,7):
        q=f.irreducible_prime(k)
        for _ in range(100):
            a=rng.randrange(1<<k);b=rng.randrange(1<<k)
            assert f.trace(a^b,q)==f.trace(a,q)^f.trace(b,q)
            linear_cases+=1
    n=80;m=1;z=2;s,ell=f.code_parameters(n,m,z)
    assert s==37 and (1<<s)>1024*(n+1)*(2*m)**6*(1<<(6*z))
    for _ in range(100):
        x=tuple(rng.randrange(2) for _ in range(n));y=tuple(rng.randrange(2) for _ in range(n))
        index=rng.getrandbits(ell)
        a=f.code_bit(f.coefficients(x,s),s,index);b=f.code_bit(f.coefficients(y,s),s,index)
        xy=tuple(a^b for a,b in zip(x,y))
        assert f.code_bit(f.coefficients(xy,s),s,index)==a^b
        linear_cases+=1
    result['trace_and_actual_local_code_linearity']=linear_cases
    result['local_code_parameters']={'n':n,'m':m,'epsilon_exponent':z,'field_degree':s,'index_bits':ell}
    # Polynomial-graph rule on a genuine small prime-degree field. This tests
    # design intersections only, not the full entropy/extraction guarantee.
    q=f.irreducible_prime(7);Q=128;sets=[]
    for a in range(4):
        for b in range(4):sets.append({alpha*Q+(f.mul(b,alpha,q)^a) for alpha in range(74)})
    graph_cases=0
    for i,j in itertools.combinations(range(len(sets)),2):
        assert len(sets[i]&sets[j])<=1;graph_cases+=1
    result['finite_linear_graph_intersections']=graph_cases
    # The common-radical rank-halving assertion, with actual low-rank spaces
    # in a larger ambient space, rather than only full-rank vacuous examples.
    spaces=0;rank_checks=0
    for trial in range(64):
        n=12;c=1+trial%3;generators=[]
        for _ in range(4):
            rows=[0]*n
            for i in range(c):
                for j in range(i+1,n):
                    if rng.randrange(2):rows[i]^=1<<j;rows[j]^=1<<i
            generators.append(tuple(rows))
        forms=[]
        for mask in range(16):
            rows=[0]*n
            for j,g in enumerate(generators):
                if mask>>j&1:rows=[a^b for a,b in zip(rows,g)]
            forms.append(tuple(rows))
        initial=max(rank(a) for a in forms);basis=tuple(1<<i for i in range(n))
        while True:
            restricted=[restrict_form(a,basis) for a in forms]
            ranks=[rank(a) for a in restricted];maximum=max(ranks)
            if maximum==0:break
            kernel=f.kernel(restricted[ranks.index(maximum)],len(basis))
            nextbasis=tuple(xor_basis(basis,v) for v in kernel)
            for a in forms:
                assert 2*rank(restrict_form(a,nextbasis))<=maximum;rank_checks+=1
            basis=nextbasis
        assert n-len(basis)<=2*initial
        assert all(rank(restrict_form(a,basis))==0 for a in forms);spaces+=1
    result['common_radical_spaces']=spaces;result['rank_halving_checks']=rank_checks
    permutations=list(itertools.permutations(range(3)));stats={};support_cases=0
    for trial in range(1000):
        k=rng.randrange(83,501);d=(2*k+4)//5+10;end=k-d;L=k-2*d+4;U=d-4
        support={0,end}
        for i in range(1,end):
            if trial%2==0 and L<=i<=U:continue
            if rng.randrange(2):support.add(i)
        ordered=sorted(support)
        gaps=[(ordered[(i+1)%len(ordered)]-x-1)%k for i,x in enumerate(ordered)]
        if max(gaps)!=d-1 or sum(g==d-1 for g in gaps)!=1:continue
        middle=sorted(i for i in support if L<=i<=U)
        if middle:columns=[0,middle[0]+1,end+2];case='middle'
        else:
            a=max(i for i in support if i<L);b=min(i for i in support if i>U);A=a+1;C=end-b+1
            if A>C:columns=[0,end+1,a+2];case='first_island'
            elif C>A:columns=[b,a+1,end+2];case='second_island'
            else:columns=[0,b+1,a+2];case='equal_islands'
        assert len(set(columns))==3
        for rho in range(k):
            matrix=[[((b-u-rho)%k) in support for b in columns] for u in range(3)]
            match=sum(all(matrix[u][p[u]] for u in range(3)) for p in permutations)
            assert match==(1 if rho==0 else 0)
        support_cases+=1;stats[case]=stats.get(case,0)+1
    assert len(stats)==4
    result['cubic_support_patterns']=support_cases;result['isolation_cases']=stats
    interval_cases=0
    for n in range(1,257):
        K=f.ceil_root(n,100)
        assert (K-1)**100<n<=K**100
        for t in (2,4,8):
            level=0
            while t**level<=n:
                intervals=[f.interval_bounds(n,t,level,i) for i in range(t**level)]
                assert intervals[0][0]==0 and intervals[-1][1]==n
                assert all(a[1]==b[0] for a,b in zip(intervals,intervals[1:]))
                lengths=[b-a for a,b in intervals]
                assert min(lengths)>=1 and max(lengths)-min(lengths)<=1
                if level:
                    for i in range(t**(level-1)):
                        parent=f.interval_bounds(n,t,level-1,i)
                        children=intervals[i*t:(i+1)*t]
                        assert children[0][0]==parent[0] and children[-1][1]==parent[1]
                interval_cases+=len(intervals);level+=1
    result['actual_interval_partition_nodes']=interval_cases
    rounding_cases=0
    for L in range(1,2049):
        B=f.ceil_log2(L+1);U=f.ceil_root(L**3,4);V=f.ceil_root(B,4)
        Z=U*V;E=5120*Z;A=f.ceil_root(L,4)*f.ceil_root(B**3,4)
        ell_log=f.ceil_root(L,2)*f.ceil_root(B,2)
        ratio=(L+B-1)//B;q=f.ceil_root(ratio,2)
        assert (U-1)**4<L**3<=U**4 and (V-1)**4<B<=V**4
        assert (1<<(B-1))<L+1<=(1<<B)
        assert (q-1)**2<ratio<=q*q and E%2==0
        assert A>=1 and ell_log>=1
        rounding_cases+=1
    for denominator in (1,2,3,10,100,1000,10000,100000):
        L=(1<<160)*denominator**8;B=f.ceil_log2(L+1)
        Z=f.ceil_root(L**3,4)*f.ceil_root(B,4)
        assert 3+5120*Z<=(L-1)//denominator
        rounding_cases+=1
    result['quarter_entropy_rounding_cases']=rounding_cases
    extension_cases=0;cycles=[]
    for degree in (1,2,3,4):
        ff,modulus,g=f.singer(3,degree);Q=8**degree
        assert f.extension_irreducible(modulus,3,ff)
        seen=set();v=1
        for _ in range(Q-1):
            seen.add(v);v=f.extension_mul(v,g,3,modulus,ff)
        assert v==1 and len(seen)==Q-1
        cycles.append([degree,Q-1])
        for _ in range(200):
            a,b,c=(rng.randrange(Q) for _ in range(3))
            assert f.extension_mul(a,b^c,3,modulus,ff)==f.extension_mul(a,b,3,modulus,ff)^f.extension_mul(a,c,3,modulus,ff)
            assert f.extension_power(a,Q,3,modulus,ff)==a
            extension_cases+=1
    result['singer_nonzero_cycles']=cycles
    result['extension_field_identity_cases']=extension_cases
    monomial_cases=0
    for degree in range(1,7):
        for h in range(1,8):
            exps=list(f.monomial_exponents(degree,h))
            assert len(exps)==f.comb(h+degree-1,degree)==len(set(exps))
            assert all(len(e)==degree and sum(e)<h for e in exps)
            monomial_cases+=1
    result['graded_monomial_enumeration_cases']=monomial_cases
    optimal_linearity=0
    for degree in (1,2,3):
        for _ in range(100):
            n=12;h=14;seed=rng.randrange(1<<(3*(degree+1)))
            x=tuple(rng.randrange(2) for _ in range(n));y=tuple(rng.randrange(2) for _ in range(n))
            xy=tuple(a^b for a,b in zip(x,y))
            a=f.optimal_extract_with_parameters(x,seed,3,degree,h,4)
            b=f.optimal_extract_with_parameters(y,seed,3,degree,h,4)
            c=f.optimal_extract_with_parameters(xy,seed,3,degree,h,4)
            assert c==tuple(v^w for v,w in zip(a,b))
            optimal_linearity+=1
    result['optimal_extractor_finite_linearity_cases']=optimal_linearity
    assert Fraction(1,8)+Fraction(1,10)==Fraction(9,40)<Fraction(49,200)
    assert Fraction(49,300)>Fraction(1,8)
    assert Fraction(1,64)>Fraction(1,100)
    assert Fraction(1,300)>Fraction(1,310)
    assert Fraction(5120,2*310)>8
    assert Fraction(5120,6)>9
    assert Fraction(9,5120)<Fraction(1,2)<Fraction(2,3)<Fraction(3,4)
    assert Fraction(2,5120)<Fraction(1,100)
    assert Fraction(1,5120)<Fraction(1,100)
    assert Fraction(37*25,10)<93
    assert Fraction(93,5120)<Fraction(1,2)
    assert 6<7<8
    assert Fraction(5,6)<1 and Fraction(22,3)<8 and Fraction(11,3)<Fraction(37,10)
    assert Fraction(99,100)-Fraction(1,5120)>Fraction(3,4)
    result['quarter_exact_margin_checks']=14
    row_geometry=0
    for bits in (400,800,1200):
        initial=1<<bits;large=1<<(bits//4)
        for rows in (2,3,large-1,large):
            assert rows**4<=initial
            rounds=f.ceil_log2(rows)
            last=initial//(4**rounds)
            assert last*8*rows*rows>=initial
            assert last**100>=initial**49
            for j in sorted(set((0,1,rounds//2,max(0,rounds-2),rounds-1))):
                if not 0<=j<rounds:continue
                count=(rows+(1<<j)-1)//(1<<j);length=initial//(4**j)
                assert count**10<=length**3
                row_geometry+=1
    result['precision_row_finite_geometry_cases']=row_geometry
    fallback=0
    for n in range(257):
        assert f.tree_parameters(n) is None
        for pattern in (0,1,2):
            x=tuple(0 if pattern==0 else 1 if pattern==1 else rng.randrange(2) for _ in range(n))
            assert f.evaluate(x)==0;fallback+=1
    result['small_total_fallback_inputs']=fallback
    result['all_passed']=True
    result['scope']='Finite component predicates and total fallback only; no asymptotically nontrivial full tree execution or universal extraction/circuit acceptance.'
    return result

if __name__=='__main__':
    request=json.loads(input())
    if request!={}:raise ValueError('empty input object required')
    print(json.dumps(run(),sort_keys=True))
