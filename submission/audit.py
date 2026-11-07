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
        Z=U*V;E=13*Z;base=f.ceil_root(L,4)*f.ceil_root(B**3,4)
        A=(27*base+24)//25
        assert 25*(A-1)<27*base<=25*A
        ell_log=f.ceil_root(L,2)*f.ceil_root(B,2)
        ratio=(L+B-1)//B;qbase=f.ceil_root(ratio,2);q=(22*qbase+4)//5
        assert (U-1)**4<L**3<=U**4 and (V-1)**4<B<=V**4
        assert (1<<(B-1))<L+1<=(1<<B)
        assert (qbase-1)**2<ratio<=qbase*qbase and 5*(q-1)<22*qbase<=5*q and E==13*Z
        pexp=(843*Z+1999)//2000
        assert 2000*(pexp-1)<843*Z<=2000*pexp
        assert A>=1 and ell_log>=1
        rounding_cases+=1
    for denominator in (1,2,3,10,100,1000,10000,100000):
        L=(1<<160)*denominator**8;B=f.ceil_log2(L+1)
        Z=f.ceil_root(L**3,4)*f.ceil_root(B,4)
        assert 3+13*Z<=(L-1)//denominator
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
    F=Fraction
    eta=F(5,13);aa=F(27,25);rho=F(247,325)
    active=F(2563,6500);hidden=F(5009,13000);final=F(1001,2600)
    gap=F(4753,1000);ordinary=F(116,1000)
    zeta=F(10553,125000);single=F(21106010553,125000000000);near=F(21106135553,125000000000)
    assert F(2000001,1000000)*zeta==single and single+F(1,1000000)==near<F(337699,2000000)
    assert F(3377,10000)>F(337699,1000000)
    assert F(6623,10000)+F(337699,1000000)==F(999999,1000000)<1
    assert F(8442399,100000000)<zeta<F(337699,2000000)
    assert F(8442399,100000000)>F(4221199,50000000)
    assert 5*F(4221199,50000000)>F(843,2000)
    assert 5*F(3377,10000)>F(423,1000)
    assert F(423,1000*13)<eta<final<hidden
    assert F(1,1000)+gap/13<F(3,8)<hidden<active
    assert rho-gap/13>active and active-ordinary/13>hidden
    assert F(1519787,500000)*F(1963,500)<12
    assert 12*aa/13<F(999,1000)
    assert F(2000001,1000000)*F(125,594)<F(421,1000)<F(843,2000)<F(423,1000)
    assert F(3377,10000)+F(10000,15099)==F(150989323,150990000)<1
    assert F(75989323,25000000)<F(1519787,500000)<8
    assert F(10000,15099)+4*F(3377,10000)==F(75989323,37747500)
    assert F(15099,10000)*F(75989323,37747500)==F(75989323,25000000)
    assert F(6623,10000)*eta>1/F(1963,500)
    assert F(463,1000)*F(25001,100000)<ordinary
    assert F(1519787,500000)*aa/F(6623,2000)<F(9914,10000)
    assert F(999999,1000000)+F(1,1000000)==1 and zeta>F(1,2000)
    result['scaled_gap_exact_margin_checks']=20
    separated_ceil=0
    for Z in range(1,1025):
        dexp=(247*Z+24)//25
        assert 25*(dexp-1)<247*Z<=25*dexp<247*Z+25
        assert 13*(5*Z)==5*(13*Z)
        separated_ceil+=1
    result['scaled_gap_dimension_component_cases']=separated_ceil
    # Finite exact witnesses for BOTH exhaustive branches of the new leaf bound.
    powers=[j**1000 for j in range(97)]
    leaf_cases=0;leaf_depth_cases=[0,0]
    for n in range(2,97):
        for t in range(2,6):
            J=0
            while t**(J+1)<=n:J+=1
            for spacing in range(1,6):
                for threshold in range(1,9):
                    for delta in (0,1):
                        kval=powers[threshold]+delta
                        if powers[n]<kval:continue
                        height=0
                        while height+spacing<=J and powers[n//t**(height+spacing)]>=kval:height+=spacing
                        assert powers[n//t**height]>=kval
                        maximum=(n+t**height-1)//t**height
                        if height+spacing<=J:
                            assert powers[n//t**(height+spacing)]<kval
                            leaf_depth_cases[0]+=1
                        else:
                            assert J+1-height<=spacing and maximum<=t**spacing
                            leaf_depth_cases[1]+=1
                        ceiling_root=threshold+(delta>0)
                        assert maximum<=t**spacing*(ceiling_root+1)+1
                        leaf_cases+=1
    assert all(leaf_depth_cases)
    result['scaled_gap_leaf_upper_bound_cases']=leaf_cases
    result['scaled_gap_leaf_depth_branch_counts']=leaf_depth_cases
    output_roots=0
    for r in (0,1,2,3,4,7,8):
        m=f.root(r**999999,1000000)
        assert m**1000000<=r**999999<(m+1)**1000000
        output_roots+=1
    result['larger_smalloutput_integer_roots']=output_roots
    reservoir_plans=0;reservoir_rounds=0
    for bits in (4000000,4200000,4400000):
        K=1<<bits;g=1<<((bits*6623)//10000);w=1<<((bits*337699)//1000000)
        R=1<<((bits*3377)//10000);M=K//4
        assert w<min(R,M)
        for initial in (1,2,3,g-1,g):
            J=f.ceil_log2(initial);count=initial;spent=0;steps=0
            assert (2*initial+J)*w<=K//4
            if initial<=3:
                while count>1:
                    spent+=count*w;count=(count+1)//2;steps+=1
            else:
                # Exact sum of ceil(initial/2^j), j=0,...,J-1.
                spent=(2*initial+J-2-(initial-1).bit_count())*w
                steps=J
                for j in (0,1,J//2,J-1):
                    power=1<<j
                    count=(initial+power-1)//power
                    assert (count+1)//2==(initial+2*power-1)//(2*power)
            assert steps==J and spent<=(2*initial+J)*w
            assert K-spent>=K//2+M
            reservoir_plans+=1;reservoir_rounds+=steps
    result['reservoir_row_count_budget_cases']=reservoir_plans
    result['reservoir_row_count_rounds']=reservoir_rounds
    kernel_fibers=0;kernel_sources=0
    for trial in range(32):
        dimension=8;observations=tuple(rng.randrange(1<<dimension) for _ in range(1+trial%5))
        kernel=f.kernel(observations+(0,)*(dimension-len(observations)),dimension)
        kernel_values={xor_basis(kernel,m) for m in range(1<<len(kernel))}
        fibers={}
        for x in range(1<<dimension):
            z=sum(((x&row).bit_count()%2)<<j for j,row in enumerate(observations))
            fibers.setdefault(z,set()).add(x)
        assert len(kernel)==dimension-rank(observations)
        assert len(fibers)*len(kernel_values)==1<<dimension
        for points in fibers.values():
            representative=min(points)
            assert {x^representative for x in points}==kernel_values
            kernel_fibers+=1
        kernel_sources+=1
    result['reservoir_affine_kernel_sources']=kernel_sources
    result['reservoir_affine_kernel_fibers']=kernel_fibers
    confidence_geometry=0
    for ell_log in list(range(1,601))+[1999999,2000000,2000001,3999999,4000000,4000001,6000000]:
        unit=(ell_log+1999999)//2000000;separation=2000001*unit
        gap=separation-ell_log
        assert gap>=unit and 2000000*(unit-1)<ell_log<=2000000*unit
        if ell_log>=2000001:assert gap<=2*ell_log
        for remaining in range(1,21):
            log1=4000002*remaining*unit;log2=log1+separation
            logD=ell_log+log2;ancestor=4000002*(remaining+1)*unit
            assert log2-log1-ell_log==gap
            assert ancestor-logD==gap
            assert 3*ell_log+log2==(4000002*remaining+2000001)*unit+3*ell_log
            confidence_geometry+=1
    result['quantized_confidence_exponent_geometry_cases']=confidence_geometry
    degree_rounding=0
    for L in range(1,2049):
        for w in sorted(set((1,2,3,L))):
            d=(15099*L+10000*w-1)//(10000*w)
            assert 10000*w*(d-1)<15099*L<=10000*w*d
            degree_rounding+=1
    result['sharpordinary_degree_rounding_cases']=degree_rounding
    fallback=0
    for n in range(257):
        assert f.tree_parameters(n) is None
        for pattern in (0,1,2):
            x=tuple(0 if pattern==0 else 1 if pattern==1 else rng.randrange(2) for _ in range(n))
            assert f.evaluate(x)==0;fallback+=1
    result['small_total_fallback_inputs']=fallback
    # Observe the ACTUAL component before its integer feasibility guard returns.
    # These powers of two distinguish tree ceil(log2 n) from component ceil(log2(n+1)).
    import sys
    component_cases=0;accepted_parameter_cases=0
    target=f.optimal_parameters.__wrapped__.__code__
    for tree_log in (2,3,4,8,16,32,64,128):
        n=1<<tree_log
        for w in sorted(set((1,max(1,tree_log//2),tree_log))):
            seen={}
            def trace_component(frame,event,arg):
                if frame.f_code is target:
                    if event=='return':
                        seen.update({k:frame.f_locals[k] for k in ('L','w','d')})
                    return trace_component
                return None
            f.optimal_parameters.cache_clear()
            sys.settrace(trace_component)
            try:
                actual=f.optimal_parameters(n,1<<w)
            finally:
                sys.settrace(None)
            local_log=f.ceil_log2(n+1)
            assert local_log==tree_log+1 and seen['L']==local_log and seen['w']==w
            expected=(15099*local_log+10000*w-1)//(10000*w)
            assert seen['d']==expected
            if actual is not None:
                assert actual[2]==expected
                accepted_parameter_cases+=1
            component_cases+=1
    result['actual_component_local_log_power_two_cases']=component_cases
    result['actual_component_local_log_accepted_cases']=accepted_parameter_cases
    variable_cases=0;variable_branches=[0,0];partition_cases=0
    for n in range(2,100):
        for sequence in ((1,2,1,3),(2,1,2,1),(3,1,1,2)):
            logs=[0];ar=[];j=0
            while True:
                a=sequence[j%len(sequence)];ar.append(a)
                if logs[-1]+a>n.bit_length()-1:break
                logs.append(logs[-1]+a);j+=1
            for level in range(len(logs)):
                pieces=[f.level_interval_bounds(n,logs,level,i) for i in range(1<<logs[level])]
                assert pieces[0][0]==0 and pieces[-1][1]==n
                assert all(pieces[i][1]==pieces[i+1][0] for i in range(len(pieces)-1))
                assert max(b-a for a,b in pieces)-min(b-a for a,b in pieces)<=1
                if level+1<len(logs):
                    for i,(start,end) in enumerate(pieces):
                        children=[f.level_interval_bounds(n,logs,level+1,(i<<ar[level])+r) for r in range(1<<ar[level])]
                        assert children[0][0]==start and children[-1][1]==end
                partition_cases+=1
            for gapvalue in range(max(ar)+1,max(ar)+4):
                indices=f.active_level_indices(logs,gapvalue)
                assert indices[0]==0 and len(indices)==len(set(indices))
                assert all(logs[indices[i+1]]-logs[indices[i]]<=gapvalue+max(ar) for i in range(len(indices)-1))
                for threshold in range(1,n+1):
                    eligible=[j for j in indices if n//(1<<logs[j])>=threshold]
                    h=eligible[-1];pos=indices.index(h);maximum=(n+(1<<logs[h])-1)//(1<<logs[h])
                    if pos+1<len(indices):
                        nxt=indices[pos+1];assert n//(1<<logs[nxt])<threshold;variable_branches[0]+=1
                    else:
                        assert logs[-1]-logs[h]<gapvalue
                        assert n<(1<<(logs[-1]+ar[-1]));variable_branches[1]+=1
                    assert maximum<=(1<<(gapvalue+max(ar)))*(threshold+1)+1
                    variable_cases+=1
    assert all(variable_branches)
    result['local_variable_partition_levels']=partition_cases
    result['local_variable_active_leaf_cases']=variable_cases
    result['local_variable_leaf_branch_counts']=variable_branches
    depth_cases=0
    # Exponent arithmetic ONLY: L is symbolic logarithmic size; no N=2^L,
    # component construction or asymptotic full tree is allocated or executed.
    for exponent in (160,192,224,256,320):
        L=1<<exponent;B=f.ceil_log2(L+1);Z=f.ceil_root(L**3,4)*f.ceil_root(B,4)
        base=f.ceil_root(L,4)*f.ceil_root(B**3,4);W=5*Z;clip=f.ceil_root(W**20,19)
        M=max(1,27*base*L//(25*clip));gapvalue=(594*Z+124)//125
        denominator=54*base*L-50*L-25*M*M
        assert denominator>0
        count=(25*L*L+denominator-1)//denominator+(clip+M-1)//M+2
        assert (1000*count)**4*B**3<=463**4*L**3
        assert 1000*(gapvalue+M)<=4753*Z
        assert W**3>=(L+1)**2 and clip**19>=W**20
        depth_cases+=1
    result['symbolic_local_depth_gap_witnesses']=depth_cases
    result['all_passed']=True
    result['scope']='Finite component predicates and total fallback only; no asymptotically nontrivial full tree execution or universal extraction/circuit acceptance.'
    return result

if __name__=='__main__':
    request=json.loads(input())
    if request!={}:raise ValueError('empty input object required')
    print(json.dumps(run(),sort_keys=True))
