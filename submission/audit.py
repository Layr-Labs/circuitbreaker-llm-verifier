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
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def table_value(table, x, y):
    return (table >> (2 * x + y)) & 1


def run():
    f = load_family()
    results = {}
    and_tables=[v for v in range(16) if v.bit_count() in (1,3)]
    points=list(itertools.product((0,1),repeat=2))
    from types import SimpleNamespace
    finite_cases=0
    for n in (0,1,2,7,16,64,127,128,129,256,1024):
        for pattern in (0,1,2):
            bits=[0 if pattern==0 else 1 if pattern==1 else j%2 for j in range(n)]
            assert f.evaluate(bits)==0
            finite_cases+=1
    results['zero_cutoff_and_early_rejection_cases']=finite_cases
    eng=f
    roots = 0
    for degree in (*range(1, 13), 16, 20, 25, 50):
        for integer in range(257):
            lo = eng.root_floor(integer, degree)
            hi = eng.ceil_root(integer, degree)
            assert lo**degree <= integer < (lo+1)**degree
            assert integer <= hi**degree
            assert hi == 0 or (hi-1)**degree < integer
            roots += 1
    results['advanced_integer_roots'] = roots
    # Retained regression: the former near-half dimension 100000 exceeded Python's
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
    # Current elimination Section1: finite semantics of promised pairs, not global coverage.
    triangle = [(0, 0), (0, 1), (1, 0)]
    pair_lines = [[(0, 0), (0, 1)], [(0, 0), (1, 0)], [(0, 1), (1, 0)]]
    common_pair_tables = 0
    for table in range(16):
        line = next(line for line in pair_lines
                    if table_value(table, *line[0]) == table_value(table, *line[1]))
        assert all(x*y == 0 for x,y in line)
        assert len(set(line)) == 2
        common_pair_tables += 1
    results['promised_common_child_tables'] = common_pair_tables
    pair_creation = 0
    for table in and_tables:
        alpha,beta,gamma = next((a,b,c) for a,b,c in itertools.product((0,1), repeat=3)
                               if all(table_value(table,x,y) == ((x^a)&(y^b))^c
                                      for x,y in points))
        for child in range(16):
            for newx,newy in triangle:
                x,y = newx^alpha,newy^beta
                gate = table_value(table,x,y)
                assert gate == gamma
                child_values = [table_value(child,gate,z) for z in (0,1)]
                assert child_values in ([0,0],[0,1],[1,0],[1,1])
                for z in (0,1):
                    assert table_value(child,gate,z) == table_value(child,gamma,z)
                    pair_creation += 1
    results['promised_pair_creation_child_assignments'] = pair_creation
    domain_maps = 0
    for q in range(3):
        # First 2q coordinates paired, then two unpaired coordinates.
        n = 2*q+2
        valid = lambda u: all(u[2*j]*u[2*j+1] == 0 for j in range(q))
        for tail in itertools.product((0,1), repeat=n-1):
            # AR changes only the last (unpaired) coordinate by an affine sum.
            for shift in (0,1):
                for w in (0,1):
                    u = list(tail)+[w]
                    mapped = list(tail)+[w ^ (sum(tail)%2) ^ shift]
                    assert valid(u) == valid(mapped)
                    domain_maps += 1
        if q:
            for tail in itertools.product((0,1), repeat=n-1):
                h = tail[0]
                # Reinsert first member p by the three legal PF maps.
                for pairvalues in ((0,h), (h,0), (h^1,h)):
                    mapped = list(pairvalues)+list(tail[1:])
                    if all(tail[2*j-1]*tail[2*j] == 0 for j in range(1,q)):
                        assert valid(mapped)
                        domain_maps += 1
    results['promised_domain_map_assignments'] = domain_maps
    promise_potential = 0
    for extra in range(9):
        # >=two children in paired cleanup; newly unpaired partner <=6.
        for released in (0,6):
            assert -4-(2+extra)+released <= 0
            assert -4-(2+extra)+released <= 0 # paired row cleanup
            promise_potential += 2
        assert -(2+extra)+4-12 <= -10 # both unpaired: add pair
        assert -(2+extra)-4-6 <= -10 # exactly one paired
        assert -(2+extra)-8 <= -10   # two distinct pairs
        promise_potential += 3
    results['promised_potential_cases'] = promise_potential
    # Current elimination Section2's five canonical sharing patterns. Gate parents below
    # are node identifiers, and unconstrained external values are opaque
    # independent bits. This checks the stated patterns, not all DAGs.
    patterns = [
        [('G','x','y'),('Hx','hx','a'),('Hy','hy','b'),('D','G','c')],
        [('G','x','y'),('Hx','hx','G'),('Hy','hy','b'),('E','Hx','c')],
        [('G','x','y'),('Hx','hx','G'),('Hy','hy','Hx'),('F','Hy','c')],
        [('G','x','y'),('B','hx','hy'),('D','G','c'),('E','B','a')],
        [('G','x','y'),('B','hx','hy'),('D','G','B'),('F','D','c')],
    ]
    four_deletion_cases = 0
    for pattern in patterns:
        for gate_table in and_tables:
            forcing = next(x for x in (0,1)
                           if table_value(gate_table,x,0) == table_value(gate_table,x,1))
            for other_tables in itertools.product(range(16), repeat=3):
                mapped = {'x':(None,forcing), 'hx':(None,0), 'hy':(None,0)}
                mapped.update({v:(v,0) for v in ('y','a','b','c')})
                removed = 0
                for (node,pa,pb),table in zip(pattern,(gate_table,)+other_tables):
                    aa,bb = mapped[pa],mapped[pb]
                    variables = {v for v,_ in (aa,bb) if v is not None}
                    assert len(variables) <= 1
                    variable = next(iter(variables),None)
                    vals = []
                    for bit in (0,1):
                        av = aa[1] if aa[0] is None else bit ^ aa[1]
                        bv = bb[1] if bb[0] is None else bit ^ bb[1]
                        vals.append(table_value(table,av,bv))
                    mapped[node] = ((None,vals[0]) if vals[0] == vals[1]
                                    else (variable,vals[0]))
                    removed += 1
                assert removed == 4
                four_deletion_cases += 1
    results['four_object_canonical_table_cases'] = four_deletion_cases
    finite_eight_cases = 0
    for extra in range(9):
        for release in (0,4):
            assert -2-(2+extra)+release <= 0
            finite_eight_cases += 1
        assert -(2+extra)+2-8 <= -8
        assert -(2+extra)-2-4 <= -8
        assert -(4+extra)-4 <= -8
        finite_eight_cases += 3
    results['finite_eight_d_potential_cases'] = finite_eight_cases
    # New polylog-branch components: these exact small cases do not
    # instantiate the immense complete asymptotic construction.
    guarded_poly = eng.PolyPlan(1 << 128)
    assert not guarded_poly.ready
    results['scalar_rejected_plan_cases'] = 1
    twelve_plan=eng.PolyPlan(1<<1024,1<<20)
    assert twelve_plan.lext.su.D==11
    assert twelve_plan.lext.error==eng.ceil_div(1024,28)+10
    results['balanced_initial_dimension_precision_plan_cases']=1
    a1_exact_cases=0
    for n in (16,256,1<<20):
        for dimension in (3,5,11):
            for precision in (1,4,10):
                par=eng.seeded_params(n,max(1,n//2),dimension,1,precision)
                rr=16*dimension+16
                base=(1<<20)*rr*par.h*((1<<precision)*8*par.m)**4
                assert (1<<par.s)**(rr-dimension)>=base**rr
                assert par.s<eng.ceil_log(base**8)
                a1_exact_cases+=1
    results['rational_a1_field_scalar_cases']=a1_exact_cases
    tight_pv_cases=0
    for ambient in range(1,17):
        for entropy in range(1,ambient+1):
            for error in (1,2,3):
                plan=eng.PVPlan(ambient,entropy,error)
                assert plan.s==eng.next_pow2(eng.ceil_log(ambient)+eng.ceil_log(entropy)+plan.a+error+3)
                assert (plan.D-1)*((1<<plan.a)-1)*plan.m*(1<<error)<=1<<plan.s
                assert plan.output<=4*entropy+2*entropy*eng.ceil_log(ambient)//plan.a+2*eng.ceil_log(ambient)+4*plan.a
                tight_pv_cases+=1
    results['tight_rank_condenser_scalar_cases']=tight_pv_cases
    entropy_guards=0
    for length in ((1<<1024),(1<<1024)+1,(1<<2048)):
        bad=eng.PolyPlan(length,1<<eng.ceil_log(length))
        assert not bad.ready and bad.entropy_log>length.bit_length()-1
        entropy_guards+=1
    results['actual_entropy_guard_cases']=entropy_guards
    dimension_cases = 0
    for dimension in (3,5,9):
        # Bounded scalar plans only: no field construction or huge input
        # string is allocated. Precision=3 keeps all integer bit sizes small.
        plan = eng.AffineSeedPlan(1<<120,1<<120,1,3,dimension=dimension)
        assert plan.valid and plan.calls == 1
        assert plan.su.D == dimension
        assert plan.data == (dimension*dimension+1)*plan.su.s
        dimension_cases += 1
    results['polylog_variable_dimension_scalar_cases'] = dimension_cases
    # New B2.2 implementation predicates; toy sizes do NOT meet extractor entropy margins.
    for n,K,M,E in [(1<<100,1<<50,1<<38,3),(1<<120,1<<60,1<<48,7)]:
        plan=eng.HalfSeedPlan(n,K,M,E,1<<50)
        assert plan.valid and plan.data==0
        assert plan.actual_seed==plan.pv.s+eng.ceil_div(2*plan.nw_s**2,plan.overlap_log)+2*plan.nw_s
        assert plan.nw_s&(plan.nw_s-1)==0
        assert K>=32*plan.budget
        assert K*plan.block>=4*plan.overlap_weight*plan.nw_s*M
        old_seed=plan.pv.s+eng.ceil_div(2*(2*plan.nw_s)**2,plan.overlap_log)+2*plan.nw_s
        assert plan.actual_seed<old_seed
    assert not eng.HalfSeedPlan(16,16,4,10,1<<30).valid
    results['direct_half_scalar_valid_plans']=2
    results['direct_half_scalar_rejected_plans']=1
    overlap_cases=0
    for log_weight in (1,2,3,4):
        size=eng.ceil_div(2*4**2,log_weight)+4
        rows=eng.weak_design(4,12,size,1<<log_weight,1)
        masks=[sum(1<<j for j in row) for row in rows]
        for i,row in enumerate(masks):
            assert sum(1<<(row&old).bit_count() for old in masks[:i])<=(1<<log_weight)*i
            overlap_cases+=1
    results['entropy_paid_overlap_design_rows']=overlap_cases
    low=eng.HalfSeedPlan(1<<100,1<<50,1<<36,3,1<<50,1)
    high=eng.HalfSeedPlan(1<<100,1<<50,1<<36,3,1<<50,5)
    assert low.valid and high.valid and high.actual_seed<low.actual_seed
    assert not eng.HalfSeedPlan(1<<100,1<<45,1<<40,3,1<<50,5).valid
    results['entropy_paid_overlap_scalar_plan_cases']=3
    class IdentityPV:
        def value(self,bits,seed):return bits
    toy=SimpleNamespace(n=8,M=2,block=2,nw_s=2,nw_a=2,universe=1600,overlap_weight=2,
                        pv=SimpleNamespace(s=0,output=8))
    direct=eng.HalfSeed(toy,IdentityPV())
    direct_linear=0
    for seed in range(32):
        columns=[eng.to_int(direct.value(eng.to_bits(1<<i,8),seed)) for i in range(8)]
        for x in range(256):
            expected=0
            for i,column in enumerate(columns):
                if x>>i&1:expected^=column
            assert eng.to_int(direct.value(eng.to_bits(x,8),seed))==expected
            direct_linear+=1
    results['direct_half_toy_linearity_assignments']=direct_linear
    results['direct_half_toy_design_rows']=len(direct.rows)
    half3=eng.HalfSeed(SimpleNamespace(n=8,M=3,block=2,nw_s=2,nw_a=2,universe=1600,overlap_weight=2,pv=SimpleNamespace(s=0,output=8)),IdentityPV())
    half4=eng.HalfSeed(SimpleNamespace(n=8,M=4,block=2,nw_s=2,nw_a=2,universe=1600,overlap_weight=2,pv=SimpleNamespace(s=0,output=8)),IdentityPV())
    partial_half=0
    for seed in range(16):
        for x in range(16):
            bits=eng.to_bits(x,8)
            assert half3.value(bits,seed)==half4.value(bits,seed)[:3]
            partial_half+=1
    results['block_half_partial_last_block_cases']=partial_half
    for n,K,M,E in [(1<<100,1<<50,1<<20,3),(1<<120,1<<60,1<<25,7)]:
        plan=eng.ShortSeedPlan(n,K,M,E,1<<50)
        assert plan.valid and plan.data==0
        assert plan.actual_seed==2*plan.hash_s
        assert K>=M+2*E+10
        assert plan.hash_s>=M+2*E+eng.ceil_log(n)+10
    assert not eng.ShortSeedPlan(16,16,4,10,1<<30).valid
    results['direct_short_scalar_valid_plans']=2
    results['direct_short_scalar_rejected_plans']=1
    direct_short=0
    for count in (1,2):
        hash_plan=SimpleNamespace(n=8,M=count,hash_s=2)
        short=eng.ShortSeed(hash_plan)
        for seed in range(16):
            columns=[eng.to_int(short.value(eng.to_bits(1<<i,8),seed)) for i in range(8)]
            for x in range(256):
                expected=0
                for i,column in enumerate(columns):
                    if x>>i&1:expected^=column
                assert eng.to_int(short.value(eng.to_bits(x,8),seed))==expected
                direct_short+=1
        for difference in range(1,256):
            roots=0
            for y in range(4):
                value=0
                for i in range(3,-1,-1):
                    value=eng.field_mul(value,y,2)^((difference>>(2*i))&3)
                roots+=value==0
            assert roots<=3
            counts=[0]*(1<<count)
            for seed in range(16):
                counts[eng.to_int(short.value(eng.to_bits(difference,8),seed))]+=1
            for output in range(1<<count):
                expected=(4-roots)*(1<<(2-count))+(4*roots if output==0 else 0)
                assert counts[output]==expected
    results['direct_short_toy_linearity_assignments']=direct_short
    results['direct_short_toy_difference_distributions']=510
    header_cases = 0
    for n in (4,5,8,9,16,17):
        length = (n-1).bit_length()
        for code in range(1<<length):
            for payload in (0,1):
                bits = [(code>>j)&1 for j in range(length)]+[payload]*(n-length)
                assert f.scale_header(bits) == (length,code+1)
                header_cases += 1
    results['scale_header_parse_cases'] = header_cases
    class DataReader:
        def __init__(self,length):self.length=length
        def bind(self,prefix):return prefix.number(self.length)
    offset_cases = 0
    toy = SimpleNamespace(header=2,data_prefix=8,prefix=10,
        seeds=[DataReader(1) for _ in range(6)],design=DataReader(2))
    for word in range(1<<10):
        for payload in (0,1):
            bits = [(word>>j)&1 for j in range(10)]+[payload]*6
            obj = eng.PolyEngine(toy,bits)
            assert [obj.lext,obj.q,obj.w,obj.one,obj.two,obj.three] == bits[2:8]
            assert obj.design == bits[8]+2*bits[9]
            assert obj.x == bits[10:]+[0]*10
            offset_cases += 1
    results['scale_header_data_offset_cases'] = offset_cases
    zero_offsets=0
    toy2=SimpleNamespace(header=2,data_prefix=2,prefix=4,
        seeds=[DataReader(i) for i in (1,0,0,0,0,0)],design=DataReader(1))
    for word in range(1<<4):
        bits=eng.to_bits(word,4)+[0,1,0,1]
        obj=eng.PolyEngine(toy2,bits)
        assert obj.lext==bits[2]
        assert obj.q==obj.w==obj.one==obj.two==obj.three==0
        assert obj.design==bits[3] and obj.x==bits[4:]+[0]*4
        zero_offsets+=1
    results['direct_half_and_short_zero_data_prefix_cases']=zero_offsets
    tower_nonzero = 0
    for degree in (1,2,4):
        top = 1 << (degree-1)
        tr = 0
        value = top
        for _ in range(degree):
            tr ^= value
            value = eng.field_mul(value,value,degree)
        assert tr == 1
        for value in range(1,1<<degree):
            assert eng.field_pow(value,(1<<degree)-1,degree) == 1
            tower_nonzero += 1
    results['polylog_tower_nonzero_elements'] = tower_nonzero
    modulus = eng.relative_modulus(2,2)
    assert modulus == [2,1,1]
    maps = [[y+4*eng.pv_value(x,y,modulus,2,1,2) for y in range(4)]
            for x in range(16)]
    expanded = 0
    for cardinality in range(1,5):
        for subset in itertools.combinations(range(16),cardinality):
            neighbors = set().union(*(set(maps[x]) for x in subset))
            assert len(neighbors) >= 2*cardinality
            expanded += 1
    linear_pv = 0
    for seed in range(4):
        for x in range(16):
            for y in range(16):
                assert eng.pv_value(x^y,seed,modulus,2,1,2) == (
                    eng.pv_value(x,seed,modulus,2,1,2) ^
                    eng.pv_value(y,seed,modulus,2,1,2))
                linear_pv += 1
    results['polylog_pv_expanded_sets'] = expanded
    results['polylog_pv_linearity_assignments'] = linear_pv
    hash_cases = 0
    for difference in range(1,256):
        counts = [0]*4
        for y in range(16):
            value = (difference & 15) ^ eng.mul(difference>>4,y,0b10011)
            for z in range(16):
                counts[eng.mul(z,value,0b10011)&3] += 1
                hash_cases += 1
        # Degree-one nonzero polynomial: at most one zero evaluation.
        assert max(counts) <= 16*16//4+16
    results['polylog_shift_hash_seed_assignments'] = hash_cases
    conditional_cases = 0
    masks = [3,12,48]
    for selected in range(128):
        remaining = 127^selected
        positions = [j for j in range(7) if remaining>>j&1]
        for slots in range(min(2,len(positions))+1):
            numerator,denominator = eng.weak_design_expectation(masks,selected,remaining,slots)
            total = 0
            subsets = 0
            for extra in itertools.combinations(positions,slots):
                completed = selected|sum(1<<j for j in extra)
                total += sum(2**((completed&prior).bit_count()) for prior in masks)
                subsets += 1
            assert numerator == total and denominator == subsets
            conditional_cases += 1
    results['polylog_conditional_expectation_components'] = conditional_cases
    larger_rows = eng.weak_design(8,2)
    assert set(larger_rows[0]).isdisjoint(larger_rows[1])
    results['polylog_conditional_large_design_rows'] = len(larger_rows)
    design_cases = 0
    for width in (1,2,3):
        rows = eng.weak_design(width,8,4*width*width)
        sets = [set(row) for row in rows]
        for i,row in enumerate(sets):
            assert len(row) == width
            assert 2*sum(2**len(row&prior) for prior in sets[:i]) <= 3*i
            design_cases += 1
    results['polylog_greedy_weak_design_rows'] = design_cases
    shared_design_rows=0
    for width in (1,2,3):
        rows=eng.weak_design(width,8,5*width**2+width,5,4)
        masks=[sum(1<<j for j in row) for row in rows]
        for i,row in enumerate(masks):
            assert 4*sum(1<<(row&old).bit_count() for old in masks[:i])<=5*i
            shared_design_rows+=1
    results['shared_multiplier_five_quarters_design_rows']=shared_design_rows
    raw_cases=0
    raw_design=eng.RecycledHash.__new__(eng.RecycledHash)
    raw_design.p=SimpleNamespace(r=8,code_s=2);raw_design.Q=7
    for x in range(256):
        for y in range(4):
            raw=0
            for index in range(3,-1,-1):
                raw=eng.mul(raw,y,7)^((x>>(2*index))&3)
            for z in range(4):
                recovered=eng.mul(z,raw,7)
                assert raw_design.word(x,y|(z<<2))==recovered
                assert direct.word(x,y|(z<<2))==recovered
                for count in (1,2):
                    assert (recovered&((1<<count)-1))==(raw_design.word(x,y|(z<<2))&((1<<count)-1))
                    raw_cases+=1
    results['shared_multiplier_raw_table_identity_cases']=raw_cases
    recycled_cases = 0
    recycled = object.__new__(eng.RecycledHash)
    recycled.p = SimpleNamespace(r=4,m=3,code_s=2,code_c=1,code_block=1)
    recycled.Q = 7
    recycled.rows = [(0,1),(2,3)]
    for seed in range(0,1<<12,257):
        for left in range(16):
            for right in range(16):
                assert recycled.value(left^right,seed) == (recycled.value(left,seed)^recycled.value(right,seed))
                recycled_cases += 1
    results['weak_hash_component_linearity_assignments'] = recycled_cases
    block_plan=eng.WeakHashPlan(10**9,2*10**8,10**8,10**9)
    assert block_plan.valid and block_plan.code_block==750*eng.ceil_log(10**9)
    assert eng.ceil_div(block_plan.m-block_plan.code_c,block_plan.code_block)<block_plan.m-block_plan.code_c
    results['weak_hash_component_scalar_rows']=eng.ceil_div(block_plan.m-block_plan.code_c,block_plan.code_block)
    assert block_plan.actual_seed==block_plan.code_s+5*block_plan.code_s**2+block_plan.code_s
    assert block_plan.actual_seed<2*block_plan.code_s+(5*(2*block_plan.code_s)**2+1)//2+2*block_plan.code_s
    blocktoy=SimpleNamespace(r=8,m=6,code_s=4,code_c=2,code_block=2,code_a=4,universe=6400)
    generator=eng.RecycledHash(blocktoy,0b10011)
    assert len(generator.rows)==2
    block_linear=0
    for seed in (i*0x010101 for i in range(16)):
        columns=[generator.value(1<<i,seed) for i in range(8)]
        for x in range(256):
            expected=0
            for i,column in enumerate(columns):
                if x>>i&1:expected^=column
            assert generator.value(x,seed)==expected
            block_linear+=1
    results['weak_hash_component_toy_linearity_assignments']=block_linear
    partialtoy=SimpleNamespace(r=8,m=5,code_s=4,code_c=2,code_block=2,code_a=4,universe=6400)
    partial=eng.RecycledHash(partialtoy,0b10011)
    partial_cases=0
    for seed in (i*0x010101 for i in range(16)):
        for x in range(256):
            assert partial.value(x,seed)==(generator.value(x,seed)&31)
            partial_cases+=1
    results['weak_hash_component_partial_last_block_cases']=partial_cases
    packing_cases=0
    for serial in range(64):
        half_seed=serial
        generator_seed=(serial*0x915b)&((1<<20)-1)
        for sample in range(16):
            x=(sample*37)&255
            hp=half3.p;common_z=half_seed&3;uword=half_seed>>2
            expected=[]
            for row in half3.rows:
                y=sum(((uword>>position)&1)<<index for index,position in enumerate(row))
                raw=0
                for index in range(3,-1,-1):raw=eng.mul(raw,y,7)^((x>>(2*index))&3)
                expected.extend(eng.to_bits(eng.mul(common_z,raw,7),min(hp.block,hp.M-len(expected))))
            assert half3.value(eng.to_bits(x,8),half_seed)==expected
            gp=generator.p
            def oracle_word(y,z):
                raw=eng.mul(x>>4,y,0b10011)^(x&15)
                return eng.mul(z,raw,0b10011)
            value=oracle_word(generator_seed&15,(generator_seed>>4)&15)&3
            common_z=(generator_seed>>8)&15;uword=generator_seed>>12;offset=gp.code_c
            for row in generator.rows:
                y=sum(((uword>>position)&1)<<index for index,position in enumerate(row))
                count=min(gp.code_block,gp.m-offset)
                value|=(oracle_word(y,common_z)&((1<<count)-1))<<offset;offset+=count
            assert generator.value(x,generator_seed)==value
            packing_cases+=2
    results['shared_multiplier_seed_packing_cases']=packing_cases
    layouts = 0
    for seed_length in (16,32,64,128,256):
        for scale in (1,4,9):
            rr,mm,bb,uu = eng.layout(seed_length,5,scale)
            assert eng.ceil_log(uu)+mm+bb <= seed_length
            rr += 1
            mm = rr//5
            bb = eng.root_floor(scale,2)*eng.ceil_log(rr)//750
            uu = 2*bb+2*eng.ceil_log(rr)+10
            assert eng.ceil_log(uu)+mm+bb > seed_length
            layouts += 1
    results['polylog_layout_maximality_cases'] = layouts
    balanced_layouts=0
    for seed_length in (32,64,128):
        for scale in (1,10**6,10**12):
            rr,mm,bb,uu=eng.layout(seed_length,scale=scale)
            assert mm==rr//125
            assert bb==eng.root_floor(scale,2)*eng.ceil_log(rr)//(30*125**2)
            assert eng.ceil_log(uu)+mm+bb<=seed_length
            balanced_layouts+=1
    results['balanced_default_layout_cases']=balanced_layouts
    thresholds = 0
    for rr,mm,bb,uu in ((3,2,1,6),(5,2,2,8),(8,3,2,8)):
        threshold = eng.biased_threshold(rr,mm,bb,uu)
        den,B,M,R = 1<<uu,1<<bb,1<<mm,1<<rr
        assert R*(den**B-(den-threshold)**B)**M <= den**(B*M)
        assert R*(den**B-(den-threshold-1)**B)**M > den**(B*M)
        thresholds += 1
    results['polylog_exact_threshold_cases'] = thresholds
    conditional_cases = 0
    for length in range(1,5):
        for events in itertools.permutations(range(1,8),length):
            survivors = set(range(8))
            for i,mask in enumerate(events):
                event = {x for x in range(8) if x&mask == mask}
                conditional = Fraction(len(survivors&event),len(survivors))
                probability = Fraction(len(event),8)
                overlap = sum((Fraction(sum((x&mask == mask) and
                    (x&other == other) for x in range(8)),8)
                    for other in events[:i] if mask&other),Fraction(0))
                assert conditional >= probability-overlap
                survivors -= event
                conditional_cases += 1
    results['polylog_janson_conditional_cases'] = conditional_cases
    class MockLinear:
        def __init__(self,tag,n,m):self.tag,self.n,self.m,self.calls=tag,n,m,0
        def value(self,bits,seed):
            self.calls += 1
            padded = list(bits)+[0]*(self.n-len(bits))
            return [sum(bit for j,bit in enumerate(padded)
                        if ((seed+self.tag*(a+1)+j*j)>>(j%3))&1)%2
                    for a in range(self.m)]
    row_cases = 0
    for value in range(64):
        for index in range(4):
            obj = eng.PolyEngine.__new__(eng.PolyEngine)
            obj.p = SimpleNamespace(ell=1,d=20,h=2,alphabet=2,blocks=2,d0=2)
            obj.x = [(value>>j)&1 for j in range(6)]
            obj.lext = MockLinear(1,6,24)
            obj.q = MockLinear(2,20,1)
            obj.w = MockLinear(3,6,1)
            obj.one = MockLinear(4,6,1)
            obj.two = MockLinear(5,24,18)
            obj.three = MockLinear(6,24,20)
            actual = obj.row_bit(index)
            assert obj.one.calls == obj.two.calls == 2
            assert obj.three.calls == 1
            Y = obj.lext.value(obj.x,index)
            current = Y[:20]
            digits = [int(x) for x in format(index,'02b')]
            for position,digit in enumerate(digits):
                R0 = obj.one.value(obj.x,eng.to_int(current))
                Q = obj.two.value(Y,eng.to_int(R0))
                S = Q[:1]
                for time in range(digit+1):
                    R = obj.w.value(obj.x,eng.to_int(S))
                    if time<digit:S = obj.q.value(Q,eng.to_int(R))
                if position == 0:current = obj.three.value(Y,eng.to_int(R))
            assert actual == R[0]
            row_cases += 1
    results['polylog_refreshed_index_mock_cases'] = row_cases
    # Exact combinatorial tree cap; no asymptotic source margin is inferred.
    def tree_bound(weights,n):
        total=sum(weights);probs=[Fraction(x,total) for x in weights];prefix={}
        for x,probability in enumerate(probs):
            for i in range(n+1):
                key=(i,x&((1<<i)-1));prefix[key]=prefix.get(key,Fraction(0))+probability
        cases=0
        for selected in range(1,1<<n):
            for entropy in range(n+1):
                cap=Fraction(1,1<<entropy);masses={}
                for x,probability in enumerate(probs):
                    if not probability:continue
                    product_probability=Fraction(1)
                    for i in range(n):
                        if selected>>i&1:
                            product_probability*=prefix[i+1,x&((1<<(i+1))-1)]/prefix[i,x&((1<<i)-1)]
                    if product_probability<=cap:
                        word=x&selected;masses[word]=masses.get(word,Fraction(0))+probability
                assert all(mass<=cap for mass in masses.values())
                cases+=1
        return cases
    flat_trees=0
    for n in (1,2,3):
        for support in range(1,1<<(1<<n)):
            flat_trees+=tree_bound([(support>>x)&1 for x in range(1<<n)],n)
    weighted_trees=0
    for weights in itertools.product(range(3),repeat=4):
        if sum(weights):weighted_trees+=tree_bound(weights,2)
    results['entropy_sampler_flat_tree_caps']=flat_trees
    results['entropy_sampler_nonuniform_tree_caps']=weighted_trees
    vandermonde=0;field4_positions=[]
    for degree,Q,length in ((2,7,4),(3,11,3)):
        sampler=eng.EntropySampler(SimpleNamespace(w=degree,k=length,length=length),Q)
        counts={}
        for seed in range(1<<(degree*length)):
            coordinate_outputs=[]
            for coordinate in range(degree):
                truth=sum(((index>>coordinate)&1)<<index for index in range(1<<degree))
                coordinate_outputs.append(sampler.value(truth,seed))
            positions=tuple(sum(output[index]<<coordinate for coordinate,output in enumerate(coordinate_outputs)) for index in range(length))
            counts[positions]=counts.get(positions,0)+1
            if degree==2:field4_positions.append(positions)
            vandermonde+=1
        assert len(counts)==1<<(degree*length) and set(counts.values())=={1}
    results['entropy_sampler_vandermonde_assignments']=vandermonde
    mean=Fraction(6,4);second=Fraction(0);fourth=Fraction(0)
    for positions in field4_positions:
        collision=sum(positions[i]==positions[j] for i in range(4) for j in range(i))
        second+=(collision-mean)**2;fourth+=(collision-mean)**4
    assert second/len(field4_positions)==Fraction(9,8)
    for t,moment in ((1,second/len(field4_positions)),(2,fourth/len(field4_positions))):
        constant=(2*t)*(2*t)**(4*t)*(1<<(2*t))
        assert moment<=constant*4**t
    results['entropy_sampler_collision_moments']=2
    giant=1<<512
    bootstrap_plan=eng.DesignPlan(giant,giant//125,1<<100)
    assert bootstrap_plan.valid and bootstrap_plan.sampler.k==1024
    assert bootstrap_plan.rounds<=1<<41
    assert bootstrap_plan.helper.valid() and bootstrap_plan.helper.m>=bootstrap_plan.weak.actual_seed
    assert bootstrap_plan.data==0
    assert bootstrap_plan.actual_seed==2*bootstrap_plan.code_s+bootstrap_plan.rounds*(2*bootstrap_plan.sampler.seed_length+bootstrap_plan.helper.d)
    assert not eng.DesignPlan(10**9,10**8,10**12).valid
    assert not eng.DesignPlan(giant,giant//125,1).valid
    results['bootstrap_positive_scalar_plans']=1
    results['bootstrap_rejected_scalar_plans']=2
    class MockHelper:
        def value(self,bits,seed,count):return eng.to_bits(eng.to_int(bits)^(seed*29),count)
    weak_parameters=SimpleNamespace(r=8,m=2,code_s=4,code_c=0,code_block=2,actual_seed=8)
    weak=eng.RecycledHash.__new__(eng.RecycledHash);weak.p=weak_parameters;weak.Q=0b10011;weak.rows=[(0,1,2,3)]
    sampler_parameters=SimpleNamespace(w=3,k=2,length=8,seed_length=6)
    sampler=eng.EntropySampler(sampler_parameters,0b1011)
    bootstrap=eng.Design.__new__(eng.Design)
    bootstrap.p=SimpleNamespace(r=8,m=4,code_s=4,code_c=1,sampler=sampler_parameters,helper=SimpleNamespace(d=3),weak=weak_parameters,rounds=2)
    bootstrap.Q=0b10011;bootstrap.sampler=sampler;bootstrap.helper=MockHelper();bootstrap.weak=weak
    packing=0
    for serial in range(64):
        seed=(serial*0x152915b)&((1<<38)-1)
        for sample in range(16):
            value=(sample*37)&255
            def sample_oracle(sample_seed):
                constant=sample_seed&7;slope=(sample_seed>>3)&7
                return sum(((value>>(eng.mul(slope,j,0b1011)^constant))&1)<<j for j in range(8))
            raw=(value&15)^eng.mul(value>>4,seed&15,0b10011)
            expected=eng.mul((seed>>4)&15,raw,0b10011)&1
            cursor=seed>>8;offset=1
            for _ in range(2):
                left=cursor&63;cursor>>=6;right=cursor&63;cursor>>=6;helper_seed=cursor&7;cursor>>=3
                first=sample_oracle(left);second_source=sample_oracle(right)
                internal=(second_source^(helper_seed*29))&255
                raw=(first&15)^eng.mul(first>>4,(internal>>4)&15,0b10011)
                block=eng.mul(internal&15,raw,0b10011)
                count=min(2,4-offset);expected|=(block&((1<<count)-1))<<offset;offset+=count
            assert bootstrap.value(value,seed)==expected
            packing+=1
    results['bootstrap_backward_seed_packing_cases']=packing
    # Exact log-field rewrite, without an exponentially large planner integer.
    exact_logs=0
    for n in (16,256,1<<20):
        for dimension in (2,4,11,1024):
            for m in (1,2,8):
                for error in (1,4,10):
                    par=eng.seeded_params(n,n//2,dimension,m,error)
                    R=16*dimension+16
                    oldbase=(1<<20)*R*par.h*((1<<error)*8*m)**4
                    assert par.s==eng.ceil_div(R*eng.ceil_log(oldbase),R-dimension)
                    exact_logs+=1
    results['logarithmic_a1_field_rewrite_cases']=exact_logs
    # Giant SCALAR plans only: no field, seed array or input is constructed.
    L=1<<4096;z=eng.ceil_log(L);scale=1<<100;C=1<<30;dimension=max(2,z)
    h=eng.next_pow2(scale*L**2*z)
    ell=eng.next_pow2(scale**3*L**2*z**2)
    d=eng.next_pow2(scale**9*L**4*z**5)
    Y=eng.next_pow2(scale**11*L**7*z**6)
    k=eng.next_pow2(scale**13*L**9*z**7)
    q0=(4*h+C*h*z*z+10)*ell
    scalar_cases=0
    for ambient,threshold,output in ((d,C*h*z*z*ell,ell),(Y,C*h*z*z*q0,q0),(Y,C*h*z*z*d,d)):
        full=eng.AffineSeedPlan(ambient,threshold,output,h,ell,dimension=dimension)
        par=eng.AffineSeedPlan(ambient,threshold,output,h,ell,dimension=dimension,companion=True)
        assert par.valid and par.calls<=2 and par.su.D==dimension
        assert par.data==(dimension+1)*par.su.s
        assert par.su==full.su and par.actual_seed==full.actual_seed
        assert full.data-par.data==(dimension**2-dimension)*par.su.s
        assert par.actual_seed==par.pv.s+par.calls*par.su.d<=ell
        assert par.su.h<=128*z and 16*dimension+16<=32*z
        assert par.su.s<5*h and par.pv.s<5*h and par.actual_seed<25*z*h
        scalar_cases+=1
    assert q0-2*h*ell>=C*h*z*z*ell and d>=q0
    assert Y-100*5*L*h*d>=C*h*z*z*max(d,q0)
    assert k-h*Y-100*5*L*h*2*ell>=10*d+C*h*z*z*ell
    results['variable_dimension_positive_scalar_plans']=scalar_cases
    streamed=0
    for dimension in (2,3,8,17):
        for ambient in (1,4,16):
            H=dimension*eng.ceil_root(ambient,dimension)
            exponents=list(eng.monomials(ambient,dimension,H))
            assert len(exponents)==len(set(exponents))==ambient
            assert all(len(t)==dimension and sum(t)<H for t in exponents)
            assert [sum(t) for t in exponents]==sorted(sum(t) for t in exponents)
            streamed+=1
    results['variable_dimension_streamed_monomial_cases']=streamed
    # This toy tests actual A1 package parsing and AffineSeed seed slicing only.
    supar=eng.SeedParams(4,2,2,4,1,1,2)
    package=[1,1,0,1,1,0,1,1,0,1]
    reader=eng.Prefix(package);su=reader.seeded(supar)
    assert reader.pos==10 and su.Q==7
    assert su.A==[[2,1],[3,2]]
    class TinyPV:
        def value(self,bits,seed):return bits[seed%4:]+bits[:seed%4]
    par=SimpleNamespace(n=4,pv=SimpleNamespace(s=2),su=supar,calls=2,M=2)
    wrapper=eng.AffineSeed(par,TinyPV(),su)
    packing=0
    for seedno in range(64):
        seed=(seedno*257)&((1<<14)-1)
        for value in range(16):
            bits=eng.to_bits(value,4);condensed=TinyPV().value(bits,seed&3)
            expected=su.value(condensed,(seed>>2)&63)+su.value(condensed,(seed>>8)&63)
            assert wrapper.value(bits,seed)==expected
            assert wrapper.value(bits,seed|(7<<14))==expected
            packing+=1
    results['charged_internal_package_and_call_packing_cases']=packing
    canonical_cases=0
    for degree,dimension in ((1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(4,1)):
        Q,matrix=eng.canonical_cyclic_matrix(degree,dimension)
        assert Q==eng.canonical_modulus(degree)
        for factor_degree in range(1,degree//2+1):
            for low in range(1<<factor_degree):
                assert eng.polynomial_remainder(Q,(1<<factor_degree)|low)
        start=(1,)+(0,)*(dimension-1);value=start;seen=set()
        for step in range((1<<(degree*dimension))-1):
            assert any(value) and value not in seen;seen.add(value)
            value=tuple(__import__('functools').reduce(int.__xor__,(eng.mul(a,b,Q) for a,b in zip(row,value)),0) for row in matrix)
        assert value==start and len(seen)==(1<<(degree*dimension))-1
        # No earlier coefficient tuple has a full companion orbit.
        coeffs=[row[-1] for row in matrix]
        number=sum(c<<(degree*j) for j,c in enumerate(coeffs));period=len(seen)
        for earlier in range(number):
            coefficients=[(earlier>>(degree*j))&((1<<degree)-1) for j in range(dimension)]
            value=list(start)
            for step in range(1,period+1):
                top=value[-1];value=[eng.mul(coefficients[j],top,Q)^(value[j-1] if j else 0) for j in range(dimension)]
                if not any(value) or value==list(start):break
            assert not(step==period and value==list(start))
        canonical_cases+=1
    results['canonical_small_field_and_companion_cycles']=canonical_cases
    # Canonical and paid wrappers have identical scalar seeds and different data.
    canonical_plans=0
    for ambient,entropy,output,error in ((1<<1024,1<<128,16,8),(1<<2048,1<<160,32,10)):
        paid=eng.AffineSeedPlan(ambient,entropy,output,error,dimension=11)
        canon=eng.AffineSeedPlan(ambient,entropy,output,error,dimension=11,canonical=True)
        assert paid.valid and canon.valid
        assert paid.actual_seed==canon.actual_seed and paid.su==canon.su
        assert paid.data==122*paid.su.s and canon.data==0
        canonical_plans+=1
    results['canonical_initial_scalar_data_savings']=canonical_plans
    assert bootstrap_plan.data==0
    assert bootstrap_plan.helper.s<=2048*bootstrap_plan.w
    results['canonical_b5_scalar_data_savings']=1
    companion_cases=0
    for dimension in (2,3,8):
        pars=eng.SeedParams(4,2,dimension,4*dimension,1,1,2)
        for coefficients in __import__('itertools').product(range(4),repeat=2):
            coefficients=list(coefficients)+[0]*(dimension-2)
            package=eng.to_bits(3,2)+[bit for coefficient in coefficients for bit in eng.to_bits(coefficient,2)]
            reader=eng.Prefix(package);obj=reader.companion_seeded(pars)
            assert reader.pos==(dimension+1)*2 and obj.Q==7
            expected=eng.companion_matrix(coefficients)
            assert obj.A==expected
            # Compare independently packed FULL matrix using the generic parser.
            paid_package=eng.to_bits(3,2)+[bit for row in expected for coefficient in row for bit in eng.to_bits(coefficient,2)]
            paid_reader=eng.Prefix(paid_package);paid=paid_reader.seeded(pars)
            assert paid_reader.pos==(dimension**2+1)*2 and paid.A==obj.A
            for seed in (0,1,3,5,13,31):
                for value in (0,1,3,7,15):
                    bits=eng.to_bits(value,4)
                    assert obj.value(bits,seed)==paid.value(bits,seed)
                    companion_cases+=1
    results['compressed_companion_package_output_cases']=companion_cases
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
