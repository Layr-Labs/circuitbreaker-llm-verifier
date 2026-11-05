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
    tower=0
    for degree in (1,2,4):
        for value in range(1,1<<degree):
            assert eng.field_pow(value,(1<<degree)-1,degree)==1
            tower+=1
    results['tower_nonzero_orders']=tower
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
    layouts = 0
    for seed_length in (16,32,64,128,256):
        for scale in (1,4,9):
            rr,mm,bb,uu = eng.layout(seed_length,5,scale)
            assert mm+bb <= seed_length
            rr += 1
            mm = rr//5
            bb = eng.root_floor(scale,2)*eng.ceil_log(rr)//750
            uu = 2*bb+2*eng.ceil_log(rr)+10
            assert mm+bb > seed_length
            layouts += 1
    results['polylog_layout_maximality_cases'] = layouts
    balanced_layouts=0
    for seed_length in (32,64,128):
        for scale in (1,10**6,10**12):
            rr,mm,bb,uu=eng.layout(seed_length,scale=scale)
            assert mm==rr//125
            assert bb==eng.root_floor(scale,2)*eng.ceil_log(rr)//(30*125**2)
            assert mm+bb<=seed_length
            balanced_layouts+=1
    results['balanced_default_layout_cases']=balanced_layouts
    def envelope(L,scale):
        z=eng.ceil_log(L);C=1<<30
        h=eng.next_pow2(scale*L**2);e=eng.next_pow2(scale*L**2*z)
        ell=eng.next_pow2(128*z*e)
        ep=eng.next_pow2(512*h*ell+4*e)
        q=(4*h+C*e*z*z+10)*ell
        d=eng.next_pow2(16*q+128*ep+128*L+1280)
        Y=eng.next_pow2(10000*L*h*d+4*C*e*z*z*d)
        return h,e,ell,ep,q,d,Y,2*h*Y
    L=1<<4096;z=eng.ceil_log(L);scale=1<<100;C=1<<30;dimension=max(2,z)
    h,error,ell,eprime,q0,d,Y,k=envelope(L,scale)
    scalar_cases=0
    for ambient,threshold,output in ((d,C*error*z*z*ell,ell),(Y,C*error*z*z*q0,q0),(Y,C*error*z*z*d,d)):
        full=eng.AffineSeedPlan(ambient,threshold,output,error,ell,dimension=dimension)
        companion=eng.AffineSeedPlan(ambient,threshold,output,error,ell,dimension=dimension,companion=True)
        par=eng.AffineSeedPlan(ambient,threshold,output,error,ell,dimension=dimension,primitive=True)
        assert par.valid and par.calls<=2 and par.su.D==dimension
        assert par.data==dimension*par.su.s-eng.primitive_package_saving(dimension*par.su.s)
        assert par.su==full.su and par.actual_seed==full.actual_seed
        assert full.data-companion.data==(dimension**2-dimension)*par.su.s
        assert companion.data-par.data==par.su.s+eng.primitive_package_saving(dimension*par.su.s)
        assert par.actual_seed==par.pv.s+par.calls*par.su.d<=ell
        assert par.su.h<=128*z and 16*dimension+16<=32*z
        assert par.su.s<5*error and par.pv.s<5*error and par.actual_seed<25*z*error
        scalar_cases+=1
    assert q0-2*h*ell>=C*error*z*z*ell and d>=q0
    assert Y-100*5*L*h*d>=C*error*z*z*max(d,q0)
    assert k-h*Y-100*5*L*h*2*ell>=10*d+C*error*z*z*ell
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
    giant=1<<512
    bootstrap_plan=eng.DesignPlan(giant,giant//125,1<<100)
    assert bootstrap_plan.valid and bootstrap_plan.rounds<=1<<41
    assert bootstrap_plan.data==0
    assert bootstrap_plan.helper.s<=2048*bootstrap_plan.w
    results['canonical_b5_scalar_data_savings']=1
    primitive_products=primitive_orbits=primitive_outputs=primitive_parsers=0
    for field_degree,dimension,R in ((1,2,7),(1,3,11),(2,2,19),(2,3,67),(3,2,67),(3,3,529)):
        t=field_degree*dimension;period=(1<<t)-1
        ring_value=1;ring_seen=set()
        for _ in range(period):
            assert ring_value and ring_value not in ring_seen;ring_seen.add(ring_value)
            ring_value=eng.mul(ring_value,2,R)
        assert ring_value==1 and len(ring_seen)==period
        pars=eng.SeedParams(4,2,dimension,4*dimension,2,1,field_degree)
        obj=eng.primitive_seeded(pars,R)
        assert isinstance(obj,eng.RelativeSeeded) and obj.field.vectors[0]==1
        for left in range(1<<field_degree):
            assert obj.multiply(left,1)==left and obj.multiply(left,0)==0
            assert obj.field.decoder.solve(obj.field.encode(left))==left
            for right in range(1<<field_degree):
                assert obj.field.encode(obj.multiply(left,right))==eng.mul(obj.field.encode(left),obj.field.encode(right),R)
                primitive_products+=1
        def embed(vector):
            out=0;alpha_power=1
            for coefficient in vector:
                out^=eng.mul(alpha_power,obj.field.encode(coefficient),R)
                alpha_power=eng.mul(alpha_power,2,R)
            return out
        vector=[1]+[0]*(dimension-1);seen=set()
        for _ in range(period):
            assert any(vector) and tuple(vector) not in seen;seen.add(tuple(vector))
            previous=embed(vector);vector=[obj.dot(row,vector) for row in obj.A]
            assert embed(vector)==eng.mul(previous,2,R)
            primitive_orbits+=1
        assert vector==[1]+[0]*(dimension-1)
        package_index=next(i for i in range(1<<(t-2)) if eng.reciprocal_package_polynomial(t,i)==R)
        reader=eng.Prefix(eng.to_bits(package_index,t-2)+[1,0,1])
        parsed=reader.primitive_seeded(pars);assert reader.pos==t-2
        assert isinstance(parsed,eng.RelativeSeeded) and parsed.Q==R
        primitive_parsers+=1
        mask=(1<<field_degree)-1
        for seed_number in range(32):
            seed=(seed_number*37)&((1<<((dimension+1)*field_degree))-1)
            for value in range(16):
                bits=eng.to_bits(value,4)
                vector=[(seed>>(j*field_degree))&mask for j in range(dimension)]
                zword=(seed>>(dimension*field_degree))&mask;expected=[]
                for _ in range(pars.m):
                    vector=[obj.dot(row,vector) for row in obj.A]
                    ring_result=0
                    for bit,exponents in zip(bits,eng.monomials(4,dimension,pars.h)):
                        if bit:
                            term=1
                            for coefficient,exponent in zip(vector,exponents):
                                term=eng.mul(term,eng.power(obj.field.encode(coefficient),exponent,R),R)
                            ring_result^=term
                    coordinate=obj.field.decoder.solve(ring_result);assert coordinate is not None
                    expected.append((coordinate&zword).bit_count()%2)
                assert obj.value(bits,seed)==expected==parsed.value(bits,seed)
                assert obj.value(bits,seed,count=1)==expected[:1]
                primitive_outputs+=1
    results['primitive_binary_field_product_cases']=primitive_products
    results['primitive_binary_relative_orbit_steps']=primitive_orbits
    results['primitive_binary_seeded_output_cases']=primitive_outputs
    results['primitive_binary_package_parser_cases']=primitive_parsers
    bad_zero=bad_total=0
    for middle in range(32):
        R=(1<<6)|1|(middle<<1);pars=eng.SeedParams(4,2,3,12,2,1,2)
        obj=eng.primitive_seeded(pars,R)
        if isinstance(obj,eng.ZeroSeeded):
            assert obj.value([1,0,1,1],37)==[0,0];bad_zero+=1
        for seed in (0,3,7,19):
            for left in range(8):
                right=(left*3)&15
                a=obj.value(eng.to_bits(left,4),seed);b=obj.value(eng.to_bits(right,4),seed)
                c=obj.value(eng.to_bits(left^right,4),seed)
                assert c==[x^y for x,y in zip(a,b)];bad_total+=1
    assert bad_zero>0
    results['primitive_binary_guarded_zero_packages']=bad_zero
    results['primitive_binary_arbitrary_package_linearity_cases']=bad_total
    reciprocal_cases=0
    for degree in range(3,14):
        length=degree-1;values=[]
        for index in range(1<<(degree-2)):
            R=eng.reciprocal_package_polynomial(degree,index)
            if R is None:continue
            middle=(R>>1)&((1<<length)-1)
            reverse=sum(((middle>>j)&1)<<(length-1-j) for j in range(length))
            assert middle<reverse
            values.append(middle);reciprocal_cases+=1
        expected=[middle for middle in range(1<<length) if middle<sum(((middle>>j)&1)<<(length-1-j) for j in range(length))]
        assert len(values)==len(set(values)) and set(values)==set(expected)
        assert len(values)==((1<<length)-(1<<((length+1)//2)))//2
    assert eng.reciprocal_package_polynomial(2,0)==7
    results['reciprocal_package_decoder_cases']=reciprocal_cases
    class GroupDesign:
        def value(self,vertex,seed):return (vertex^seed)&1
    grouped_cases=0;grouped_values=set()
    for encoded in range(1<<12):
        words=[(encoded>>(3*j))&7 for j in range(4)]
        obj=eng.PolyEngine.__new__(eng.PolyEngine)
        obj.p=SimpleNamespace(r=2,m=1,seed_b=1,u=3)
        obj.design=GroupDesign();calls=[]
        def row_word(index,count):
            calls.append((index,count));return eng.to_bits(words[index],count)
        obj.row_word=row_word
        actual=obj.evaluate();threshold=eng.biased_threshold(2,1,1,3)
        groups=[word>=threshold for word in words]
        expected=int(all(any(all(groups[z*2+(((vertex^z)&1)^shift)] for z in range(2)) for shift in range(2)) for vertex in range(4)))
        assert actual==expected and calls==[(j,3) for j in range(4)]
        grouped_values.add(actual);grouped_cases+=1
    assert grouped_values=={0,1}
    results['whole_row_group_packing_cases']=grouped_cases
    unit_codec_cases=unit_invalid=0
    for degree in (32,64,512):
        saving=eng.primitive_package_saving(degree)
        modulus,components=eng.small_factor_code(degree.bit_length()-3)
        d=modulus.bit_length()-1;phi=1
        for factor,radix,partial,inverse in components:
            phi*=radix
            assert eng.polynomial_remainder(eng.binary_polynomial_product(partial,inverse),factor)==1
        assert phi<=1<<(d-saving) and d<=degree//2
        indices=range(1<<(d-saving)) if degree==32 else (0,1,2,phi-1,phi)
        for quotient_low in (0,1,7):
            for index in indices:
                code=(index<<(degree-d))|quotient_low
                R=eng.decode_internal_polynomial(degree,code)
                if index>=phi:
                    assert R is None;unit_invalid+=1;continue
                assert R.bit_length()==degree+1 and R&1 and R.bit_count()%2==1
                quotient,remainder=eng.binary_polynomial_divrem(R,modulus)
                assert quotient==(1<<(degree-d))|quotient_low
                encoded=0;weight=1
                for factor,radix,partial,inverse in components:
                    residue=eng.polynomial_remainder(R,factor)
                    assert 0<residue<=radix
                    encoded+=(residue-1)*weight;weight*=radix
                assert encoded==index;unit_codec_cases+=1
    results['small_factor_unit_codec_cases']=unit_codec_cases
    results['small_factor_unit_invalid_codes']=unit_invalid
    # A positive degree32 primitive fixture, certified by its complete order.
    degree=32;R=0x1000000af;period=(1<<degree)-1
    prime_divisors=(3,5,17,257,65537)
    assert __import__('math').prod(prime_divisors)==period
    for prime in prime_divisors:
        assert all(prime%divisor for divisor in range(2,eng.root_floor(prime,2)+1))
    assert eng.power(2,period,R)==1
    assert all(eng.power(2,period//prime,R)!=1 for prime in prime_divisors)
    modulus,components=eng.small_factor_code(degree.bit_length()-3)
    quotient,remainder=eng.binary_polynomial_divrem(R,modulus);d=modulus.bit_length()-1
    index=0;weight=1
    for factor,radix,partial,inverse in components:
        residue=eng.polynomial_remainder(R,factor);assert residue
        index+=(residue-1)*weight;weight*=radix
    code=(index<<(degree-d))|(quotient^ (1<<(degree-d)))
    assert eng.decode_internal_polynomial(degree,code)==R
    params=eng.SeedParams(16,8,8,16,2,1,4)
    reader=eng.Prefix(eng.to_bits(code,degree-eng.primitive_package_saving(degree)))
    compiled=reader.primitive_seeded(params)
    assert reader.pos==degree-eng.primitive_package_saving(degree)
    assert isinstance(compiled,eng.RelativeSeeded)
    products=0
    for left in range(16):
        for right in range(16):
            assert compiled.field.encode(compiled.multiply(left,right))==eng.mul(compiled.field.encode(left),compiled.field.encode(right),R)
            products+=1
    results['small_factor_positive_primitive32_products']=products
    # Actual scalar sharing/revalidation, with no field or matrix allocation.
    L=1<<4096;z=eng.ceil_log(L);scale=1<<100;C=1<<30;dimension=max(2,z)
    h,error,ell,eprime,q0,d,Y,k=envelope(L,scale)
    shared_plans=[]
    for ambient,threshold,output in ((d,C*error*z*z*ell,ell),(Y,C*error*z*z*q0,q0),(Y,C*error*z*z*d,d)):
        plan=eng.AffineSeedPlan(ambient,threshold,output,error,ell,dimension=dimension,primitive=True)
        assert plan.valid;shared_plans.append(plan)
    old_data=sum(plan.data for plan in shared_plans)
    common_degree=max(plan.su.s for plan in shared_plans)
    old_seed_params=[plan.su for plan in shared_plans]
    old_seed_lengths=[plan.actual_seed for plan in shared_plans]
    one_package=dimension*common_degree-eng.primitive_package_saving(dimension*common_degree)
    for plan,old,old_seed in zip(shared_plans,old_seed_params,old_seed_lengths):
        plan.use_shared_field(common_degree)
        assert plan.valid and plan.shared and plan.data==0
        assert plan.su.s==common_degree and plan.su.n==old.n and plan.su.m==old.m and plan.su.h==old.h
        assert plan.actual_seed>=old_seed and plan.actual_seed<=ell
    assert one_package<old_data
    results['common_field_positive_scalar_plans']=len(shared_plans)
    # Shared C3 state is reparameterized for each source/output plan.
    master=eng.SeedParams(3,2,2,4,1,1,2)
    template=eng.primitive_seeded(master,19)
    assert isinstance(template,eng.RelativeSeeded)
    class IdentityPlan:
        s=2
        def bind(self):return self
        def value(self,bits,seed):return bits[seed%len(bits):]+bits[:seed%len(bits)]
    wrappers=[];shared_outputs=0
    for ambient,output,calls in ((4,1,1),(5,2,2),(6,3,2),(4,2,1)):
        params=eng.SeedParams(ambient,ambient,2,4,output,1,2)
        plan=eng.AffineSeedPlan.__new__(eng.AffineSeedPlan)
        plan.n=ambient;plan.M=calls*output-(calls>1);plan.calls=calls;plan.su=params
        plan.pv=IdentityPlan();plan.shared=True;plan.data=0
        bound=plan.bind(eng.Prefix([]),template)
        independent=eng.AffineSeed(plan,plan.pv.bind(),eng.primitive_seeded(params,19))
        assert bound.su.p==params and bound.su.field is template.field and bound.su.A is template.A
        for seed in range(32):
            for value in range(1<<ambient):
                bits=eng.to_bits(value,ambient)
                assert bound.value(bits,seed)==independent.value(bits,seed)
                shared_outputs+=1
        zero=plan.bind(eng.Prefix([]),eng.ZeroSeeded(master))
        assert zero.su.p==params and zero.value([1]*ambient,13)==[0]*plan.M
        wrappers.append(plan)
    results['common_field_reparameterized_output_cases']=shared_outputs
    class NoData:
        shared=False
        def bind(self,prefix):return 0
    # One input code is consumed ONCE, even though four maps use it.
    shared_cursor=0
    for code in (0,2):
        for payload in range(16):
            toy=SimpleNamespace(header=2,prefix=4,data_prefix=2,shared_parameters=master,
                seeds=[NoData(),wrappers[0],wrappers[1],NoData(),wrappers[2],wrappers[3]],design=NoData())
            bits=[0,0]+eng.to_bits(code,2)+eng.to_bits(payload,4)
            engine=eng.PolyEngine(toy,bits)
            assert engine.x==bits[4:]+[0]*4
            for index,obj in enumerate((engine.q,engine.w,engine.two,engine.three)):
                assert obj.su.p==wrappers[index].su
                if code==0:assert isinstance(obj.su,eng.RelativeSeeded)
                else:assert isinstance(obj.su,eng.ZeroSeeded)
            if code==0:assert engine.q.su.field is engine.w.su.field is engine.two.su.field is engine.three.su.field
            shared_cursor+=1
    results['common_field_single_prefix_consumption_cases']=shared_cursor
    def gamma(scale):
        word=[int(bit) for bit in bin(scale)[2:]]
        return [0]*(len(word)-1)+word
    gamma_cases=0
    for n in (0,1,2,3,7,16,31,64,127,128,129,256,1024):
        L=eng.ceil_log(n) if n else 0
        assert eng.scale_header([0]*n) is None
        for scale in range(1,min(1<<L,80)+1):
            code=gamma(scale);h=len(code)
            if h>n:continue
            for fill in (0,1):
                bits=code+[fill]*(n-h)
                assert eng.scale_header(bits)==(h,scale)
                gamma_cases+=1
        for zeros in range(min(L,8)+1):
            truncated=[0]*zeros+[1]+[0]*max(0,zeros-1)
            if zeros:assert eng.scale_header(truncated) is None
        if n>=2*L+1:
            too_large=gamma((1<<L)+1)
            if len(too_large)<=n:
                assert eng.scale_header(too_large+[0]*(n-len(too_large))) is None
    results['gamma_parse_append_and_rejection_cases']=gamma_cases
    variable_cursor=0
    for scale in (1,2,3,4,7,16):
        code=gamma(scale);h=len(code)
        for payload in range(8):
            toy=SimpleNamespace(header=h,prefix=h+2,data_prefix=2,shared_parameters=master,
                seeds=[NoData(),wrappers[0],wrappers[1],NoData(),wrappers[2],wrappers[3]],design=NoData())
            bits=code+eng.to_bits(0,2)+eng.to_bits(payload,4)
            assert eng.scale_header(bits)==(h,scale)
            engine=eng.PolyEngine(toy,bits)
            assert engine.x==bits[h+2:]+[0]*(h+2)
            variable_cursor+=1
    results['gamma_variable_prefix_cursor_cases']=variable_cursor
    envelope_cases=0
    for logL in (32,256,1024,4096):
        L=1<<logL;z=eng.ceil_log(L)
        for logS in (10,40,100):
            scale=1<<logS;C=1<<30
            h,e,ell,ep,q,d,Y,k=envelope(L,scale)
            old=2*eng.next_pow2(scale*L**2)*eng.next_pow2(scale**11*L**7*z**5)
            assert k==1<<(1+eng.ceil_log(h)+eng.ceil_log(Y))
            assert k-h*Y-100*5*L*h*2*ell>=10*d+C*e*z*z*ell
            assert q-2*h*ell>=C*e*z*z*ell and d>=q
            assert Y>=1000*5*L*h*d
            assert Y-100*5*L*h*d>=C*e*z*z*max(d,q)
            assert ep>=100*h*2*ell+2*e and ell>=128*z*e
            short=eng.ShortSeedPlan(20*d,10*d,ell,ep,d)
            assert short.valid and short.actual_seed<=d
            assert 2*eng.next_pow2(8*(ell+ep+L+10))<=d
            assert k<old
            if logL>=256:
                assert k*scale**8<=2**65*old
                assert ell<=512*scale*L**2*z*z
                assert ep<=2**21*scale**2*L**4*z*z
                assert q<=2**11*C*scale**2*L**4*z**5
                assert d<=2**17*C*scale**2*L**4*z**5
                assert Y<=2**34*C*scale**3*L**7*z**5
                assert k<=2**36*C*scale**4*L**9*z**5
            envelope_cases+=1
    results['explicit_resource_envelope_scalar_cases']=envelope_cases
    return {'all_passed': True, 'checks': results,
            'scope': 'Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
