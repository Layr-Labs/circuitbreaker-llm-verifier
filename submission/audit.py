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
        assert not bad.ready
        L=eng.ceil_log(length)
        smallest=min(eng.split_threshold_resource_candidate(L,bad.h,bad.E,bad.ell,step,eng.ceil_div(L,step))[0] for step in range(1,eng.ceil_log(bad.h)+1))
        assert eng.ceil_log(smallest)>length.bit_length()-1
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
    def round_cap(L,step):return eng.ceil_div(357*L,100*step)
    def envelope(L,scale):
        z=eng.ceil_log(L);C=1<<20
        u,v,bias=eng.unpack_scales(scale)
        h=eng.next_pow2(u*L**2);e=eng.next_pow2(v*L**2*z)
        ell=eng.next_pow2(16*z*e)
        a=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))));step=eng.ceil_log(a);cap=round_cap(L,step)
        ep=100*h*a*ell+e+2
        q=(2*h*(a-1)+C*e*z*z)*ell
        d=max(q,(2*h*a+10)*ell,2*eng.next_pow2(8*(ell+ep+L+10)))
        Y=100*cap*h*d+C*e*z*z*d
        return h,e,ell,ep,q,d,Y,h*Y+100*cap*h*a*ell+10*d+C*e*z*z*ell
    L=1<<4096;z=eng.ceil_log(L);scale=1<<100;C=1<<20;dimension=max(2,z)
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
    assert q0-2*h*(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))-1)*ell>=C*error*z*z*ell and d>=q0
    assert Y-100*round_cap(L,eng.ceil_log(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))))*h*d>=C*error*z*z*max(d,q0)
    assert k-h*Y-100*round_cap(L,eng.ceil_log(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))))*h*eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))*ell>=10*d+C*error*z*z*ell
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
    L=1<<4096;z=eng.ceil_log(L);scale=1<<100;C=1<<20;dimension=max(2,z)
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
            scale=1<<logS;C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,scale)
            old=2*eng.next_pow2(scale*L**2)*eng.next_pow2(scale**11*L**7*z**5)
            assert k<=1<<eng.ceil_log(k) and k>1<<(eng.ceil_log(k)-1)
            assert k-h*Y-100*round_cap(L,eng.ceil_log(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))))*h*eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))*ell>=10*d+C*e*z*z*ell
            assert q-2*h*(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))-1)*ell>=C*e*z*z*ell and d>=q
            assert Y-100*round_cap(L,eng.ceil_log(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))))*h*d>=C*e*z*z*max(d,q)
            assert ep>=100*h*eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))*ell+e+2 and ell>=16*z*e
            short=eng.ShortSeedPlan(20*d,10*d,ell,ep,d)
            assert short.valid and short.actual_seed<=d
            assert 2*eng.next_pow2(8*(ell+ep+L+10))<=d
            assert k<old
            step=eng.ceil_log(eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z)))))
            old_q=(4*h+C*e*z*z+10)*ell
            old_ep=eng.next_pow2(512*h*ell+4*e)
            old_d=eng.next_pow2(16*old_q+128*old_ep+128*L+1280)
            old_Y=eng.next_pow2(10000*L*h*old_d+4*C*e*z*z*old_d)
            old_envelope_k=2*h*old_Y
            assert k<=old_envelope_k
            if logL>=256:
                assert k<old_envelope_k
                assert k*step<=2**12*old_envelope_k
                assert Y*step<=2**35*C*scale**3*L**7*z**5
                assert k*step<=2**37*C*scale**4*L**9*z**5
                assert k*scale**8<=2**65*old
                assert ell<=512*scale*L**2*z*z
                assert ep<=2**25*scale**2*L**4*z**5
                assert q<=2**11*C*scale**2*L**4*z**5
                assert d<=2**17*C*scale**2*L**4*z**5
                assert Y<=2**34*C*scale**3*L**7*z**5
                assert k<=2**36*C*scale**4*L**9*z**5
            envelope_cases+=1
    results['growing_alphabet_resource_scalar_cases']=envelope_cases
    ordered_words=0
    for d0 in range(1,11):
        for alphabet in (4,8,16,64):
            step=eng.ceil_log(alphabet);blocks=eng.ceil_div(d0,step)
            previous=None
            for index in range(1<<d0):
                shifted=index<<(blocks*step-d0)
                digits=tuple((shifted>>((blocks-1-j)*step))&(alphabet-1) for j in range(blocks))
                assert all(0<=digit<alphabet for digit in digits)
                assert previous is None or previous<digits
                assert sum(digit<<(step*(blocks-1-j)) for j,digit in enumerate(digits))==shifted
                previous=digits;ordered_words+=1
    results['growing_alphabet_ordered_index_words']=ordered_words
    non_dyadic=0
    for requested in (3,5,7,9):
        plan=eng.AffineSeedPlan(1<<64,1<<50,requested,1,dimension=11)
        assert plan.valid and plan.calls==2 and plan.su.m<=requested<2*plan.su.m
        non_dyadic+=1
    results['exact_budget_non_dyadic_scalar_outputs']=non_dyadic
    exact_cases=0
    for logL in (256,1024,4096):
        for logS in (10,40,100):
            L=1<<logL;scale=1<<logS;z=eng.ceil_log(L);C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,scale)
            a=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))));cap=round_cap(L,eng.ceil_log(a))
            assert d>=q and d>=2*eng.next_pow2(8*(ell+ep+L+10))
            assert eng.next_pow2(Y)&(eng.next_pow2(Y)-1)==0
            assert k==h*Y+100*cap*h*a*ell+10*d+C*e*z*z*ell
            old_d=eng.next_pow2(d)
            old_Y=eng.next_pow2(10000*cap*h*old_d+4*C*e*z*z*old_d)
            old_k=2*h*old_Y
            assert d<=old_d and Y<=old_Y and 4*k<=3*old_k
            exact_cases+=1
    results['prefix_initial_budget_scalar_cases']=exact_cases
    initial_prefix_cases=0
    class SourcePrefix:
        def __init__(self,word):self.word=word
        def value(self,bits,seed):return self.word
    class CheckY:
        def __init__(self,expected):self.expected=expected;self.calls=0
        def value(self,bits,seed):
            assert bits==self.expected;self.calls+=1
            return [1,0]
    class First:
        def value(self,bits,seed):return [0,1]
    for word in range(16):
        full=eng.to_bits(word,4)
        for index in range(8):
            obj=eng.PolyEngine.__new__(eng.PolyEngine)
            obj.p=SimpleNamespace(Y=3,d=2,ell=2,alphabet=4,blocks=2,d0=3)
            obj.x=[0]*6;obj.lext=SourcePrefix(full);obj.one=First()
            obj.two=CheckY(full[:3]);obj.three=CheckY(full[:3])
            obj.lookahead=lambda Q,t:[t&1,(t>>1)&1]
            assert len(obj.row_word(index,2))==2
            assert obj.two.calls==2 and obj.three.calls==1
            initial_prefix_cases+=1
    for requested in (3,5,7,9):
        plan=eng.AffineSeedPlan(1<<64,1<<50,eng.next_pow2(requested),1,dimension=11)
        assert plan.valid and plan.calls==1 and plan.M>=requested
        initial_prefix_cases+=1
    results['dyadic_initial_prefix_output_cases']=initial_prefix_cases
    prefix_compare=0
    for logL in (256,1024,4096):
        for logS in (10,40,100):
            L=1<<logL;z=eng.ceil_log(L);scale=1<<logS;C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,scale)
            old_a=eng.next_pow2(z*z);old_cap=round_cap(L,eng.ceil_log(old_a))
            old_ep=eng.next_pow2(512*h*old_a*ell+4*e)
            old_q=(2*h*old_a+C*e*z*z+10)*ell
            old_d=16*old_q+128*old_ep+128*L+1280
            old_Y=eng.next_pow2(10000*old_cap*h*old_d+4*C*e*z*z*old_d)
            old_k=h*old_Y+100*old_cap*h*old_a*ell+10*old_d+C*e*z*z*ell
            assert k<old_k and 10*k<=9*old_k
            initial=eng.next_pow2(Y)
            assert Y<=initial<2*Y and initial&(initial-1)==0
            prefix_compare+=1
    results['larger_alphabet_loss_comparison_scalar_cases']=prefix_compare
    def pack_scales(row,error,bias):
        values=[row-1,error-1,bias-1];code=0
        for position in range(max(value.bit_length() for value in values)):
            for lane,value in enumerate(values):
                code|=((value>>position)&1)<<(3*position+lane)
        return code+1
    codec_cases=0
    for code in range(1,4097):
        triple=eng.unpack_scales(code)
        assert pack_scales(*triple)==code and all(1<=value<=code for value in triple)
        codec_cases+=1
    for row in range(1,17):
        for error in range(1,17):
            for bias in range(1,17):
                assert eng.unpack_scales(pack_scales(row,error,bias))==(row,error,bias)
                codec_cases+=1
    results['packed_independent_scale_roundtrip_cases']=codec_cases
    comparisons=0
    R0=1<<20;K0=1<<30
    for logL in (256,1024,4096):
        for logt in (40,60,80):
            L=1<<logL;t=1<<logt
            now=envelope(L,pack_scales(R0,K0*t*t,t*t))
            old=envelope(L,pack_scales(t*t,t**3,t*t))
            assert now[0]*t*t==old[0]*R0
            assert now[1]*t==old[1]*K0 and now[2]*t==old[2]*K0
            numerator=R0*R0*K0*K0;denominator=t**6
            assert now[-1]*denominator<=2*old[-1]*numerator
            assert 2*now[-1]*denominator>=old[-1]*numerator
            assert now[-1]<old[-1]
            comparisons+=1
    results['packed_independent_scale_resource_cases']=comparisons
    exact_y_cases=0
    for logL in (256,1024,4096):
        for logS in (40,60,100):
            L=1<<logL;code=1<<logS;z=eng.ceil_log(L);C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,code)
            alphabet=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))
            cap=round_cap(L,eng.ceil_log(alphabet))
            assert Y==100*cap*h*d+C*e*z*z*d
            for rounds in (1,max(1,cap//2),cap):
                assert Y-100*rounds*h*d>=C*e*z*z*max(d,q)
            assert Y-100*cap*h*d==C*e*z*z*d
            old_Y=10000*cap*h*d+4*C*e*z*z*d
            old_k=h*old_Y+100*cap*h*alphabet*ell+10*d+C*e*z*z*ell
            assert 50*k<=old_k and 200*k>=old_k
            exact_y_cases+=1
    results['exact_lookahead_entropy_scalar_cases']=exact_y_cases
    maxd_cases=0
    for logL in (256,1024,4096):
        for logS in (40,60,100):
            L=1<<logL;code=1<<logS;z=eng.ceil_log(L);C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,code)
            alphabet=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))
            cap=round_cap(L,eng.ceil_log(alphabet))
            seed=2*eng.next_pow2(8*(ell+ep+L+10))
            assert q-2*h*(alphabet-1)*ell==C*e*z*z*ell
            assert d==max(q,(2*h*alphabet+10)*ell,seed)
            assert d>=q and d>=(2*h*alphabet+10)*ell and d>=seed
            assert 10*d>=ell+2*ep+10 and d==q
            old_q=(2*h*alphabet+C*e*z*z+10)*ell
            old_d=16*old_q+128*ep+128*L+1280
            F=100*cap*h+C*e*z*z
            old_k=h*F*old_d+100*cap*h*alphabet*ell+10*old_d+C*e*z*z*ell
            assert d<=old_d and old_k-k==(h*F+10)*(old_d-d)>0
            assert 8*k<=old_k and 32*k>=old_k
            maxd_cases+=1
    results['exact_q_and_max_seed_budget_scalar_cases']=maxd_cases
    coefficient_cases=0
    assert (1<<20)//2-327680==196608
    for logL in (256,1024,4096):
        for logS in (40,60,100):
            L=1<<logL;z=eng.ceil_log(L);C=1<<20;old_C=1<<30
            h,e,ell,ep,q,d,Y,k=envelope(L,1<<logS)
            a=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))
            cap=round_cap(L,eng.ceil_log(a))
            old_q=(2*h*(a-1)+old_C*e*z*z)*ell
            old_d=max(old_q,(2*h*a+10)*ell,2*eng.next_pow2(8*(ell+ep+L+10)))
            old_Y=100*cap*h*old_d+old_C*e*z*z*old_d
            old_k=h*old_Y+100*cap*h*a*ell+10*old_d+old_C*e*z*z*ell
            assert d==q and old_d==old_q
            assert old_q-q==(old_C-C)*e*z*z*ell>0
            assert old_k-k==h*(old_Y-Y)+10*(old_d-d)+(old_C-C)*e*z*z*ell>0
            assert 512*k<=old_k and 2048*k>=old_k
            for ambient,output in ((d,ell),(Y,q),(Y,d)):
                plan=eng.AffineSeedPlan(ambient,C*e*z*z*output,output,e,ell,dimension=max(2,z),primitive=True)
                assert plan.valid and plan.calls<=2 and plan.actual_seed<=ell
            coefficient_cases+=1
    results['reduced_entropy_coefficient_scalar_cases']=coefficient_cases
    seed_envelope_cases=0
    for logL in (256,1024,4096):
        for logS in (40,60,100):
            L=1<<logL;z=eng.ceil_log(L);C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,1<<logS)
            a=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))
            cap=round_cap(L,eng.ceil_log(a))
            old_ell=eng.next_pow2(128*z*e)
            assert old_ell==8*ell and (10*z+15)*e<=12*z*e<ell
            old_ep=100*h*a*old_ell+e+2
            old_q=(2*h*(a-1)+C*e*z*z)*old_ell
            old_d=max(old_q,(2*h*a+10)*old_ell,2*eng.next_pow2(8*(old_ell+old_ep+L+10)))
            old_Y=100*cap*h*old_d+C*e*z*z*old_d
            old_k=h*old_Y+100*cap*h*a*old_ell+10*old_d+C*e*z*z*old_ell
            assert d==q and old_d==old_q==8*q
            assert old_Y==8*Y and old_k==8*k
            plans=[eng.AffineSeedPlan(ambient,C*e*z*z*output,output,e,ell,dimension=max(2,z),primitive=True) for ambient,output in ((d,ell),(Y,q),(Y,d))]
            degree=max(plan.su.s for plan in plans)
            for plan in plans:
                plan.use_shared_field(degree)
                assert plan.valid and plan.calls<=2 and plan.actual_seed<12*z*e<=ell
            seed_envelope_cases+=1
    results['reduced_seed_envelope_scalar_cases']=seed_envelope_cases
    initial_precision_cases=0
    from fractions import Fraction
    assert Fraction(357,100)-Fraction(646,181)==Fraction(17,18100)
    assert Fraction(1,35)-Fraction(357,12500)==Fraction(1,87500)
    for logL in (1024,2048,4096):
        for code in (1,1<<20,1<<40):
            L=35*(1<<logL)//37;z=eng.ceil_log(L);C=1<<20
            h,e,ell,ep,q,d,Y,k=envelope(L,code)
            a=eng.next_pow2(max(2,z**3//max(1,eng.ceil_log(z))))
            step=eng.ceil_log(a);cap=round_cap(L,step)
            old_cap=eng.ceil_div(5*L,step)
            old_Y=100*old_cap*h*d+C*e*z*z*d
            old_k=h*old_Y+100*old_cap*h*a*ell+10*d+C*e*z*z*ell
            assert old_k-k==100*(old_cap-cap)*h*(h*d+a*ell)>0
            assert 7*old_k<=10*k and 10*k<=8*old_k
            E0=eng.ceil_div(L,35)+10;logK=eng.ceil_log(k)
            pv_a=logK+E0+20
            pv_s=eng.next_pow2(L+2*logK+2*E0+33)
            nprime=eng.ceil_div(k,pv_a)*pv_s
            su=eng.seeded_params(nprime,k//2,11,eng.next_pow2(Y),E0+logK+10)
            assert su.valid() and eng.next_pow2(Y)<=k//4
            seed=pv_s+su.d
            assert 100*seed<357*L and eng.ceil_div(seed,step)<=cap
            initial_precision_cases+=1
    results['initial_precision_and_round_cap_scalar_cases']=initial_precision_cases
    optimized_cases=0
    for logL in (128,256,512):
        for code in (1,1<<20,1<<40):
            L=1<<logL;z=eng.ceil_log(L);U,V,T=eng.unpack_scales(code)
            h=eng.next_pow2(U*L*L);e=eng.next_pow2(V*L*L*z);ell=eng.next_pow2(16*z*e);C=1<<20
            selected=eng.optimized_alphabet_resources(L,h,e,ell)
            k,step,ep,q,d,cap,Y=selected;a=1<<step
            assert 2<=a<=h and cap==round_cap(L,step)
            assert q-2*h*(a-1)*ell==C*e*z*z*ell
            assert ep==100*h*a*ell+e+2
            assert d==max(q,(2*h*a+10)*ell,2*eng.next_pow2(8*(ell+ep+L+10)))
            assert Y==100*cap*h*d+C*e*z*z*d
            assert k==h*Y+100*cap*h*a*ell+10*d+C*e*z*z*ell
            old_h,old_e,old_ell,old_ep,old_q,old_d,old_Y,old_k=envelope(L,code)
            assert k<old_k and 2*k>old_k and 10*k<8*old_k
            plans=[eng.AffineSeedPlan(ambient,C*e*z*z*output,output,e,ell,dimension=max(2,z),primitive=True) for ambient,output in ((d,ell),(Y,q),(Y,d))]
            degree=max(plan.su.s for plan in plans)
            for plan in plans:
                plan.use_shared_field(degree)
                assert plan.valid and plan.calls<=2 and plan.actual_seed<12*z*e
            optimized_cases+=1
    results['optimized_alphabet_actual_resource_cases']=optimized_cases
    joint_cases=0
    for N in (128,256,512,1024):
        L=eng.ceil_log(N);z=eng.ceil_log(L)
        h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
        assert eng.joint_initial_resources(N,L,h,e,ell) is None
        assert not eng.PolyPlan(N,1).ready
        joint_cases+=1
    results['joint_initial_finite_rejection_cases']=joint_cases
    feasible_cases=0
    assert Fraction(337,100)-Fraction(10346,3071)==Fraction(327,307100)
    assert Fraction(1,37)-Fraction(337,12500)==Fraction(31,462500)
    for logL in (1024,2048,4096):
        for code in (1,1<<20,1<<40):
            L=37*(1<<logL)//39;z=eng.ceil_log(L);U,V,T=eng.unpack_scales(code)
            h=eng.next_pow2(U*L*L);e=eng.next_pow2(V*L*L*z);ell=eng.next_pow2(16*z*e)
            parent=eng.optimized_alphabet_resources(L,h,e,ell)
            old_k,step,_,_,_,old_cap,_=parent
            cap=eng.ceil_div(337*L,100*step)
            k,_,ep,q,d,_,Y=eng.alphabet_resource_candidate(L,h,e,ell,step,cap,coefficient=1<<20)
            assert cap<=eng.ceil_div(5*L,step) and k<old_k
            E0=eng.ceil_div(L,37)+10;logK=eng.ceil_log(k)
            pv_a=logK+E0+20;pv_s=eng.next_pow2(L+2*logK+2*E0+33)
            nprime=eng.ceil_div(k,pv_a)*pv_s
            su=eng.seeded_params(nprime,k//2,10,eng.next_pow2(Y),E0+logK+10)
            assert su.valid() and eng.next_pow2(Y)<=k//4
            seed=pv_s+su.d
            assert 100*seed<337*L and eng.ceil_div(seed,step)<=cap
            assert 125*E0-seed>=125*eng.ceil_root(L,2)
            assert E0<=L and 100*k<95*old_k
            feasible_cases+=1
    results['joint_initial_feasible_scalar_witness_cases']=feasible_cases
    prefix_cases=0
    Q,matrix=eng.canonical_cyclic_matrix(2,2)
    for keep in (3,5):
        params=eng.SeedParams(3,3,2,2,eng.next_pow2(keep),1,2)
        full=eng.Seeded(params,Q,matrix);view=eng.OutputPrefix(full,keep)
        for seed in range(64):
            basis=[eng.to_int(view.value(eng.to_bits(1<<i,3),seed)) for i in range(3)]
            for x in range(8):
                projected=eng.to_int(view.value(eng.to_bits(x,3),seed))
                assert len(view.value(eng.to_bits(x,3),seed))==keep
                assert projected==eng.to_int(full.value(eng.to_bits(x,3),seed))&((1<<keep)-1)
                assert projected==__import__('functools').reduce(int.__xor__,(basis[i] for i in range(3) if (x>>i)&1),0)
                prefix_cases+=1
    results['dyadic_prefix_value_and_linearity_cases']=prefix_cases
    reduced_full_cases=0
    assert 5*72000-359996==4 and 144000*25==150000*24
    # Current neighbor/initial check; it does not claim internal validity at this logL.
    L=1<<4096;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));j=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*j);E0=eng.ceil_div(L,37)+10
    old=eng.full_request_resource_candidate(L,h,e,ell,j,cap,coefficient=150000)
    def seed_metadata(source,output):
        logK=eng.ceil_log(source);pv_a=logK+E0+20
        pv_s=eng.next_pow2(L+2*logK+2*E0+33)
        pn=eng.ceil_div(source,pv_a)*pv_s
        su=eng.seeded_params(pn,source//2,10,eng.next_pow2(output),E0+logK+10)
        assert su.valid() and eng.next_pow2(output)<=source//4
        return pv_s+su.d
    old_seed=seed_metadata(old[0],old[6]);newcap=eng.ceil_div(old_seed,j-1)
    k,step,ep,q,d,_,Y=eng.full_request_resource_candidate(L,h,e,ell,j-1,newcap,coefficient=144000)
    assert 3*q>=old[3] and 25*q<=24*old[3] and 3*d>=old[4] and 25*d<=24*old[4]
    assert 40*k>=old[0] and 50*k<49*old[0] and Y<old[6]
    seed=seed_metadata(k,Y)
    assert seed<=old_seed<=5*L and eng.ceil_div(seed,step)<=newcap
    assert eng.ceil_div(L,step)<=newcap<=eng.ceil_div(5*L,step)
    assert 125*E0-seed>=125*eng.ceil_root(L,2)
    assert Y-100*newcap*h*d==144000*e*z*z*eng.next_pow2(d)
    reduced_full_cases+=1
    results['reduced_full_request_neighbor_initial_case']=reduced_full_cases
    # Current actual shared-port check at a larger scalar length; no initial D10 here.
    L=1<<8192;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));step=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*step)
    k,_,ep,q,d,_,Y=eng.full_request_resource_candidate(L,h,e,ell,step,cap,coefficient=144000);C=144000
    plans=[eng.DyadicPrefixPlan(n,C*e*z*z*eng.next_pow2(m),m,e,ell,dimension=z,primitive=True) for n,m in ((d,ell),(Y,q),(Y,d))]
    Kw=C*e*z*z*ell;logK=eng.ceil_log(Kw);pv=eng.next_pow2(L+2*logK+2*e+33)
    pn=eng.ceil_div(Kw,logK+e+20)*pv;wsu=eng.seeded_params(pn,Kw//2,z,ell,e+logK+10)
    degree=max([wsu.s]+[p.su.s for p in plans]);wsu=eng.SeedParams(wsu.n,wsu.K,wsu.D,wsu.h,wsu.m,wsu.errlog,degree)
    assert wsu.valid() and pv+wsu.d<6*z*e<=ell
    for p in plans:
        p.use_shared_field(degree)
        assert p.valid and p.calls==1 and p.M==eng.next_pow2(p.prefix_output)
        assert p.su.h<=65*z and 10*(16*z+16)<=161*z and 10*degree<43*e
        assert p.actual_seed<6*z*e<=ell
    results['reduced_full_request_actual_shared_port_case']=1
    seed_budget_cases=0
    assert 10**29>2**96 and 73000*144==144000*73
    # Current neighbor/initial check; it does not claim internal validity at this logL.
    L=1<<4096;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));j=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*j);E0=eng.ceil_div(L,37)+10
    old=eng.full_request_resource_candidate(L,h,e,ell,j,cap,coefficient=144000)
    def seed_metadata(source,output):
        logK=eng.ceil_log(source);pv_a=logK+E0+20
        pv_s=eng.next_pow2(L+2*logK+2*E0+33)
        pn=eng.ceil_div(source,pv_a)*pv_s
        su=eng.seeded_params(pn,source//2,10,eng.next_pow2(output),E0+logK+10)
        assert su.valid() and eng.next_pow2(output)<=source//4
        return pv_s+su.d
    old_seed=seed_metadata(old[0],old[6]);newcap=eng.ceil_div(old_seed,j-1)
    k,step,ep,q,d,_,Y=eng.full_request_resource_candidate(L,h,e,ell,j-1,newcap)
    assert 3*q>=old[3] and 144*q<=73*old[3] and 3*d>=old[4] and 144*d<=73*old[4]
    assert 40*k>=old[0] and 4*k<3*old[0] and Y<old[6]
    seed=seed_metadata(k,Y)
    assert seed<=old_seed<=5*L and eng.ceil_div(seed,step)<=newcap
    assert eng.ceil_div(L,step)<=newcap<=eng.ceil_div(5*L,step)
    assert 125*E0-seed>=125*eng.ceil_root(L,2)
    assert Y-100*newcap*h*d==73000*e*z*z*eng.next_pow2(d)
    seed_budget_cases+=1
    results['seed_budget_dimension_neighbor_initial_case']=seed_budget_cases
    # Current actual shared-port check at a larger scalar length; no initial D10 here.
    L=1<<8192;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));step=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*step)
    k,_,ep,q,d,_,Y=eng.full_request_resource_candidate(L,h,e,ell,step,cap);C=73000
    plans=[eng.DyadicPrefixPlan(n,C*e*z*z*eng.next_pow2(m),m,e,ell,dimension=eng.ceil_div(29*z,16),primitive=True) for n,m in ((d,ell),(Y,q),(Y,d))]
    Kw=C*e*z*z*ell;logK=eng.ceil_log(Kw);pv=eng.next_pow2(L+2*logK+2*e+33)
    pn=eng.ceil_div(Kw,logK+e+20)*pv;wsu=eng.seeded_params(pn,Kw//2,eng.ceil_div(29*z,16),ell,e+logK+10)
    degree=max([wsu.s]+[p.su.s for p in plans]);wsu=eng.SeedParams(wsu.n,wsu.K,wsu.D,wsu.h,wsu.m,wsu.errlog,degree)
    assert wsu.valid() and 10*(pv+wsu.d)<79*z*e and pv+wsu.d<=ell
    for p in plans:
        p.use_shared_field(degree)
        assert p.valid and p.calls==1 and p.M==eng.next_pow2(p.prefix_output)
        assert p.su.D==eng.ceil_div(29*z,16) and p.su.h<=10*p.su.D and 100*degree<427*e
        assert 10*p.actual_seed<79*z*e and p.actual_seed<=ell
    results['seed_budget_dimension_actual_shared_port_case']=1
    split_cases=0
    assert 5**29>2**64 and 36500*2==73000
    # Current neighbor/initial check; it does not claim internal validity at this logL.
    L=1<<4096;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));j=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*j);E0=eng.ceil_div(L,37)+10
    old=eng.full_request_resource_candidate(L,h,e,ell,j,cap,coefficient=73000)
    def seed_metadata(source,output):
        logK=eng.ceil_log(source);pv_a=logK+E0+20
        pv_s=eng.next_pow2(L+2*logK+2*E0+33)
        pn=eng.ceil_div(source,pv_a)*pv_s
        su=eng.seeded_params(pn,source//2,10,eng.next_pow2(output),E0+logK+10)
        assert su.valid() and eng.next_pow2(output)<=source//4
        return pv_s+su.d
    old_seed=seed_metadata(old[0],old[6]);newcap=eng.ceil_div(old_seed,j-1)
    k,step,ep,q,d,_,Y=eng.split_threshold_resource_candidate(L,h,e,ell,j-1,newcap,lookahead_coefficient=36500,row_coefficient=73000,row_loss=100)
    assert 3*q>=old[3] and 2*q<=old[3] and 3*d>=old[4] and 2*d<=old[4]
    assert 40*k>=old[0] and 4*k<3*old[0] and Y<old[6]
    seed=seed_metadata(k,Y)
    assert seed<=old_seed<=5*L and eng.ceil_div(seed,step)<=newcap
    assert eng.ceil_div(L,step)<=newcap<=eng.ceil_div(5*L,step)
    assert 125*E0-seed>=125*eng.ceil_root(L,2)
    assert Y-100*newcap*h*d==73000*e*z*z*eng.next_pow2(d)
    assert k-h*Y-100*newcap*h*(a//2)*ell==10*d+36500*e*z*z*ell
    split_cases+=1
    results['split_threshold_neighbor_initial_case']=split_cases
    lookahead_dimension_cases=0
    assert 7**10>2**28 and 19**10>2**42 and 31500*73==36500*63
    # Current neighbor/initial check; it does not claim internal validity at this logL.
    L=1<<4096;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));j=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*j);E0=eng.ceil_div(L,37)+10
    old=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,lookahead_coefficient=36500,row_coefficient=73000,row_loss=100)
    def seed_metadata(source,output):
        logK=eng.ceil_log(source);pv_a=logK+E0+20
        pv_s=eng.next_pow2(L+2*logK+2*E0+33)
        pn=eng.ceil_div(source,pv_a)*pv_s
        su=eng.seeded_params(pn,source//2,10,eng.next_pow2(output),E0+logK+10)
        assert su.valid() and eng.next_pow2(output)<=source//4
        return pv_s+su.d
    old_seed=seed_metadata(old[0],old[6]);newcap=eng.ceil_div(old_seed,j-1)
    k,step,ep,q,d,_,Y=eng.split_threshold_resource_candidate(L,h,e,ell,j-1,newcap,row_loss=100)
    assert 3*q>=old[3] and 73*q<=63*old[3] and 3*d>=old[4] and 73*d<=63*old[4]
    assert 40*k>=old[0] and 73*k<68*old[0] and Y<old[6]
    seed=seed_metadata(k,Y)
    assert seed<=old_seed<=5*L and eng.ceil_div(seed,step)<=newcap
    assert eng.ceil_div(L,step)<=newcap<=eng.ceil_div(5*L,step)
    assert 125*E0-seed>=125*eng.ceil_root(L,2)
    assert Y-100*newcap*h*d==86000*e*z*z*eng.next_pow2(d)
    assert k-h*Y-100*newcap*h*(a//2)*ell==10*d+31500*e*z*z*ell
    lookahead_dimension_cases+=1
    results['lookahead_dimension_neighbor_initial_case']=lookahead_dimension_cases
    # Current same-alphabet Y-budget witness; exact initial comparison only.
    parent=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,row_loss=100)
    current=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,row_loss=32)
    k,step,ep,q,d,_,Y=current
    assert current[1:6]==parent[1:6] and 25*k>=8*parent[0] and 2*k<parent[0]
    assert 25*Y>=8*parent[6] and Y<parent[6]
    ps=seed_metadata(parent[0],parent[6]);ns=seed_metadata(k,Y)
    assert ns<=ps<=5*L and eng.ceil_div(ns,step)<=cap
    assert 125*E0-ns>=125*eng.ceil_root(L,2)
    assert Y-32*cap*h*d==86000*e*z*z*eng.next_pow2(d)
    assert k-h*Y-100*cap*h*(1<<step)*ell==10*d+31500*e*z*z*ell
    # Explicit ledger bound, including equality in the d guard.
    for hh in (2,4,16,128):
        for aa in (2,hh):
            ll=16;dd=(2*hh*aa+10)*ll
            assert 4*hh*dd+32*hh*hh*aa*ll+8*hh*ll<=24*hh*dd<32*hh*dd
    results['round_budget_same_alphabet_initial_and_ledger_case']=1
    # Current uninflated ledger: compare with actual32hd parent at the same a,B,E0.
    parent=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,row_loss=32)
    current=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,row_loss=8)
    k,step,ep,q,d,_,Y=current
    assert current[1:6]==parent[1:6] and 4*k>=parent[0] and 2*k<parent[0]
    assert 4*Y>=parent[6] and Y<parent[6]
    assert seed_metadata(k,Y)<=seed_metadata(parent[0],parent[6])<=5*L
    assert eng.ceil_div(seed_metadata(k,Y),step)<=cap
    assert 125*E0-seed_metadata(k,Y)>=125*eng.ceil_root(L,2)
    assert Y-8*cap*h*d==86000*e*z*z*eng.next_pow2(d)
    for hh in (2,4,16,128):
        for aa in (2,hh):
            ll=16;dd=(2*hh*aa+10)*ll
            assert 4*hh*dd+8*hh*aa*ll+8*hh*ll<=6*hh*dd<8*hh*dd
    results['uninflated_row_same_alphabet_initial_and_ledger_case']=1
    # Current three-list witness, preserving the original8-budget comparison above.
    parent=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap,row_loss=8)
    k,step,ep,q,d,_,Y=eng.split_threshold_resource_candidate(L,h,e,ell,j,cap)
    assert 2*k>=parent[0] and 4*k<3*parent[0] and 2*Y>=parent[6] and Y<parent[6]
    ns=seed_metadata(k,Y);ps=seed_metadata(parent[0],parent[6])
    assert ns<=ps<=5*L and eng.ceil_div(ns,step)<=cap
    assert 125*E0-ns>=125*eng.ceil_root(L,2)
    assert Y-4*cap*h*d==86000*e*z*z*eng.next_pow2(d)
    for hh in (4,8,16,128):
        for aa in (2,hh):
            ll=16;dd=(2*hh*aa+10)*ll
            assert 3*hh*dd+8*hh*aa*ll+8*hh*ll<=4*hh*dd
    results['three_list_same_alphabet_initial_and_ledger_case']=1
    # Current actual shared-port check at a larger scalar length; no initial D10 here.
    L=1<<16384;z=eng.ceil_log(L);h=eng.next_pow2(L*L);e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
    a=eng.next_pow2(z**3//eng.ceil_log(z));step=eng.ceil_log(a);cap=eng.ceil_div(337*L,100*step)
    k,_,ep,q,d,_,Y=eng.split_threshold_resource_candidate(L,h,e,ell,step,cap);CW,CY=31500,86000
    plans=[eng.DyadicPrefixPlan(n,C*e*z*z*eng.next_pow2(m),m,e,ell,dimension=eng.ceil_div(10*z,7),primitive=True) for n,C,m in ((d,CW,ell),(Y,CY,q),(Y,CY,d))]
    Kw=CW*e*z*z*ell;logK=eng.ceil_log(Kw);pv=eng.next_pow2(L+2*logK+2*e+33)
    pn=eng.ceil_div(Kw,logK+e+20)*pv;wsu=eng.seeded_params(pn,Kw//2,eng.ceil_div(10*z,7),ell,e+logK+10)
    degree=max([wsu.s]+[p.su.s for p in plans]);wsu=eng.SeedParams(wsu.n,wsu.K,wsu.D,wsu.h,wsu.m,wsu.errlog,degree)
    assert wsu.valid() and pv+wsu.d<7*z*e and pv+wsu.d<=ell
    for p in plans:
        p.use_shared_field(degree)
        assert p.valid and p.calls==1 and p.M==eng.next_pow2(p.prefix_output)
        assert p.su.D==eng.ceil_div(10*z,7) and p.su.h<=19*p.su.D and 1000*degree<4267*e
        assert p.actual_seed<7*z*e and p.actual_seed<=ell
    assert plans[0].su.h<=7*plans[0].su.D and wsu.h<=7*wsu.D
    assert Y-4*cap*h*d==CY*e*z*z*eng.next_pow2(d)
    assert q-2*h*(a-1)*ell==CW*e*z*z*ell
    assert k-h*Y-100*cap*h*a*ell==10*d+CW*e*z*z*ell
    results['three_list_actual_shared_port_case']=1
    return {'all_passed': True, 'checks': results,
            'scope': 'C23 dyadic predicates and older fixed/optimized-alphabet predicates remain construction-class checks; historical C24 tight-coefficient predicates;historical C25 full-request predicates;historical C26 reduced-full-request predicates;historical C27 dimension predicates;historical C28 split-threshold neighbor;retired old C28 actual shared-port repetition;historical C29 dimension initial-neighbor;historical32hd same-alphabet initial comparison;historical8hd same-alphabet initial comparison;current4hd same-alphabet initial comparison,three-list ledger arithmetic and actual shared-port predicates check actual full-output shared plans; current negative joint search and dyadic-prefix ports have scoped rejection,linear-output and exact positive scalar plan predicates without ambient allocation;parent feasibility comparisons explicitly retain their old coefficient. Exact finite audit only. Universal coverage, uniformity, and asymptotics are proved in the supplied Markdown files.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared audit takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
