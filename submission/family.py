"""Complete fixed total scale-header language; small lengths output zero."""

def scale_header(bits):
    """Bounded standard Elias gamma code for positive S; invalid returns None."""
    n=len(bits);L=(n-1).bit_length() if n else 0
    zeros=0
    while zeros<n and bits[zeros]==0:
        zeros+=1
        if zeros>L:return None
    length=2*zeros+1
    if zeros==n or length>n:return None
    scale=1
    for bit in bits[zeros+1:length]:scale=(scale<<1)|bit
    if scale>1<<L:return None
    return length,scale

def evaluate(bits):
    """One fixed total language, with every construction datum charged."""
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("the input alphabet is {0,1}")
    if len(bits).bit_length() < 8:
        return 0
    parsed = scale_header(bits)
    if parsed is None:return 0
    header, scale = parsed
    plan = PolyPlan(len(bits), scale)
    if not plan.ready or plan.header!=header:
        return 0
    return PolyEngine(plan, bits).evaluate()

# Complete large-length evaluator, kept in this same module.
"""Exact prefix-data evaluator for the supplied subpolynomial-entropy construction."""
from dataclasses import dataclass
from math import comb
from itertools import islice, combinations

def root_floor(a,d):
    if a<2:return a
    lo,hi=0,1<<((a.bit_length()+d-1)//d)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**d<=a:lo=m
        else:hi=m
    return hi if hi**d<=a else lo

def ceil_root(a,d):
    b=root_floor(a,d)
    return b+(b**d<a)

def ceil_log(a):return (a-1).bit_length()
def ceil_div(a,b):return (a+b-1)//b

def mul(a,b,Q):
    s=Q.bit_length()-1;hi=1<<s;out=0
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a&hi:a^=Q
    return out

def power(a,e,Q):
    z=1
    while e:
        if e&1:z=mul(z,a,Q)
        e>>=1
        if e:a=mul(a,a,Q)
    return z

@dataclass(frozen=True)
class SeedParams:
    n:int
    K:int
    D:int
    h:int
    m:int
    errlog:int
    s:int
    @property
    def d(self):return (self.D+1)*self.s
    @property
    def data(self):return (self.D*self.D+1)*self.s
    def valid(self):
        R=16*self.D+16;rho_log=self.errlog+3+ceil_log(self.m)
        return (self.m>0 and self.m&(self.m-1)==0 and
                self.K>=16*self.m*R*self.h*self.s+4*rho_log+10)

def seeded_params(n,K,D,m,errlog):
    h=D*ceil_root(n,D);R=16*D+16
    # rho^-1=2^errlog*8*m, with m a power of two.
    base_log=32+4*errlog+4*ceil_log(m)+ceil_log(R*h)
    s=ceil_div(R*base_log,R-D)
    return SeedParams(n,K,D,h,m,errlog,s)

def compositions(total,D):
    # Iterative lexicographic enumeration; D can exceed recursion limits.
    exps=[0]*D
    exps[-1]=total
    while True:
        yield tuple(exps)
        available=exps[-1]
        for i in range(D-2,-1,-1):
            if available:
                exps[i]+=1
                exps[i+1:]=[0]*(D-i-1)
                exps[-1]=available-1
                break
            available+=exps[i]
        else:
            return


def monomials(n,D,h):
    remain=n
    for degree in range(h):
        for exps in compositions(degree,D):
            yield exps;remain-=1
            if not remain:return
    if remain:raise ValueError('insufficient monomials')

class Seeded:
    def __init__(self,p,Q,A):self.p,self.Q,self.A=p,Q,A
    def value(self,bits,seed,count=None):
        p=self.p;Q=self.Q;mask=(1<<p.s)-1
        v=[(seed>>(j*p.s))&mask for j in range(p.D)]
        z=(seed>>(p.D*p.s))&mask
        count=p.m if count is None else count
        out=[]
        for _ in range(count):
            v=[self.dot(row,v) for row in self.A]
            value=0
            for bit,exps in zip(bits,monomials(p.n,p.D,p.h)):
                if bit:
                    a=1
                    for x,e in zip(v,exps):
                        if e:a=self.multiply(a,self.exponentiate(x,e))
                    value^=a
            out.append((value&z).bit_count()%2)
        return out
    def dot(self,a,b):
        z=0
        for x,y in zip(a,b):z^=self.multiply(x,y)
        return z

    def multiply(self,a,b):return mul(a,b,self.Q)
    def exponentiate(self,a,e):return power(a,e,self.Q)

class Prefix:
    def __init__(self,bits):self.bits=bits;self.pos=0
    def number(self,n):
        z=sum(b<<j for j,b in enumerate(self.bits[self.pos:self.pos+n]))
        self.pos+=n;return z
    def seeded(self,p):
        Q=(1<<p.s)|self.number(p.s)
        A=[[self.number(p.s) for _ in range(p.D)] for _ in range(p.D)]
        return Seeded(p,Q,A)

    def companion_seeded(self,p):
        Q=(1<<p.s)|self.number(p.s)
        coefficients=[self.number(p.s) for _ in range(p.D)]
        return Seeded(p,Q,companion_matrix(coefficients))

    def primitive_seeded(self,p):
        degree=p.D*p.s
        R=decode_internal_polynomial(degree,self.number(degree-primitive_package_saving(degree)))
        return ZeroSeeded(p) if R is None else primitive_seeded(p,R)

def companion_matrix(coefficients):
    D=len(coefficients);matrix=[[0]*D for _ in range(D)]
    for j in range(D):
        if j:matrix[j][j-1]=1
        matrix[j][-1]=coefficients[j]
    return matrix

# Deterministic field and extraction components for the large-length branch.

from dataclasses import dataclass
from math import comb


def next_pow2(n):return 1<<ceil_log(n)


def to_bits(x,n):return [(x>>j)&1 for j in range(n)]


def to_int(bits):return sum(x<<j for j,x in enumerate(bits))


from itertools import combinations


def field_mul(a,b,m):
    # Explicit postorder evaluation avoids an unbounded Python call stack.
    key=lambda x,y,d:(min(x,y),max(x,y),d)
    target=key(a,b,m);memo={};stack=[(target,0)]
    while stack:
        t,stage=stack.pop()
        if t in memo:continue
        x,y,d=t
        if x==0 or y==0:memo[t]=0;continue
        if d==1:memo[t]=x&y;continue
        h=d//2;mask=(1<<h)-1
        u,v=x&mask,x>>h;r,s=y&mask,y>>h
        ur,us,vr,vs=[key(c,e,h) for c,e in [(u,r),(u,s),(v,r),(v,s)]]
        if stage==0:
            stack.append((t,1))
            for z in [ur,us,vr,vs]:
                if z not in memo:stack.append((z,0))
        elif stage==1:
            avs=key(1<<(h-1),memo[vs],h)
            stack.append((t,2))
            if avs not in memo:stack.append((avs,0))
        else:
            avs=key(1<<(h-1),memo[vs],h)
            memo[t]=(memo[ur]^memo[avs])|((memo[us]^memo[vr]^memo[vs])<<h)
    return memo[target]


def field_pow(a,e,m):
    z=1
    while e:
        if e&1:z=field_mul(z,a,m)
        a=field_mul(a,a,m);e>>=1
    return z


def relative_modulus(s,D):
    assert s>0 and s&(s-1)==0 and D>=2 and D&(D-1)==0
    m=s*D;theta=1<<(m//2);a=theta;coeff=[1]
    for _ in range(D):
        out=[0]*(len(coeff)+1)
        for j,c in enumerate(coeff):out[j]^=field_mul(c,a,m);out[j+1]^=c
        coeff=out;a=field_pow(a,1<<s,m)
    assert a==theta and all(c<(1<<s) for c in coeff)
    return coeff


def quotient_mul(a,b,E,s):
    D=len(E)-1;out=[0]*(2*D-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]^=field_mul(c,d,s)
    for i in range(2*D-2,D-1,-1):
        c=out[i]
        for j,d in enumerate(E[:-1]):out[i-D+j]^=field_mul(c,d,s)
    return out[:D]


def pv_value(x,y,E,s,a,m):
    D=len(E)-1;mask=(1<<s)-1
    f=[(x>>(s*i))&mask for i in range(D)];out=0
    for j in range(m):
        v=0
        for c in reversed(f):v=field_mul(v,y,s)^c
        out|=v<<(j*s)
        for _ in range(a):f=quotient_mul(f,f,E,s)
    return out


class PVPlan:
    def __init__(self,n,K,E):
        self.n,self.K,self.error=n,K,E
        self.a=ceil_log(K)+E+10
        self.s=next_pow2(ceil_log(n)+ceil_log(K)+self.a+E+3)
        self.D=next_pow2(max(2,ceil_div(n,self.s)))
        self.m=ceil_div(K,self.a)
        self.output=self.m*self.s
    def bind(self):return PV(self)


class PV:
    def __init__(self,p):self.p=p;self.E=None
    def value(self,bits,seed):
        p=self.p
        if self.E is None:self.E=relative_modulus(p.s,p.D)
        x=to_int(bits)
        return to_bits(pv_value(x,seed,self.E,p.s,p.a,p.m),p.output)


class AffineSeedPlan:
    def __init__(self,n,K,M,E,seed_length=None,dimension=11,canonical=False,companion=False,primitive=False):
        self.n,self.K,self.M,self.error=n,K,M,E
        self.valid=False;self.data=0;self.canonical=canonical;self.companion=companion;self.primitive=primitive;self.shared=False
        if not(0<K<=n and 0<M<=K//2):return
        self.pv=PVPlan(n,K,E+10)
        candidate=1<<(min(M,K//4).bit_length()-1)
        while candidate:
            su=seeded_params(self.pv.output,K//2,dimension,candidate,E+ceil_log(K)+10)
            if su.valid():break
            candidate//=2
        if not candidate:return
        self.su=su;self.calls=ceil_div(M,candidate)
        self.actual_seed=self.pv.s+self.calls*su.d
        self.seed_length=self.actual_seed if seed_length is None else seed_length
        self.data=0 if canonical else (su.D*su.s-primitive_package_saving(su.D*su.s) if primitive else ((su.D+1)*su.s if companion else su.data))
        self.valid=self.actual_seed<=self.seed_length
    def use_shared_field(self,degree):
        old=self.su
        if degree<old.s:
            self.valid=False;return
        self.su=SeedParams(old.n,old.K,old.D,old.h,old.m,old.errlog,degree)
        self.actual_seed=self.pv.s+self.calls*self.su.d
        self.valid=self.su.valid() and self.actual_seed<=self.seed_length
        self.shared=True;self.data=0
    def bind(self,prefix,shared=None):
        if self.shared:
            su=reparameterize_shared(shared,self.su)
        elif self.canonical:su=canonical_seeded(self.su)
        elif self.primitive:su=prefix.primitive_seeded(self.su)
        elif self.companion:su=prefix.companion_seeded(self.su)
        else:su=prefix.seeded(self.su)
        return AffineSeed(self,self.pv.bind(),su)


class AffineSeed:
    def __init__(self,p,pv,su):self.p,self.pv,self.su=p,pv,su
    def value(self,bits,seed):
        p=self.p;assert len(bits)<=p.n
        bits=list(bits)+[0]*(p.n-len(bits))
        y=seed&((1<<p.pv.s)-1);seed>>=p.pv.s
        condensed=self.pv.value(bits,y);out=[]
        mask=(1<<p.su.d)-1
        for _ in range(p.calls):
            out.extend(self.su.value(condensed,seed&mask));seed>>=p.su.d
        return out[:p.M]


def weak_design_expectation(masks,selected,remaining,k):
    t=remaining.bit_count();num=0
    for previous in masks:
        h=(selected&previous).bit_count();l=(remaining&previous).bit_count()
        num+=(1<<h)*sum(comb(l,j)*comb(t-j,k-j) for j in range(min(k,l)+1))
    return num,comb(t,k)

def weak_design(a,count,universe=None,weight_num=3,weight_den=2):
    universe=((5*a*a+1)//2+a) if universe is None else universe
    masks=[];rows=[]
    for i in range(count):
        selected=0;remaining=(1<<universe)-1;k=a
        for pos in range(universe):
            t=universe-pos
            if k==0:break
            if k==t:selected|=remaining;break
            bit=1<<pos;remaining^=bit
            inc=weak_design_expectation(masks,selected|bit,remaining,k-1)
            exc=weak_design_expectation(masks,selected,remaining,k)
            if inc[0]*exc[1]<=exc[0]*inc[1]:selected|=bit;k-=1
        assert selected.bit_count()==a
        assert weight_den*sum(1<<(selected&prior).bit_count() for prior in masks)<=weight_num*i
        masks.append(selected)
        row=[];bits=selected
        while bits:
            bit=bits&-bits;row.append(bit.bit_length()-1);bits^=bit
        rows.append(tuple(row))
    return rows


class ShortSeedPlan:
    """Strong short-output polynomial hash, with deterministic field data."""
    def __init__(self,n,K,M,E,seed_length):
        self.n,self.K,self.M,self.error=n,K,M,E
        self.valid=False;self.data=0
        if not(0<K<=n and 0<M<=K//2):return
        self.hash_s=next_pow2(8*(M+E+ceil_log(n)+10))
        self.actual_seed=2*self.hash_s
        self.seed_length=seed_length
        self.valid=(K>=M+2*E+10 and self.actual_seed<=seed_length)
    def bind(self,prefix):return ShortSeed(self)


class ShortSeed:
    def __init__(self,p):self.p=p
    def value(self,bits,seed):
        p=self.p;assert len(bits)<=p.n
        x=to_int(bits);s=p.hash_s;mask=(1<<s)-1
        y,z=seed&mask,(seed>>s)&mask;value=0
        for i in range(ceil_div(p.n,s)-1,-1,-1):
            value=field_mul(value,y,s)^((x>>(i*s))&mask)
        value=field_mul(z,value,s)
        return to_bits(value,p.M)


class WeakHashPlan:
    """General-source NW hash; its derived seed excludes the zero-length prefix."""
    def __init__(self,n,K,M,outer_r):
        self.r=n;self.K=K;self.m=M;self.valid=False;self.data=0
        w=ceil_log(outer_r)
        if w<1 or not(0<K<=n and 0<M<=n and n<=outer_r):return
        self.code_s=1000*w;self.code_c=0;self.code_block=750*w
        self.code_a=self.code_s;self.universe=5*self.code_s**2+self.code_s
        self.actual_seed=self.code_s+self.universe
        self.valid=(3*(K-150*w-self.code_s)>=5*M)


class EntropySamplerPlan:
    """A fixed 1024-wise polynomial sampler; all bounds use scalar integers."""
    def __init__(self,r):
        self.valid=False;self.w=ceil_log(r);self.k=1024;self.t=256;self.b=17
        if self.w<1:return
        self.n=1<<self.w;self.length=self.n//(64*2000*self.b)
        self.entropy=self.length//16000;self.seed_length=self.k*self.w
        if self.length<self.t or not self.entropy:return
        # The tail bound is compared logarithmically; never allocate2**n.
        if self.n<16000*152*self.w:return
        t,b,L=self.t,self.b,self.length
        if 4*(128*t*b*b*2000**2)**t*r**150>L**t:return
        moment=(2*t)*(2*t)**(4*t)*(1<<(2*t))
        if 4*moment*self.n**t*(128*b*2000)**(2*t)*r**150>(7*L)**(2*t):return
        self.valid=True


class EntropySampler:
    def __init__(self,p,Q):self.p,self.Q=p,Q
    def value(self,v,seed):
        p=self.p;mask=(1<<p.w)-1
        coefficients=[(seed>>(i*p.w))&mask for i in range(p.k)]
        out=[]
        for argument in range(p.length):
            position=0
            for coefficient in reversed(coefficients):
                position=mul(position,argument,self.Q)^coefficient
            out.append((v>>position)&1)
        return out


class RecycledHash:
    def __init__(self,p,Q):
        self.p,self.Q=p,Q
        self.rows=weak_design(p.code_a,ceil_div(p.m-p.code_c,p.code_block),p.universe,5,4)
    def word(self,v,seed):
        p=self.p;mask=(1<<p.code_s)-1
        y=seed&mask;z=(seed>>p.code_s)&mask
        value=0
        for i in range(ceil_div(p.r,p.code_s)-1,-1,-1):
            value=mul(value,y,self.Q)^((v>>(i*p.code_s))&mask)
        return mul(z,value,self.Q)
    def value(self,v,seed):
        p=self.p
        out=self.word(v,seed)&((1<<p.code_c)-1)
        seed>>=2*p.code_s
        shared_z=seed&((1<<p.code_s)-1);seed>>=p.code_s
        offset=p.code_c
        for row in self.rows:
            local=sum(((seed>>j)&1)<<k for k,j in enumerate(row))
            count=min(p.code_block,p.m-offset)
            out|=(self.word(v,local|(shared_z<<p.code_s))&((1<<count)-1))<<offset
            offset+=count
        return out



class DesignPlan:
    """Strong general-source bootstrap; only scalar planning happens here."""
    def __init__(self,r,m,b):
        self.r,self.m,self.b=r,m,b;self.valid=False;self.data=0
        w=ceil_log(r);self.w=w;self.code_s=1000*w;self.code_c=250*w
        if not(0<2*m<=r and self.code_c<m):return
        self.sampler=EntropySamplerPlan(r)
        if not self.sampler.valid:return
        ell,kappa=self.sampler.length,self.sampler.entropy
        K0=m//2
        if 1000*K0<self.sampler.n:return
        if 2000*(K0-ell-150*w)<self.sampler.n:return
        if 4*K0+600*w>3*m or m<4*(self.code_c+150*w):return
        output=kappa//8
        if not output:return
        self.weak=WeakHashPlan(ell,kappa,output,r)
        if not self.weak.valid:return
        self.helper=seeded_params(ell,kappa,4,next_pow2(self.weak.actual_seed),160*w+10)
        if not self.helper.valid() or self.helper.m<self.weak.actual_seed:return
        self.rounds=ceil_div(m-self.code_c,output)
        self.round_seed=2*self.sampler.seed_length+self.helper.d
        self.actual_seed=2*self.code_s+self.rounds*self.round_seed
        self.data=0
        self.valid=self.actual_seed<=b
    def bind(self,prefix):
        Q=canonical_modulus(self.code_s)
        sampler_Q=canonical_modulus(self.sampler.w)
        helper=canonical_seeded(self.helper)
        return Design(self,Q,sampler_Q,helper)


class Design:
    word=RecycledHash.word
    def __init__(self,p,Q,sampler_Q,helper):
        self.p,self.Q,self.helper=p,Q,helper
        self.sampler=EntropySampler(p.sampler,sampler_Q)
        self.weak=RecycledHash(p.weak,Q)
    def value(self,v,seed):
        p=self.p
        out=self.word(v,seed)&((1<<p.code_c)-1)
        seed>>=2*p.code_s;offset=p.code_c
        sample_mask=(1<<p.sampler.seed_length)-1;helper_mask=(1<<p.helper.d)-1
        for _ in range(p.rounds):
            first=seed&sample_mask;seed>>=p.sampler.seed_length
            second=seed&sample_mask;seed>>=p.sampler.seed_length
            helper_seed=seed&helper_mask;seed>>=p.helper.d
            Y1=self.sampler.value(v,first);Y2=self.sampler.value(v,second)
            internal=to_int(self.helper.value(Y2,helper_seed,count=p.weak.actual_seed))
            block=self.weak.value(to_int(Y1),internal<<(2*p.weak.code_s))
            count=min(p.weak.m,p.m-offset)
            out|=(block&((1<<count)-1))<<offset;offset+=count
        assert offset==p.m
        return out


def layout(d0,den=125,scale=1):
    def at(r):
        m=r//den;b=root_floor(scale,2)*ceil_log(r)//(30*den*den)
        u=2*b+2*ceil_log(r)+10
        return m,b,u
    lo,hi=1,(d0+1)*den
    while lo+1<hi:
        mid=(lo+hi)//2;m,b,u=at(mid)
        if m+b<=d0:lo=mid
        else:hi=mid
    m,b,u=at(lo)
    return lo,m,b,u


def biased_threshold(r,m,b,u):
    # Exact comparison; integer bit lengths are polynomial in group input size.
    R=1<<r;B=1<<b;M=1<<m;den=1<<u
    denB=1<<(u*B);denBM=1<<(u*B*M)
    lo,hi=0,den
    while lo+1<hi:
        j=(lo+hi)//2
        if R*pow(denB-pow(den-j,B),M)<=denBM:lo=j
        else:hi=j
    return lo


class PolyPlan:
    def __init__(self,N,scale=1):
        self.N=N;self.ready=False
        if N.bit_length()<8:return
        L=ceil_log(N);z=ceil_log(L);self.header=2*(scale.bit_length()-1)+1;self.scale=scale
        if not 1<=scale<=1<<L:return
        self.entropy_unit=z*z
        A=self.entropy_unit
        internal_dimension=max(2,z)
        entropy_coefficient=1<<30
        self.alphabet=next_pow2(max(2,z**3//max(1,ceil_log(z))))
        self.index_step=ceil_log(self.alphabet)
        self.block_cap=ceil_div(5*L,self.index_step)
        self.h=next_pow2(scale*L**2)
        if self.alphabet>self.h:return
        self.E=next_pow2(scale*L**2*z)
        self.ell=next_pow2(128*z*self.E)
        self.Eprime=next_pow2(512*self.h*self.alphabet*self.ell+4*self.E)
        h,ell=self.h,self.ell
        q0=(2*h*self.alphabet+entropy_coefficient*self.E*A+10)*ell
        self.d=16*q0+128*self.Eprime+128*L+1280
        d=self.d
        self.Y=10000*self.block_cap*h*d+4*entropy_coefficient*self.E*A*d
        self.initial_output=next_pow2(self.Y)
        k=h*self.Y+100*self.block_cap*h*self.alphabet*ell+10*d+entropy_coefficient*self.E*A*ell
        self.entropy_log=ceil_log(k)
        if self.entropy_log>N.bit_length()-1:return
        self.k=k
        if d>self.Y or q0>d:return
        self.initial_error=ceil_div(L,28)+10
        self.lext=AffineSeedPlan(N,k,self.initial_output,self.initial_error,dimension=11,canonical=True)
        if not self.lext.valid:return
        self.d0=self.lext.seed_length
        if self.d0>5*L:return
        self.blocks=ceil_div(self.d0,ceil_log(self.alphabet))
        b=self.blocks
        self.overlap_log=max(1,L.bit_length()-1)
        self.overlap_weight=1<<self.overlap_log
        rho=self.overlap_weight;t=self.overlap_log
        self.q=AffineSeedPlan(d,entropy_coefficient*self.E*A*ell,ell,self.E,ell,dimension=internal_dimension,primitive=True)
        self.w=AffineSeedPlan(N,entropy_coefficient*self.E*A*ell,ell,self.E,ell,dimension=internal_dimension,primitive=True)
        self.one=ShortSeedPlan(N,10*d,ell,self.Eprime,d)
        self.two=AffineSeedPlan(self.Y,entropy_coefficient*self.E*A*q0,q0,self.E,ell,dimension=internal_dimension,primitive=True)
        self.three=AffineSeedPlan(self.Y,entropy_coefficient*self.E*A*d,d,self.E,ell,dimension=internal_dimension,primitive=True)
        self.seeds=[self.lext,self.q,self.w,self.one,self.two,self.three]
        if not all(x.valid for x in self.seeds):return
        internal=(self.q,self.w,self.two,self.three)
        self.individual_internal_data=sum(plan.data for plan in internal)
        common_degree=max(plan.su.s for plan in internal)
        for plan in internal:plan.use_shared_field(common_degree)
        if not all(plan.valid for plan in internal):return
        self.shared_parameters=self.q.su
        shared_t=self.shared_parameters.D*common_degree
        self.shared_data=shared_t-primitive_package_saving(shared_t)
        self.r,self.m,self.seed_b,self.u=layout(self.d0,scale=scale)
        self.design=DesignPlan(self.r,self.m,self.seed_b)
        if not self.design.valid:return
        self.data_prefix=self.shared_data+sum(x.data for x in self.seeds)+self.design.data
        self.prefix=self.header+self.data_prefix
        self.ready=(self.prefix+k<=N and
            k-h*self.Y-100*b*h*self.alphabet*ell>=10*d+entropy_coefficient*self.E*A*ell and
            self.Y-100*b*h*d>=entropy_coefficient*self.E*A*max(d,q0) and
            q0-2*h*(self.alphabet-1)*ell>=entropy_coefficient*self.E*A*ell and
            self.Y>=1000*b*h*d and
            self.Eprime>=100*h*self.alphabet*ell+2*self.E and
            ell>=self.u and h>=self.r**2 and self.E>=100*self.r**2*self.u)


class PolyEngine:
    def __init__(self,p,bits):
        self.p=p;prefix=Prefix(bits[p.header:])
        shared=prefix.primitive_seeded(p.shared_parameters) if hasattr(p,'shared_parameters') else None
        self.lext,self.q,self.w,self.one,self.two,self.three=[x.bind(prefix,shared) if getattr(x,'shared',False) else x.bind(prefix) for x in p.seeds]
        self.design=p.design.bind(prefix)
        assert prefix.pos==p.data_prefix
        self.x=list(bits[p.prefix:])+[0]*p.prefix
    def lookahead(self,Q,t):
        p=self.p;S=Q[:p.ell]
        for step in range(t):
            R=self.w.value(self.x,to_int(S))
            if step+1<t:S=self.q.value(Q,to_int(R))
        return R
    def row_word(self,i,count):
        p=self.p;assert 0<count<=p.ell
        Y=self.lext.value(self.x,i)[:p.Y];seed=Y[:p.d]
        step=ceil_log(p.alphabet);digits=i<<(p.blocks*step-p.d0)
        for j in range(p.blocks):
            index=(digits>>((p.blocks-1-j)*step))&(p.alphabet-1)
            R0=self.one.value(self.x,to_int(seed))
            Q=self.two.value(Y,to_int(R0))
            R=self.lookahead(Q,index+1)
            if j+1<p.blocks:seed=self.three.value(Y,to_int(R))
            else:return R[:count]
    def row_bit(self,i):return self.row_word(i,1)[0]
    def evaluate(self):
        p=self.p;B=1<<p.seed_b;M=1<<p.m;S=B*M
        threshold=biased_threshold(p.r,p.m,p.seed_b,p.u)
        groups=[]
        for i in range(S):
            word=to_int(self.row_word(i,p.u))
            groups.append(word>=threshold)
        for v in range(1<<p.r):
            neighbors=[self.design.value(v,z) for z in range(B)]
            passed=False
            for shift in range(M):
                if all(groups[z*M+(neighbors[z]^shift)] for z in range(B)):
                    passed=True;break
            if not passed:return 0
        return 1



# Canonical searches are used only at fixed dimension and bounded degree.
from functools import lru_cache

def polynomial_remainder(a,b):
    while a and a.bit_length()>=b.bit_length():
        a^=b<<(a.bit_length()-b.bit_length())
    return a

@lru_cache(maxsize=16)
def canonical_modulus(s):
    for low in range(1<<s):
        Q=(1<<s)|low
        reducible=False
        for degree in range(1,s//2+1):
            for lower in range(1<<degree):
                if polynomial_remainder(Q,(1<<degree)|lower)==0:
                    reducible=True;break
            if reducible:break
        if not reducible:return Q
    raise AssertionError('finite-field existence')

@lru_cache(maxsize=8)
def canonical_cyclic_matrix(s,D):
    Q=canonical_modulus(s);q=1<<s;period=(1<<(s*D))-1
    start=[1]+[0]*(D-1);mask=q-1
    for number in range(1<<(s*D)):
        coefficients=[(number>>(s*j))&mask for j in range(D)]
        # Companion action avoids building a full matrix for each trial.
        v=start[:]
        for step in range(1,period+1):
            top=v[-1]
            v=[mul(coefficients[j],top,Q)^(v[j-1] if j else 0) for j in range(D)]
            if not any(v):break
            if v==start:break
        if step==period and v==start:
            matrix=[[0]*D for _ in range(D)]
            for j in range(D):
                if j:matrix[j][j-1]=1
                matrix[j][-1]=coefficients[j]
            return Q,tuple(tuple(row) for row in matrix)
    raise AssertionError('primitive-polynomial existence')

def canonical_seeded(parameters):
    Q,A=canonical_cyclic_matrix(parameters.s,parameters.D)
    return Seeded(parameters,Q,A)


class BinaryColumns:
    """Exact XOR column elimination, retaining original coefficient vectors."""
    def __init__(self,columns,ambient):
        self.basis=[None]*ambient;self.kernel=[];self.rank=0
        for index,column in enumerate(columns):
            value,coefficient=column,1<<index
            while value:
                pivot=value.bit_length()-1;entry=self.basis[pivot]
                if entry is None:
                    self.basis[pivot]=(value,coefficient);self.rank+=1;break
                value^=entry[0];coefficient^=entry[1]
            if not value:self.kernel.append(coefficient)
    def solve(self,value):
        coefficient=0
        while value:
            pivot=value.bit_length()-1
            if pivot>=len(self.basis) or self.basis[pivot] is None:return None
            entry=self.basis[pivot];value^=entry[0];coefficient^=entry[1]
        return coefficient

class Subfield:
    def __init__(self,R,s,vectors):
        self.R,self.s,self.vectors=R,s,tuple(vectors)
        self.decoder=BinaryColumns(self.vectors,R.bit_length()-1)
    def encode(self,value):
        result=0
        for index,vector in enumerate(self.vectors):
            if (value>>index)&1:result^=vector
        return result
    def multiply(self,a,b):
        value=mul(self.encode(a),self.encode(b),self.R)
        result=self.decoder.solve(value)
        return 0 if result is None else result

class RelativeSeeded(Seeded):
    def __init__(self,p,R,A,field):
        super().__init__(p,R,A);self.field=field
    def multiply(self,a,b):return self.field.multiply(a,b)
    def exponentiate(self,a,e):
        result=1
        while e:
            if e&1:result=self.multiply(result,a)
            e>>=1
            if e:a=self.multiply(a,a)
        return result

class ZeroSeeded:
    def __init__(self,p):self.p=p
    def value(self,bits,seed,count=None):return [0]*(self.p.m if count is None else count)

def primitive_seeded(p,R):
    s,D=p.s,p.D;t=s*D
    columns=[]
    for index in range(t):
        original=1<<index;value=original
        for _ in range(s):value=mul(value,value,R)
        columns.append(value^original)
    fixed=BinaryColumns(columns,t)
    if len(fixed.kernel)!=s:return ZeroSeeded(p)
    # Coordinate 1 is the field identity, as required by monomial evaluation.
    vectors=[1];span=BinaryColumns(vectors,t)
    for vector in sorted(fixed.kernel):
        if span.solve(vector) is None:
            vectors.append(vector);span=BinaryColumns(vectors,t)
    if len(vectors)!=s:return ZeroSeeded(p)
    field=Subfield(R,s,vectors)
    alpha=polynomial_remainder(2,R);alpha_power=1;module=[]
    for _ in range(D):
        module.extend(mul(alpha_power,vector,R) for vector in vectors)
        alpha_power=mul(alpha_power,alpha,R)
    coordinates=BinaryColumns(module,t)
    if coordinates.rank!=t:return ZeroSeeded(p)
    expansion=coordinates.solve(alpha_power)
    if expansion is None:return ZeroSeeded(p)
    coefficients=[(expansion>>(j*s))&((1<<s)-1) for j in range(D)]
    return RelativeSeeded(p,R,companion_matrix(coefficients),field)


def reciprocal_package_polynomial(degree,index):
    """Decode non-palindromic reciprocal representatives without searching."""
    if degree==2:return 7 if index==0 else None
    if degree<2 or index<0:return None
    length=degree-1
    for unequal in range(length//2):
        inner_length=length-2*unequal-2
        block_size=1<<(length-unequal-2)
        if index>=block_size:
            index-=block_size;continue
        inner=index&((1<<inner_length)-1)
        equal=index>>inner_length
        middle=(inner<<(unequal+1))|(1<<unequal)
        for position in range(unequal):
            if (equal>>position)&1:
                middle|=(1<<position)|(1<<(length-1-position))
        return (1<<degree)|1|(middle<<1)
    return None


def binary_polynomial_product(a,b):
    result=0
    while b:
        if b&1:result^=a
        b>>=1;a<<=1
    return result

def binary_polynomial_divrem(a,b):
    quotient=0
    while a and a.bit_length()>=b.bit_length():
        shift=a.bit_length()-b.bit_length()
        quotient^=1<<shift;a^=b<<shift
    return quotient,a

def primitive_package_saving(degree):
    if degree<32:return min(2,degree)
    cutoff=degree.bit_length()-3
    return (cutoff+1).bit_length()-1

@lru_cache(maxsize=4)
def small_factor_code(cutoff):
    factors=[]
    for degree in range(1,cutoff+1):
        for low in range(1<<degree):
            candidate=(1<<degree)|low
            if all(polynomial_remainder(candidate,factor) for factor in factors if factor.bit_length()-1<=degree//2):
                factors.append(candidate)
    modulus=1
    for factor in factors:modulus=binary_polynomial_product(modulus,factor)
    components=[]
    for factor in factors:
        quotient,remainder=binary_polynomial_divrem(modulus,factor)
        assert remainder==0
        degree=factor.bit_length()-1
        inverse=power(polynomial_remainder(quotient,factor),(1<<degree)-2,factor)
        components.append((factor,(1<<degree)-1,quotient,inverse))
    return modulus,tuple(components)

def decode_internal_polynomial(degree,code):
    saving=primitive_package_saving(degree)
    if code<0 or code.bit_length()>degree-saving:return None
    if degree==1:return 3 if code==0 else None
    if degree<32:return reciprocal_package_polynomial(degree,code)
    modulus,components=small_factor_code(degree.bit_length()-3)
    d=modulus.bit_length()-1;quotient_bits=degree-d
    quotient=(1<<quotient_bits)|(code&((1<<quotient_bits)-1))
    index=code>>quotient_bits;residues=[]
    for factor,radix,partial,inverse in components:
        residues.append(index%radix+1);index//=radix
    if index:return None
    remainder=0
    for residue,(factor,radix,partial,inverse) in zip(residues,components):
        coefficient=mul(residue,inverse,factor)
        remainder^=binary_polynomial_product(partial,coefficient)
    return binary_polynomial_product(quotient,modulus)^remainder


def reparameterize_shared(template,parameters):
    if not isinstance(template,RelativeSeeded):return ZeroSeeded(parameters)
    if template.p.s!=parameters.s or template.p.D!=parameters.D:return ZeroSeeded(parameters)
    return RelativeSeeded(parameters,template.Q,template.A,template.field)
