"""Complete fixed total scale-header language; small lengths output zero."""

def scale_header(bits):
    n = len(bits)
    length = (n-1).bit_length() if n else 0
    return length, 1 + sum(bit << j for j, bit in enumerate(bits[:length]))

def evaluate(bits):
    """One fixed total language, with every construction datum charged."""
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("the input alphabet is {0,1}")
    if len(bits).bit_length() < 8:
        return 0
    L, scale = scale_header(bits)
    plan = PolyPlan(len(bits), scale)
    if not plan.ready:
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
    base=(1<<20)*R*h*((1<<errlog)*8*m)**4
    s=ceil_log(base**8)
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
                        if e:a=mul(a,power(x,e,Q),Q)
                    value^=a
            out.append((value&z).bit_count()%2)
        return out
    def dot(self,a,b):
        z=0
        for x,y in zip(a,b):z^=mul(x,y,self.Q)
        return z

class Prefix:
    def __init__(self,bits):self.bits=bits;self.pos=0
    def number(self,n):
        z=sum(b<<j for j,b in enumerate(self.bits[self.pos:self.pos+n]))
        self.pos+=n;return z
    def seeded(self,p):
        Q=(1<<p.s)|self.number(p.s)
        A=[[self.number(p.s) for _ in range(p.D)] for _ in range(p.D)]
        return Seeded(p,Q,A)

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
        self.a=ceil_log(n)+ceil_log(K)+E+10
        self.s=next_pow2(2*self.a)
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
    def __init__(self,n,K,M,E,seed_length=None,dimension=11):
        self.n,self.K,self.M,self.error=n,K,M,E
        self.valid=False;self.data=0
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
        self.data=su.data
        self.valid=self.actual_seed<=self.seed_length
    def bind(self,prefix):return AffineSeed(self,self.pv.bind(),prefix.seeded(self.su))


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


class HalfSeedPlan:
    """Linear half-output extraction by condensation and a weak-design hybrid."""
    def __init__(self,n,K,M,E,seed_length,overlap_log=1):
        self.n,self.K,self.M,self.error=n,K,M,E
        self.valid=False;self.data=0
        if not isinstance(overlap_log,int) or overlap_log<1:return
        self.overlap_log=overlap_log;self.overlap_weight=1<<overlap_log
        if not(0<K<=n and 0<M and 4*self.overlap_weight*M<=K):return
        self.pv=PVPlan(n,K,E+10)
        self.budget=E+ceil_log(M)+10
        self.block=self.budget;self.row_count=ceil_div(M,self.block)
        self.nw_s=next_pow2(8*(E+ceil_log(self.pv.output)+ceil_log(M)+10))
        self.nw_a=2*self.nw_s
        self.universe=ceil_div(2*self.nw_a**2,self.overlap_log)+self.nw_a
        self.actual_seed=self.pv.s+self.universe
        self.seed_length=seed_length
        self.valid=(K>=32*self.budget and self.actual_seed<=seed_length)
    def bind(self,prefix):return HalfSeed(self,self.pv.bind())


class HalfSeed:
    def __init__(self,p,pv):
        self.p,self.pv=p,pv
        self.rows=weak_design(p.nw_a,ceil_div(p.M,p.block),p.universe,p.overlap_weight,1)
    def word(self,x,seed):
        p=self.p;s=p.nw_s;mask=(1<<s)-1
        y,z=seed&mask,(seed>>s)&mask
        value=0
        for i in range(ceil_div(p.pv.output,s)-1,-1,-1):
            value=field_mul(value,y,s)^((x>>(i*s))&mask)
        return field_mul(z,value,s)
    def value(self,bits,seed):
        p=self.p
        assert len(bits)<=p.n
        bits=list(bits)+[0]*(p.n-len(bits))
        y=seed&((1<<p.pv.s)-1);seed>>=p.pv.s
        x=to_int(self.pv.value(bits,y))
        out=[]
        for row in self.rows:
            local=sum(((seed>>j)&1)<<k for k,j in enumerate(row))
            count=min(p.block,p.M-len(out))
            out.extend(to_bits(self.word(x,local),count))
        return out


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


class DesignPlan:
    def __init__(self,r,m,b):
        self.r,self.m,self.b=r,m,b;self.valid=False;self.data=0
        self.code_s=1000*ceil_log(r);self.code_c=ceil_div(self.code_s,4)
        self.code_block=100*ceil_log(r)
        if not(0<2*m<=r and self.code_c<m):return
        self.code_a=2*self.code_s
        self.universe=(5*self.code_a**2+1)//2+self.code_a
        self.actual_seed=2*self.code_s+self.universe
        self.data=self.code_s
        self.valid=self.actual_seed<=b
    def bind(self,prefix):
        Q=(1<<self.code_s)|prefix.number(self.code_s)
        return Design(self,Q)


class Design:
    def __init__(self,p,Q):
        self.p,self.Q=p,Q
        self.rows=weak_design(p.code_a,ceil_div(p.m-p.code_c,p.code_block),p.universe)
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
        offset=p.code_c
        for row in self.rows:
            local=sum(((seed>>j)&1)<<k for k,j in enumerate(row))
            count=min(p.code_block,p.m-offset)
            out|=(self.word(v,local)&((1<<count)-1))<<offset
            offset+=count
        return out


def layout(d0,den=1000,scale=1):
    def at(r):
        m=r//den;b=root_floor(scale,2)*ceil_log(r)**2//(30*den*den)
        u=2*b+2*ceil_log(r)+10
        return m,b,u
    lo,hi=1,(d0+1)*den
    while lo+1<hi:
        mid=(lo+hi)//2;m,b,u=at(mid)
        if ceil_log(u)+m+b<=d0:lo=mid
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
        L=ceil_log(N);z=ceil_log(L);self.header=L;self.scale=scale
        if not 1<=scale<=1<<L:return
        self.entropy_log=ceil_log(scale**13*L**11*z**9)
        if self.entropy_log>N.bit_length()-1:return
        k=1<<self.entropy_log;self.k=k
        self.alphabet=2
        self.h=next_pow2(scale*L**2*z**2)
        self.ell=next_pow2(scale**3*L**4*z**3)
        self.d=next_pow2(scale**9*L**6*z**5)
        self.Y=next_pow2(scale**11*L**9*z**7)
        h,ell,d=self.h,self.ell,self.d
        q0=(2*h*self.alphabet+10)*ell
        if d>self.Y or q0>d:return
        self.E=self.h;self.Eprime=next_pow2(scale**5*L**6*z**5)
        self.initial_error=ceil_div(L,100)+10
        self.lext=AffineSeedPlan(N,k,self.Y,self.initial_error,dimension=11)
        if not self.lext.valid:return
        self.d0=self.lext.seed_length
        if self.d0>10*L:return
        self.blocks=ceil_div(self.d0,ceil_log(self.alphabet))
        b=self.blocks
        self.overlap_log=max(1,L.bit_length()-1)
        self.overlap_weight=1<<self.overlap_log
        rho=self.overlap_weight;t=self.overlap_log
        self.q=HalfSeedPlan(d,4*rho*ell,ell,self.E,ell,t)
        self.w=HalfSeedPlan(N,4*rho*ell,ell,self.E,ell,t)
        self.one=ShortSeedPlan(N,10*d,ell,self.Eprime,d)
        self.two=HalfSeedPlan(self.Y,4*rho*q0,q0,self.E,ell,t)
        self.three=HalfSeedPlan(self.Y,4*rho*d,d,self.E,ell,t)
        self.seeds=[self.lext,self.q,self.w,self.one,self.two,self.three]
        if not all(x.valid for x in self.seeds):return
        self.r,self.m,self.seed_b,self.u=layout(self.d0,scale=scale)
        self.design=DesignPlan(self.r,self.m,self.seed_b)
        if not self.design.valid:return
        self.data_prefix=sum(x.data for x in self.seeds)+self.design.data
        self.prefix=self.header+self.data_prefix
        self.ready=(self.prefix+k<=N and
            k-h*self.Y-100*b*h*self.alphabet*ell>=10*d+4*rho*ell and
            self.Y-100*b*h*d>=4*rho*max(d,q0) and
            q0-2*h*(self.alphabet-1)*ell>=4*rho*ell and
            self.Y>=1000*b*h*d and
            self.Eprime>=100*h*self.alphabet*ell+2*self.E and
            h>=self.r**2*self.u and self.E>=100*self.r**2*self.u)


class PolyEngine:
    def __init__(self,p,bits):
        self.p=p;prefix=Prefix(bits[p.header:])
        self.lext,self.q,self.w,self.one,self.two,self.three=[x.bind(prefix) for x in p.seeds]
        self.design=p.design.bind(prefix)
        assert prefix.pos==p.data_prefix
        self.x=list(bits[p.prefix:])+[0]*p.prefix
    def lookahead(self,Q,t):
        p=self.p;S=Q[:p.ell]
        for step in range(t):
            R=self.w.value(self.x,to_int(S))
            if step+1<t:S=self.q.value(Q,to_int(R))
        return R
    def row_bit(self,i):
        p=self.p;Y=self.lext.value(self.x,i);seed=Y[:p.d]
        step=ceil_log(p.alphabet);digits=i<<(p.blocks*step-p.d0)
        for j in range(p.blocks):
            index=(digits>>((p.blocks-1-j)*step))&(p.alphabet-1)
            R0=self.one.value(self.x,to_int(seed))
            Q=self.two.value(Y,to_int(R0))
            R=self.lookahead(Q,index+1)
            if j+1<p.blocks:seed=self.three.value(Y,to_int(R))
            else:return R[0]
    def evaluate(self):
        p=self.p;B=1<<p.seed_b;M=1<<p.m;S=B*M
        threshold=biased_threshold(p.r,p.m,p.seed_b,p.u)
        groups=[]
        for i in range(S):
            word=sum(self.row_bit(i*p.u+j)<<j for j in range(p.u))
            groups.append(word>=threshold)
        for v in range(1<<p.r):
            neighbors=[self.design.value(v,z) for z in range(B)]
            passed=False
            for shift in range(M):
                if all(groups[z*M+(neighbors[z]^shift)] for z in range(B)):
                    passed=True;break
            if not passed:return 0
        return 1

