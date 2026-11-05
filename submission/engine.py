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

@dataclass(frozen=True)
class Round:
    rows:int
    width:int
    slice:int
    degree:int
    near:SeedParams
    calls:int
    @property
    def seed(self):return self.calls*self.near.d
    @property
    def data(self):return self.degree+self.near.data
    def valid(self):
        return self.near.valid() and self.seed<=root_floor(self.degree**49,100)

class Plan:
    def __init__(self,N):
        self.N=N;L=ceil_log(N);self.L=L
        z=ceil_log(L)
        scale=ceil_root(L**3*z,4)
        self.entropy_log=300*scale
        # Compare exponents before allocating the value 2^entropy_log.
        if self.entropy_log>N.bit_length()-1:
            self.ready=False
            return
        self.k=1<<self.entropy_log
        self.t=1<<ceil_div(3*ceil_root(L*z**3,4),50)
        self.ell=1<<5*ceil_root(L*z,2)
        self.pmax=1<<ceil_div(27*scale,10)
        self.q=400*root_floor(L//z,2)
        logt=ceil_log(self.t);H=ceil_div(L,logt)
        self.n=self.t**H
        self.h=0;leaf=self.n
        while leaf>=root_floor(self.k**2,5):
            self.h+=1;leaf//=self.t
        self.hA=self.h//self.q
        K=1<<(ceil_log(self.k)//2);r=1<<(ceil_log(K)//2)
        D=ceil_div(100*ceil_log(self.n),ceil_log(K))
        self.E0=seeded_params(self.n,K,D,r,3)
        self.E1=seeded_params(self.n,K,D,r,root_floor(K,10))
        self.r=r;self.T=1<<((49*ceil_log(r))//100)
        logT=ceil_log(self.T)
        self.b=1+ceil_div(self.E0.d,max(1,logT))
        self.rounds=[]
        rows,width=self.T,r
        while rows>1:
            w=root_floor(width**49,100);p=2*w+1
            m0=1<<((7*ceil_log(width))//8)
            near=seeded_params(rows*width,width//4,100,m0,root_floor(width,10))
            calls=ceil_div(width//2,m0)
            self.rounds.append(Round(rows,width,w,p,near,calls))
            rows//=2;width//=2
        self.prefix=self.E0.data+self.E1.data+sum(x.data for x in self.rounds)
        # No field/matrix data are allocated during plan construction.
        self.ready=(self.T>=2 and self.t>=self.b and self.hA>=1 and self.prefix+self.k<=N and self.E0.valid() and
                    self.E1.valid() and self.E0.d<=1000000*ceil_log(self.n) and
                    self.ell*self.pmax<=r and self.E1.d<=width and
                    self.h*self.t**self.b<=self.ell and
                    self.ell**(6*self.hA+6)<=self.pmax and
                    all(x.valid() for x in self.rounds))

class Engine:
    def __init__(self,plan,bits):
        self.p=plan;prefix=Prefix(bits)
        self.e0=prefix.seeded(plan.E0);self.e1=prefix.seeded(plan.E1)
        self.rounds=[]
        for row in plan.rounds:
            Q=(1<<row.degree)|prefix.number(row.degree)
            near=prefix.seeded(row.near)
            self.rounds.append((row,Q,near))
        assert prefix.pos==plan.prefix
        self.x=list(bits[prefix.pos:])+[0]*(plan.n-(len(bits)-prefix.pos))
        self.ch={};self.child={}
    def base(self,bits,Q,length):
        p=Q.bit_length()-1;x=sum(b<<j for j,b in enumerate(bits))
        x=power(x,7,Q);out=[]
        for j in range(length):
            a=mul(1<<j,x,Q);tr=0
            for _ in range(p):tr^=a;a=mul(a,a,Q)
            out.append(tr&1)
        return out
    def sr(self,rows):
        for spec,Q,near in self.rounds:
            sliced=[x[:spec.slice] for x in rows];flat=sum(rows,[]);new=[]
            for i in range(0,len(rows),2):
                seedbits=self.base(sliced[i]+sliced[i+1],Q,spec.seed)
                out=[]
                for j in range(spec.calls):
                    sb=seedbits[j*near.p.d:(j+1)*near.p.d]
                    seed=sum(b<<a for a,b in enumerate(sb))
                    out.extend(near.value(flat,seed))
                new.append(out[:spec.width//2])
            rows=new
        return rows[0]
    def be(self,blocks):
        rows=[self.e0.value(blocks[0],seed) for seed in range(1<<self.e0.p.d)]
        for x in blocks[1:]:
            if len(rows)==1:break
            out=[]
            for i in range(0,len(rows),self.p.T):
                chunk=rows[i:i+self.p.T]
                chunk+= [chunk[-1]]*(self.p.T-len(chunk))
                sb=self.sr(chunk)[:self.e1.p.d]
                seed=sum(b<<a for a,b in enumerate(sb))
                out.append(self.e1.value(x,seed))
            rows=out
        assert len(rows)==1
        return rows[0][:self.p.pmax]
    def source(self,node):
        a,b=node
        return self.x[a:b]+[0]*(self.p.n-(b-a))
    def children(self,node):
        a,b=node;w=(b-a)//self.p.t
        return [(a+j*w,a+(j+1)*w) for j in range(self.p.t)]
    def cr(self,node,target,p):
        source=self.source(node);rows=self.ch[target]
        for seed in range(1<<self.e0.p.d):
            response=self.e0.value(source,seed,self.p.ell*self.p.pmax)
            if all(response[j*self.p.pmax:j*self.p.pmax+p]==
                   (rows[j][:p] if j<len(rows) else [0]*p)
                   for j in range(self.p.ell)):
                return 0 # fail
        return 1 # pass
    def challenges(self):
        # Bottom-up enumeration avoids recursion at unbounded tree depth.
        for depth in range(self.p.h,-1,-1):
            count=self.p.t**depth
            width=self.p.n//count
            for index in range(count):
                node=(index*width,(index+1)*width)
                if depth==self.p.h:
                    self.ch[node]=[]
                    continue
                children=self.children(node)
                added=[self.be([self.source(children[j]) for j in choice])
                       for choice in combinations(range(self.p.t),self.p.b)]
                chosen=next((v for v in children
                             if self.cr(node,v,self.p.pmax)),children[0])
                self.child[node]=chosen
                self.ch[node]=self.ch[chosen]+added
    def evaluate(self):
        root=(0,self.p.n);self.challenges()
        path=[root]
        for _ in range(self.p.h):path.append(self.child[path[-1]])
        for u in range(self.p.hA):
            node=path[u*self.p.q];target=path[(u+1)*self.p.q]
            p1=self.p.ell**(6*(self.p.hA-u)+3)
            if self.cr(node,target,p1)==0:
                return self.cr(node,target,self.p.ell**(6*(self.p.hA-u)+6))
        return 0


# Supplied polylogarithmic branch; older routines remain finite-audit fixtures.

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
        self.s=next_pow2(4*self.a)
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
    def __init__(self,n,K,M,E,seed_length=None,dimension=10000):
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

def weak_design(a,count,universe=None):
    universe=100*a*a if universe is None else universe
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
        assert 2*sum(1<<(selected&prior).bit_count() for prior in masks)<=3*i
        masks.append(selected)
        row=[];bits=selected
        while bits:
            bit=bits&-bits;row.append(bit.bit_length()-1);bits^=bit
        rows.append(tuple(row))
    return rows


class HalfSeedPlan:
    """Linear half-output extraction by condensation and a weak-design hybrid."""
    def __init__(self,n,K,M,E,seed_length):
        self.n,self.K,self.M,self.error=n,K,M,E
        self.valid=False;self.data=0
        if not(0<K<=n and 0<M<=K//2):return
        self.pv=PVPlan(n,K,E+10)
        self.budget=E+ceil_log(M)+10
        self.nw_s=next_pow2(8*(E+ceil_log(self.pv.output)+ceil_log(M)+10))
        self.nw_a=2*self.nw_s;self.universe=100*self.nw_a**2
        self.actual_seed=self.pv.s+self.universe
        self.seed_length=seed_length
        self.valid=(K>=32*self.budget and self.actual_seed<=seed_length)
    def bind(self,prefix):return HalfSeed(self,self.pv.bind())


class HalfSeed:
    def __init__(self,p,pv):
        self.p,self.pv=p,pv
        self.rows=weak_design(p.nw_a,p.M,p.universe)
    def one_bit(self,x,seed):
        p=self.p;s=p.nw_s;mask=(1<<s)-1
        y,z=seed&mask,(seed>>s)&mask
        value=0
        for i in range(ceil_div(p.pv.output,s)-1,-1,-1):
            value=field_mul(value,y,s)^((x>>(i*s))&mask)
        return field_mul(z,value,s).bit_count()&1
    def value(self,bits,seed):
        p=self.p
        assert len(bits)<=p.n
        bits=list(bits)+[0]*(p.n-len(bits))
        y=seed&((1<<p.pv.s)-1);seed>>=p.pv.s
        x=to_int(self.pv.value(bits,y))
        return [self.one_bit(x,sum(((seed>>j)&1)<<k for k,j in enumerate(row)))
                for row in self.rows]


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
        if not(0<2*m<=r and self.code_c<m):return
        self.code_a=2*self.code_s
        self.universe=100*self.code_a**2
        self.actual_seed=2*self.code_s+self.universe
        self.data=self.code_s
        self.valid=self.actual_seed<=b
    def bind(self,prefix):
        Q=(1<<self.code_s)|prefix.number(self.code_s)
        return Design(self,Q)


class Design:
    def __init__(self,p,Q):
        self.p,self.Q=p,Q
        self.rows=weak_design(p.code_a,p.m-p.code_c,p.universe)
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
        for i,row in enumerate(self.rows):
            local=sum(((seed>>j)&1)<<k for k,j in enumerate(row))
            out|=(self.word(v,local).bit_count()&1)<<(p.code_c+i)
        return out


def layout(d0,den=10**12,scale=1):
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
        if N.bit_length()<1024:return
        L=ceil_log(N);z=ceil_log(L);self.header=L;self.scale=scale
        if not 1<=scale<=1<<L:return
        self.entropy_log=ceil_log(scale**13*L**11*z**10)
        if self.entropy_log>N.bit_length()-1:return
        k=1<<self.entropy_log;self.k=k
        self.alphabet=2
        self.h=next_pow2(scale*L**2*z**2)
        self.ell=next_pow2(scale**3*L**4*z**4)
        self.d=next_pow2(scale**9*L**6*z**6)
        self.Y=next_pow2(scale**11*L**9*z**8)
        h,ell,d=self.h,self.ell,self.d
        q0=(2*h*self.alphabet+10)*ell
        if d>self.Y or q0>d:return
        self.E=self.h;self.Eprime=next_pow2(scale**5*L**6*z**6)
        self.lext=AffineSeedPlan(N,k,self.Y,4*L+10)
        if not self.lext.valid:return
        self.d0=self.lext.seed_length
        if self.d0>1000000000*L:return
        self.blocks=ceil_div(self.d0,ceil_log(self.alphabet))
        b=self.blocks
        self.q=HalfSeedPlan(d,10*ell,ell,self.E,ell)
        self.w=HalfSeedPlan(N,10*ell,ell,self.E,ell)
        self.one=ShortSeedPlan(N,10*d,ell,self.Eprime,d)
        self.two=HalfSeedPlan(self.Y,2*q0,q0,self.E,ell)
        self.three=HalfSeedPlan(self.Y,2*d,d,self.E,ell)
        self.seeds=[self.lext,self.q,self.w,self.one,self.two,self.three]
        if not all(x.valid for x in self.seeds):return
        self.r,self.m,self.seed_b,self.u=layout(self.d0,scale=scale)
        self.design=DesignPlan(self.r,self.m,self.seed_b)
        if not self.design.valid:return
        self.data_prefix=sum(x.data for x in self.seeds)+self.design.data
        self.prefix=self.header+self.data_prefix
        self.ready=(self.prefix+k<=N and
            k-h*self.Y-100*b*h*self.alphabet*ell>=10*d and
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

