"""Exact prefix-data evaluator for the supplied subpolynomial-entropy construction."""
from dataclasses import dataclass
from itertools import islice, combinations
from math import isqrt

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

def isprime(n):
    return n>=2 and (n==2 or n%2 and all(n%d for d in range(3,isqrt(n)+1,2)))

def prime_at_least(n):
    while not isprime(n):n+=1
    return n

@dataclass(frozen=True)
class Round:
    rows:int
    width:int
    slice:int
    prime:int
    near:SeedParams
    calls:int
    @property
    def seed(self):return self.calls*self.near.d
    @property
    def data(self):return self.prime+self.near.data
    def valid(self):
        return self.near.valid() and self.seed<=root_floor(self.prime,8)

class Plan:
    def __init__(self,N):
        self.N=N;L=ceil_log(N);self.L=L
        self.k=1<<ceil_root(L**9,10)
        self.t=1<<ceil_root(L**3,10)
        self.ell=1<<ceil_root(L,2)
        self.pmax=1<<ceil_root(L**4,5)
        self.q=ceil_root(L,2)
        logt=ceil_log(self.t);H=ceil_div(L,logt)
        self.n=self.t**H
        self.h=0;leaf=self.n
        while leaf>=root_floor(self.k**2,5):
            self.h+=1;leaf//=self.t
        self.hA=self.h//self.q
        K=1<<(ceil_log(self.k)//2);r=1<<(ceil_log(K)//2)
        D=ceil_div(100*ceil_log(self.n),ceil_log(K))
        self.E0=seeded_params(self.n,K,D,r,3)
        self.E1=seeded_params(self.n,K,D,r,root_floor(K,1000))
        self.r=r;self.T=1<<(ceil_log(r)//10)
        logT=ceil_log(self.T)
        self.b=1+ceil_div(self.E0.d,max(1,logT))
        self.rounds=[]
        rows,width=self.T,r
        while rows>1:
            w=root_floor(width,10);p=prime_at_least(2*w)
            m0=1<<((999*ceil_log(width))//1000)
            near=seeded_params(rows*width,width//4,100000,m0,root_floor(width,100000))
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
            Q=(1<<row.prime)|prefix.number(row.prime)
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
