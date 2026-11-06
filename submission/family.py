"""Total tree disperser evaluator. Proofs in the supplied evidence."""
from math import isqrt, comb
from functools import lru_cache
from itertools import combinations

def is_prime(n):
    return n>=2 and (n==2 or (n%2 and all(n%d for d in range(3,isqrt(n)+1,2))))

@lru_cache(None)
def prime_at_least(n):
    n=max(3,n)
    if n%2==0:n+=1
    while not is_prime(n):n+=2
    return n

def ceil_log2(n):
    return (n-1).bit_length() if n>1 else 0

def root(n,d):
    if n<0 or d<1:raise ValueError("nonnegative radicand, positive root")
    if n<2:return n
    lo=0;hi=1<<((n.bit_length()+d-1)//d)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**d<=n:lo=mid
        else:hi=mid
    return lo

def ceil_root(n,d):
    a=root(n,d)
    return a+(a**d<n)

def pack(bits):
    return sum(b<<j for j,b in enumerate(bits))

def power(a,e,f):
    v=1
    while e:
        if e&1:v=mul(v,a,f)
        a=mul(a,a,f);e>>=1
    return v

def trace(a,f):
    v=0
    for _ in range(f.bit_length()-1):
        v^=a;a=mul(a,a,f)
    assert v in (0,1)
    return v

# Deterministic prime-degree field construction proved in cubic.md Appendix A.
def mod(a,b):
    while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a

def gcd(a,b):
    while b:a,b=b,mod(a,b)
    return a

def quotient(a,b):
    q=0
    while a and a.bit_length()>=b.bit_length():
        s=a.bit_length()-b.bit_length();q^=1<<s;a^=b<<s
    assert not a
    return q

def kernel(rows,n):
    rows=list(rows);pivots=[];rank=0
    for col in range(n):
        pivot=next((i for i in range(rank,n) if (rows[i]>>col)&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(n):
            if i!=rank and (rows[i]>>col)&1:rows[i]^=rows[rank]
        pivots.append(col);rank+=1
    basis=[]
    for free in range(n):
        if free in pivots:continue
        v=1<<free
        for row,pivot in enumerate(pivots):
            if (rows[row]>>free)&1:v^=1<<pivot
        basis.append(v)
    return basis

def factor_one(f):
    while True:
        d=f.bit_length()-1
        cols=[mod(1<<(2*j),f)^(1<<j) for j in range(d)]
        rows=[sum(((cols[j]>>i)&1)<<j for j in range(d)) for i in range(d)]
        basis=kernel(rows,d)
        if len(basis)==1:return f
        b=next(b for b in basis if b not in (0,1))
        g=gcd(f,b);h=quotient(f,g)
        assert g not in (1,f) and h not in (1,f)
        f=min(g,h)

def mul(a,b,f):
    out=0;degree=f.bit_length()-1
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if (a>>degree)&1:a^=f
    return out

@lru_cache(None)
def irreducible_prime(k):
    assert k>=3 and all(k%d for d in range(2,isqrt(k)+1))
    residue=2%k;m=1
    while residue!=1:residue=(2*residue)%k;m+=1
    order=(1<<m)-1;valuation=0
    while order%k==0:order//=k;valuation+=1
    f=factor_one((1<<k)-1)
    assert f.bit_length()-1==m
    for _ in range(1,valuation):
        lifted=sum(((f>>j)&1)<<(j*k) for j in range(m+1))
        f=factor_one(lifted)
        assert f.bit_length()-1==m
    a=2
    def kmul(x,y):return mul(x,y,f)
    zero=(0,)*k;one=(1,)+(0,)*(k-1)
    beta=(0,1)+(0,)*(k-2)
    def add(x,y):return tuple(a^b for a,b in zip(x,y))
    def square(x):
        out=[0]*k
        for j,c in enumerate(x):
            if not c:continue
            c=kmul(c,c)
            if 2*j>=k:c=kmul(c,a)
            out[(2*j)%k]^=c
        return tuple(out)
    def emul(x,y):
        tmp=[0]*(2*k-1)
        for i,c in enumerate(x):
            if c:
                for j,d in enumerate(y):
                    if d:tmp[i+j]^=kmul(c,d)
        for j in range(2*k-2,k-1,-1):
            if tmp[j]:tmp[j-k]^=kmul(tmp[j],a)
        return tuple(tmp[:k])
    current=beta;gamma=zero
    for _ in range(m):
        gamma=add(gamma,current)
        for _ in range(k):current=square(current)
    assert current==beta and any(gamma[1:])
    poly=[one];root=gamma
    for _ in range(k):
        out=[zero]*(len(poly)+1)
        for j,c in enumerate(poly):
            out[j+1]=add(out[j+1],c)
            out[j]=add(out[j],emul(c,root))
        poly=out;root=square(root)
    assert root==gamma
    assert all(not any(c[1:]) and c[0] in (0,1) for c in poly)
    return sum(c[0]<<j for j,c in enumerate(poly))


@lru_cache(None)
def code_parameters(n,m,z):
    # epsilon=2^-z; rho=epsilon/(2m).
    if n<1 or m<1 or z<0:return None
    bound=1024*(n+1)*(2*m)**6*(1<<(6*z))
    s=prime_at_least(bound.bit_length())
    return s,2*s

def coefficients(bits,s):
    return tuple(pack(bits[j:j+s]) for j in range(0,len(bits),s))

def code_bit(coeff,s,index):
    f=irreducible_prime(s);mask=(1<<s)-1
    a=index&mask;b=(index>>s)&mask
    value=0
    for c in reversed(coeff):value=mul(value,a,f)^c
    return trace(mul(b,value,f),f)

@lru_cache(None)
def constant_parameters(n,k):
    if not 1<=k<=n or k**5<n:return None
    m=root(k,8);s,ell=code_parameters(n,m,2)
    overlap=(k.bit_length()-1)//4;t=overlap+1
    d=(32*ell*ell+t-1)//t
    need=d+m*(1<<overlap)+s+ceil_log2(16*m)
    if t>ell or d<2*ell or m>=4**t or need>k:return None
    return m,s,ell,overlap,d

@lru_cache(None)
def constant_design(n,k):
    p=constant_parameters(n,k)
    if p is None:return None
    m,s,ell,a,d=p;sets=[];masks=[]
    for indices in combinations(range(d),ell):
        mask=sum(1<<j for j in indices)
        if all((mask&old).bit_count()<=a for old in masks):
            sets.append(indices);masks.append(mask)
            if len(sets)==m:return tuple(sets)
    raise AssertionError('proved greedy-design existence failed')

def restrict_seed(seed,indices):
    return sum(((seed>>j)&1)<<i for i,j in enumerate(indices))

def constant_extract(bits,seed,k,take):
    p=constant_parameters(len(bits),k)
    if p is None or not 0<=take<=p[0]:raise ValueError('invalid constant extractor')
    coeff=coefficients(bits,p[1]);design=constant_design(len(bits),k)
    return tuple(code_bit(coeff,p[1],restrict_seed(seed,design[i])) for i in range(take))

@lru_cache(None)
def small_parameters(n,r):
    if not 1<=r<=n or r**5<n:return None
    m=root(r**9,10);z=root(r,1000)
    if z<2:return None
    s,ell=code_parameters(n,m,z)
    p=prime_at_least(max(3,ceil_log2(ell)));Q=1<<p
    count=1;A=0
    while count<m:count*=Q;A+=1
    A=max(1,A)
    d=Q*Q;need=d+m*(1<<(A-1))+s+z+ceil_log2(4*m)
    if ell>Q or A>ell or 2*need>r:return None
    return m,z,s,ell,p,Q,A,d

@lru_cache(None)
def small_set(n,r,i):
    pp=small_parameters(n,r)
    if pp is None:raise ValueError('invalid small-error design')
    m,z,s,ell,p,Q,A,d=pp
    if not 0<=i<m:raise ValueError('design row')
    f=irreducible_prime(p);coeff=[];value=i
    for _ in range(A):coeff.append(value%Q);value//=Q
    assert value==0
    result=[]
    for alpha in range(ell):
        beta=0
        for c in reversed(coeff):beta=mul(beta,alpha,f)^c
        result.append(alpha*Q+beta)
    return tuple(result)

def small_extract(bits,seed,r,take):
    p=small_parameters(len(bits),r)
    if p is None or not 0<=take<=p[0]:raise ValueError('invalid small-error extractor')
    coeff=coefficients(bits,p[2])
    return tuple(code_bit(coeff,p[2],restrict_seed(seed,small_set(len(bits),r,i))) for i in range(take))

@lru_cache(None)
def near_parameters(n,r):
    p=small_parameters(n,r);M=r//4
    if p is None or M<1:return None
    J=(M+p[0]-1)//p[0]
    return M,J,J*p[-1]

def near_extract(bits,seed,r):
    pp=near_parameters(len(bits),r)
    if pp is None:raise ValueError('invalid iterated extractor')
    M,J,dtotal=pp;p=small_parameters(len(bits),r);d=p[-1];mask=(1<<d)-1;out=[]
    for j in range(J):
        take=min(p[0],M-len(out))
        out.extend(small_extract(bits,(seed>>(j*d))&mask,r,take))
    return tuple(out)

@lru_cache(None)
def cubic_parameters(w):
    if w<1:return None
    k=prime_at_least(2*w);out=isqrt(w)//100
    d0=(2*k+4)//5+10
    if out<1 or 5*k>11*w or 10*(w-d0)<w:return None
    return k,out

def cubic_extract(bits,take):
    if len(bits)%2:raise ValueError('two equal slices')
    w=len(bits)//2;p=cubic_parameters(w)
    if p is None or not 0<=take<=p[1]:raise ValueError('invalid cubic seed')
    f=irreducible_prime(p[0]);value=power(pack(bits),7,f)
    return tuple((value>>i)&1 for i in range(take))

@lru_cache(None)
def somewhere_parameters(t,r):
    if t<1 or r<1 or t**100>r:return None
    initial_r=r;plan=[]
    while t>1:
        if t**50>r:return None
        w=root(r*r,5);near=near_parameters(t*r,r);cubic=cubic_parameters(w)
        if near is None or cubic is None or cubic[1]<near[2]:return None
        plan.append((t,r,w,near[0],near[2]))
        t=(t+1)//2;r=near[0]
    take=root(initial_r**97,100)
    if take>r:return None
    return tuple(plan),take

def somewhere_extract(rows):
    t=len(rows);r=len(rows[0]) if rows else 0
    if any(len(row)!=r for row in rows):raise ValueError('ragged source')
    p=somewhere_parameters(t,r)
    if p is None:raise ValueError('invalid somewhere-random parameters')
    rows=tuple(tuple(row) for row in rows)
    for t,r,w,M,seedlen in p[0]:
        flat=tuple(b for row in rows for b in row);out=[]
        for i in range(0,t,2):
            slices=rows[i][:w]+rows[min(i+1,t-1)][:w]
            seed=pack(cubic_extract(slices,seedlen))
            out.append(near_extract(flat,seed,r))
        rows=tuple(out)
    return rows[0][:p[1]]

@lru_cache(None)
def block_parameters(n,K):
    cp=constant_parameters(n,K);sp=small_parameters(n,K);R=root(K,10)
    if cp is None or sp is None or R<1 or R>min(cp[0],sp[0]):return None
    g=root(R,100)
    if g<2:return None
    sr=somewhere_parameters(g,R)
    if sr is None or sr[1]<sp[-1]:return None
    v=1<<cp[-1];groups=1;b=1
    while groups<v:groups*=g;b+=1
    return R,g,b,sp[-1]

def block_extract(blocks,K):
    n=len(blocks[0]) if blocks else 0;p=block_parameters(n,K)
    if p is None or len(blocks)!=p[2] or any(len(x)!=n for x in blocks):
        raise ValueError('invalid block-source parameters')
    R,g,b,seedlen=p;cp=constant_parameters(n,K)
    rows=tuple(constant_extract(blocks[0],seed,K,R) for seed in range(1<<cp[-1]))
    for x in blocks[1:]:
        new=[]
        for i in range(0,len(rows),g):
            group=rows[i:i+g]
            group+=((0,)*R,)*(g-len(group))
            seed=pack(somewhere_extract(group)[:seedlen])
            new.append(small_extract(x,seed,K,R))
        rows=tuple(new)
    assert len(rows)==1
    return rows[0]

@lru_cache(None)
def tree_parameters(N):
    if N<2:return None
    L0=ceil_log2(N);A=ceil_root(L0**3,10);t=1<<A;J=0;n=1
    while n<N:n*=t;J+=1
    L=A*J;k=isqrt(n);K=isqrt(k);q=isqrt(L)+(isqrt(L)**2<L)
    ell=1<<q;P=1<<ceil_root(L**4,5)
    cp=constant_parameters(n,K)
    if cp is None or ell*P>cp[0]:return None
    bp=block_parameters(n,K)
    if bp is None or bp[0]<P:return None
    H=0
    while H+q<=J and (n//(t**(H+q)))**4>=k:H+=q
    if H<q:return None
    hp=H//q;b=bp[2]
    if t<b or H*comb(t,b)>ell:return None
    confidences=tuple((ell**(6*(hp-u)+3),ell**(6*(hp-u)+6)) for u in range(hp))
    if any(p2>P for p1,p2 in confidences):return None
    return n,t,k,K,q,ell,P,H,b,confidences

def evaluate(bits):
    if any(type(b) is not int or b not in (0,1) for b in bits):raise ValueError('binary input')
    p=tree_parameters(len(bits))
    if p is None:return 0
    n,t,k,K,q,ell,P,H,b,confidences=p;data=tuple(bits)+(0,)*(n-len(bits))
    challenges={};choices={}
    def ambient(level,index):
        length=n//t**level;start=index*length
        return (0,)*start+data[start:start+length]+(0,)*(n-start-length)
    def test(level,index,challenge,confidence):
        x=ambient(level,index);cp=constant_parameters(n,K)
        target=challenge+((0,)*P,)*(ell-len(challenge))
        for seed in range(1<<cp[-1]):
            response=constant_extract(x,seed,K,ell*P)
            if all(response[j*P:j*P+confidence]==target[j][:confidence] for j in range(ell)):
                return 0 # FAIL
        return 1 # PASS
    for index in range(t**H):challenges[H,index]=()
    for level in range(H-1,-1,-1):
        child_length=n//t**(level+1)
        for index in range(t**level):
            selected=0
            for j in range(t):
                if test(level,index,challenges[level+1,index*t+j],P):selected=j;break
            choices[level,index]=selected
            rows=list(challenges[level+1,index*t+selected])
            for indices in combinations(range(t),b):
                blocks=[]
                for j in indices:
                    start=(index*t+j)*child_length
                    blocks.append(data[start:start+child_length]+(0,)*(n-child_length))
                rows.append(block_extract(tuple(blocks),K)[:P])
            challenges[level,index]=tuple(rows)
    path=[0]
    for level in range(H):path.append(path[-1]*t+choices[level,path[-1]])
    for u,(p1,p2) in enumerate(confidences):
        level=u*q;nextlevel=(u+1)*q;challenge=challenges[nextlevel,path[nextlevel]]
        if not test(level,path[level],challenge,p1):return test(level,path[level],challenge,p2)
    return 0
