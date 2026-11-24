"""Directed centered sign moments, central-coordinate interpolation and certified GL32."""
from functools import lru_cache
from math import isqrt, comb
from fractions import Fraction as F
from itertools import combinations
import json
S=2**72
def ceil(p,q): return -((-p)//q)
@lru_cache(16384)
def rat(p,q=1):
    f=F(p)/q;return f.numerator*S//f.denominator,ceil(f.numerator*S,f.denominator)
ZERO,ONE,TWO=rat(0),rat(1),rat(2)
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return -x[1],-x[0]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
    if x[0]>=0 and y[0]>=0:return x[0]*y[0]//S,ceil(x[1]*y[1],S)
    p=(x[0]*y[0],x[1]*y[0],x[0]*y[1],x[1]*y[1])
    return min(p)//S,ceil(max(p),S)
def sq(x):
    a,b=x;return (0 if a<=0<=b else min(a*a,b*b)//S, ceil(max(a*a,b*b),S))
def div(x,y):
    if y[0]<=0:raise ArithmeticError('nonpositive divisor')
    return (x[0]*S//(y[0] if x[0]<0 else y[1]),ceil(x[1]*S,y[0] if x[1]>=0 else y[1]))
def root(x):
    if x[0]<0:raise ArithmeticError('negative square root')
    a=isqrt(x[0]*S);b=isqrt(x[1]*S);return a,b+int(b*b<x[1]*S)
def intdiv(x,d):return x[0]//d,ceil(x[1],d)
def scale(x,n):return (x[0]*n,x[1]*n) if n>=0 else neg(scale(x,-n))
def clamp(x,a=-S,b=S):
    y=max(a,x[0]),min(b,x[1])
    if y[0]>y[1]:raise ArithmeticError('correlation inconsistent')
    return y
@lru_cache(16384)
def atan(x):
    if not 0<=x[0]<=x[1]<=S//2:raise ArithmeticError('atan range')
    x2=sq(x);k=16;acc=rat(1,2*k+1)
    for j in range(k-1,-1,-1):acc=sub(rat(1,2*j+1),mul(x2,acc))
    a=mul(x,acc);p=x
    for _ in range(17):p=mul(p,x2)
    # Taylor remainder after index 16 has negative sign; upper suffices at even index
    return a[0]-intdiv(p,35)[1],a[1]
PI=sub(scale(atan(rat(1,5)),16),scale(atan(rat(1,239)),4));K=div(TWO,PI)
def apos(x):
    w=div(x,add(ONE,root(sub(ONE,sq(x)))))
    z=div(w,add(ONE,root(add(ONE,sq(w)))))
    return scale(atan(z),4)
@lru_cache(16384)
def asin(x):
    a,b=clamp(x)
    if a>=0:return apos((a,b))
    if b<=0:return neg(apos((-b,-a)))
    return -apos((-a,-a))[1],apos((b,b))[1]
N=32
ROOTS=(228126501309447433224, 682249549093029178565, 1130002619273779262979, 1567205164133916262516, 1989775143913665128883, 2393767140111235560609, 2775409193448129022088, 3131138022810122571336, 3457632296895626724748, 3751843649262056999671, 4011025150725659818630, 4232756983876586986725, 4414969117639565221176, 4555960939724158464571, 4654418767251863546368, 4709445435775136009011)
def legsign(x):
    a,b=1,x
    for k in range(1,N):a,b=b,(2*k+1)*x*b-k*k*S*S*a
    return (b>0)-(b<0)
def leg_interval(n,x):
    a,b=ONE,x
    for k in range(1,n):a,b=b,intdiv(sub(scale(mul(x,b),2*k+1),scale(a,k)),k+1)
    return b
def quad():
    qs=[];prior=-S
    for l,h in [(-x-1,-x) for x in reversed(ROOTS)]+[(x,x+1) for x in ROOTS]:
        if not prior<l<h<S or legsign(l)*legsign(h)!=-1:raise ArithmeticError('all 32 root brackets')
        prior=h;x=l,h;arg=intdiv(add(ONE,x),2)
        w=div(sub(ONE,sq(x)),scale(sq(leg_interval(N-1,x)),N*N))
        if w[0]<=0 or not (0<arg[0]<=arg[1]<S):raise ArithmeticError('quadrature positive')
        qs.append((arg,sq(arg),w))
    if not sum(w[0] for _,_,w in qs)<=S<=sum(w[1] for _,_,w in qs):raise ArithmeticError('weight sum')
    return qs
Q=quad();err=rat(F(4*512)*F(2,3)**(2*N)/(1-F(2,3)))[1]
def cov(d):
    d=abs(d);a,t=F(28,57),F(41,57)
    return (F(1) if d==0 else t**d/2+a*t**(d-1)+(F(d-1,2)*a*a*t**(d-2) if d>1 else 0))
def ellipse():
    m=6;G=[[cov(i-j) for j in range(m)]+[F(i==j) for j in range(m)] for i in range(m)]
    for i in range(m):
        p=G[i][i]
        if p<=0:raise ArithmeticError('SPD')
        G[i]=[a/p for a in G[i]]
        for j in range(m):
            if i!=j:
                f=G[j][i];G[j]=[u-f*v for u,v in zip(G[j],G[i])]
    if any(G[i][m+i]>F(20,3) for i in range(m)):raise ArithmeticError('analytic domain')
ellipse()
def moment(A):
    n=len(A)
    if n==2:return mul(K,asin(A[0][1]))
    # Disconnect the higher of the two central positions, independent of state.
    center=n//2;idx=[center]+[i for i in range(n) if i!=center];A=[[A[i][j] for j in idx] for i in idx]
    total=ZERO
    for j in range(1,n):
        b=A[0][j];bs=sq(b);ids=[k for k in range(1,n) if k!=j]
        v={k:sub(A[0][k],mul(b,A[j][k])) for k in ids}
        c={(k,l):sub(A[k][l],mul(A[k][j],A[j][l])) for ik,k in enumerate(ids) for l in ids[ik:]}
        vv={(k,l):mul(v[k],v[l]) for ik,k in enumerate(ids) for l in ids[ik:]}
        acc=ZERO
        for t,tt,w in Q:
            d=sub(ONE,mul(tt,bs));f=div(tt,d);B={p:sub(val,mul(f,vv[p])) for p,val in c.items()}
            def cr(k,l):return div(B[k,l],root(mul(B[k,k],B[l,l])))
            if n==4:m=mul(K,asin(cr(*ids)))
            else:m=moment([[ONE if k==l else cr(*sorted((k,l))) for l in ids] for k in ids])
            acc=add(acc,mul(w,div(m,root(d))))
        total=add(total,mul(b,acc))
    y=mul(K,total);return clamp((y[0]-err,y[1]+err))
def key(s):
    x=tuple(s[j+1]-s[j] for j in range(len(s)-1));return min(x,x[::-1])
def compute():
    rows={():ONE};alt={}
    for n in range(2,7):
        acc=ZERO
        for d in range(0,n+1,2):
            for ss in combinations(range(n),d):
                k=key(ss)
                if k not in rows:rows[k]=moment([[rat(cov(i-j)) for j in ss] for i in ss])
                q=rows[k];acc=add(acc,neg(q) if sum(ss)%2 else q)
        alt[n-1]=intdiv(acc,1<<(n-1))
    u=ZERO
    for k in range(1,6):u=add(u,alt[k] if k%2 else neg(alt[k]))
    return rows,u
if __name__=='__main__':
    rows,u=compute()
    if u[1]*2**39>73670242301*S:raise ArithmeticError('U5 cap missed')
    print(json.dumps({'N':N,'S':S,'U5':u,'rows':{''.join(map(str,k)):v for k,v in sorted(rows.items())},'method':'Gauss-Legendre32 middle-coordinate'}))
