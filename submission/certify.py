#!/usr/bin/env python3
"""Outward-rounded certificate of gamma5 <= U < 1341/10000. Standard library only."""
import heapq
import json
import math
import sys
import time
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

def dn(x): return math.nextafter(x, -math.inf)
def up(x): return math.nextafter(x, math.inf)
def decimal_bound(x, upper=False):
    with localcontext() as context:
        context.prec=100
        return str(Decimal.from_float(x).quantize(Decimal("1e-18"),rounding=ROUND_CEILING if upper else ROUND_FLOOR))
def pair(a): return [decimal_bound(a.l),decimal_bound(a.h,True)]

class DomainError(ArithmeticError): pass

class I:
    __slots__ = ('l', 'h')
    def __init__(self, l=0.0, h=None):
        if isinstance(l, I): self.l,self.h=l.l,l.h; return
        l=float(l); h=l if h is None else float(h)
        if not (math.isfinite(l) and math.isfinite(h) and l<=h): raise DomainError('invalid interval')
        self.l,self.h=l,h
    def __add__(a,b):
        if isinstance(b,J): return NotImplemented
        b=as_i(b); return I(dn(a.l+b.l),up(a.h+b.h))
    __radd__=__add__
    def __neg__(a): return I(-a.h,-a.l)
    def __sub__(a,b): return a+-b if isinstance(b,(I,J)) else a+I(-b)
    def __rsub__(a,b): return as_i(b)+-a
    def __mul__(a,b):
        if isinstance(b,J): return NotImplemented
        b=as_i(b); p=(a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h)
        if not all(math.isfinite(x) for x in p): raise DomainError('nonfinite product')
        return I(dn(min(p)),up(max(p)))
    __rmul__=__mul__
    def __truediv__(a,b):
        if isinstance(b,J): return NotImplemented
        return a*inverse(as_i(b))
    def __rtruediv__(a,b): return as_i(b)*inverse(a)

def as_i(x): return x if isinstance(x,I) else I(x)
def mag(x): return max(abs(x.l),abs(x.h))
def inverse(a):
    if isinstance(a,J):
        r=inverse(a.v); r2=square(r)
        return J(r,-a.d*r2,2*square(a.d)*r2*r-a.dd*r2)
    if a.l<=0<=a.h: raise DomainError('zero divisor')
    return I(dn(1/a.h),up(1/a.l))
def square(a):
    if isinstance(a,J): return J(square(a.v),2*a.v*a.d,2*square(a.d)+2*a.v*a.dd)
    if a.l<=0<=a.h: return I(0,up(max(a.l*a.l,a.h*a.h)))
    return I(dn(min(a.l*a.l,a.h*a.h)),up(max(a.l*a.l,a.h*a.h)))
def root(a):
    if isinstance(a,J):
        v=root(a.v); r=1/(2*v)
        return J(v,r*a.d,r*a.dd-square(a.d)/(4*a.v*v))
    if a.l<=0: raise DomainError('sqrt domain')
    l=math.sqrt(a.l); h=math.sqrt(a.h)
    while up(l*l)>a.l: l=dn(l)
    while dn(h*h)<a.h: h=up(h)
    return I(l,h)
def atan_small(z):
    zz=square(z); term=z; result=z
    for k in range(1,25):
        term=term*-zz; result=result+term/(2*k+1)
    e=(I(mag(term))*I(mag(zz))/51).h
    return result+I(-e,e)
def atan_i(z):
    if z.h<0: return -atan_i(-z)
    if z.l<0: return I(atan_i(I(z.l)).l,atan_i(I(z.h)).h)
    k=0
    while z.h>.2:
        z=z/(1+root(1+square(z))); k+=1
    return atan_small(z)*math.ldexp(1.,k)
def asin(a):
    if isinstance(a,J):
        z=1-square(a.v); r=1/root(z)
        return J(asin(a.v),r*a.d,r*a.dd+a.v*r/z*square(a.d))
    if a.l<=-1 or a.h>=1: raise DomainError('asin domain')
    return atan_i(a/root(1-square(a)))

class J:
    __slots__=('v','d','dd')
    def __init__(self,v=0,d=0,dd=0): self.v,self.d,self.dd=as_i(v),as_i(d),as_i(dd)
    def __add__(a,b):
        b=as_j(b); return J(a.v+b.v,a.d+b.d,a.dd+b.dd)
    __radd__=__add__
    def __neg__(a): return J(-a.v,-a.d,-a.dd)
    def __sub__(a,b): return a+-as_j(b)
    def __rsub__(a,b): return as_j(b)+-a
    def __mul__(a,b):
        b=as_j(b)
        return J(a.v*b.v,a.d*b.v+a.v*b.d,a.dd*b.v+2*a.d*b.d+a.v*b.dd)
    __rmul__=__mul__
    def __truediv__(a,b): return a*inverse(as_j(b))
    def __rtruediv__(a,b): return as_j(b)*inverse(a)
def as_j(x): return x if isinstance(x,J) else J(x)

PI=16*atan_small(I(1)/5)-4*atan_small(I(1)/239)
def corr(n):
    r=1/root(I(2)); a=I(1)/4+r; cov=[I(1)]
    rd=I(1)
    for d in range(1,n):
        rd=rd*r; cov.append((a+I(d)/4)*rd*(-1 if d%2 else 1))
    return [[cov[abs(i-j)] for j in range(n)] for i in range(n)]
def conditioned(a,i,j):
    keep=[k for k in range(len(a)) if k!=i and k!=j]
    den=a[i][i]*a[j][j]-square(a[i][j])
    return [[a[x][y]-(a[j][j]*a[x][i]*a[y][i]+a[i][i]*a[x][j]*a[y][j]-a[i][j]*(a[x][i]*a[y][j]+a[x][j]*a[y][i]))/den for y in keep] for x in keep]
def smallorth(a):
    n=len(a)
    if n==0:return I(1)
    if n==1:return I(.5)
    result=I(.25 if n==2 else .125)
    for i in range(n):
        for j in range(i+1,n):
            result=result+asin(a[i][j]/root(a[i][i]*a[j][j]))/((2 if n==2 else 4)*PI)
    return result
def derivorth(a,t):
    n=len(a); b=[[a[i][j] if i==j else a[i][j]*t for j in range(n)] for i in range(n)]
    result=I(0)
    for i in range(n):
        for j in range(i+1,n):
            k=a[i][j]/(2*PI*root(a[i][i]*a[j][j]-square(b[i][j])))
            result=result+k*smallorth(conditioned(b,i,j))
    return result

RUNS=[]
def integrate(a,tolerance):
    sequence=0
    def evaluate(lo,hi,depth):
        nonlocal sequence
        # At depth <= 50 all grid endpoints, widths, and midpoints are exact binary64 values.
        if depth>50: raise ArithmeticError('dyadic depth limit')
        sequence+=1
        try:
            allv=derivorth(a,J(I(lo,hi),1,0)); mid=derivorth(a,I((lo+hi)/2))
            err=(square(I(hi-lo))*I(mag(allv.dd))/24).h
            estimate=(mid+I(-err,err))*(hi-lo); natural=allv.v*(hi-lo)
            low=max(estimate.l,natural.l); high=min(estimate.h,natural.h)
            if not (math.isfinite(low) and math.isfinite(high) and low<=high):raise DomainError('bad enclosure')
        except DomainError:low,high=-math.inf,math.inf
        return (-(high-low),sequence,lo,hi,depth,low,high)
    queue=[evaluate(0.,1.,0)]
    def split(cell):
        _,_,lo,hi,depth,_,_=cell; m=(lo+hi)/2
        if not lo<m<hi: raise ArithmeticError('degenerate subdivision')
        return evaluate(lo,m,depth+1),evaluate(m,hi,depth+1)
    while not math.isfinite(queue[0][5]):
        for child in split(heapq.heappop(queue)):heapq.heappush(queue,child)
    low=sum(x[5] for x in queue); high=sum(x[6] for x in queue)
    while high-low>tolerance:
        cell=heapq.heappop(queue);low-=cell[5];high-=cell[6]
        for child in split(cell):
            if not math.isfinite(child[5]): raise ArithmeticError('lost finite enclosure')
            low+=child[5];high+=child[6];heapq.heappush(queue,child)
        if len(queue)>100000:raise ArithmeticError('cell budget')
    result=I(0)
    for cell in sorted(queue):result=result+I(cell[5],cell[6])
    RUNS.append({'dimension':len(a),'leaves':len(queue),'max_depth':max(x[4] for x in queue),'integral':pair(result)})
    return result
def orthant(a,tolerance):return smallorth(a) if len(a)<=3 else I(math.ldexp(1.,-len(a)))+integrate(a,tolerance)
def exp_small(x):
    term=I(1);result=I(1)
    for k in range(1,25):term=term*x/k;result=result+term
    nxt=term*x/25;err=(I(mag(nxt))/(1-I(mag(x))/26)).h
    return result+I(-err,err)
def cdf_small(x,phi0):
    term=x;result=x
    for k in range(1,21):term=term*-square(x)/(2*k);result=result+term/(2*k+1)
    nxt=term*-square(x)/(42*43)
    return I(.5)+phi0*(result+I(-mag(nxt),mag(nxt)))
def ldl_check(a):
    n=len(a);lower=[[I(0) for _ in range(n)] for _ in range(n)];d=[]
    for i in range(n):
        pivot=a[i][i]
        for k in range(i):pivot=pivot-square(lower[i][k])*d[k]
        if not pivot.l>0:raise ArithmeticError('SPD failure')
        d.append(pivot);lower[i][i]=I(1)
        for j in range(i+1,n):
            v=a[j][i]
            for k in range(i):v=v-lower[j][k]*lower[i][k]*d[k]
            lower[j][i]=v/pivot
    return min(x.l for x in d)
def main():
    start=time.monotonic();tolerance=3e-6
    a6=corr(6);pivot=ldl_check(a6);a=corr(5)
    p1=2*smallorth(corr(2));p2=2*smallorth(corr(3));p3=2*orthant(corr(4),tolerance);q=orthant(a,tolerance)
    r=a6[5][:5];aug=[a[i][:]+[r[i]] for i in range(5)]
    for k in range(5):
        pivot_i=aug[k][k]
        if not pivot_i.l>0:raise ArithmeticError('linear-solve pivot')
        for j in range(k,6):aug[k][j]=aug[k][j]/pivot_i
        for i in range(5):
            if i==k:continue
            factor=aug[i][k]
            for j in range(k,6):aug[i][j]=aug[i][j]-factor*aug[k][j]
    variance=I(1)
    for i in range(5):variance=variance-r[i]*aug[i][5]
    sigma=root(variance);phi0=1/root(2*PI);b=[None]*5
    for i in range(3):
        keep=[j for j in range(5) if j!=i]
        face=[[a[j][k]-a[j][i]*a[k][i] for k in keep] for j in keep]
        b[i]=b[4-i]=phi0*orthant(face,tolerance)
    h=[[I(0) for _ in range(5)] for _ in range(5)]
    for i in range(5):
        for j in range(i+1,5):h[i][j]=h[j][i]=smallorth(conditioned(a,i,j))/(2*PI*root(1-square(a[i][j])))
    for i in range(5):
        for j in range(5):
            if i!=j:h[i][i]=h[i][i]-a[i][j]*h[i][j]
    rb=I(0);rhr=I(0)
    for i in range(5):
        rb=rb+r[i]*b[i]
        for j in range(5):rhr=rhr+r[i]*r[j]*h[i][j]
    center=I(-1)/2;kappa=phi0*exp_small(I(-1)/2)/2
    cdf=cdf_small(center,phi0);density=phi0*exp_small(-square(center)/2)
    aq=2*(cdf-center*density+kappa*square(center)+kappa*(1-variance)/variance-1)
    br=2*(density-2*kappa*center)/sigma;cr=2*kappa/variance
    majorant=p1-p2+p3+aq*q+br*rb+cr*rhr
    target=I(1341)/10000;passed=majorant.h<target.l
    result={'certified':passed,'target':{'numerator':1341,'denominator':10000},'gamma5_majorant':pair(majorant),'gamma5_upper_bound':decimal_bound(majorant.h,True),'p1':pair(p1),'p2':pair(p2),'p3':pair(p3),'q':pair(q),'rb':pair(rb),'rHr':pair(rhr),'conditional_variance':pair(variance),'ldl_pivot_lower_bound':decimal_bound(pivot),'quadrature':RUNS,'seconds':time.monotonic()-start}
    print(json.dumps(result,allow_nan=False,separators=(',',':')))
    return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
