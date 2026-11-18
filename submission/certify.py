#!/usr/bin/env python3
"""Certified alternating orthants of the critical cubic vertex Gaussian wave.
Directed Decimal initialization, exact fixed-point contractions, and analytic
localized Gaussian-moment Gauss-Legendre error bounds and triple tails.
All quadrature, truncation and rounding errors are enclosed explicitly.
"""
import decimal as dec, math, json, sys, time
from operator import mul
P=100
dec.getcontext().prec=P
lo=dec.Context(prec=P,rounding=dec.ROUND_FLOOR)
hi=dec.Context(prec=P,rounding=dec.ROUND_CEILING)
mid=dec.Context(prec=P,rounding=dec.ROUND_HALF_EVEN)
D=dec.Decimal
# Positive fixed-point intervals represent [a/SCALE, b/SCALE].
# Python integers have unlimited precision. All contraction sums and products
# below are exact; floor/ceiling are applied only at the stated divisions.
FP=128
SCALE=1<<FP
SCALE2=SCALE*SCALE

def fixed(a,b):
 assert 0<=a<=b
 # Decimal.as_integer_ratio is exact, including tiny positive endpoints.
 na,da=a.as_integer_ratio();nb,db=b.as_integer_ratio()
 return na*SCALE//da, (nb*SCALE+db-1)//db

def decimal_fixed(a,b):
 # Decimal(int) is exact; the two divisions are explicitly directed.
 return I(lo.divide(D(a),D(SCALE)),hi.divide(D(b),D(SCALE)))
Z=D(0); ONE=D(1)
class I:
 def __init__(self,a=0,b=None):
  if isinstance(a,I): self.a,self.b=a.a,a.b
  else: self.a=D(a);self.b=self.a if b is None else D(b)
 def __add__(a,b):
  b=I(b);return I(lo.add(a.a,b.a),hi.add(a.b,b.b))
 __radd__=__add__
 def __neg__(a):return I(a.b.copy_negate(),a.a.copy_negate())
 def __sub__(a,b):return a+-I(b)
 def __rsub__(a,b):return I(b)+-a
 def __mul__(a,b):
  b=I(b); ll=[lo.multiply(x,y) for x in (a.a,a.b) for y in (b.a,b.b)]; hh=[hi.multiply(x,y) for x in (a.a,a.b) for y in (b.a,b.b)];return I(min(ll),max(hh))
 __rmul__=__mul__
 def __truediv__(a,b):
  b=I(b);assert not b.a<=0<=b.b;return a*I(lo.divide(ONE,b.b),hi.divide(ONE,b.a))
 def __rtruediv__(a,b):return I(b)/a
 def sqrt(a):
  assert a.a>=0
  return I(lo.next_minus(mid.sqrt(a.a)),hi.next_plus(mid.sqrt(a.b)))
 def exp(a):return I(lo.next_minus(mid.exp(a.a)),hi.next_plus(mid.exp(a.b)))
 def square(a):
  if a.a>=0:return I(lo.multiply(a.a,a.a),hi.multiply(a.b,a.b))
  if a.b<=0:return I(lo.multiply(a.b,a.b),hi.multiply(a.a,a.a))
  return I(0,max(hi.multiply(a.a,a.a),hi.multiply(a.b,a.b)))
 def bounds(a):return [str(a.a),str(a.b)]
 def center(a):return mid.divide(mid.add(a.a,a.b),D(2))

def pi_interval():
 def atan(q):
  x=I(1)/q;xx=x*x;term=x;out=term
  for k in range(1,100):term=-term*xx;out=out+term/(2*k+1)
  err=(term*xx/201);radius=max(err.a.copy_abs(),err.b.copy_abs());return out+I(radius.copy_negate(),radius)
 return 16*atan(5)-4*atan(239)

def legendre(n,x):
 a=I(1);b=x
 if not n:return a,I(0)
 for k in range(2,n+1):a,b=b,((2*k-1)*x*b-(k-1)*a)/k
 return b,n*(x*b-a)/(x*x-1)

def nodes_weights(n,L):
 nodes=[];weights=[];brackets=[]
 for k in range(1,n+1):
  x=D(str(math.cos(math.pi*(k-.25)/(n+.5))))
  for it in range(8):
   p,dp=legendre(n,I(x));x=mid.subtract(x,mid.divide(p.center(),dp.center()))
  eps=D('1e-65');br=I(lo.subtract(x,eps),hi.add(x,eps))
  pa,_=legendre(n,I(br.a));pb,_=legendre(n,I(br.b))
  assert (pa.b<0 and pb.a>0) or (pb.b<0 and pa.a>0),(k,pa.bounds(),pb.bounds())
  assert -1<br.a<br.b<1
  _,dp=legendre(n,br)
  assert not dp.a<=0<=dp.b
  w=2/((1-br*br)*dp*dp)
  nodes.append((br+1)*L/2);weights.append(w*L/2);brackets.append(br)
 assert all(brackets[i].a>brackets[i+1].b for i in range(n-1))
 return nodes,weights

def main():
 t=time.monotonic();N=80;J=65;Ln=3;Ld=1;L=I(Ln)/Ld
 assert N>=1 and J>=1 and J%2==1
 pi=pi_interval();r2=I(2).sqrt();c2=3/(2*pi);c1=(6/pi).sqrt()
 xs,ws=nodes_weights(N,L)
 print('nodes certified',N,'seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 state=[]
 for x,wx in zip(xs,ws):
  row=[]
  for y,wy in zip(xs,ws):
   f=c2*(-I(9)/2*(x*x+y*y)-6*r2*x*y).exp()*wx*wy
   row.append(fixed(f.a,f.b))
  state.append(row)
 # Exact factorization of exp[-6(z+sqrt(2)y+x/2)^2].
 # Only O(N^2) transcendental evaluations; all remaining products are positive.
 ax=[(-I(3)/2*x.square()).exp() for x in xs]
 by=[(-12*y.square()).exp() for y in xs]
 cz=[(-6*z.square()).exp() for z in xs]
 pair_xy=[];pair_xz=[];pair_yz=[]
 for i,x in enumerate(xs):
  pair_xy.append([ax[i]*(-6*r2*x*y).exp() for y in xs])
  pair_xz.append([(-6*x*z).exp() for z in xs])
 for j,y in enumerate(xs):
  pair_yz.append([c1*ws[k]*by[j]*cz[k]*(-12*r2*y*z).exp() for k,z in enumerate(xs)])
 # K[j,k,i] = pair_xy[i,j] * pair_xz[i,k] * pair_yz[j,k].
 # Keep only the quadratic factors. Each update first weights state by
 # pair_xy, contracts against pair_xz, then multiplies by 4*pair_yz.
 # At path length j, state encloses 4^(j-1) times the unscaled tensor state.
 # All factors are nonnegative. Quantize their already certified Decimal
 # bounds outward, then enclose this same finite tensor rule by exact integers.
 # The rounding argument is given beside the contractions below.
 A=[[fixed(x.a,x.b) for x in row] for row in pair_xy]
 Bcols=[[fixed(pair_xz[i][k].a,pair_xz[i][k].b) for i in range(N)] for k in range(N)]
 C=[[fixed(x.a,x.b) for x in row] for row in pair_yz]
 Blo=[[a for a,b in row] for row in Bcols]
 Bhi=[[b for a,b in row] for row in Bcols]
 print('kernel factors certified','seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 # Positive quadrature norm used by the localized derivative proof.
 qhalf=sum((w*c1*c for w,c in zip(ws,cz)),I(0))
 h=I(max(D('0.5'),qhalf.b))
 print('positive quadrature norm',qhalf.bounds(),'h',h.bounds(),file=sys.stderr,flush=True)
 # Exact two-step positive completion, with directed alternating atan series.
 x=I(1)/r2; xx=x*x;term=x;atan=term
 for k in range(1,220):
  term=-term*xx;atan=atan+term/(2*k+1)
 err=term*xx/441;rad=max(abs(err.a),abs(err.b))
 kappa=(atan+I(-rad,rad))/(2*pi)
 print('two-step bound',kappa.bounds(),file=sys.stderr,flush=True)
 x=I(1)/(2*r2);xx=x*x;term=x;atan=term
 for k in range(1,120):
  term=-term*xx;atan=atan+term/(2*k+1)
 err=term*xx/241;rad=max(abs(err.a),abs(err.b))
 exact_p1=(atan+I(-rad,rad))/pi
 endpoint_tail=c2*c1*(-6*L*L).exp()/(864*r2*L*L*L)
 middle_tail=c2*c1*(-I(33)*L*L/2).exp()/(9504*L*L*L)
 print('triple outlier bounds',endpoint_tail.bounds(),middle_tail.bounds(),file=sys.stderr,flush=True)
 prefactor=(I(Ln**(2*N+1))/Ld**(2*N+1))*math.factorial(N)**3/((2*N+1)*2**N*math.factorial(2*N)**2)*c2
 results=[];gamma=I(0)
 for j in range(1,J+1):
  sa=sum(a for row in state for a,b in row)
  sb=sum(b for row in state for a,b in row)
  s=decimal_fixed(sa,sb)/4**(j-1)
  # Localized Gaussian even-moment error. Reversibly anchor at a pair
  # containing coordinate i. Its derivative touches only s_i neighbors,
  # leaving n-1-s_i positive completion operators with norm at most h.
  n=j+1
  qs=[9,9] if n==2 else ([12,33,12] if n==3 else [12,36]+[39]*(n-4)+[36,12])
  qsum=I(0)
  for ii,q in enumerate(qs):
   neighbors=min(2,ii)+min(2,n-1-ii)
   term=I(q**N)
   for _ in range(neighbors-1):term=term*c1
   for _ in range(neighbors):term=term*L
   for _ in range(n-1-neighbors):term=term*h
   qsum=qsum+term
  quad=prefactor*qsum
  # For n>=3, center an anchor triple on each interior outlier;
  # at the two endpoints use the nearest triple. Its exterior arms
  # contain arbitrary positive integral operators, stripped in blocks
  # of two with mass kappa and one leftover step with mass 1/2.
  tail=I(0)
  if n>=3:
   for ii in range(n):
    start=0 if ii==0 else n-3 if ii==n-1 else ii-1
    left=start;right=n-start-3
    term=endpoint_tail if ii in (0,n-1) else middle_tail
    for _ in range(left//2+right//2):term=term*kappa
    term=term/(2**(left%2+right%2))
    tail=tail+term
   p=2*(s+I(-quad.b,quad.b)+I(0,tail.b))
  else:
   p=exact_p1;quad=I(0);tail=I(0)
  coeff=I(3)*I(2**max(j-2,0))/(2 if j==1 else 1)
  gamma=gamma+((-1 if j%2==0 else 1)*coeff)*p
  results.append({'j':j,'probability':p.bounds(),'quadrature_error':str(quad.b),'tail_error':str(tail.b),'partial_gamma':gamma.bounds()})
  print('j',j,'p',p.bounds(),'gamma',gamma.bounds(),flush=True,file=sys.stderr)
  if j==J:break
  # If m/S and a/S enclose M,A, the product has denominator S^2.
  # Flooring/ceiling m*a/S restores endpoints with denominator S.
  scaled=[]
  for jj in range(N):
   col=[]
   for ii in range(N):
    sa,sb=state[ii][jj];a,b=A[ii][jj]
    col.append(((sa*a)>>FP,(sb*b+SCALE-1)>>FP))
   scaled.append(col)
  Slo=[[a for a,b in row] for row in scaled]
  Shi=[[b for a,b in row] for row in scaled]
  nxt=[]
  for jj in range(N):
   row=[]
   for kk in range(N):
    a=sum(map(mul,Slo[jj],Blo[kk]))
    b=sum(map(mul,Shi[jj],Bhi[kk]))
    ca,cb=C[jj][kk]
    # The exact dot products a,b have denominator S^2. Multiplication
    # by 4*C contributes another S and the exact state-rescaling factor 4.
    # Divide by S^2 to restore the represented denominator S.
    # No floating-point multiply, addition, or matrix library is used.
    row.append(((4*a*ca)>>(2*FP),(4*b*cb+SCALE2-1)>>(2*FP)))
   nxt.append(row)
  state=nxt
 target=I(130647566423167664341749381)/1000000000000000000000000000;passed=gamma.b<target.a
 print(json.dumps({'certified':passed,'target':{'numerator':130647566423167664341749381,'denominator':1000000000000000000000000000},'gamma_upper_bound':str(gamma.b),'N':N,'J':J,'L':{'numerator':Ln,'denominator':Ld},'positive_quadrature_norm':qhalf.bounds(),'positive_operator_bound':h.bounds(),'completion_constant':kappa.bounds(),'exact_p1':exact_p1.bounds(),'endpoint_tail':endpoint_tail.bounds(),'middle_tail':middle_tail.bounds(),'state_rescaling_base':4,'error_bound':'localized-gaussian-moments-triple-outliers-two-step-completions','precision':P,'fixed_point_bits':FP,'intervals':results,'seconds':time.monotonic()-t},indent=2))
 return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
