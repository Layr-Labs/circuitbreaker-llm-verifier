#!/usr/bin/env python3
"""Certified alternating orthants of the critical cubic vertex Gaussian wave.
A standalone Decimal interval computation and tensor Gauss-Legendre error bound.
All quadrature, truncation and rounding errors are enclosed explicitly.
"""
import decimal as dec, math, json, sys, time
P=100
dec.getcontext().prec=P
lo=dec.Context(prec=P,rounding=dec.ROUND_FLOOR)
hi=dec.Context(prec=P,rounding=dec.ROUND_CEILING)
mid=dec.Context(prec=P,rounding=dec.ROUND_HALF_EVEN)
D=dec.Decimal
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
  eps=D('1e-60');br=I(lo.subtract(x,eps),hi.add(x,eps))
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
 t=time.monotonic();N=int(sys.argv[1]) if len(sys.argv)>1 else 72;J=int(sys.argv[2]) if len(sys.argv)>2 else 13;L=3
 assert N>=1 and J>=1 and J%2==1
 pi=pi_interval();r2=I(2).sqrt();c2=3/(2*pi);c1=(6/pi).sqrt()
 xs,ws=nodes_weights(N,L)
 print('nodes certified',N,'seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 state=[]
 for x,wx in zip(xs,ws):
  row=[]
  for y,wy in zip(xs,ws):
   f=c2*(-I(9)/2*(x*x+y*y)-6*r2*x*y).exp()*wx*wy
   row.append((f.a,f.b))
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
 # K[j][k][i] permits contiguous contraction of previous index i.
 K=[]
 for j in range(N):
  layer=[]
  for k in range(N):
   row=[];c=pair_yz[j][k]
   for i in range(N):
    a=pair_xy[i][j];b=pair_xz[i][k]
    row.append((lo.multiply(lo.multiply(a.a,b.a),c.a),hi.multiply(hi.multiply(a.b,b.b),c.b)))
   layer.append(row)
  K.append(layer)
 print('kernel certified','seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 results=[];gamma=I(0)
 for j in range(1,J+1):
  s=I(0)
  for row in state:
   for a,b in row:s=s+I(a,b)
  # Uniform Cauchy bound: f_n(0)<2^n, max precision diagonal <=39.
  # exp(N)<3^N and r^2=2N/39. Tensor rule has positive weights and norm L.
  n=j+1
  numer=n*L**n*2**n*math.factorial(N)**4*(L*L*117)**N
  denom=(2*N+1)*math.factorial(2*N)**2*(2*N)**N
  quad=I(numer)/denom
  # Tail union bound using an adjacent negative-correlated pair:
  # integral_{x>L,y>0} f2 <= 3/(2pi)*exp(-9L²/2)/(9L)*sqrt(pi/18)
  tail=n*c2*(-I(9)*L*L/2).exp()/(9*L)*(pi/18).sqrt()
  p=2*(s+I(-quad.b,quad.b)+I(0,tail.b))
  coeff=I(3)*I(2**max(j-2,0))/(2 if j==1 else 1)
  gamma=gamma+((-1 if j%2==0 else 1)*coeff)*p
  results.append({'j':j,'probability':p.bounds(),'quadrature_error':str(quad.b),'tail_error':str(tail.b),'partial_gamma':gamma.bounds()})
  print('j',j,'p',p.bounds(),'gamma',gamma.bounds(),flush=True,file=sys.stderr)
  if j==J:break
  nxt=[]
  for jj in range(N):
   row=[]
   for kk in range(N):
    a=Z;b=Z
    for ii,(ka,kb) in enumerate(K[jj][kk]):
     sa,sb=state[ii][jj];a=lo.fma(sa,ka,a);b=hi.fma(sb,kb,b)
    row.append((a,b))
   nxt.append(row)
  state=nxt
 target=I(2613)/20000;passed=gamma.b<target.a
 print(json.dumps({'certified':passed,'target':{'numerator':2613,'denominator':20000},'gamma_upper_bound':str(gamma.b),'N':N,'J':J,'precision':P,'intervals':results,'seconds':time.monotonic()-t},indent=2))
 return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
