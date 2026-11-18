#!/usr/bin/env python3
"""Certified alternating orthants of the critical cubic vertex Gaussian wave.
Directed Decimal initialization, exact fixed-point contractions, and analytic
localized Gaussian-moment Gauss-Legendre error bounds and triple tails.
All quadrature, truncation and rounding errors are enclosed explicitly.
Exact buffered-event arithmetic also supports the structural graph improvement.
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

def certify_buffered_curvature():
    """Verify the rational curvature constants; return compact exact results."""
    from fractions import Fraction as Q

    scale = 10**30
    exp_n = 65536

    def exp_lower(n, d):
        # Lower bound for exp(-n/d), in units of 1/scale.
        assert 0 <= n < exp_n*d
        x = scale - (scale*n + exp_n*d - 1)//(exp_n*d)
        for unused in range(16):
            x = x*x//scale
        return x

    def density(n, d):
        # (69/50)*(1-6*(n/d)^2/65536)^65536, rounded down.
        return 69*exp_lower(6*n*n, d*d)//50

    def buffer(n, d=200000):
        # Integrated cubic lower bound, evaluated at min(n/d,1/2).
        assert n >= 0
        if 2*n > d:
            n, d = 1, 2
        a = 69*(35*n*d**6 - 70*n**3*d**4 + 126*n**5*d**2 - 180*n**7)
        b = 25*35*d**7
        assert a >= 0
        return scale*a//b

    def one_root_lower(m, step):
        # All coordinate units are 1/400; a<=1/2 and b<=4/5.
        assert (200-m) % step == 0 and (320-m) % step == 0
        total = 0
        for a in range(m, 200, step):
            A = a + step
            inner = 0
            for b in range(m, 320, step):
                B = b + step
                f = density(200*B+283*A, 80000)
                g = buffer(707*b+250*a-500*m)
                inner += f*g//scale
            inner = inner*step//400
            term = density(A, 400)*(inner*inner//scale)//scale
            total += term*step//400
        return total

    # Rational sqrt(2) brackets and sqrt(6/pi)>69/50, using pi<22/7.
    assert 707**2 < 2*500**2 and 283**2 > 2*200**2
    assert Q(69, 50)**2 < Q(21, 11)
    # The residual coefficient l1 bound (13/4+2sqrt(2))^2 is <37.
    assert Q(297, 16) + 13*Q(283, 200) < 37

    buffered = one_root_lower(20, 1)
    unbuffered = one_root_lower(0, 10)
    assert buffered == 4024054722670646630717842689
    assert unbuffered == 7239337205339735179200724088
    assert 250*buffered > scale
    assert 20000*unbuffered**2 > scale**2

    # Conditional innovation covariance and coupling bounds for zeta<=1e-5.
    zmax = Q(1, 100000)
    assert 10*(zmax+40*zmax*zmax) <= Q(101, 10)*zmax
    assert Q(1, 12)-Q(101, 10)*zmax > Q(7, 25)**2
    assert Q(101, 10)/Q(14, 25) < 19
    assert Q(21, 4)*19**2 < 1896
    assert 37 < Q(61, 10)**2 and 1896 < 44**2
    assert zmax < Q(1, 316)**2
    assert 2*Q(61, 10)*44/Q(316)+1896*zmax < 3

    q, theta, zeta = Q(3, 100000), Q(1, 100), Q(1, 700000)
    eta, d0 = Q(99, 100)*zeta, Q(159, 1000)
    h = Q(1, 20)-theta/2
    assert h*h/(80*zeta) == Q(567, 32)
    # exp(567/32)>(79/48)^35*(39/32)>40,000,000.
    assert 39*79**35 > 40000000*32*48**35
    assert Q(1, 62500)-Q(1, 1000000) == q/2
    assert zeta/100+eta == zeta

    assert (1-theta*theta)*Q(7, 44) > d0
    assert Q(99, 70)**2 > 2
    assert (1-Q(2, 3)*Q(99, 70))/2 == Q(1, 35)
    assert 1-theta*theta/2 > Q(999, 1000)
    A = eta*q*d0
    B = 3*q*Q(1, 35)*d0*Q(999, 1000)*theta*theta/2
    alpha, delta = Q(1, 6370000), Q(1, 10**18)
    assert A == Q(47223, 7000000000000000)
    assert B == Q(1429569, 70000000000000000)
    assert A*alpha > delta and B > delta

    return {
        'certified': True,
        'scale': str(scale),
        'buffered_one_root_lower_numerator': str(buffered),
        'unbuffered_one_root_lower_numerator': str(unbuffered),
        'q': str(q),
        'theta': str(theta),
        'zeta_c': str(zeta),
        'initial_covariance_error_max': str(zeta/100),
        'eta': str(eta),
        'd0_lower': str(d0),
        'A_lower': str(A),
        'B_star_lower': str(B),
        'A_times_alpha_lower': str(A*alpha),
        'structural_delta': str(delta),
    }

def certify_buffered_gain():
    """Exact arithmetic for the buffered gain and approved structural contract.

    It uses a local standard-library import, has no output or other side effects,
    and returns a JSON-compatible certificate. The accompanying analytic
    subsection proves the density, coupling, and graph inequalities whose
    rational arithmetic is checked here.
    """
    from fractions import Fraction as F

    # Critical root/innovation box and all fourteen sign margins.
    slo, shi = F(7, 5), F(99, 70)
    assert slo*slo < 2 < shi*shi
    vlo, vhi = F(1, 5), F(21, 100)
    wlo, whi = F(103, 100), F(53, 50)
    xiv = F(1, 100)
    xiwlo, xiwhi = -F(33, 20), -F(8, 5)
    xiab, xicd, margin = F(9, 50), F(7, 25), F(1, 5)
    abhi = shi*vhi-wlo/2+xiv
    chi = shi*whi-vlo/2+xiwhi
    dlo = slo*wlo-vhi/2-xiwhi
    assert abhi == -F(26, 125) < -margin
    assert chi == -F(703, 3500) < -margin
    assert dlo == F(2937, 1000) > margin
    assert vlo >= margin and wlo > margin
    assert slo*abhi-vlo/2+xiab == -F(132, 625) < -margin
    assert slo*chi-wlo/2+xicd == -F(2581, 5000) < -margin
    assert slo*dlo-whi/2-xicd == F(16509, 5000) > margin

    # pi<10/3 follows by integrating 1/(1+x^2)<=1-x^2/2.
    assert F(3, 2)/F(10, 3) == F(9, 20)
    assert 6/F(10, 3) > F(4, 3)**2
    root_q = vhi*vhi/2+F(9, 2)*(whi-F(14, 15)*vlo)**2
    assert whi-F(14, 15)*vlo == F(131, 150)
    assert root_q == F(13817, 4000) < F(7, 2)
    total_q = F(7, 2)+6*(xiv*xiv+xiwlo*xiwlo+2*xiab**2+2*xicd**2)
    assert total_q == F(52913, 2500) < F(43, 2)
    pref = ((vhi-vlo)*(whi-wlo)*F(9, 20)*(2*xiv)*(xiwhi-xiwlo)
            *(2*xiab)**2*(2*xicd)**2*F(4, 3)**6)
    assert pref == F(4704, 152587890625)
    p = F(1, 10**17)

    # 8/3<e<11/4 and sqrt(e)<5/3, from the exponential series.
    assert F(8, 3)+F(1, 24)/(1-F(1, 5)) == F(87, 32) < F(11, 4)
    assert F(11, 4) < F(5, 3)**2
    assert F(11, 4)**21*F(5, 3) < 3*10**9
    box_probability_lower = pref/(3*10**9)
    assert box_probability_lower > p

    # Singular square-root coupling, fourteen coordinates, error margin 19/100.
    theta, zeta = F(1, 100), F(1, 35000)
    h = margin-theta
    assert h == F(19, 100)
    assert h*h/(2*14*zeta) == F(361, 8) > 45
    coupling_failure_upper = 28*F(3, 8)**45
    assert coupling_failure_upper < p/2
    assert zeta < (1-F(2, 3)*shi)/4

    # Diameter-five Hellinger bound, orientation multiplicity, and tradeoff.
    for distance in range(1, 6):
        assert distance*distance*2**(distance-1) <= 20**2
    assert 4*2*5 < 7**2
    assert 2+4+8 == 14 and 1+4+8 == 13 and 2+8+16 == 26
    tau = zeta/21
    alpha = 3*tau/26
    assert 20*tau+zeta/21 == zeta
    assert tau == F(1, 735000)
    assert alpha == F(1, 6370000)

    # Approved companion curvature contract and primary structural loss.
    q, zeta_c, d0_lower = F(3, 100000), F(1, 700000), F(159, 1000)
    initial_error = zeta_c/100
    eta = 99*zeta_c/100
    assert eta == F(99, 70000000)
    assert initial_error+eta == zeta_c
    assert (1-theta*theta)*F(7, 44) > d0_lower
    assert (1-F(2, 3)*shi)/2 == F(1, 35)
    assert 1-theta*theta/2 > F(999, 1000)
    A_lower = eta*q*d0_lower
    B_star_lower = F(3, 2)*q*F(1, 35)*d0_lower*F(999, 1000)*theta*theta
    epsilon0 = F(1, 10**18)
    assert A_lower == F(47223, 7000000000000000)
    assert A_lower*alpha == F(47223, 44590000000000000000000)
    assert B_star_lower == F(1429569, 70000000000000000)
    assert min(A_lower*alpha, p/8, B_star_lower, F(1, 100)) > epsilon0

    return {
        'certified': True,
        'arithmetic': 'exact integers and fractions.Fraction',
        'p': str(p),
        'critical_sign_margin': str(margin),
        'box_probability_lower': str(box_probability_lower),
        'theta': str(theta),
        'zeta_g': str(zeta),
        'coupling_failure_upper': str(coupling_failure_upper),
        'window_vertices': 14,
        'window_edges': 13,
        'window_diameter': 5,
        'hellinger_coefficient': 20,
        'oriented_window_multiplicity': 26,
        'tau': str(tau),
        'alpha': str(alpha),
        'initial_curvature_error': str(initial_error),
        'eta': str(eta),
        'A_lower': str(A_lower),
        'A_times_alpha_lower': str(A_lower*alpha),
        'gain_loss_lower': str(p/8),
        'B_star_lower': str(B_star_lower),
        'epsilon0': str(epsilon0),
    }

def certify_structural_constants():
 from fractions import Fraction as F
 curvature=certify_buffered_curvature()
 gain=certify_buffered_gain()
 base=F(130647566423167664341749381,10**27)
 delta=F(1,10**18)
 gamma=base-delta
 assert gamma==F(130647566423167663341749381,10**27)>0
 c0=1+1/(2*gamma)
 target=1+1/(2*(gamma+F(1,10**27)))
 previous=F(315323783211583832170874691,65323783211583832170874691)
 assert c0==F(630647566423167663341749381,130647566423167663341749381)
 assert target==F(315323783211583831670874691,65323783211583831670874691)
 assert c0>target>previous
 assert c0-target==F(250000000000000000000000000,8534393306148003315042418524038500659446754377816271)
 assert F(curvature['A_lower'])==F(gain['A_lower'])
 assert F(curvature['B_star_lower'])==F(gain['B_star_lower'])
 assert F(curvature['structural_delta'])==F(gain['epsilon0'])==delta
 return {'certified':True,'curvature':curvature,'gain':gain,
         'gamma_base':str(base),'structural_delta':str(delta),
         'gamma_new':str(gamma),'intermediate_coefficient':str(c0),
         'eventual_coefficient':str(target),'positive_slack':str(c0-target),
         'gain_over_previous':str(target-previous)}

def main():
 t=time.monotonic();structural=certify_structural_constants();N=80;J=65;Ln=3;Ld=1;L=I(Ln)/Ld
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
 print(json.dumps({'certified':passed,'structural_certificate':structural,'target':{'numerator':130647566423167664341749381,'denominator':1000000000000000000000000000},'gamma_upper_bound':str(gamma.b),'N':N,'J':J,'L':{'numerator':Ln,'denominator':Ld},'positive_quadrature_norm':qhalf.bounds(),'positive_operator_bound':h.bounds(),'completion_constant':kappa.bounds(),'exact_p1':exact_p1.bounds(),'endpoint_tail':endpoint_tail.bounds(),'middle_tail':middle_tail.bounds(),'state_rescaling_base':4,'error_bound':'localized-gaussian-moments-triple-outliers-two-step-completions','precision':P,'fixed_point_bits':FP,'intervals':results,'seconds':time.monotonic()-t},indent=2))
 return 0 if passed else 1
if __name__=='__main__':sys.exit(main())
