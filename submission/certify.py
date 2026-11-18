#!/usr/bin/env python3
"""Certified alternating orthants of the critical cubic vertex Gaussian wave.
Directed Decimal initialization, exact fixed-point contractions, and analytic
Gaussian-moment Gauss-Legendre error bounds.
All quadrature, truncation and rounding errors are enclosed explicitly.
"""
import decimal as dec, math, json, sys, time
P=100
dec.getcontext().prec=P
lo=dec.Context(prec=P,rounding=dec.ROUND_FLOOR)
hi=dec.Context(prec=P,rounding=dec.ROUND_CEILING)
mid=dec.Context(prec=P,rounding=dec.ROUND_HALF_EVEN)
D=dec.Decimal
# Positive fixed-point intervals represent [a/SCALE, b/SCALE].
# Python integers have unlimited precision. All contraction sums and products
# below are exact; floor/ceiling are applied only at the stated divisions.
FP=224
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

def exact_legendre(n, x):
 # At x=a/b store k! b**k P_k(x) as an integer.
 a,b=x.as_integer_ratio()
 previous, value = 1, a
 for k in range(2,n+1):
  previous,value=value,(2*k-1)*a*value-(k-1)**2*b*b*previous
 return previous, value, a, b

def nodes_weights(n,L):
 nodes=[];weights=[];brackets=[]
 for k in range(1,n+1):
  x=D(str(math.cos(math.pi*(k-.25)/(n+.5))))
  for iteration in range(7):
   previous, value = D(1), x
   for degree in range(2,n+1):
    previous, value = value, ((2*degree-1)*x*value-(degree-1)*previous)/degree
   derivative = n*(x*value-previous)/(x*x-1)
   x -= value/derivative
  eps=D('1e-72');br=I(lo.subtract(x,eps),hi.add(x,eps))
  left=exact_legendre(n,br.a)[1]
  right=exact_legendre(n,br.b)[1]
  assert left*right<0 and -1<br.a<br.b<1
  previous,value,a,b=exact_legendre(n,x)
  # P_n'(x)=n(x P_n(x)-P_(n-1)(x))/(x*x-1).
  dp=I(n*(a*value-n*b*b*previous)*b)/(math.factorial(n)*b**n*(a*a-b*b))
  width=max(hi.subtract(x,br.a),hi.subtract(br.b,x))
  error=hi.multiply(width,D(n**4))
  dp=dp+I(error.copy_negate(),error)
  assert not dp.a<=0<=dp.b
  w=2/((1-br*br)*dp*dp)
  nodes.append((br+1)*L/2);weights.append(w*L/2);brackets.append(br)
 assert all(brackets[i].a>brackets[i+1].b for i in range(n-1))
 return nodes,weights

def pivotal_certificate():
    from fractions import Fraction as F
    from math import factorial


    def half_density(coefficient, rate, value):
        # exp(-z) >= (1-z/64)^64 on 0 <= z <= 64.
        base = 1 - rate * value * value / 64
        assert 0 <= base <= 1
        return coefficient * base**64


    def half_d_cdf(value):
        # Integrate the odd Taylor lower polynomial for exp(-3*x^2/2).
        value = min(value, F(5, 4))
        assert value > 0
        rate = F(-3, 2)
        integral = sum(rate**i * value**(2*i+1) / (factorial(i)*(2*i+1))
                       for i in range(8))
        assert integral > 0
        return F(11, 8) * integral


    width = F(1, 100)
    one_side = F(0)
    one_prescribed_sign = F(0)
    for ai in range(4, 60):
        a = ai * width
        one_neighbor = F(0)
        for ki in range(4, 80):
            k = ki * width
            largest_contrast = 2*(k+width) + F(2829, 1000)*(a+width)
            smallest_child_sum = F(2828, 1000)*k + a - F(8, 100)
            density = half_density(F(11, 8), F(3, 2), largest_contrast)
            one_neighbor += 2*width*density*half_d_cdf(smallest_child_sum)
        amplitude_weight = width*half_density(F(11, 4), F(6), a+width)
        one_side += amplitude_weight*one_neighbor**2
        one_prescribed_sign += amplitude_weight*one_neighbor

    four_private = one_side**2
    two_private = one_prescribed_sign**2
    pivotal_floor = (two_private*10**12).__floor__()
    assert pivotal_floor == 18079200775 and two_private > F(18, 1000)
    return pivotal_floor


def check_good_centers():
    """Integer lower sum for the good-center flip at eta=.000425, |t|<=.67."""
    from fractions import Fraction as F


    def good_trigger_certificate():
        scale = 10**9

        def exp_lower(x):
            # For all boxes 0 <= x <= 23. Return a fixed-point lower bound.
            q = 1 - x/4096
            u = q.numerator*scale//q.denominator
            for _ in range(12):
                u = u*u//scale
            return u

        def rectangle(t):
            # A lower bound on the probability that |N(0,1)| <= t/scale.
            t = min(scale, max(0, t))
            return (797*t*(6*scale**2-t*t))//(6000*scale**2)

        wx = [exp_lower(F((670+5*(i+1))**2, 2*10**6)) for i in range(40)]
        wy = [exp_lower(F(1035, 2000)*F((j+1)**2, 40**2)) for j in range(140)]
        wz = [exp_lower(F(1030, 2000)*F((k+1)**2, 40**2)) for k in range(260)]
        # Upper bounds on the conditional means, in units of 1/10000.
        means = [(9433, -1610, 172), (8358, -2196, 368),
                 (8358, 4730, -2810), (7148, 5055, -3920)]
        # Defect tolerances r and lower bounds on each two-sided normal mass.
        tolerances = [(17_050_000, 748), (23_700_000, 889), (31_850_000, 967)]
        total = 0
        for i in range(40):
            x = 5_000_000*(i+1)
            for j in range(140):
                y = 25_000_000*j
                for k in range(260):
                    z = 25_000_000*k
                    A, B, C, D = [-(p*x + q*(y+25_000_000*(q>0))
                                    + s*(z+25_000_000*(s>0)))//10_000
                                  for p, q, s in means]
                    if min(A, B, C, D) <= 0:
                        continue
                    best = 0
                    for r, m in tolerances:
                        if 1000*C<=972*r or 1000*D<=2043*r:
                            continue
                        third = (r+2)//3
                        half = (r+1)//2
                        first_bound = min(A-third, (B-half)*1000//1415-third)
                        d3 = D-(2043*r+999)//1000
                        if first_bound<=0:
                            continue
                        for num in (10, 16):
                            d1 = first_bound*num//20
                            d2 = B-half-(1415*(d1+third)+999)//1000
                            term = rectangle(346*d1//100)*rectangle(346*d2//100)
                            best = max(best, term*rectangle(346*d3//100)*m**4)
                    total += wx[i]*wy[j]*wz[k]*best
        # Root and both tilt-coordinate densities, then the uniform cell volume.
        return F(total, scale**6*1000**4)*F(615, 10000)*F(1, 200*40*40)
    value = good_trigger_certificate()
    assert value > F('0.00000002') * 15
    return (10**15 * value).__floor__()

def check_edge_credit():
    from fractions import Fraction as F
    from functools import lru_cache
    from math import factorial, isqrt

    UNIT=10**18
    H=F(1,100)
    ERROR=F(1,10000)


    def down(x):
        return max(0,x.numerator*UNIT//x.denominator)


    @lru_cache(None)
    def exp_lower(x):
        # The odd Taylor polynomial underestimates exp(-z), including on z>=0.
        z=x/64
        y=down(sum((-z)**k/factorial(k) for k in range(8)))
        for _ in range(6):
            y=y*y//UNIT
        return y


    @lru_cache(None)
    def normal_integral(x,rate):
        # Lower bound for integral_0^x exp(-rate*t^2).
        if x<=0:
            return 0
        return sum((-rate)**i*x**(2*i+1)/(factorial(i)*(2*i+1))
                   for i in range(20))


    @lru_cache(None)
    def half_cdf(x):
        return down(F(11,8)*normal_integral(min(x,F(5,4)),F(3,2)))


    def tail_upper(m,s,sd):
        a=(m-ERROR)**2*(1-s)/(s*sd**2)
        t=F(isqrt((a.numerator*1000**2)//a.denominator),1000)
        return F(1,2)-F(3988,10000)*normal_integral(min(t,F(3)),F(1,2))


    def sides(m1,m2,m3):
        total=0
        for ai in range(m1,65):
            x=H*ai
            a=11*exp_lower(6*(x+H)**2)//4
            subtotal=0
            for ki in range(m2,85):
                k=H*ki
                arg=2*(k+H)+F(2829,1000)*(x+H)
                g=11*exp_lower(F(3,2)*arg**2)//8
                v=half_cdf(F(2828,1000)*k+x-2*m3*H)
                subtotal+=2*g*v//(100*UNIT)
            total+=a*subtotal//(100*UNIT)
        return F(total,UNIT)**2


    # Each tuple specifies an interpolation endpoint in units of .00001,
    # and wave score margins at depths 1,2,3 in units of .01.
    cells=[(10,2,2,4),(25,2,4,4),(50,2,6,6),(75,4,6,6),
           (100,4,6,8),(150,4,8,8),(200,4,8,10),(300,4,10,12),
           (400,6,10,12),(600,6,14,16),(800,8,14,18),(1000,8,16,20),
           (1500,8,20,24),(2000,10,20,26),(3000,12,24,32)]
    tilted=True
    if tilted:
        cells=[(50,2,2,2),(100,2,2,4),(200,2,4,4),(300,2,4,4),
            (500,2,4,6),(800,2,6,6),(1200,2,6,8),(1800,2,8,10),
            (2500,4,8,10),(3500,4,10,12),(5000,4,10,14),(7500,4,14,18),
            (10000,6,14,18),(14000,6,16,22),(18000,6,20,26),
            (24000,8,22,30),(32000,10,24,34),(42000,12,30,42)]
    last=0
    integral=F(0)
    for end,m1,m2,m3 in cells:
        s=F(end,100000)
        p=sides(m1,m2,m3)-F(1,10**9)
        stds=(86,208,336) if tilted else (501,816,1032)
        safe=1-4*sum(tail_upper(x*H,s,F(sd,1000))
                     for x,sd in zip((m1,m2,m3),stds))
        if tilted:
            safe*=(1-s)**3
        assert p>0 and safe>0
        integral+=(s-last)*p*safe
        last=s

    # The lower pair density .47 is justified by limiting rho=2sqrt(2)/3.
    # Integrate two sign-reflected certificates. If both hold, mixed difference
    # is -2; otherwise a single successful certificate implies at most -1.
    heavy=2*F(47,100)*integral
    edge_floor=down(heavy)
    assert edge_floor == 672905379119700
    return edge_floor

def main():
 t=time.monotonic();N=int(sys.argv[1]) if len(sys.argv)>1 else 93;J=int(sys.argv[2]) if len(sys.argv)>2 else 47;L=3
 assert N>=1 and J>=3 and J%2==1
 pi=pi_interval();r2=I(2).sqrt();c2=3/(2*pi);c1=(6/pi).sqrt()
 xs,ws=nodes_weights(N,L)
 remote_mass = sum((weight * (-I(39)/2 * x.square()).exp() for x, weight in zip(xs, ws)), I(0))
 assert remote_mass.b < D('0.25')
 print('nodes certified',N,'seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 # The Gaussian exponents use two symmetric products. Reuse them in both
 # densities, and square the sqrt(2) cross factor for the y,z coefficient.
 def symmetric_exponential(rate):
  values = [[None]*N for _ in xs]
  for i,x in enumerate(xs):
   for k in range(i+1):
    entry = (-rate*x*xs[k]).exp()
    values[i][k] = values[k][i] = entry
  return values

 cross_root = symmetric_exponential(6*r2)
 cross_plain = symmetric_exponential(I(6))
 pair_diagonal = [(-I(9)/2*x.square()).exp()*w for x,w in zip(xs,ws)]
 x_factor = [(-I(3)/2*x.square()).exp() for x in xs]
 y_factor = [(-12*x.square()).exp() for x in xs]
 z_factor = [c1*w*(-6*x.square()).exp() for x,w in zip(xs,ws)]

 def to_fixed(interval):
  return fixed(interval.a, interval.b)

 state = [[to_fixed(c2*left*right*cross_root[i][j])
           for j,right in enumerate(pair_diagonal)]
          for i,left in enumerate(pair_diagonal)]
 A = [[to_fixed(x_factor[i]*cross_root[i][j]) for j in range(N)] for i in range(N)]
 Bcols = [[to_fixed(cross_plain[i][k]) for i in range(N)] for k in range(N)]
 assert all(hi_b<=SCALE for column in Bcols for _,hi_b in column)
 C = [[to_fixed(y_factor[j]*z_factor[k]*cross_root[j][k].square())
       for k in range(N)] for j in range(N)]
 print('kernel factors certified','seconds',time.monotonic()-t,flush=True,file=sys.stderr)
 results=[];gamma=I(0)
 for j in range(1,J+1):
  sa=sum(a for row in state for a,b in row)
  sb=sum(b for row in state for a,b in row)
  s=decimal_fixed(sa,sb)
  # Exact density at zero and the Gaussian even-moment derivative bound.
  # Precision diagonals are [9,9] at n=2, [12,33,12] at n=3,
  # and [12,36,39,...,39,36,12] at n>=4. The tensor-error proof
  # sums their Nth powers, rather than using an empirical derivative estimate.
  n=j+1
  if n==2:qsum=2*9**N
  elif n==3:qsum=2*12**N+33**N
  else:qsum=2*12**N+2*36**N+(n-4)*39**N
  density=c2
  for _ in range(n-2):density=density*c1
  numer=L**(n+2*N)*math.factorial(N)**3*qsum
  denom=(2*N+1)*2**N*math.factorial(2*N)**2
  # Coordinates outside the two ends and the local precision band have
  # envelope exp(-39*x*x/2), with mass <1/4 for both positive operators.
  quad=(I(numer)/denom)*density / (4*L)**max(0,n-9)
  # Anchor at an adjacent positive pair containing an outlying coordinate.
  # Every other positive forward or reverse Markov step has mass <=1/2.
  # Retaining the pair-density cross term bounds its one-coordinate tail
  # by c2*exp(-9*L^2/2)/(54*sqrt(2)*L^2). A union bound over n
  # coordinates and the n-2 completion steps gives the expression below.
  # Reverse and forward arms have worst two-step completion mass < 1/10.
  def arm_denominator(length):
   return 10**(length//2)*2**(length%2)
  completion=I(0)
  for index in range(n):
   anchor=min(index,n-2)
   completion=completion+I(1)/(arm_denominator(anchor)*arm_denominator(n-2-anchor))
  tail=completion*c2*(-I(9)*L*L/2).exp()/(54*r2*L*L)
  p=2*(s+I(-quad.b,quad.b)+I(0,tail.b))
  coeff=I(3)*I(2**max(j-2,0))/(2 if j==1 else 1)
  gamma=gamma+((-1 if j%2==0 else 1)*coeff)*p
  results.append({'j': j, 'probability': p.bounds()})
  if j == J-1: previous_even = gamma.a
  if j%10==0: print('completed term',j,flush=True,file=sys.stderr)
  if j==J:break
  # If m/S and a/S enclose M,A, the product has denominator S^2.
  # Flooring/ceiling m*a/S restores endpoints with denominator S.
  nxt=[]
  for jj in range(N):
   active=[]
   skipped=0
   total_upper=0
   for ii in range(N):
    sl,su=state[ii][jj]; al,au=A[ii][jj]
    upper=(su*au+SCALE-1)//SCALE
    total_upper+=upper
    if upper<=8:
     skipped+=upper
    else:
     active.append((ii,sl*al//SCALE,upper))
   row=[]
   for kk in range(N):
    cl,cu=C[jj][kk]
    if total_upper*cu<=SCALE:
     row.append((0,1))
     continue
    # All B entries lie in [0,1]. Lost lower contributions are nonnegative;
    # each skipped scaled upper endpoint contributes at most itself times S.
    lower_sum=0
    upper_sum=skipped*SCALE
    column=Bcols[kk]
    for index,sl,su in active:
     bl,bu=column[index]
     lower_sum+=sl*bl
     upper_sum+=su*bu
    row.append((lower_sum*cl//SCALE2,(upper_sum*cu+SCALE2-1)//SCALE2))
   nxt.append(row)
  state=nxt
 lower_target = D('0.1306475664231676642')
 upper_target = D('0.130647566423167665')
 pivotal_floor = pivotal_certificate()
 good_floor = check_good_centers()
 edge_floor = check_edge_credit()
 from fractions import Fraction as F
 q = F(good_floor, 10**15)
 K = F(edge_floor, 10**18)
 assert F(4,5)**7 > F(1,5)
 assert min(q/5, F(99,100)*F(17,40000)*K/6) > F(47,10**9)
 passed = previous_even > lower_target and gamma.b < upper_target
 print(json.dumps({
  'certified': passed, 'N': N, 'J': J, 'L': L, 'precision': P, 'fixed_point_bits': FP,
  'gamma_odd_lower': str(gamma.a), 'gamma_odd_upper': str(gamma.b),
  'gamma_even_lower': str(previous_even), 'target_series_upper': str(upper_target),
  'structural_saving_proved_in_manuscript': '0.000000047', 'good_floor_scaled_1e15': good_floor, 'integrated_edge_floor_scaled_1e18': edge_floor, 'pivotal_floor_scaled_1e12': pivotal_floor,
  'intervals': results
 }))
 return 0 if passed else 1

if __name__ == '__main__':
 sys.exit(main())
