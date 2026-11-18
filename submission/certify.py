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


from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt

UNIT = 10**18
spatial = 200
time_steps = 80
def down(q):
    return max(0, q.numerator*UNIT//q.denominator)
def mul(a,b):
    return a*b//UNIT
def exp_lo(x):
    # The seventh Taylor polynomial at x/64 is below exp(-x/64).
    z=x/64
    ret=down(sum((-z)**j/factorial(j) for j in range(8)))
    for _ in range(6): ret=mul(ret,ret)
    return ret
@lru_cache(None)
def central_normal(i):
    # Lower bound on Pr(|N(0,1)| <= i/1000).
    t=F(min(i,3300),1000)
    integral=sum(F(-1,2)**j*t**(2*j+1)/(factorial(j)*(2*j+1))
                 for j in range(24))
    return down(F(7978,10000)*integral)
def strip_at(j, s, st, divisor=1):
    # h*j/(divisor*st*sqrt(r)), rounded downward to .001.
    a=F(j*1000,spatial*divisor)**2*(1-s)/(st**2*s)
    return central_normal(isqrt(a.numerator//a.denominator))


def coefficients():
    nx=(7*spatial)//10
    nk=(9*spatial)//10
    nl=nx+3*nk
    # Lower coefficients follow from pi < 3.1416.
    assert F(31416,10000)*F(27638,10000)**2 < 24
    assert F(31416,10000)*F(13819,10000)**2 < 6
    assert F(31416,10000)*F(7978,10000)**2 < 2
    az=[27638*exp_lo(6*F(i+1,spatial)**2)//10000 for i in range(nx)]
    dz=[13819*exp_lo(F(3,2)*F(i,spatial)**2)//10000 for i in range(nl)]
    cells=[]
    # Density is decreasing. Quantize the upper argument upwards to 1/1000
    # for efficient shared exponential evaluation.
    density={}
    for i in range(nx):
        row=[]
        for j in range(nk):
            dh=(200000*(j+1)+282843*(i+1))
            den=100000*spatial
            index=(1000*dh+den-1)//den
            if index not in density:
                density[index]=13819*exp_lo(F(3,2)*F(index,1000)**2)//10000
            lower_l=(282842*j+100000*i)//100000
            # Two exterior cell widths and the k change-of-variable factor.
            coeff=2*mul(az[i],density[index])//(spatial**2)
            row.append((lower_l,coeff))
        cells.append(row)
    return nx,nk,nl,dz,cells


def pivotal_integral(stds=(F(86,1000),F(208,1000),F(336,1000))):
    nx,nk,nl,dz,cells=coefficients()
    old=F(0); integrated=F(0)
    for step in range(1,time_steps):
        s=F(step,time_steps)**2
        sx=[strip_at(i,s,stds[0]) for i in range(nx)]
        sk=[strip_at(j,s,stds[1]) for j in range(nk)]
        sz=[strip_at(t,s,stds[2],2) for t in range(nl)]
        # For l bins, the z upper coordinate is q/spatial; the remaining
        # next-generation amplitude is at least (l-q)/(2*spatial).
        conv=[sum(dz[q]*sz[l-q] for q in range(1,l+1))//(UNIT*spatial)
              for l in range(nl)]
        side=0
        for i,row in enumerate(cells):
            row_value=sum(mul(mul(c,sk[j]),conv[l]) for j,(l,c) in enumerate(row))
            side+=mul(sx[i],row_value)
        integrated+=(s-old)*F(side,UNIT)**2*(1-s)**3
        old=s
    return integrated



def sqrt_lower(x):
    return F(isqrt(x.numerator*UNIT**2//x.denominator), UNIT)

def strips(j,s,coef):
    # |N(0,1)| threshold j/spatial*coef/sqrt(r).
    r=F(s,1-s)
    arg=(F(j*1000,spatial)*coef)**2/r
    return central_normal(isqrt(arg.numerator//arg.denominator))

def double_path():
    nx=(8*spatial)//10
    # |N(0,1/33)| and |N(0,1/3)|.
    assert F(45834,10000)**2*F(31416,10000)<66
    assert F(13819,10000)**2*F(31416,10000)<6
    xa=[45834*exp_lo(F(33,2)*F(i+1,spatial)**2)//10000 for i in range(nx)]
    dp=[13819*exp_lo(F(3,2)*F(j,spatial)**2)//10000 for j in range(nx)]
    integral=F(0)
    old=0
    for step in range(1,time_steps):
        s=F(step,time_steps)**2
        # c=3-2sqrt(2) < .17158; lambda > 2.82842.
        sx=[strips(i,s,F(282842,17158)) for i in range(nx)]
        sd=[strips(i,s,F(100000,17158)) for i in range(nx)]
        total=0
        for i in range(nx):
            conv=sum(dp[j]*sd[i-j-1] for j in range(i))//(UNIT*spatial)
            total+=mul(mul(xa[i],sx[i]),mul(conv,conv))//spatial
        r=s/(1-s)
        # Variance of the one root coordinate is <3.91; contrasts at most 1.
        determinant=sqrt_lower(1/(1+F(391,100)*r))*(1-s)
        integral+=(s-old)*F(total,UNIT)*determinant
        old=s
    return integral


def all_center_trigger():
    S = 10**9

    def exp_lower(x):
        v = 1-x/4096
        v = v.numerator*S//v.denominator
        for _ in range(12):
            v = v*v//S
        return v

    def normal_strip(v):
        v = min(16*S//10, max(0, v))
        poly = 1680*S**6-280*v**2*S**4+42*v**4*S**2-5*v**6
        return 797*v*poly//(1_680_000*S**6)

    wx = [exp_lower((F(11,200)+F(i+1,250))**2/2) for i in range(50)]
    wy = [exp_lower(F(1038,2000)*F((j+1)**2,50**2)) for j in range(180)]
    wz = [exp_lower(F(1033,2000)*F((k+1)**2,50**2)) for k in range(340)]
    rows = [(94329,-16039,1921,315), (83603,-21780,4113,917),
            (83603,47316,-28012,917), (71567,50606,-39010,1611),
            (-94280,-32705,522,315), (-83332,-45350,-28012,917)]
    rs = [(18_000_000,748), (25_040_000,889), (33_650_000,967)]
    # Tabulate the normal strip on a 1/10000 grid. Monotonicity allows
    # rounding its nonnegative argument down before indexing.
    strip_table = [normal_strip(i*100_000) for i in range(16_001)]

    def strip_at(value):
        return strip_table[min(16_000, max(0, value//100_000))]

    total = 0
    for i in range(50):
        x = 4_000_000*i
        for j in range(180):
            y = 20_000_000*j
            partial = [(p*(x+4_000_000*(p>0))+q*(y+20_000_000*(q>0))+h*S,c)
                       for p,q,c,h in rows]
            # For fixed x,y each worst-endpoint margin numerator is
            # affine in k. Restrict to the intersection of six positive
            # integer-margin ranges before evaluating the rectangles.
            lines = [(-base-c*20_000_000*(c>0), -c*20_000_000) for base,c in partial]
            begin, end = 0, 339
            for intercept, slope in lines:
                if slope > 0:
                    begin = max(begin, -((intercept-100_000)//slope))
                elif slope < 0:
                    end = min(end, (intercept-100_000)//(-slope))
                elif intercept < 100_000:
                    end = -1
            for k in range(begin, end+1):
                A,B,C,D,M,N = [(intercept+slope*k)//100_000 for intercept,slope in lines]
                best = 0
                for r,m in rs:
                    if min(C,N)*1000<=972*r or 3*M<=r or 1000*D<=2043*r:
                        continue
                    third,half = (r+2)//3,(r+1)//2
                    d1max = min(A-third,(B-half)*1000//1415-third)
                    d3 = D-(2043*r+999)//1000
                    if d1max<=0:
                        continue
                    common = strip_at(346*d3//100)*m**4
                    for f in (5,8):
                        d1 = d1max*f//10
                        d2 = B-half-(1415*(d1+third)+999)//1000
                        prod = common*strip_at(346*d1//100)*strip_at(346*d2//100)
                        best = max(best,prod)
                total += wx[i]*wy[j]*wz[k]*best
    return F(total,S**6*1000**4)*F(621,10_000)*F(1,250*50*50)



def main():
 t=time.monotonic();N=int(sys.argv[1]) if len(sys.argv)>1 else 75;J=int(sys.argv[2]) if len(sys.argv)>2 else 27;L=3
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
 lower_target = D('0.130647566419')
 upper_target = D('0.130647566425')
 pivotal_floor = pivotal_certificate()
 q = all_center_trigger()
 first = F(954,1000)*pivotal_integral()
 assert F(2878,10000)**2*11*F(31416,10000)**2 < 9
 second = F(2878,10000)*double_path()
 K = first + F(565684,100000)*second
 assert int(q*10**15)==224201878 and int(6*q*10**15)==1345211272
 assert down(first)==3607233186504135 and down(second)==2211584760808344
 assert K >= F('0.016')
 delta = F('0.00000125')
 assert min(6*q,K*F(19,40000)/6)>delta
 outer = 3*F('.39')*F('.056')*F('1.1')*F('.0179')*F('.055')**2/2
 assert outer > delta
 structural_gamma = F(upper_target)-delta
 circuit_limit = 1+1/(2*structural_gamma)
 assert F((circuit_limit*10**32).__floor__(),10**32) == F('96542533104947432504697194718477/20000000000000000000000000000000') < circuit_limit
 passed = previous_even > lower_target and gamma.b < upper_target
 print(json.dumps({
  'certified': passed, 'N': N, 'J': J, 'L': L, 'precision': P, 'fixed_point_bits': FP,
  'gamma_odd_lower': str(gamma.a), 'gamma_odd_upper': str(gamma.b),
  'gamma_even_lower': str(previous_even), 'target_series_upper': str(upper_target),
  'structural_saving_proved_in_manuscript':'0.00000125', 'pivotal_floor_scaled_1e12':pivotal_floor,
  'one_orientation_floor_scaled_1e15':int(q*10**15), 'six_orientation_floor_scaled_1e15':int(6*q*10**15),
  'adjacent_credit_floor_scaled_1e18':down(first), 'distance_two_credit_floor_scaled_1e18':down(second),
  'combined_credit_floor_scaled_1e18':down(K), 'outer_floor_scaled_1e15':int(outer*10**15),
  'intervals': results
 }))
 return 0 if passed else 1

if __name__ == '__main__':
 sys.exit(main())
