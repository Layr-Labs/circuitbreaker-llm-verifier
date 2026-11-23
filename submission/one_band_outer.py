"""Distance-seven and outer threshold arithmetic."""
from types import SimpleNamespace


"""Standalone exact D7 contraction with regenerated size-zero/two arm inputs."""
from fractions import Fraction as F
from functools import lru_cache
from bisect import bisect_right
from hashlib import sha256
from math import factorial,isqrt,lcm
from itertools import permutations,product
import json

U,N,SPATIAL=10**18,36,80
STEPS=tuple(range(1,24))
FAC7=factorial(7)
C_HI,Q_HI=F(8579,50000),F(70711,100000)
SQRT2_LO,SQRT2_HI=F(141421,100000),F(70711,50000)
ARMS_DIGEST="8e900c99d144cf8f46deb165bf2844b6d47735fd5806b27632014ce8df9b6ae9"


def digest(value):
  return sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def exp_lower(x):
  assert 0<=x<=64
  n,d=x.numerator,64*x.denominator
  value=0
  for j in range(7,-1,-1):
    value=value*(-n)+(FAC7//factorial(j))*d**(7-j)
  value=max(0,U*value//(FAC7*d**7))
  for _ in range(6):value=value*value//U
  return value


CENTRAL_DEN=lcm(*[2**j*factorial(j)*(2*j+1)for j in range(24)])
CENTRAL_COEFF=[(-1)**j*(CENTRAL_DEN//(2**j*factorial(j)*(2*j+1)))*1000**(2*(23-j))for j in range(24)]
CENTRAL_DIV=10000*CENTRAL_DEN*1000**47


@lru_cache(None)
def central_strip(index):
  t=min(index,3300);assert t>=0;square=t*t;value=0
  for c in reversed(CENTRAL_COEFF):value=value*square+c
  return max(0,7978*U*t*value//CENTRAL_DIV)


def strip(index,s,sd,divisor=1):
  value=index*index*1000**2*(s.denominator-s.numerator)*sd.denominator**2//((SPATIAL*divisor)**2*s.numerator*sd.numerator**2)
  return central_strip(isqrt(value))


def scalar_controls():
  def atan(x,last):
    return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
  pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
  pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
  assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
  assert SQRT2_LO**2<2<SQRT2_HI**2 and Q_HI**2>F(1,2)
  assert C_HI==3-2*SQRT2_LO
  assert pi_hi*F(27638,10000)**2<24
  assert pi_hi*F(13819,10000)**2<6
  assert pi_hi*F(7978,10000)**2<2
  assert 0<exp_lower(F(1))<3*U//8
  assert central_strip(0)==0 and 0<central_strip(1000)<central_strip(2000)<U
  return pi_hi


def arm_inputs():
  """Recompute all 23 pairs of length-64 arrays, without a receipt input."""
  cells,central=48,64
  def fd(x):return 13819*exp_lower(F(3,2)*x*x)//10000
  @lru_cache(None)
  def fa(index):return 27638*exp_lower(6*F(index,1000)**2)//10000
  def ceil_index(q):return (1000*q.numerator+q.denominator-1)//q.denominator
  roots=[[2*fd(F(x+1+2*(a+1),SPATIAL))//SPATIAL for a in range(cells)]for x in range(central)]
  first=[[fa(ceil_index((F(b+1)+SQRT2_HI*(a+1))/SPATIAL))//SPATIAL for b in range(cells)]for a in range(cells)]
  lower_m=[[int(2*SQRT2_LO*b+a)for b in range(cells)]for a in range(cells)]
  zmax=max(central-1,max(map(max,lower_m)))
  dz=[0]+[fd(F(z,SPATIAL))//SPATIAL for z in range(1,zmax+1)]
  ell=[F(0)]+[C_HI/2*sum(((j+1)*Q_HI**j for j in range(d)),F(0))for d in range(1,4)]
  rows=[]
  for step in STEPS:
    s=F(step,24)**2
    root_strip=[strip(a,s,ell[1])for a in range(cells)]
    child_strip=[strip(b,s,ell[2])for b in range(cells)]
    terminal_strip=[strip(z,s,ell[3],2)for z in range(zmax+1)]
    base_strip=[strip(z,s,C_HI)for z in range(zmax+1)]
    terminal=[sum(dz[z]*terminal_strip[m-z]for z in range(1,m+1))//U for m in range(zmax+1)]
    baseline=[sum(dz[z]*base_strip[x-z]for z in range(1,x+1))//U for x in range(central)]
    inner=[sum((first[a][b]*child_strip[b]//U)*terminal[lower_m[a][b]]for b in range(cells))//U for a in range(cells)]
    arms=[sum((roots[x][a]*root_strip[a]//U)*inner[a]for a in range(cells))//U for x in range(central)]
    assert all(a<=b for a,b in zip(baseline,baseline[1:]))
    assert all(a>=b for a,b in zip(arms,arms[1:]))
    assert all(0<=v<=U for v in baseline+arms)
    rows.append({'step':step,'s':[s.numerator,s.denominator],'baseline_scaled':baseline,'size2_arm_scaled':arms})
  assert digest(rows)==ARMS_DIGEST
  return rows


EXP_COEFF=[FAC7//factorial(j)*64000**(7-j)for j in range(8)]
EXP_DEN=FAC7*64000**7


@lru_cache(None)
def exp_index(index):
  assert 0<=index<=64000
  value=0
  for c in reversed(EXP_COEFF):value=value*(-index)+c
  value=max(0,U*value//EXP_DEN)
  for _ in range(6):value=value*value//U
  return value


def margin(a,b,c=0):return (282842*b+100000*(a+c))//100000


def geometry():
  initial=[]
  for i in range(N):
    row=[];I=2*i+1
    for j in range(N):
      J=2*j+1;n=1800000*I*I+1950000*J*J+2545596*I*J
      row.append((exp_index((n+1250000+2559999)//2560000),margin(j,i)))
    initial.append(row)
  inner=[];last=[];maximum=0
  for b in range(N):
    B=2*b+1
    for c in range(N):
      C=2*c+1;middle=[];final=[]
      for a in range(N):
        A=2*a+1;m=margin(a,b,c);maximum=max(maximum,m)
        n=1950000*C*C+2545596*B*C+600000*A*C
        middle.append((exp_index((n+650000+2559999)//2560000),m))
        n=1800000*C*C+2545596*B*C+600000*A*C
        final.append((exp_index((n+600000+2559999)//2560000),m))
      inner.append(middle);last.append(final)
  return initial,inner,last,maximum


from math import comb
from operator import mul

PARAMETERS = ((7, F(1997, 50), F(333, 2000), F(1129, 1000), F(1, 2500000), 532), (8, F(281, 5), F(327, 2000), F(553, 500), F(19, 250000000), 1040), (9, F(783, 10), F(101, 625), F(1091, 1000), F(7, 500000000), 1968), (10, F(546, 5), F(803, 5000), F(537, 500), F(53, 20000000000), 3620))
PATHS = PARAMETERS[1:]
PANELS = (2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 23)

def translate(polynomial, center):
  return tuple(sum(a * comb(i, j) * center**(i-j)
                   for i, a in enumerate(polynomial) if i >= j)
               for j in range(len(polynomial)))

def longer_controls():
  pi_upper = scalar_controls()
  p8, p9, p10 = (24, -104, 70, -15, 1), (112, -200, 96, -17, 1), (-48, 320, -340, 126, -19, 1)
  assert translate(p8, 8) == (88, 184, 94, 17, 1)
  assert translate(p9, 8) == (48, 120, 72, 15, 1)
  assert translate(p10, 8) == (208, 640, 508, 158, 21, 1)
  assert 4 * 94 > 17**2 and 4 * 72 > 15**2
  assert translate((508, 158, 21, 1), -8) == (76, 14, -3, 1) and 4 * 14 > 3**2
  center = F(9, 2)
  negative_side_bound = 3 * sum(abs(a)*center**i for i, a in enumerate(translate(p9, center)))
  assert negative_side_bound < 1920*2*SQRT2_LO
  assert 384 * F(11, 48) == 88 and 1536 * F(13, 96) == 208
  assert 768*F(1, 8) == 48*2
  assert 736 * SQRT2_LO > PATHS[0][-1]
  assert 48 + 120 * 16 == PATHS[1][-1]
  assert 2560 * SQRT2_LO > PATHS[2][-1]
  scale = 10**12
  sqrt10_lower = F(isqrt(10 * scale**2), scale)
  assert sqrt10_lower**2 < 10
  theta_upper = 38 - 12 * sqrt10_lower
  determinants = (F(263, 288), 1-F(11, 48)**2, F(31, 32), 1-F(13, 96)**2)
  correlation_upper = (5*SQRT2_HI/24, F(11, 48), SQRT2_HI/8, F(13, 96))
  assert 7*8**3-65*8**2+144*8-44 == PARAMETERS[0][-1]
  assert 280*2*SQRT2_LO-708 > 0
  for params, endpoints_det, w in zip(PARAMETERS, determinants, correlation_upper, strict=True):
    distance, normalizer, density, mean, _, multiplier = params
    middle = distance-1
    middle_det = 9 * 12**middle * endpoints_det
    assert normalizer**2 * (2*pi_upper)**middle < 4 * middle_det
    assert (2*pi_upper*density)**2 * endpoints_det < 1
    assert 1 + F(2*(middle//2), middle)*theta_upper + 2*w/(middle*(1+w)) < mean


def incidence_time_moment(left, right, modes, mean):
  a, b = 1-left, 1-right
  degree = (modes+4)//2
  den = 1+(mean-1)*right
  if modes % 2:
    lo = F(isqrt(a.numerator*U**2//a.denominator), U)
    hi = F(isqrt(b.numerator*U**2//b.denominator)+1, U)
    moment = (a**degree*lo-b**degree*hi)/F(2*degree+1, 2)
    reciprocal = 1/den**modes
    determinant = F(isqrt(reciprocal.numerator*U**2//reciprocal.denominator), U)
  else:
    moment = (a**degree-b**degree)/degree
    determinant = 1/den**(modes//2)
  tangent = 1+F(modes, 4)*(mean-1)*(right-left)/den
  assert moment > 0
  return moment*determinant*tangent

def advance(state, weights):
  return [[sum(map(mul, column, coefficients))//U for coefficients in row]
          for column, row in zip(zip(*state, strict=True), weights, strict=True)]

def panel_masses(s, arm, shape):
  initial, interior, last, maximum = shape
  sides = [b+(s.denominator-s.numerator)*a//s.denominator
           for b, a in zip(arm["baseline_scaled"][:N], arm["size2_arm_scaled"][:N], strict=True)]
  strips = [strip(m, s, C_HI) for m in range(maximum+1)]
  inner, end = (
      [[[ex*strips[m]//U//SPATIAL for ex, m in row] for row in factors[i:i+N]]
       for i in range(0, N*N, N)]
      for factors in (interior, last)
  )
  state = [[ex*sides[i]//U*strips[m]//U//SPATIAL**2 for ex, m in row]
           for i, row in enumerate(initial)]
  for _ in range(3):
    state = advance(state, inner)
  tail = [[strips[margin(i, j)] for j in range(N)] for i in range(N)]
  for distance in range(7, 11):
    final = advance(state, end)
    mass = sum(sum(v*y//U*z//U for v, y, z in zip(row, sides, t, strict=True))
               for row, t in zip(final, tail, strict=True))
    yield mass
    if distance < 10:
      state = advance(state, inner)

def certificate():
  longer_controls()
  arms, shape = arm_inputs(), geometry()
  totals = [F(0)]*len(PARAMETERS)
  previous = F(0)
  for step in PANELS:
    s = F(step, 24)**2
    for i, (p, mass) in enumerate(zip(PARAMETERS, panel_masses(s, arms[step-1], shape), strict=True)):
      distance, normalizer, density, mean, _, _ = p
      totals[i] += density*normalizer*F(mass, U)*incidence_time_moment(previous, s, distance-1, mean)
    previous = s
  rows = []
  for (distance, _, _, _, floor, multiplier), total in zip(PARAMETERS, totals, strict=True):
    assert total > floor
    rows.append({"distance": distance, "exact_raw_price_sum": str(total),
                 "usable_raw_price_floor": str(floor), "K_increment": str(floor*multiplier)})
  return {"scope": "exact D7–D10 shared Price contraction", "middle_cells": N, "time_steps": list(PANELS),
          "rows": rows,
          "global_K_increment": str(sum(F(row["K_increment"]) for row in rows))}


d7 = SimpleNamespace(certificate=certificate, scalar_controls=scalar_controls, C_HI=C_HI, Q_HI=Q_HI, SQRT2_HI=SQRT2_HI, SQRT2_LO=SQRT2_LO, U=U, exp_index=exp_index, exp_lower=exp_lower, margin=margin, strip=strip)


STEP, T = F(1, 100), F(9, 400)
LO, HI, ROOT3 = SQRT2_LO, SQRT2_HI, F(173205, 100000)
PHI = F(797, 2000)


def moving_controls():
  def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
  pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
  pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
  assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
  assert LO**2<2<HI**2 and ROOT3**2<3
  assert F(229,50)**2*pi_hi<66 and F(13819,10000)**2*pi_hi<6
  assert F(27638,10000)**2*pi_hi<24 and F(7978,10000)**2*pi_hi<2
  assert F(81027,10000)*pi_hi<18*LO
  endpoint,child,closing=F(49,10000),F(1,500),F(39,10000)
  assert 0<2-14*HI/11 and T*(2-14*LO/11)<endpoint
  assert 0<3-2*HI and T*(3-2*LO)/2<child
  assert T*(3-2*LO)<closing
  assert max(T*(8*HI/11-1),T/(22*LO),endpoint)<STEP
  assert (1+LO)*STEP>child and (1+2*LO)*STEP>closing
  bias=HI/2-F(2,3);c=3-2*LO
  assert 0<LO/2-F(2,3) and bias*T<STEP
  assert T*(F(1,3)-LO/6)<STEP
  assert T*(bias+c)<endpoint
  assert F(627,250)**2>2*pi_hi and (1-T*T/2)/F(627,250)>PHI
  assert F(33,200)<F(1,6)<F(17,100)
  assert F(54,100)**2<F(11,36)<F(56,100)**2
  assert (1-(F(17,100)*T/F(54,100))**2/2)/(F(56,100)*F(627,250))>F(71,100)
  assert F(73,250)<1-HI/2<1-LO/2<F(147,500)
  assert F(707,1000)<LO/2<HI/2<F(177,250)
  assert (1-(F(147,500)*T)**2/(2*F(707,1000)**2))/(F(177,250)*F(627,250))>F(14,25)
  return pi_hi


def closing_strip(xlow,ylow):
  numerator=282842*ylow+100000*xlow-39000
  assert numerator>0
  return central_strip(173205*numerator//1000000000)


BASELINE_FLOORS = ('241/5000', '419/100000', '89/100000', '89/500000')
BASELINE_FLOORS = tuple(map(F, BASELINE_FLOORS))

def strip_box(left,right,offset,squared=False):
  a,b=left-offset,right-offset
  assert 0<=a<b<=F(7,10)
  shape=((1,F(1)),(3,F(-1,2)),(5,F(9,40)),(7,F(-9,112)))
  terms=[(i+j,c*d) for i,c in shape for j,d in shape] if squared else shape
  moments=[sum(c*(b**(k+j+1)-a**(k+j+1))/(k+j+1)for k,c in terms)for j in range(3)]
  m0,m1,m2=moments;assert m0>0
  return F(13819,10000)**(2 if squared else 1)*m0,m1/m0+offset,m2/m0+2*offset*m1/m0+offset**2

def baseline_controls():
  def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
  pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
  pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
  assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
  assert LO**2<2<HI**2 and F(13819,10000)**2*pi_hi<6
  assert F(627,250)**2>2*pi_hi
  assert (1-F(3,100)**2/2)/F(627,250)>F(797,2000)
  normalizers=F(229,50),F(81027,10000),F(128,10),F(19)
  norm2,norm3,norm4,norm5=normalizers
  assert norm2**2*pi_hi<66
  assert norm3*pi_hi<18*LO
  assert norm4**2*(2*pi_hi)**3<4*10260
  assert norm5**2*(2*pi_hi)**4<4*7*12**4
  return pi_hi,normalizers

def weighted_distances(normalizers, threshold):
  budgets, left = [], F(1,100)
  shifts = (2-14*LO/11, F(7,3)-3*LO/2, 2-24*LO/19, (61-39*LO)/28)
  assert 0 < threshold <= F(3,100) and all(q > 0 for q in shifts)
  assert threshold*max(8*HI/11-1,F(1,11),HI/2-F(2,3),
                       F(1,3)-LO/6) < left
  for distance,normalizer,usable in zip(range(2,6),normalizers[:4],BASELINE_FLOORS,strict=True):
    single=distance==2;n=distance-1
    cells = 96 if single else (30 if distance == 3 else 12)
    step = F(69,100*cells) if distance <= 3 else F(1,25)
    delta = threshold*shifts[distance-2]
    common = (1-2*delta**2)**2
    sides, mids, ratios = [], [], []
    for i in range(cells):
      a, b = left+i*step, left+(i+1)*step
      assert 0 < delta < a < b and b+delta < F(71,100)
      assert F(3,2) < 2*(1-b*b/2)
      sides.append(strip_box(a,b,F(0),squared=single))
      mids.append((step,(a+b)/2,(a*a+a*b+b*b)/3))
      ratios.append(delta/(a-a**3/2))
        # For odd paths the two reflected strip displacements can oppose.
    factors = [[(U*(common-(x*y if distance%2 else 0))).__floor__()
                for y in ratios] for x in ratios]
    assert 0 < min(map(min,factors)) <= max(map(max,factors)) <= U
        # Scale the exact quadratic exponents to one denominator before
        # traversing the boxes. The final directed division is unchanged.
    from math import lcm
    total = 0
    if single:
      for i, (_, _, mean_square) in enumerate(sides):
        exponent = F(33,2)*mean_square
        index = ceil(1000*exponent.numerator, exponent.denominator)
        weight = normalizer*sides[i][0]*F(factors[i][i],U)*exp_index(index)
        total += weight.numerator//weight.denominator
    else:
      axes = (sides,)+(mids,)*(n-2)+(sides,)
      diagonal = [[(18 if h in (0,n-1) else F(39,2))*r[2]
                   for r in row] for h,row in enumerate(axes)]
      near = [[[18*HI*a[1]*b[1] for b in axes[h+1]] for a in axes[h]]
              for h in range(n-1)]
      far = [[[6*a[1]*b[1] for b in axes[h+2]] for a in axes[h]]
             for h in range(n-2)]
      denom = lcm(*(x.denominator for row in diagonal for x in row),
                  *(x.denominator for table in near+far for row in table for x in row))
      def ticks(row):
        return tuple(int(denom*x) for x in row)
      diag = tuple(ticks(row) for row in diagonal)
      nn,ff = (tuple(tuple(ticks(row) for row in table) for table in source)
               for source in (near,far))
      scaled = [[normalizer*a[0]*b[0]*step**(n-2)*F(factors[i][j],U)
                 for j,b in enumerate(sides)] for i,a in enumerate(sides)]
      weights = [[(w.numerator,w.denominator) for w in row] for row in scaled]
      for indices in product(range(cells),repeat=n):
        rev = indices[::-1]
        if indices > rev:
          continue
        exponent = sum(row[i] for row,i in zip(diag,indices))
        exponent += sum(table[i][j] for table,i,j in zip(nn,indices,indices[1:]))
        exponent += sum(table[i][j] for table,i,j in zip(ff,indices,indices[2:]))
        index = ceil(1000*exponent,denom)
        numerator,denominator = weights[indices[0]][indices[-1]]
        total += (numerator*exp_index(index)//denominator)*(1 if indices==rev else 2)
    exact=F(total,U)
    floor = F((exact*10**8).__floor__()-1,10**8)
    assert exact>floor>usable
    usable = floor
    budgets.append({'distance':distance,'exact_probability_lower_sum':str(exact),'usable_probability':str(usable),
                    'cells':cells**n,'cell_width':str(step),'maximum_shift':str(delta),
                    'threshold_abs_max':str(threshold)})
  return budgets


CELLS = 60
SHORT_FLOORS = tuple(F(i, 10**7) for i in (552889, 313397, 114197, 75179, 142689, 32000))
LONG_FLOORS = tuple(F(i, 10**9) for i in (4439266, 874902, 166833))
LONG_PARAMETERS = (
    (F(64, 5), F(95, 144), F(52, 125), F(49, 100)),
    (F(96, 5), F(7, 9), F(66, 125), F(9, 20)),
    (F(28), F(55, 64), F(78, 125), F(429, 1000)),
)


def ceil(numerator, denominator):
  return (numerator+denominator-1)//denominator


def child_weight(previous, current, following):
    # Input indices locate upper cell endpoints, with zero denoting an absent
    # predecessor. These are the doubled mean and twelve times the variance of
    # 50000*y + 70711*x + 25000*p in cell-index coordinates.
  mean2 = 70711*(2*current-1)+50000*(2*following-1)
  variance12 = 70711**2+50000**2
  if previous:
    mean2 += 25000*(2*previous-1)
    variance12 += 25000**2
  index = ceil(3*mean2*mean2+variance12, 50_000_000_000)
  return 27638*exp_index(index)//1_000_000


def root_weight(middle_upper, first):
    # Use the worst middle endpoint, averaging only the arm coordinate.
  mean2 = 100*middle_upper+200*first-100
  return 69*exp_index(ceil(3*mean2*mean2+10000, 200000))//2500


def dot(left, right):
  return sum(a*b for a, b in zip(left, right, strict=True))//U


def moving_arm_inputs():
  endpoints = range(2, CELLS+2)
  closing = [[closing_strip(a-1, b-1) for b in endpoints] for a in endpoints]
  first = [[child_weight(0, a, b) for b in endpoints] for a in endpoints]
  two = [dot(row, strip) for row, strip in zip(first, closing, strict=True)]
  continuation = [[[child_weight(a, b, c) for c in endpoints] for b in endpoints]
                  for a in endpoints]
  tail = closing
  for _ in range(2):
    tail = [[dot(weights, tail[j]) for j, weights in enumerate(matrix)]
            for matrix in continuation]
  four = [dot(row, strip) for row, strip in zip(first, tail, strict=True)]
  for _ in range(4):
    for _ in range(2):
      tail = [[dot(weights, tail[j]) for j, weights in enumerate(matrix)]
              for matrix in continuation]
    for i, (row, strip) in enumerate(zip(first, tail, strict=True)):
      four[i] += dot(row, strip)
  two_sided, baseline, four_sided = [], [], []
  for x in range(2, CELLS+1):
    root = [root_weight(x, a) for a in endpoints]
    two_sided.append(dot(root, two))
    baseline.append(central_strip(173205*(100*(x-1)-49)//1000000))
    four_sided.append(dot(root, four))
  return two_sided, baseline, four_sided


def short_paths(two, baseline, four):
  sums = [0]*6
  earlier = [a+b for a, b in zip(two, baseline, strict=True)]
  for i, x in enumerate(range(2, CELLS+1)):
    density2 = 229*exp_index(ceil(11*(3*x*x-3*x+1), 20))//5000
    sums[0] += density2*2*two[i]*baseline[i]//U**2
    sums[1] += density2*two[i]**2//U**2
    sums[4] += density2*(2*earlier[i]*four[i]+four[i]**2)//U**2
    for j in range(i+1):
      y = j+2
      squares = 3*(x*x+y*y-x-y)+2
      product = (2*x-1)*(2*y-1)
      exponent = ceil(18*(200000*squares+3*70711*product), 6_000_000)
      density3 = 81027*exp_index(exponent)//100_000_000
      mult = 1 if i == j else 2
      sums[2] += mult*(density3*(two[i]*baseline[j]+baseline[i]*two[j])//U**2)
      sums[3] += mult*(density3*two[i]*two[j]//U**2)
      extra = earlier[i]*four[j]+four[i]*earlier[j]+four[i]*four[j]
      sums[5] += mult*(density3*extra//U**2)
  return tuple(F(n, U) for n in sums)


def longer_paths(two, baseline, four, *, last_distance=6):
    # Same central moment factors as the earlier long path contraction.
  start = [a+b for a, b in zip(two, four, strict=True)]
  end = [a+2*b for a, b in zip(start, baseline, strict=True)]
  centers = tuple(2*i+3 for i in range(len(start)))

  def weight(previous, current, following, terminal=False):
    diagonal = 900000 if terminal else 975000
    numerator = (diagonal*(3*following*following+1)
                 +54*70711*current*following+900000*previous*following)
    return exp_index(ceil(numerator, 6_000_000))//100

  first = [a*(exp_index(ceil(9*x*x+3, 20))//100)//U
           for a, x in zip(start, centers, strict=True)]
  initial = [[a*weight(0, x, y)//U for y in centers] for a, x in zip(first, centers, strict=True)]
  terminal = [[dot((weight(x, y, z, True) for z in centers), end) for y in centers]
              for x in centers]
  continuation = [[[weight(x, y, z) for z in centers] for y in centers] for x in centers]
  sums = []
  for distance, (norm, _, _, _) in enumerate(LONG_PARAMETERS[:last_distance-3], 4):
    total = sum(a*b for row, tail in zip(initial, terminal, strict=True)
                for a, b in zip(row, tail, strict=True))//U
    sums.append(norm*F(total, U))
    if distance < last_distance:
      terminal = [[dot(weights, terminal[j]) for j, weights in enumerate(matrix)]
                  for matrix in continuation]
  return tuple(sums)


def certificate(*, last_distance=6):
  assert last_distance in (4, 5, 6)
  pi_hi = moving_controls()
  arms = moving_arm_inputs()
  penalty = 2*(F(3,2)*F(49,10000)**2+11*6*F(1,500)**2)
  paired_factor = F(99939,100000)
  assert 0 < paired_factor < 1-penalty
  short = tuple(paired_factor*value for value in short_paths(*arms))
  long = tuple(paired_factor*value for value in longer_paths(*arms, last_distance=last_distance))
  long_floors = LONG_FLOORS[:last_distance-3]
  long_parameters = LONG_PARAMETERS[:last_distance-3]
  assert all(exact > floor > 0 for exact, floor in zip(short, SHORT_FLOORS, strict=True))
  assert all(exact > floor > 0 for exact, floor in zip(long, long_floors, strict=True))
  short_floors = tuple(F((value*10**9).__floor__()-1,10**9) for value in short)
  long_floors = tuple(F((value*10**9).__floor__()-1,10**9) for value in long)
  assert all(a>b for a,b in zip(short_floors, SHORT_FLOORS, strict=True))
  assert all(a>b for a,b in zip(long_floors, LONG_FLOORS[:last_distance-3], strict=True))
  for distance, ((norm, variance, drift, density), floor) in enumerate(
          zip(long_parameters, long_floors, strict=True), 4):
    assert norm**2*(2*pi_hi)**(distance-1) < 36*12**(distance-1)*variance
  return {
      "threshold_abs_max": str(T), "common_normal_density_floor": str(PHI),
      "short_floors": list(map(str, short_floors)), "long_floors": list(map(str, long_floors)),
  }


associated_arm_certificate = certificate
EXTENT = F(9, 400)
REQUIRED_CREDIT = F('30015793279999/1000000000000000000')

ADJACENT_UNIT = 10**18
ADJACENT_WIDTH = F(1,200)
ADJACENT_EXTENT = EXTENT


def packed_chords(prefix, strips, convolve):
  n = len(strips)-1
  lower = convolve(prefix, strips[:n], n)
  upper = convolve(prefix, strips[2:]+[0], n)
  return [(3*a+b, a+3*b) for a,b in zip(lower, upper, strict=True)]


def chord_cell(i, m0, m1, j0, j1, width_squared):
  dm, dj = m1-m0, j1-j0
  assert dm >= 0 and dj >= 0
  c0 = m0*j0*j0
  c1 = dm*j0*j0+2*m0*j0*dj
  c2 = 2*dm*j0*dj+m0*dj*dj
  c3 = dm*dj*dj
  moment = 60*c0+30*c1+20*c2+15*c3
  tilt = (30*i+15)*c1+(30*i+16)*c2+(27*i+15)*c3
  value = 6*width_squared*moment-33*tilt
  assert value >= 0
  return value//(360*width_squared*U**2)


def chord_joint_rows(densities, middle, old, high, width, gain):
  divisor = 4*U*width
  old_row = joint_row = 0
  for i, (density, (a,b), (c,d)) in enumerate(zip(densities, old, high, strict=True)):
    old0, old1 = a//divisor, b//divisor
    j0 = min(a*gain//(divisor*U), c//divisor)
    j1 = min(b*gain//(divisor*U), d//divisor)
    assert j0 >= old0 and j1 >= old1
    args = (i, middle[i], middle[i+1])
    old_row += density*chord_cell(*args, old0, old1, width**2)//(U*width)
    joint_row += density*chord_cell(*args, j0, j1, width**2)//(U*width)
  return old_row, joint_row


def chord_clean_row(densities, middle, chords, width):
  divisor = 4*U*width
  return sum(
      density*chord_cell(i, middle[i], middle[i+1], a//divisor, b//divisor, width**2)//(U*width)
      for i, (density, (a,b)) in enumerate(zip(densities, chords, strict=True))
  )


HALF_CDF_TERMS = tuple(F(-3, 2)**i/(factorial(i)*(2*i+1)) for i in range(24))
HALF_CDF_DENOMINATOR = lcm(*(term.denominator for term in HALF_CDF_TERMS))
HALF_CDF_COEFFICIENTS = tuple(int(term*HALF_CDF_DENOMINATOR) for term in HALF_CDF_TERMS)


@lru_cache(maxsize=None)
def half_cdf(value):
  assert 0 < value <= 2
  n, d = value.numerator, value.denominator
  n2, d2 = n*n, d*d
  polynomial = 0
  denominator_power = 1
  for coefficient in reversed(HALF_CDF_COEFFICIENTS):
    polynomial = polynomial*n2+coefficient*denominator_power
    denominator_power *= d2
  assert polynomial > 0
  return ADJACENT_UNIT*13819*n*polynomial//(10000*HALF_CDF_DENOMINATOR*d**47)


def adjacent_margins(extent):
  def wave(distance, root):
    value = (1+F(distance,3))/2**((distance+1)//2)
    return value*root if distance%2 else value

  margins = []
  for depth in (1,2,3):
    lower = (wave(depth,LO)+wave(depth+1,LO))/(1+wave(1,HI))
    upper = (wave(depth,HI)+wave(depth+1,HI))/(1+wave(1,LO))
    required = extent*max(abs(1-lower),abs(1-upper))
    margins.append(F((required*10**6).__floor__()+1,10**6))
  return margins


def adjacent_certificate(extent=ADJACENT_EXTENT):
  assert 0 < extent <= F(3,100)
  margins = adjacent_margins(extent)
  first_margin, second_margin, closing_margin = margins
  prefactor = F(13819,10000)
  pi_upper = F(3927,1250)
  assert pi_upper*prefactor**2 < 6
  space = 10**6
  first, second, closing = [int(v*space) for v in margins]
  width = int(ADJACENT_WIDTH*space)
  ar = [first]+list(range((first//width+1)*width,space+1,width))
  kr = [second]+list(range((second//width+1)*width,6*space//5+1,width))
  rn,rd = (2*HI).as_integer_ratio()
  ln,ld = (2*LO).as_integer_ratio()
  den = 8*rd**2*space**2
  strips = [0]+[half_cdf(F(i,1000)) for i in range(1,2001)]
  one_side = cells = 0
  for a,b in zip(ar,ar[1:]):
    wa,ma = b-a,a+b
    amplitude = 27638*exp_index(ceil(1000*(3*ma*ma+wa*wa),2*space**2))//10000
    neighbor = 0
    for k,l in zip(kr,kr[1:]):
      wk = l-k
      radius = (ln*k+ld*(a-2*closing))*1000//(ld*space)
      if radius <= 0:continue
      q = 3*(2*rd*(k+l)+rn*ma)**2+4*rd**2*wk**2+rn**2*wa**2
      density = 13819*exp_index(ceil(1000*q,den))//10000
      neighbor += 2*wk*density*strips[min(radius,2000)]//(space*U)
      cells += 1
    one_side += wa*amplitude*neighbor//(space*U)

  probability = F(one_side, ADJACENT_UNIT)**2
  usable = F((probability*10**8).__floor__()-1,10**8)
  assert probability > usable

  drift_lo, drift_hi = F(56, 1000), F(58, 1000)
  sd_lo, sd_hi = F(32, 100), F(35, 100)
  density_lower = F(11, 10)
  assert (1-drift_hi)**2 < F(8, 9) < (1-drift_lo)**2
  assert sd_lo**2 < F(1, 9) < sd_hi**2
  exponent = (drift_hi*extent)**2/(2*sd_lo**2)
  assert exponent < 1
  assert 2*pi_upper*(density_lower*sd_hi)**2 < (1-exponent)**2
  return {"scope": "exact adjacent outer probability and scalar rate",
          "pivotal_cells": cells, "exact_single_orientation": str(probability),
          "usable_single_orientation": str(usable),
          "usable_probability": str(2*usable),
          "threshold_abs_max": str(extent),
          "rate_without_normal_density": str(3*drift_lo*density_lower*2*usable)}


def long_baseline(distance, left, offset, normalizer, floor, *, reflected=None):
  step = F(1, 40)
  cells = range(12)
  sides, middle = [], []
  for i in cells:
    a, b = left+i*step, left+(i+1)*step
    sides.append(strip_box(a, b, offset))
    middle.append((step, (a+b)/2, (a*a+a*b+b*b)/3))

  def exponential(q):
    return exp_index(ceil(1000*q.numerator, q.denominator))

  adjacent = 18*HI
  states = [[(a[0]*b[0]*exponential(18*a[2]+F(39, 2)*b[2]+adjacent*a[1]*b[1])).__floor__()
             for b in middle] for a in sides]
  coordinates = (sides,)+(middle,)*(distance-3)+(sides,)
  for index in range(2, distance-1):
    A, B, C = coordinates[index-2:index+1]
    diagonal = 18 if index == distance-2 else F(39, 2)
    kernel = [[[(c[0]*exponential(diagonal*c[2]+adjacent*b[1]*c[1]+6*a[1]*c[1])).__floor__()
                for c in C] for b in B] for a in A]
    states = [[sum(states[i][j]*kernel[i][j][k] for i in cells)//U
               for k in cells] for j in cells]
  exact = normalizer*F(sum(map(sum, states)), U)
  if reflected is not None:
    assert distance%2 == 0 and offset == 0 and 0 < reflected < left
    right = left+12*step
    assert right+reflected < F(61,100)
    assert F(3,2) < 2*(1-right*right/2)
    exact *= (1-2*reflected**2)**2
    usable = F((exact*10**9).__floor__()-1,10**9)
    assert usable > floor
    floor = usable
  assert exact > floor
  return {"distance": distance, "exact_probability_lower_sum": str(exact),
          "usable_probability": str(floor), "cells": 12**(distance-1),
          "cell_width": str(step), "offset": str(offset)}


def additional_baselines(pi_hi, join):
  assert 0 < join <= F(3,100)
  T = join
  assert T*max(8*HI/11-1,F(1, 33),F(4, 33),1-20*LO/33,F(5, 33)) < F(1, 100)
  delta = T*(2-14*LO/11)
  assert 0 < 2-14*HI/11 and delta < F(1,100)
  T = EXTENT
  for a, b, den in ((-470,339,526),(-63,46,263),(-155,109,263),
                    (-48,10,263),(-430,243,526),(-135,61,263)):
    assert max(abs(a+b*LO), abs(a+b*HI))*T/den < F(3, 500)
  assert 0 < (1108-713*HI)/526 and T*(1108-713*LO)/526 < F(43, 10000)
  for w_lo, w_hi, variance, drift_lo, drift_hi, sd_lo, sd_hi, density in (
      (F(3, 8), F(3, 8), F(55, 64), F(78, 125), F(627, 1000),
       F(927, 1000), F(929, 1000), F(429, 1000)),
      (5*LO/24, 5*HI/24, F(263, 288), F(141, 200), F(353, 500),
       F(191, 200), F(239, 250), F(417, 1000)),
  ):
    assert drift_lo < 1-w_hi <= 1-w_lo < drift_hi
    assert sd_lo**2 < variance < sd_hi**2
    assert 2*pi_hi < F(627, 250)**2
    assert (1-(drift_hi*T)**2/(2*sd_lo**2))/(sd_hi*F(627, 250)) > density
  norm6, norm7 = F(28), F(399, 10)
  assert norm6**2*(2*pi_hi)**5 < 4*1924560
  assert norm7**2*(2*pi_hi)**6 < 4*24541056
  return [long_baseline(6, F(1,100), F(0), norm6, F(272,10**7), reflected=delta),
          long_baseline(7, F(3, 500), F(43, 10000), norm7, F(59, 10**7))]


def extended_even_baselines(pi_hi, join):
  def wave(distance, root):
    value = (1+F(distance,3))/2**((distance+1)//2)
    return value*root if distance%2 else value

  def magnitude(lower, upper):
    assert lower <= upper
    return max(abs(lower),abs(upper))

  rows = []
  left = F(1,100)
  for distance in (8,10):
    w = wave(distance,LO)
    bias = ([(F(0),F(0))]
            + [((wave(i,LO)+wave(distance-i,LO))/(1+w)-1,
                (wave(i,HI)+wave(distance-i,HI))/(1+w)-1)
               for i in range(1,distance)]
            + [(F(0),F(0))])
    for i in range(1,distance):
      lower,upper = bias[i]
      assert join*magnitude(lower,upper) < left
      products = (2*LO*lower,2*LO*upper,2*HI*lower,2*HI*upper)
      third_lo = min(products)-bias[i-1][1]-bias[i+1][1]-(3-2*LO)
      third_hi = max(products)-bias[i-1][0]-bias[i+1][0]-(3-2*HI)
      assert join*magnitude(third_lo,third_hi) < 2*LO*left
    delta = join*magnitude(bias[1][0]+3-2*HI,bias[1][1]+3-2*LO)
    normalizer_square = 36*12**(distance-1)*(1-w*w)/(2*pi_hi)**(distance-1)
    norm = F(isqrt((normalizer_square*10**6).__floor__()),1000)
    assert 0 < norm**2 < normalizer_square
    rows.append(long_baseline(distance,left,F(0),norm,F(0),reflected=delta))
  return rows


def outer_time_moment(threshold, used, ratio):
  assert 0 < used <= threshold < 1 and 0 < ratio < 1
  moment = used**2/2-(F(1, 8)+ratio/24)*used**4+ratio*used**6/72
  moment += (threshold-used-(threshold**3-used**3)/6)*(used-ratio*used**3/6)
  assert moment > 0
  return moment


def short_adjacent_credit(full, short, adjacent_pair, pi_hi, threshold):
  assert adjacent_pair["distance"] == 1
  w = F(adjacent_pair["cover_correlation_upper"])
  ratio = (1-w)/(1+w)
  root = F(isqrt((ratio*10**24).__floor__()),10**12)
  assert 0 < root and root**2 <= ratio
  used = F(short["threshold_abs_max"])
  assert 0 < used <= F(full["threshold_abs_max"]) and used <= threshold
  gain = F(short["usable_probability"])-F(full["usable_probability"])
  assert gain > 0
  value = 3*gain*root*outer_time_moment(threshold, used, ratio)/(2*pi_hi)
  return F((value*U).__floor__(),U)


def integrated_pair_credits(baselines, arm, adjacent, pi_hi, threshold):
  assert 0 < threshold < 1
  used = min(EXTENT, threshold)
  baseline_floors = {1:F(adjacent["usable_probability"])}
  baseline_floors.update((row["distance"],F(row["usable_probability"])) for row in baselines)
  probabilities = dict(baseline_floors)
  b2, bb2, b3, bb3, c2, c3 = map(F, arm["short_floors"])
  probabilities[2] += b2+bb2+c2
  probabilities[3] += b3+bb3+c3
  for distance, floor in enumerate(arm["long_floors"], 4):
    probabilities[distance] += F(floor)

  credits = []
  for distance, probability in probabilities.items():
    critical_upper = (1+F(distance, 3))/2**((distance+1)//2)
    if distance % 2:
      critical_upper *= HI
    w = critical_upper+F(1, 10**9)
    assert 0 < w < 1 and 0 < EXTENT < 1
    ratio = (1-w)/(1+w)
    root = F(isqrt((ratio*10**24).__floor__()), 10**12)
    assert 0 < root <= 1 and root**2 <= ratio
        # Only the reflected baseline extends to the join.
    moment = outer_time_moment(threshold, used, ratio)
    extent = threshold if distance != 7 else used
    extra = outer_time_moment(threshold, extent, ratio)-moment
    assert extra >= 0
    mass = probability*moment+baseline_floors[distance]*extra
    value = 3*2**(distance-1)*root*mass/(2*pi_hi)
    credit = F((value*U).__floor__(), U)
    assert credit > 0
    credits.append({"distance": distance, "negative_mixed_response_floor": str(probability),
                    "baseline_extent": str(extent),
                    "cover_correlation_upper": str(w), "credit": str(credit)})
  return credits


def certificate():
  pi_hi, norms = baseline_controls()
  s_max = F(3, 100)
  variance = F(95, 144)
  assert 1 - F(7, 12) > F(52, 125)
  assert 2*pi_hi*variance*F(49, 100)**2 < (1-s_max**2/(2*variance))**2
  assert 0 < 14*LO/19-1 < 14*HI/19-1 < F(1, 15)
  assert s_max*max(F(2, 19), 1-12*LO/19) < F(1, 100)
  assert s_max*(2-24*LO/19) < F(13, 2000)
  variance = F(7, 9)
  assert 1-HI/3 > F(66, 125)
  assert 2*pi_hi*variance*F(9, 20)**2 < (1-s_max**2/(2*variance))**2
  assert s_max*max((17*HI-23)/28,(4*HI-5)/14,(3-LO)/14,(19-11*LO)/28) < F(1, 100)
  assert s_max*(61-39*LO)/28 < F(13,2000)
  join = F('277/10000')
  baselines = (weighted_distances(norms,join)+additional_baselines(pi_hi,join)
               +extended_even_baselines(pi_hi,join))
  adjacent = adjacent_certificate(join)
  assert EXTENT <= F(adjacent["threshold_abs_max"])
  arm = associated_arm_certificate(last_distance=6)
  assert [row["distance"] for row in baselines] == [2,3,4,5,6,7,8,10]
  assert EXTENT <= F(arm["threshold_abs_max"])
  assert F(arm["common_normal_density_floor"]) == PHI
  credits = integrated_pair_credits(baselines, arm, adjacent, pi_hi, join)
  levels = [adjacent] + [adjacent_certificate(F(j,1000)) for j in (25,20,15,10,5)]
  extra = sum(short_adjacent_credit(a,b,credits[0],pi_hi,join)
              for a,b in zip(levels,levels[1:]))
  saving = sum(F(row["credit"]) for row in credits)+extra
  assert saving > REQUIRED_CREDIT
  return {"scope": "outer with adjacent and extended even baselines through the join",
          "threshold_abs_max_of_central_band": str(join),
          "baselines": baselines,
          "pair_credits": credits, "adjacent": adjacent,
          "adjacent_levels": levels, "additional_adjacent_credit": str(extra),
          "arms": arm,
          "outer_graph_credit": str(saving),
          "slack_above_uniform_saving": str(saving-REQUIRED_CREDIT)}


outer = SimpleNamespace(certificate=certificate)


SCORE_UNIT = 1000000000
PROBABILITY_UNIT = 1000000000000
MARGIN_STEP_END = 500000

def midpoint_scales(eta, star_factor):
  root_lo, root_hi = F(141421356,10**8), F(141421357,10**8)
  assert root_lo**2 < 2 < root_hi**2
  assert 0 < eta <= F('57/50000')
  variances = [factor*(3-2*root_lo)*eta for factor in star_factor]
  sds = [F(isqrt((value*10**12).__floor__())+1,10**6) for value in variances]
  assert all(sd**2 > value for sd,value in zip(sds,variances,strict=True))
  half = F(289,1000)
  assert half**2 > F(1,12)
  midpoint_sd = (sds[0]+root_hi*sds[1]+sds[2])/2
  correlation = (root_hi/4)*eta*(star_factor[0]+star_factor[1]+star_factor[2]/2)/(midpoint_sd*half)
  assert 0 < correlation < 1
  return midpoint_sd,half,1-correlation**2


def positive_defect_bins(sigma, step):
  normal = F(79788,10**5)
  pi_upper = F(355,113)
  atan5 = sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(9)),F(0))
  assert 16*atan5-4*(F(1,239)-F(1,3*239**3)) < pi_upper
  assert normal**2*pi_upper < 2
  scale = PROBABILITY_UNIT
  bins = []
  for lo in range(0,120000000,step):
    hi = lo+step
    exponent = F(lo*lo+lo*hi+hi*hi,6*SCORE_UNIT**2)/sigma**2
    assert exponent < 4096
    value = (scale*(1-exponent/4096)).__floor__()
    for _ in range(12):
      value = value*value//scale
    mass = (F(step,SCORE_UNIT)*normal*value/sigma).__floor__()
    bins.append((hi,mass))
  return tuple(bins)


def midpoint_strip_table(eta, strip, star_factor):
  sigma,half,factor = midpoint_scales(eta, star_factor)
  bins = positive_defect_bins(sigma,1000000)
  counts = []
  for margin in range(0,300000001,MARGIN_STEP_END):
    total = 0
    numerator,denominator = 10000*half.denominator,SCORE_UNIT*half.numerator
    previous = strip[min(30000,margin*numerator//denominator)]
    for hi,mass in bins:
      available = margin-hi
      if available < 0:
        break
      endpoint = strip[min(30000,available*numerator//denominator)]
      total += mass*(previous+endpoint)
      previous = endpoint
    counts.append(total//(2*10**6))

    # Panel endpoints for sqrt(1-rho**2), rounded in opposite
    # directions for the mass multiplier and the radius increase.
  unit = 10**9
  bounds = [isqrt(((1-(1-factor)*F(i,8))*unit**2).__floor__()) for i in range(9)]
  counts = [min(bounds[i+1]*counts[min(len(counts)-1,j*unit//(bounds[i]+(i>0)))]//unit
                for i in range(8)) for j in range(len(counts))]
  assert counts == sorted(counts) and 0 < counts[-1] <= 10**15
  return tuple(counts)


ROOT_TICKS = tuple(range(6000000,60000001,4000000))
MID_STEP = 1000000
MID_MAX_INDEX = 90
CORE_STEP = 1000000
CORE_SCALE = 1000000000000
SCORE_TICKS = 1000000000
CORR_PANELS = 48

def adaptive_core_table(eta, strip, star_factor):
  lo2,hi2 = F(141421356,10**8), F(141421357,10**8)
  assert lo2**2 < 2 < hi2**2
  x,y = [factor*eta for factor in star_factor]
  c = 3-2*lo2
  assert 0 < 2*min(x,y) <= 2*max(x,y) < c
  vx,vy = c*x-x*x,c*y-y*y/2
  sx,sy = (F(isqrt((v*10**12).__floor__())+1,10**6) for v in (vx,vy))
  assert sx**2 > vx and sy**2 > vy
  bins = positive_defect_bins(sx,CORE_STEP)[:ROOT_TICKS[-1]//CORE_STEP]
  assert len(bins)*CORE_STEP == ROOT_TICKS[-1]
  table = [[CORE_SCALE]*(MID_MAX_INDEX+1) for _ in ROOT_TICKS]
  edges = [radius//CORE_STEP for radius in ROOT_TICKS]
  assert all(i*CORE_STEP == radius for i,radius in zip(edges,ROOT_TICKS,strict=True))

  for panel in range(CORR_PANELS//2,CORR_PANELS):
    left = F(2*panel-CORR_PANELS,CORR_PANELS)
    right = F(2*panel+2-CORR_PANELS,CORR_PANELS)
    assert left*right >= 0
    nearest = min(abs(left),abs(right))
    sd_square = sy**2*(1-nearest**2)
    sd_ticks = isqrt((sd_square*SCORE_TICKS**2).__floor__())+1
    multiplier = max(abs(alpha+rho*sy/sx) for alpha in (2*lo2/3,2*hi2/3)
                     for rho in (left,right))
    mean_caps = [((hi*multiplier).__ceil__(),weight) for hi,weight in bins]
    for j in range(MID_MAX_INDEX+1):
      width = 2*j*MID_STEP
      partials = [0]
      total = 0
      previous = 2*strip[min(30000,10000*width//sd_ticks)]
      for mean,weight in mean_caps:
        if width >= mean:
          inner = min(30000,10000*(width-mean)//sd_ticks)
          outer = min(30000,10000*(width+mean)//sd_ticks)
          endpoint = strip[inner]+strip[outer]
          total += weight*(previous+endpoint)
          previous = endpoint
        partials.append(total)
      for row,i in zip(table,edges,strict=True):
        row[j] = min(row[j],partials[i]//(4*10**9))

  assert all(row == sorted(row) for row in table)
  assert all(all(a <= b for a,b in zip(lo,hi,strict=True)) for lo,hi in zip(table,table[1:]))
  return table


ONE_LEAF_BONUS = F('7/200')


def sharp_scalar_dyad():
    # All endpoint arguments are in units 1/2^20; sqrt enclosures in 1/2^40.
    N=1<<20; S=N*N
    h2=isqrt(2*S*S)
    def ratios(a,b):
        # Bounds for D_phi and D_A, with lengths scaled by S.
        lp=isqrt(a*S*S//N)+h2*(N+a)//(3*N)
        up=isqrt(b*S*S//N)+1+(h2+1)*(N+b)//(3*N)+1
        ua=isqrt(b*(5*b-2*N)*S*S//(2*N*N))+1+(h2+1)*(7*b-2*N)//(6*N)+1
        return lp,up,ua
    def lohi(a,b):
        g=lambda x:(x-N//2)*(2*N-x)
        k=max(a,min(5*N//4,b))
        return min(g(a),g(b)),g(k)
    # Far corner p >= 61/32, hence q,r <= 19/32. The ratio is D(p)/D(5/2-p).
    assert 5*ratios(61*N//32,61*N//32)[0]>=9*ratios(19*N//32,19*N//32)[1]
    qmin=N//2;maxp=61*N//32
    stack=[(N, maxp,qmin,N)]; pass_count=out_count=0
    def slot(a,b):
        nlo,nhi=lohi(a,b)
        lp,_,ua=ratios(a,b)
        return (50*nlo*S//ua)-(34*nhi*S+lp-1)//lp
    while stack:
        a,b,c,d=stack.pop()
        l=max(qmin,3*N-b-d);h=min(2*N,3*N-a-c)
        if l>h or c>b:out_count+=1; continue
        extra=5*lohi(l,h)[0]*S//ratios(l,h)[1]
        if slot(a,b)+slot(c,d)+extra>=0:pass_count+=1;continue
        if b-a>=d-c:
            m=(a+b)//2;assert a<m<b;stack.extend(((a,m,c,d),(m,b,c,d)))
        else:
            m=(c+d)//2;assert c<m<d;stack.extend(((a,b,c,m),(a,b,m,d)))
    return pass_count,out_count


def medium_payment_slacks(rho, shared, payment, anchor, positive10, epsilon, t, z02, k):
  """Check the medium and run budgets in units of the Price allowance."""
  sigma = F(2021,1250)
  a0,kappa = sigma-F(3,2),F(59599,10**6)
  c,d = a0/4,(sigma-1)/8
  g = F(2,3)*(shared-c-rho/4)
  B = c+g/2
  nu,nu2 = (F(1,4)+rho)*a0,shared-F(1,8)
  assert k == F(33,100) and 0 <= t and 0 <= z02
  assert F(2,5) <= rho <= F(1,2)
  assert F(1,8)+rho/4 <= shared <= F(1,4)+rho/4
  far_reserve = k*(sigma-F(1,2))/(4*(1+k))
  anchored_four_reserve = k*(sigma-F(1,2))/2
  mixed_basic = (2*kappa+F(83,2500))/3+F(9,20)*a0
  short_long_reserve = mixed_basic+far_reserve
  one_long_reserve = kappa/2+(sigma-1)/4-a0/2
  two_long_reserve = one_long_reserve+far_reserve
  shallow_reserve = kappa+F(9,80)
  one_leaf_bonus = ONE_LEAF_BONUS
  assert one_leaf_bonus < (sigma-1)/4
  assert kappa+far_reserve > one_long_reserve
  assert kappa/2+(sigma-1)/2-a0/2 > two_long_reserve
  assert one_long_reserve+F(9,20)*a0 > short_long_reserve
  assert F(187,5000)+one_leaf_bonus < mixed_basic
  delta = max(F(0), payment-B)
  connector_rate = k*(sigma-1)/8
  assert delta <= connector_rate
  assert c > connector_rate and F(9,20)*a0 > 2*connector_rate
  assert 2*d-2*c > 2*connector_rate
  assert anchored_four_reserve-far_reserve > 4*connector_rate
  closed_four = F(2,5)*a0+F(7,10)*k
  closed_six = a0+k/2
  end_group = F(2,5)*a0+k/2
  assert F(9,40)+a0/2 > closed_four
  assert F(9,40)+end_group > max(closed_six,2*end_group)
  assert F(9,20) > max(closed_four,closed_six,2*end_group)
  share = F(47,400)
  assert share+k/4 <= F(1,5) and share <= F(1,8)
  retained = 1-share/F(9,20)
  assert share <= (1-retained)/2
  assert F(9,40)-share/2 > c and F(9,10)-share >= F(1,2)
  assert retained*end_group-c > 5*connector_rate
  assert end_group-2*c > 6*connector_rate
  assert 2*end_group-closed_six > 4*connector_rate
  no_far_base = (1-2*share)*kappa
  no_far_increment = share*k/2
  tau = share/(1+k)
  mixed_total = 1-2*tau*(1-k)
  mixed_weight = (F(9,20)-(F(4,5)-F(9,5)*k)*tau)/F(27,20)
  assert 0 < mixed_weight < mixed_total and tau*k < F(9,20)
  mixed_long = ((mixed_total-mixed_weight)*kappa+mixed_weight*F(187,3750)
                +(F(9,20)+tau*k)*a0+tau*k/2)
  assert mixed_long == F(13604918317,113400000000)
  assert F(7,20)-a0 > max(mixed_long,no_far_base+2*no_far_increment)
  assert shared+F(83,2500) > max(mixed_long,no_far_base+2*no_far_increment)
  full_far = end_group-4*delta
  short_high10_pool = F(9,20)*(k+a0)
  assert k > 2*a0 and k < 1
  assert F(9,40)+F(9,20)*a0 > short_high10_pool
  one_three = kappa/10+F(9,20)*(sigma-1)-F(9,10)*a0
  two_three = one_three+F(9,20)*(sigma-1)
  full_mixed_far = min(one_three+F(9,20)*a0,mixed_basic+full_far)
  return {
      "type01_pair": g/4-t,
      "type02_one_leaf": g/2-t-z02,
      "type02_two_leaves": kappa/2-2*t-z02,
      "type11_high_anchor": anchor-2*delta,
      "type11_one_H1_anchor": payment+t-2*delta,
      "type10_no_far": min(nu,nu2)-positive10-delta,
      "short_high10": shared-payment-positive10-t,
      "long_high10": 2*B+short_high10_pool-positive10-3*payment-t,
      "long_high10_tail": 2*B+nu+anchored_four_reserve-2*delta-positive10-3*payment-t,
      "medium_anchor_N2": 2*d+B-c-2*payment-t,
      "medium_anchor_N4": 3*B+closed_four-c-4*payment-t,
      "medium_anchor_tail": 3*B+closed_six-c-4*payment-t-2*delta,
      "closed_N2": F(9,40)*(sigma-1)-payment-t,
      "closed_N4": B+closed_four/2-2*payment-t,
      "closed_N6": 2*B+closed_six/2-3*payment-t,
      "closed_N8": 3*B+end_group-4*payment-t-delta,
      "zero_H1": F(479,2500)-2*one_leaf_bonus+2*anchor-epsilon,
      "one_high10": F(187,5000)+one_leaf_bonus+2*anchor-epsilon-positive10,
      "one_short": F(187,5000)+one_leaf_bonus+2*anchor-epsilon-payment-t,
      "two_high10": shallow_reserve+2*anchor-epsilon-2*positive10,
      "short_high10_center": shallow_reserve+2*anchor-epsilon-positive10-payment-t,
      "short_long_no_far": (min(mixed_basic-delta,mixed_long-2*delta)
                           +2*anchor-epsilon-payment-t),
      "short_long": full_mixed_far+2*anchor+2*B-epsilon-4*payment-2*t,
      "long_high10_center": full_mixed_far+2*anchor+2*B-epsilon-positive10-3*payment-t,
      "two_long": (min(two_three,one_three+full_far,kappa+2*full_far)+2*anchor+4*B-epsilon-6*payment-2*t),
      "two_long_one_far": (min(one_three,kappa+full_far)+2*anchor+2*B-epsilon-3*payment-t-2*delta),
      "two_long_no_far": (min(kappa-2*delta,no_far_base+no_far_increment-3*delta,
                              no_far_base+2*no_far_increment-4*delta)+2*anchor-epsilon),
      "long_high10_no_far": (min(mixed_basic-delta,mixed_long-2*delta)
                             +2*anchor-epsilon-positive10),
      "double_short": shallow_reserve+2*anchor-epsilon-2*payment-2*t,
  }


@lru_cache(maxsize=1)
def root_policy(private_choice,core_table,step):
  root_debits = tuple(((r+2)//3, (1415*((r+2)//3)+999)//1000) for r in ROOT_TICKS)

  first_debits = tuple(first for first,_ in root_debits)

  @lru_cache(maxsize=None)
  def reflected_probabilities(narrow, lower, upper):
    return tuple((row[lower]+(core_table[min(i,narrow)][upper] if narrow >= 0 else 0))//2
                 for i,row in enumerate(core_table))

    # The cached calculation uses these margins only through integer floors
    # and the strict root feasibility gates.
  first_cuts = tuple(sorted({0, *(df%step for df,_ in root_debits),
                              *((df+1)%step for df,_ in root_debits)}))
  second_cuts = tuple(sorted({0, *(dg%step for _,dg in root_debits)}))

  from array import array

  def shifts(cuts):
    values = array("i")
    for a,b in zip(cuts,(*cuts[1:],step),strict=True):
      values.extend(range(0,a-b,-1))
    assert len(values) == step
    return values

  @lru_cache(maxsize=None)
  def first_bins(first):
    return first//step,tuple((first-df)//step for df,_ in root_debits if df < first)

  @lru_cache(maxsize=None)
  def second_bins(second):
    return second//step,tuple((second-dg)//step for _,dg in root_debits)

  @lru_cache(maxsize=16384)
  def root_payments(first, second, lower, upper, narrow):
    first_bin,first_with_root = first_bins(first)
    second_bin,second_with_root = second_bins(second)
    if narrow < 0:
      best = 0
      for i,(a,b) in enumerate(zip(first_with_root,second_with_root)):
        if a <= 0 or b <= 0:
          break
        positive = core_table[i][lower]
        if positive:
          best = max(best,private_choice(a,b)*positive)
      return best//2,0,0
        # If every radius is eligible and the endpoints agree, each
        # reflected coefficient equals its original positive coefficient.
    if narrow == len(core_table)-1 and lower == upper:
      best_positive = 0
      for i,(a,b) in enumerate(zip(first_with_root,second_with_root)):
        if a <= 0 or b <= 0:
          break
        positive = core_table[i][lower]
        if positive:
          payment = private_choice(a,b)*positive
          if payment > best_positive:
            best_positive = payment
      pay = core_table[narrow][lower]
      negative = pay*private_choice(first_bin,second_bin) if pay else 0
      return max(best_positive,(best_positive+negative)//2),best_positive,negative
    reflected = reflected_probabilities(narrow,lower,upper)
    best_reflected = best_positive = best_both = 0
    for i,(a,b) in enumerate(zip(first_with_root,second_with_root)):
      if a <= 0 or b <= 0:
        break
      positive = core_table[i][lower]
      if not (positive or reflected[i]):
        if i <= narrow:
          best_both = best_positive
        continue
      private_with_root = private_choice(a,b)
      reflected_payment = private_with_root*reflected[i]
      if reflected_payment > best_reflected:
        best_reflected = reflected_payment
      positive_payment = private_with_root*positive
      if positive_payment > best_positive:
        best_positive = positive_payment
      if i <= narrow:
        best_both = best_positive
    pay = core_table[narrow][lower] if narrow >= 0 else 0
    negative = pay*private_choice(first_bin,second_bin) if pay else 0
    best_root = max(best_reflected,(best_positive+negative)//2)
    return best_root, best_both, negative

  return first_debits,shifts(first_cuts),shifts(second_cuts),root_payments


def folded_private_bins(eta, cap, sd, *, shifted):
  lo2,hi2 = F(1414213562373,10**12),F(1414213562374,10**12)
  assert lo2**2 < 2 < hi2**2
  assert 2*F(707107,10**6)**2 > 1 and F(1415,1000)**2 > 2
  defect = 3-2*lo2
  a,b = F(sd,10**6),F(289,1000)
  p = F(290,1000) if shifted else b
  incident = cap*eta
  assert a*a > defect*incident and b*b > F(1,12)
  if shifted:
    assert p*p > F(1,12)+(hi2+defect)*incident/8
    assert F(1,4)+hi2*defect/4 < F(1,2)
  cross = (1+F(561,50))*eta/4
  if shifted:
    cross += incident/8
  x,y,z = incident/(2*a*p),incident/(2*lo2*a*b),cross/(p*b)
  assert max(x+y,x+z,y+z) < 1
  assert p*p*(1-(x*x+z*z)/(1-y)) > F(1543,20000)
  assert b*b*(1-(y*y+z*z)/(1-x)) > F(1543,20000)
  variance = a*a*(1-(x*x+y*y)/(1-z))/4
  assert variance > 0
  pi = F(355,113)
  norm = F(379,50) if shifted else F(38,5)
  normal = F(79788,10**5)
  assert norm*pi*p*b < 2 and normal**2*pi < 2

  edges = (*range(0,32000000,1000000), *range(32000000,52000001,2000000))
  rows = []
  for left,right in zip(edges,edges[1:]):
    q = F(left**2+left*right+right**2,6*SCORE_UNIT**2)/variance
    assert 0 <= q < 4096
    value = (10**18*(1-q/4096)).__floor__()
    for _ in range(12):
      value = value*value//10**18
    mass = (normal*F(right-left,SCORE_UNIT)*2/a*F(value,10**12)).__floor__()
    assert mass >= 0
    rows.append((right,mass))
  return tuple(rows)


def _fiber():
  from functools import cache
  from fractions import Fraction as F
  from itertools import combinations
  from math import isqrt

  SCALE = 1000000000000000000000000

  def sqrt_bounds(x):
    assert x >= 0
    n = isqrt((x*SCALE*SCALE).__floor__())
    lo = F(n, SCALE)
    return lo, lo if lo*lo == x else F(n+1, SCALE)


  def distances(order, edges):
    adjacent = [[] for _ in range(order)]
    for left, right in edges:
      adjacent[left].append(right)
      adjacent[right].append(left)
    rows = []
    for root in range(order):
      row = [-1]*order
      row[root], queue = 0, [root]
      for vertex in queue:
        for neighbor in adjacent[vertex]:
          if row[neighbor] < 0:
            row[neighbor] = row[vertex]+1
            queue.append(neighbor)
      assert min(row) >= 0
      rows.append(row)
    return rows


  CORE = (("o", "u"), ("u", "w"), ("u", "w_prime"))
    # Row, score vertex, sign, attachment-tree edges, weights in hundredths.
  ROWS = (
      (0, "a", 1, (*CORE, ("o", "a")), (215, 215, 110, 160)),
      (1, "a1", 1, (*CORE, ("o", "a"), ("a", "a1")), (343, 148, 263, 364, 218)),
      (2, "w", 1, CORE, (152, 488, 50)),
      (3, "w1", 1, (*CORE, ("w", "w1")), (120, 570, 239, 188)),
      (4, "u", -1, CORE, (214, 36, 114)),
  )
  PADDING = F(1, 50)


  def critical_covariance(depth):
    lo, hi = sqrt_bounds(F(2))
    factor = (1+F(depth, 3))*F(2)**(-((depth+1)//2))
    return (factor*lo, factor*hi) if depth % 2 else (factor, factor)


  def interval_bound(distances, origin, left, right, terms):
    a = [distances[origin][left]-row[left] for row in distances]
    b = [distances[origin][right]-row[right] for row in distances]
    root_lo, root_hi = sqrt_bounds(F(2))

    def inner_coefficient(exponent):
      factor = F(2)**(exponent//2)
      return (factor*root_lo, factor*root_hi) if exponent % 2 else (factor, factor)

    data = []
    for (i, j), low, high in terms:
      if i == j:
        continue
      first, second = a[i]+b[j], a[j]+b[i]
      assert (first+second) % 2 == 0
      numerator = F(2)**first+F(2)**second-2*F(2)**((first+second)//2)
      assert numerator >= 0
      if numerator:
        data.append((i, j, low*numerator, high*numerator,
                     inner_coefficient(a[i]+a[j]), inner_coefficient(b[i]+b[j])))

    def affine(lo, hi, x, y):
      first, second = x+lo*(y-x), x+hi*(y-x)
      return min(first, second), max(first, second)

    def check(lo, hi):
      norms = [affine(lo, hi, F(2)**x, F(2)**y) for x, y in zip(a, b, strict=True)]
      lower = F(0)
      for i, j, amount_lo, amount_hi, inner_a, inner_b in data:
        assert amount_lo <= amount_hi
        if amount_lo >= 0:
          upper_den = (sqrt_bounds(norms[i][1]*norms[j][1])[1]
                       + affine(lo, hi, inner_a[1], inner_b[1])[1])
          lower += amount_lo/upper_den
        else:
          lower_den = (sqrt_bounds(norms[i][0]*norms[j][0])[0]
                       + affine(lo, hi, inner_a[0], inner_b[0])[0])
          lower += amount_lo/lower_den
      return lower

    return check


  def check_row(row, vertex, sign, edges, ticks):
    vertices = tuple(sorted({v for edge in edges for v in edge}))
    index = {v: i for i, v in enumerate(vertices)}
    ds = distances(len(index), tuple((index[a], index[b]) for a, b in edges))
    depth = ds[index[vertex]][index["o"]]
    cov_lo, cov_hi = critical_covariance(depth)
    weights = dict(zip(edges, (F(t, 100) for t in ticks), strict=True))
        # Fix the squared norms at the score and the endpoint of the first
        # negative geometric mean. The other negative edge is paid via u.
    paid_child = "w" if sign > 0 else "w_prime"
    child_weight = weights["u", paid_child]+int(vertex == "u")
    covariance_square = (1+F(depth, 3))**2 * F(2)**(-depth)
    ratio = F(5, 2)*covariance_square*(1/weights["o", "u"]**2+1/child_weight**2)
    assert ratio < 1

        # Energy minus s (D_vw' - D_vw - c (D_ow' - D_ow)).
    terms = [(edge, w, w) for edge, w in weights.items()]
    terms.extend((((vertex, "w_prime"), F(-sign), F(-sign)),
                  ((vertex, "w"), F(sign), F(sign)),
                  (("o", "w_prime"), min(sign*cov_lo, sign*cov_hi), max(sign*cov_lo, sign*cov_hi)),
                  (("o", "w"), min(-sign*cov_lo, -sign*cov_hi), max(-sign*cov_lo, -sign*cov_hi))))
    terms = [((index[a], index[b]), lo, hi) for (a, b), lo, hi in terms]
    covers = []
    for left, right in combinations(range(len(index)), 2):
      lower = interval_bound(ds, index[vertex], left, right, terms)
      pending = [(F(0), F(1), 0)]
      width = F(0)
      cells = maximum_depth = 0
      while pending:
        a, b, depth_now = pending.pop()
        if lower(a, b) >= 0:
          width += b-a
          cells += 1
          maximum_depth = max(maximum_depth, depth_now)
          continue
        assert depth_now < 25, (row, vertices[left], vertices[right], a, b)
        mid = (a+b)/2
        pending.extend(((a, mid, depth_now+1), (mid, b, depth_now+1)))
      assert width == 1
      covers.append({"attachments": [vertices[left], vertices[right]], "cells": cells,
                     "maximum_depth": maximum_depth})
    return {"row": row, "score": vertex, "sign": sign,
            "edge_weights": [[a, b, str(w)] for (a, b), w in weights.items()],
            "critical_cost": str(sum(weights.values())),
            "paid_cost": str(sum(weights.values())+PADDING),
            "negative_edge_hessian_ratio": str(ratio), "pair_covers": covers}


  @cache
  def paid_costs():
    costs = tuple(F(check_row(*row)['paid_cost']) for row in ROWS)
    return (*costs, costs[2])


  Q_ROLES = (None, 0, 4, 0, 1, 1, 2, 5, 1, 1, 3, 3, None, None)


  def revised_q_coefficient(v, j, depth, relaxed, critical, theta, paths, lo3, hi3, reflect,
                            ):
    role = Q_ROLES[v if j == 1 else reflect[v]]
    if role is None:
      return None

    cost = paid_costs()[role]
    path_cost = paths[depth]
    if role == 3 and relaxed:
      cost += F(188,100)
      path_cost *= 1+F(1, depth)
    critical *= hi3 if critical >= 0 else lo3
    return critical+hi3*(cost*theta+F(367,100)*path_cost*theta**2)

  return revised_q_coefficient, paid_costs


revised_q_coefficient, paid_costs = _fiber()
