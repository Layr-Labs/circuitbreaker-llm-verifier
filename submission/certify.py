"""Verify all numerical premises of the one-band circuit lower bound."""
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import sys
from bisect import bisect_right

script = Path(__file__).resolve()
here = script.parent / "files" if script == Path("/snapshot/entrypoint.py") else script.parent
sys.path.insert(0, str(here))
import one_band_price as price
from one_band_outer import (d7, outer, midpoint_strip_table, MARGIN_STEP_END,
                            adaptive_core_table, ROOT_TICKS, MID_STEP, MID_MAX_INDEX)


"""Standalone directed quadrature for the 47-term crossing-path bound."""
from decimal import Context as DecimalContext, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction
from math import cos, factorial, pi
from typing import NamedTuple


Q_LOW = DecimalContext(prec=58, rounding=ROUND_FLOOR)
Q_HIGH = DecimalContext(prec=58, rounding=ROUND_CEILING)
Q_MID = DecimalContext(prec=58, rounding=ROUND_HALF_EVEN)
GRID_BITS = 104
GRID = 1 << GRID_BITS
NODE_COUNT = 60
WINDOW = 3


class Enclosure(NamedTuple):
    low: Decimal
    high: Decimal

    @staticmethod
    def exact(value):
        if isinstance(value, Enclosure):
            return value
        if isinstance(value, Fraction):
            n, d = Decimal(value.numerator), Decimal(value.denominator)
            return Enclosure(Q_LOW.divide(n,d), Q_HIGH.divide(n,d))
        point = Decimal(value)
        return Enclosure(point, point)

    def __add__(self, other):
        b = self.exact(other)
        return Enclosure(Q_LOW.add(self.low,b.low), Q_HIGH.add(self.high,b.high))

    __radd__ = __add__

    def __neg__(self):
        return Enclosure(self.high.copy_negate(), self.low.copy_negate())

    def __sub__(self, other):
        return self + (-self.exact(other))

    def __rsub__(self, other):
        return self.exact(other) + (-self)

    def __mul__(self, other):
        b = self.exact(other)
        if self.low >= 0 and b.low >= 0:
            return Enclosure(Q_LOW.multiply(self.low,b.low), Q_HIGH.multiply(self.high,b.high))
        endpoints = tuple((x,y) for x in self for y in b)
        return Enclosure(min(Q_LOW.multiply(x,y) for x,y in endpoints),
                         max(Q_HIGH.multiply(x,y) for x,y in endpoints))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.exact(other)
        assert b.high < 0 or b.low > 0
        return self * Enclosure(Q_LOW.divide(1,b.high), Q_HIGH.divide(1,b.low))

    def __rtruediv__(self, other):
        return self.exact(other) / self

    def square(self):
        magnitudes = tuple(x.copy_abs() for x in self)
        minimum = Decimal(0) if self.low <= 0 <= self.high else min(magnitudes)
        maximum = max(magnitudes)
        return Enclosure(Q_LOW.multiply(minimum,minimum), Q_HIGH.multiply(maximum,maximum))

    def sqrt(self):
        assert self.low > 0
        return Enclosure(Q_LOW.next_minus(Q_MID.sqrt(self.low)),
                         Q_HIGH.next_plus(Q_MID.sqrt(self.high)))

    def exp(self):
        return Enclosure(Q_LOW.next_minus(Q_MID.exp(self.low)),
                         Q_HIGH.next_plus(Q_MID.exp(self.high)))

    def fixed(self):
        numerator_low, denominator_low = self.low.as_integer_ratio()
        numerator_high, denominator_high = self.high.as_integer_ratio()
        return (numerator_low*GRID//denominator_low,
                -((-numerator_high*GRID)//denominator_high))


def pi_enclosure():
    def atan(inv):
        values = [Fraction((-1)**k, (2*k+1)*inv**(2*k+1)) for k in range(80)]
        lower = sum(values)
        upper = lower + Fraction(1, 161*inv**161)
        return Enclosure(Enclosure.exact(lower).low, Enclosure.exact(upper).high)

    return 16*atan(5) - 4*atan(239)


def legendre_value_at_fraction(a, b):
    """Return the last two scaled Legendre values, x=a/b."""
    older, newer = 1, a
    bb = b*b
    for k in range(2, NODE_COUNT+1):
        older, newer = newer, (2*k-1)*a*newer - ((k-1)**2)*bb*older
    return older, newer


def approximate_root(index):
    from decimal import localcontext
    with localcontext(Q_MID):
        q = NODE_COUNT
        x = Decimal(str(cos(pi*(index+0.75)/(q+0.5))))
        for unused in range(7):
            p, next_p = Decimal(1), x
            for k in range(2, q+1):
                p, next_p = next_p, ((2*k-1)*x*next_p-(k-1)*p)/k
            x -= next_p*(x*x-1)/(q*(x*next_p-p))
        return x


def gauss_rule():
    rule = []
    last_left = Decimal(1)
    q = NODE_COUNT
    radius = Decimal('1e-45')
    for index in range(q):
        center = approximate_root(index)
        left = Q_LOW.subtract(center,radius)
        right = Q_HIGH.add(center,radius)
        assert -1 < left < right < last_left
        signs = [legendre_value_at_fraction(*bound.as_integer_ratio())[1]
                 for bound in (left,right)]
        assert signs[0]*signs[1] < 0
        last_left = left
        a, b = center.as_integer_ratio()
        previous, current = legendre_value_at_fraction(a,b)
        derivative = Enclosure.exact(Fraction(q*b*(a*current-q*b*b*previous),
                                               factorial(q)*b**q*(a*a-b*b)))
        width = max(Q_HIGH.subtract(center,left),Q_HIGH.subtract(right,center))
        error = Q_HIGH.multiply(width, q**4)
        derivative += Enclosure(error.copy_negate(),error)
        root = Enclosure(left,right)
        weight = WINDOW/((1-root.square())*derivative.square())
        rule.append(((root+1)*WINDOW/2, weight))
    return rule


def positive_arms(length):
    return 10**(length//2)*2**(length%2)


def tensor_error(dimension, at_zero):
    q = NODE_COUNT
    if dimension == 2:
        precision_moment = 2*9**q
    elif dimension == 3:
        precision_moment = 2*12**q + 33**q
    else:
        precision_moment = 2*12**q + 2*36**q + (dimension-4)*39**q
    weight = Fraction(WINDOW**(dimension+2*q)*factorial(q)**3*precision_moment,
                      (2*q+1)*2**q*factorial(2*q)**2*(4*WINDOW)**max(0,dimension-9))
    return (at_zero*weight).high


def path_series_enclosure():
    pi_range = pi_enclosure()
    root2 = Enclosure.exact(2).sqrt()
    density_1 = (6/pi_range).sqrt()
    density_2 = 3/(2*pi_range)
    rule = gauss_rule()
    nodes = [x for x,w in rule]
    weights = [w for x,w in rule]
    x2 = [x.square() for x in nodes]
    assert sum(w*(-s*Fraction(39,2)).exp() for w,s in zip(weights,x2)).high < Decimal('0.25')
    z_unary = [w*density_1*(-6*s).exp() for w,s in zip(weights,x2)]
    assert sum(z_unary).high < 1

    q = NODE_COUNT
    plain = [[None]*q for x in nodes]
    rooted = [[None]*q for x in nodes]
    for i in range(q):
        for k in range(i+1):
            xy = nodes[i]*nodes[k]
            plain[i][k] = plain[k][i] = (-6*xy).exp()
            rooted[i][k] = rooted[k][i] = (-6*root2*xy).exp()
    unary_a = [(-s*Fraction(3,2)).exp() for s in x2]
    unary_c = [(-12*s).exp() for s in x2]
    pair = [w*(-s*Fraction(9,2)).exp() for w,s in zip(weights,x2)]

    def tabulate(make):
        lower = []
        widest = 0
        for i in range(q):
            row = []
            for k in range(q):
                a,b = make(i,k).fixed()
                assert 0 <= a <= b
                row.append(a)
                widest = max(widest,b-a)
            lower.append(row)
        return lower, widest

    A, debit_a = tabulate(lambda i,k: unary_a[i]*rooted[i][k])
    B, debit_b = tabulate(lambda i,k: plain[i][k])
    C, debit_c = tabulate(lambda i,k: unary_c[i]*z_unary[k]*rooted[i][k].square())
    state = []
    debit = 0
    for i in range(q):
        row = []
        for k in range(q):
            a,b = (density_2*pair[i]*pair[k]*rooted[i][k]).fixed()
            assert 0 <= a <= b
            row.append(a)
            debit += b-a
        state.append(row)
    assert sum(map(sum,state)) + debit <= GRID
    step_debit = debit_a+debit_b+q*debit_c+2*q*q

    total = Enclosure.exact(0)
    pair_tail = density_2*Enclosure.exact(Fraction(-9*WINDOW**2,2)).exp()/(54*root2*WINDOW**2)
    zero = density_2
    for j in range(1,26):
        n = j+1
        mass = sum(map(sum,state))
        low = Q_LOW.divide(mass,GRID)
        high = Q_HIGH.divide(mass+debit,GRID)
        quad = tensor_error(n,zero)
        arms = sum(Fraction(1,positive_arms(min(i,n-2))*positive_arms(n-2-min(i,n-2)))
                   for i in range(n))
        error_high = Q_HIGH.add(quad,(pair_tail*arms).high)
        orthant = Enclosure(Q_LOW.subtract(low,quad),Q_HIGH.add(high,error_high))
        total += orthant * (3*2**(j-1)*((-1)**(j+1)))
        if j == 24:
            even_lower = total.low
        if j == 25:
            break

        next_state = []
        for middle in range(q):
            # B is symmetric, so B[last] indexes the old coordinate.
            active = [(i,(state[i][middle]*A[i][middle]) >> GRID_BITS)
                      for i in range(q) if state[i][middle]]
            active = [(i,a) for i,a in active if a]
            row = []
            for last, kernel in enumerate(B):
                contraction = sum(a*kernel[i] for i,a in active)
                row.append((C[middle][last]*contraction) >> (2*GRID_BITS))
            next_state.append(row)
        state = next_state
        debit += step_debit
        zero *= density_1

    assert total.high < Decimal('0.13064756648')
    return {"gamma_upper": str(total.high), "gamma_even_lower": str(even_lower),
            "initial_and_accumulated_debit_grid_units": debit,
            "per_step_debit_grid_units": step_debit}


from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial, isqrt, lcm, prod

UNIT = 1000000000000000000
EDGES = ((0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (4, 8), (4, 9), (6, 10), (6, 11))

def base_components(progress=None):
    path = path_series_enclosure()
    assert F(path["gamma_upper"]) < F('1633094581/12500000000')
    if progress is not None:
        progress("path", path)
    return {"path": path}


base = SimpleNamespace(base_components=base_components)


"""Exact radical bounds."""
from math import isqrt
INTERVAL_UNIT=1<<96
INTERVAL_ROWS=(('root_Q', '011213', (148, 0, 218), 1, 3, (2,), 367), ('support_Q', '01122324', (140, 349, 0, 294), 1, 4, (3,), 785), ('private_R', '01122324', (330, 860, 585, 0), 2, 3, (1, 4), 1780), ('private_Q', '0112233435', (221, 333, 507, 0, 417), 1, 5, (4,), 1480), ('ending_Q', '011223', (350, 350, 350), 1, 3, (1,), 1053), ('opposite_ending', '011223344546', (286, 565, 415, 726, 0, 605), 1, 6, (5,), 2600))

def ceiling(n, d):
  return -(-n // d)

def power(e):
  assert abs(e) < 96
  return INTERVAL_UNIT << e if e >= 0 else INTERVAL_UNIT >> -e

def half_power(e):
  q, odd = divmod(e, 2)
  if not odd:
    v = power(q)
    return v, v
  r = isqrt(2 * INTERVAL_UNIT * INTERVAL_UNIT)
  if q >= 0:
    return r << q, (r + 1) << q
  return r >> -q, ceiling(r + 1, 1 << -q)

def distances(edges):
  n = 1 + max(map(max, edges))
  rows = [[0 if i == j else n for j in range(n)] for i in range(n)]
  for i, j in edges:
    rows[i][j] = rows[j][i] = 1
  for k in range(n):
    for i in range(n):
      for j in range(n):
        rows[i][j] = min(rows[i][j], rows[i][k] + rows[k][j])
  assert max(map(max, rows)) < n
  return rows

def cover(ds, edges, weights, factor, target, subtract, left, right):
  n = len(ds)
  a = [ds[0][left] - ds[i][left] for i in range(n)]
  b = [ds[0][right] - ds[i][right] for i in range(n)]
  norms = [(power(x), power(y)) for x, y in zip(a, b)]
  pairs = [(0, target), *edges, *((0, j) for j in subtract)]
  coeffs = [100 * factor, *weights, *([100] * len(subtract))]
  data = []
  for (i, j), weight in zip(pairs, coeffs):
    x, y = a[i] + b[j], a[j] + b[i]
    assert (x + y) % 2 == 0
    determinant = power(x) + power(y) - 2 * power((x + y) // 2)
    assert determinant >= 0
    lo0, hi0 = half_power(a[i] + a[j])
    lo1, hi1 = half_power(b[i] + b[j])
    data.append((i, j, weight * determinant, lo0, lo1, hi0, hi1))
  if data[0][2] == 0:
    return 1, 0, 0
  pending = [(0, 0)]
  width = count = depth_max = maximum = 0
  while pending:
    index, depth = pending.pop()
    scale = 1 << depth
    bounds = []
    for x, y in norms:
      assert (y - x) % scale == 0
      u = x + (y - x) * index // scale
      v = u + (y - x) // scale
      bounds.append((min(u, v), max(u, v)))
    values = []
    for pos, (i, j, numerator, l0, l1, h0, h1) in enumerate(data):
      if not numerator:
        values.append(0)
        continue
      if pos == 0:
        inner = min(l0 + (l1-l0)*index//scale, l0 + (l1-l0)*(index+1)//scale)
        denominator = isqrt(bounds[i][0] * bounds[j][0]) + inner
        values.append(ceiling(numerator * INTERVAL_UNIT, denominator))
      else:
        inner = max(h0 + ceiling((h1-h0)*index, scale),
              h0 + ceiling((h1-h0)*(index+1), scale))
        denominator = isqrt(bounds[i][1] * bounds[j][1]) + 1 + inner
        values.append(numerator * INTERVAL_UNIT // denominator)
      assert denominator > 0
    energy = sum(values[1:1+len(edges)])
    upper = max(0, values[0] - sum(values[1+len(edges):]))
    assert energy > 0
    if upper < energy:
      ratio = ceiling(upper * INTERVAL_UNIT, energy)
      assert ratio < INTERVAL_UNIT
      maximum = max(maximum, ratio)
      count += 1
      width += 1 << (30-depth)
      depth_max = max(depth_max, depth)
    else:
      assert depth < 30
      pending.extend(((2*index, depth+1), (2*index+1, depth+1)))
  assert width == 1 << 30
  return count, depth_max, maximum

def parent_panels():
  endpoints = (0, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128,
        192, 256, 384, 512, 768, 1024, 1280, 1536, 1792, 2048)
  maximum = 0
  for l, r in zip(endpoints, endpoints[1:]):
    l, r = l * INTERVAL_UNIT // 2048, r * INTERVAL_UNIT // 2048
    cl, ch = isqrt((INTERVAL_UNIT+3*l)*INTERVAL_UNIT), isqrt((INTERVAL_UNIT+3*r)*INTERVAL_UNIT)+1
    dl, dh = isqrt((INTERVAL_UNIT+15*l)*INTERVAL_UNIT), isqrt((INTERVAL_UNIT+15*r)*INTERVAL_UNIT)+1
    jl = isqrt((10*INTERVAL_UNIT-6*r)*INTERVAL_UNIT)
    bracket = ceiling(27*INTERVAL_UNIT*INTERVAL_UNIT, jl+cl) - 6*INTERVAL_UNIT*INTERVAL_UNIT//(dh+2*ch)
    assert bracket > 0
    first = ceiling(2*(ch+INTERVAL_UNIT+r)*INTERVAL_UNIT, dl+cl)
    upper = ceiling(first*bracket, INTERVAL_UNIT)
    assert upper < 9*INTERVAL_UNIT
    maximum = max(maximum, upper)
  return 22, str(maximum)

def replay():
  rows = []
  for name, text, weights, factor, target, subtract, paid in INTERVAL_ROWS:
    edges = tuple((int(text[i]), int(text[i+1])) for i in range(0, len(text), 2))
    assert sum(weights) < paid
    if name == "opposite_ending":
      assert sum(weights) + weights[0] < 2890
    ds = distances(edges)
    controls = [cover(ds, edges, weights, factor, target, subtract, j, k)
          for j in range(len(ds)) for k in range(j+1, len(ds))]
    rows.append((name, len(controls), sum(c[0] for c in controls),
          max(c[1] for c in controls), str(max(c[2] for c in controls))))
  paths = []
  for d, critical, paid in ((2, 218, 218), (3, 373, 374), (4, 585, 586)):
    edges = tuple((i, i+1) for i in range(d))
    count, depth, maximum = cover(distances(edges), edges, (critical,)*d, 1, d, (), 0, d)
    assert maximum < INTERVAL_UNIT and critical <= paid
    paths.append((d, critical, paid, count, depth, str(maximum)))
  return {"unit": str(INTERVAL_UNIT), "rows": rows, "paths": paths, "parent": parent_panels()}
def certificate():
    return replay()


source_intervals = SimpleNamespace(certificate=certificate)


"""Exact signed-box probabilities with computed choices and no witness file."""
from math import factorial, isqrt, lcm
import json

S = 1000000000
ROWS = ((94371, -15711, 1764, 315), (83726, -21817, 2787, 917), (83726, 47472, -27708, 917), (71721, 50914, -38773, 1611), (-94280, -32378, 696, 315), (-83333, -44970, -27708, 917))
ETA = F(9,10000)
PRIVATE = ((9000000, 596995), (10500000, 670759), (11750000, 725075), (12520000, 755308), (14000000, 806692), (15375000, 846885), (16825000, 882028), (18250000, 910066))
DENSITY = (1007848001, 1000000001, 17235544)
SIGNED_NORMALIZER = F(6173823,100000000)
SIGNED_UNIT = SIGNED_NORMALIZER*F(1,4*10**85)
ROOT_EDGES = {-1: (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 95, 98, 101, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 149, 154, 159, 164, 169, 174, 180, 186, 192, 198, 204, 210, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270, 276, 282, 291, 300, 310, 320, 332, 344, 356, 368, 382, 396, 412, 428, 446, 464, 484, 504, 527, 550, 575, 600), 1: (0, 8, 16, 24, 32, 40, 48, 58, 68, 78, 88, 98, 108, 120, 132, 146, 162, 178, 194, 210, 230, 252, 278, 308, 344, 384, 430, 484, 544, 600)}
MARGIN_STEP = 4000000

def weight(q,scale=S):
 v=(scale*(1-q/4096)).__floor__()
 for _ in range(12):v=v*v//scale
 return v

@lru_cache(maxsize=1)
def strip_table():
 den=lcm(*[2**j*factorial(j)*(2*j+1)for j in range(16)])
 coeff=[(-1)**j*(den//(2**j*factorial(j)*(2*j+1)))*10000**(2*(15-j))for j in range(16)]
 divisor=100000*den*10000**31;best=0;result=[]
 for j in range(0,30001,5):
  p=0
  for c in reversed(coeff):p=p*j*j+c
  v=79788*S*j*p//divisor
  jj=min(j,16000);pol=1680*10000**6-280*jj*jj*10000**4+42*jj**4*10000**2-5*jj**6
  coarse=797*jj*pol*S//(10000*1680000*10000**6)
  best=max(best,v,coarse);assert 0<=best<S;result.append(best)
 return tuple(result[j//5] for j in range(30001))

def diamond(first,second,third,half):
 c=min(300000000,second-half-(1415*third+999)//1000)
 if c<=0:return 0
 u=min(300000000,first-third,1000*c//1415)
 if u<=0:return 0
 v=c-(1415*u+999)//1000
 moment=u*u*(c+3*v)+(c+v)*(c*c+v*v)
 den=23439*10**18;pen=(12500*u*moment+den-1)//den
 return max(0,(u*(c+v)//2-pen)*84//11)

def shifted_diamond(first, second, third, half):
    """Lower probability at scale 10^18 for one shifted private trapezoid."""
    cap = 300000000
    closing = min(cap, second-(1415*third+999)//1000)
    if closing <= 0:
        return 0
    root_guard = first-third-(707107*half+999999)//1000000
    width = min(cap, root_guard, 1000*closing//1415)
    if width <= 0:
        return 0
    edge = closing-(1415*width+999)//1000
    moment = width**2*(closing+3*edge)+(closing+edge)*(closing**2+edge**2)
    denominator = 2337*10**18
    penalty = (1250*width*moment+denominator-1)//denominator
    return max(0, (width*(closing+edge)//2-penalty)*70000//9207)


def shifted_private_covariance(eta):
    lo2, hi2, hi3 = F(141421, 100000), F(70711, 50000), F(1733, 1000)
    assert lo2**2 < 2 < hi2**2 and hi3**2 > 3
    k2, k3 = F(109, 25), F(561, 50)
    # The unconditioned private defect covariances with R and Q have
    # magnitudes at most 4D and sqrt(3)*(7sqrt(2)/4)*D.
    root_regression = 3*hi2*(2+hi3*F(367, 100))*eta
    driver_row = 3*(4+7*hi3*hi2/4)+root_regression
    assert driver_row < 25
    assert 25+25*hi2/4 < 34

    root_o, root_w = 2*eta, k2*eta/2
    driver_inverse = F(50,47)
    raw_l = F(1, 12)-(k2/2+3*hi2/4)*eta
    raw_u = F(1, 12)+(3*hi2/8+3*(3-2*lo2)/8)*eta
    variance_l = raw_l-root_o**2-(34*eta)**2*driver_inverse
    variance_w_l = F(1, 12)-k2*eta/2-root_w**2-(59*eta)**2*driver_inverse
    assert variance_l > F(77,1000) and variance_w_l > F(76,1000)
    off = ((1+k3)/4+F(3, 8))*eta+root_o*root_w+34*59*eta**2*driver_inverse
    assert off < F(6,1000)
    schur = min(variance_l-off**2/variance_w_l, variance_w_l-off**2/variance_l)
    assert schur > F(779,10000)
    deviation_product = F(837, 10000)
    assert raw_u/12 < deviation_product**2
    assert F(70000, 9207)*deviation_product == F(7, 11)
    assert F(1250,2337) == 1/(24*F(779,10000))
    assert 2*F(707107, 10**6)**2 > 1 and F(1415, 1000)**2 > 2


def private_mass(first,second,code,strips):
 third=0
 if code >= 32:
  half,mass=PRIVATE[code-32]
  return mass*shifted_diamond(first,second,third,half)
 idx,kind=divmod(code,4);half,mass=PRIVATE[idx]
 if kind==0:return mass*diamond(first,second,third,half)
 capacity=min(first-third,(second-half)*1000//1415-third)
 if capacity<=0:return 0
 balance=((second-half)*1000//1415-third)//2
 d1=(capacity//2,8*capacity//10,min(capacity-1,balance))[kind-1]
 d2=second-half-(1415*(d1+third)+999)//1000
 if min(d1,d2)<=0:return 0
 return mass*strips[min(30000,346*d1//10**7)]*strips[min(30000,346*d2//10**7)]

def prepare():
 strips=strip_table();joints={};a,b,c=DENSITY;den=4096*960000*10**9;U=10**18
 for dy,dz in ((24, 20), (30, 24)):
  joint=[]
  for j in range(0,1440,dy):
   row=[]
   for k in range(0,680,dz):
    N=a*(3*j*j+3*j*dy+dy*dy)+16*b*(3*k*k+3*k*dz+dz*dz)+6*c*(2*j+dy)*(2*k+dz)
    value=U-(N*U+den-1)//den
    for _ in range(12):value=value*value//U
    row.append(value*dy*dz)
   joint.append(row)
  joints[dy,dz]=joint
 return strips,joints


def policy_sum(side, widths, threshold, prepared, private_choice, midpoint, core_table):
    assert side in (-1, 1)
    strips, joints = prepared
    dy, dz = widths
    biases = [(threshold*h/F(11, 200)).__ceil__() if (i >= 4) == (side > 0) else 0
              for i, (*_, h) in enumerate(ROWS)]
    required = (1,1,MID_STEP,1,1,MID_STEP)
    edges = ROOT_EDGES[side]
    assert edges[0] == 0 and all(a < b for a, b in zip(edges, edges[1:]))
    shift = side*threshold
    bits = 144

    root_debits = tuple(((r+2)//3, (1415*((r+2)//3)+999)//1000) for r in ROOT_TICKS)

    first_debits = tuple(first for first,_ in root_debits)

    def face_mass(first, second, middle, terminal, root_guard, middle_guard):
        lower = min(MID_MAX_INDEX, min(middle, middle_guard)//MID_STEP)
        narrow = bisect_right(first_debits, root_guard)-1
        terminal_pair = midpoint[min(len(midpoint)-1, terminal//MARGIN_STEP_END)]
        sibling_only = 10**6*strips[min(30000, 346*terminal//10**7)]
        no_root_debit = private_choice(first//MARGIN_STEP, second//MARGIN_STEP)

        best_positive = best_two_sign = 0
        for i, ((first_debit, second_debit), row) in enumerate(
                zip(root_debits, core_table, strict=True)):
            if first_debit >= first:
                break
            full_debit = private_choice((first-first_debit)//MARGIN_STEP,
                                        (second-second_debit)//MARGIN_STEP)
            best_positive = max(best_positive, full_debit*row[lower])
            if i <= narrow:
                two_sign = row[lower]*(full_debit*sibling_only+no_root_debit*terminal_pair)//2
                best_two_sign = max(best_two_sign, two_sign)

        negative = core_table[narrow][lower]*no_root_debit if narrow >= 0 else 0
        split_root = (best_positive+negative)//2*terminal_pair
        return max(split_root, best_two_sign) >> bits

    def face_layer(x):
        values = []
        for y in range(0, 1441, dy):
            offsets, slopes = [], []
            low, high = 0, 680
            for (a, b, c, _), bias, req in zip(ROWS, biases, required, strict=True):
                offset = -(5*a*x+25*b*y+10000*bias)
                slope = -100*c
                if slope > 0:
                    low = max(low, -((offset-req)//slope))
                else:
                    high = min(high, (offset-req)//(-slope))
                offsets.append(offset)
                slopes.append(slope)
            start = -(-low//dz)*dz
            row = []
            f,g,m,t,u,v = (a+b*start for a,b in zip(offsets,slopes,strict=True))
            df,dg,dm,dt,du,dv = (value*dz for value in slopes)
            for _ in range(start//dz,high//dz+1):
                row.append(face_mass(f,g,m,t,u,v))
                f+=df;g+=dg;m+=dm;t+=dt;u+=du;v+=dv
            values.append((start//dz, row))
        products = []
        for (a0,row0),(a1,row1) in zip(values,values[1:]):
            start = max(a0,a1)
            end = min(a0+len(row0),a1+len(row1))-1
            products.append((start,[row0[k-a0]*row0[k+1-a0]*row1[k-a1]*row1[k+1-a1]
                                    for k in range(start,end)]))
        return products

    total = count = 0
    lower = face_layer(edges[0])
    for xlo, xhi in zip(edges, edges[1:]):
        left, right = F(xlo, 2000), F(xhi, 2000)
        q = max((left*left+left*right+right*right)/3,
                ((left+shift)**2+(left+shift)*(right+shift)+(right+shift)**2)/3)/2
        xweight = weight(q)
        upper = face_layer(xhi)
        root_sum = 0
        for y, joint_row in enumerate(joints[dy,dz]):
            a0, row0 = lower[y]
            a1, row1 = upper[y]
            start = max(a0,a1)
            end = min(a0+len(row0),a1+len(row1))
            for k in range(start,end):
                product = row0[k-a0]*row1[k-a1]
                root_sum += joint_row[k]*isqrt(isqrt(isqrt(product)))
            count += max(0,end-start)
        total += (xhi-xlo)*xweight*root_sum
        lower = upper
    return SIGNED_UNIT*(total << bits)/2, count


def signed_cutoff_controls(eta, rows, density, normalizer, private):
    lo2, hi2 = F(1414213562373, 10**12), F(1414213562374, 10**12)
    lo3, hi3 = F(1732050807568, 10**12), F(1732050807569, 10**12)
    assert lo2**2 < 2 < hi2**2 and lo3**2 < 3 < hi3**2
    k2, k3 = F(109, 25), F(561, 50)
    paths = (F(0), F(1), k2, k3)
    lower = (F(0), 2*lo2/3, F(5, 6), lo2/2)
    upper = (F(0), 2*hi2/3, F(5, 6), hi2/2)
    r_abs, q_abs = 2*eta, hi3*F(367, 100)*eta
    error_y = (2*k2, F(89, 5), 2, 2*k2, 2*k2, 2*k3)
    error_z = (F(157, 20), F(74, 5), k2, F(1053, 100), 1, k2)
    ideal_y = (-F(1, 6), -lo2/6, hi2/3, F(1, 2), -F(1, 3), -lo2/3)
    ideal_z = (F(0), F(0), -1/(2*hi3), -lo2/(2*hi3), F(0), -1/(2*hi3))
    for i, d in enumerate((1, 2, 2, 3, 1, 2)):
        covariance = upper[d]+paths[d]*eta
        assert covariance < 1
        bounds = (covariance if i < 4 else -lower[d],
                  ideal_y[i]+error_y[i]*eta+covariance*r_abs,
                  ideal_z[i]+hi3*error_z[i]*eta+covariance*q_abs,
                  (1-lower[d])*F(11, 200))
        assert all(F(tick, 10**5) > bound for tick, bound in zip(rows[i], bounds, strict=True))

    off = hi3*F(1105, 100)*eta+r_abs*q_abs
    r_lo, q_lo = 1-8*k2*eta-r_abs**2, 1-6*k2*eta-q_abs**2
    assert min(r_lo, q_lo)-off > F(47,50)
    assert all(F(tick, 10**9) > bound for tick, bound in
               zip(density, (1+2*k2*eta, F(1), off), strict=True))
    pi_upper = F(355, 113)
    atan5 = sum((F((-1)**j, (2*j+1)*5**(2*j+1)) for j in range(9)), F(0))
    assert 16*atan5-4*(F(1, 239)-F(1, 3*239**3)) < pi_upper
    assert normalizer**2*(2*pi_upper)**3 < r_lo*q_lo-off**2

    sd, prefactor = F(5381,250000), F(79788, 10**5)
    assert sd**2 > 3*(3-2*lo2)*eta and prefactor**2*pi_upper < 2
    for half,mass in private:
        z = F(2*half, 10**9)/sd
        probability = prefactor*sum((F((-1)**j, 2**j*factorial(j)*(2*j+1))*z**(2*j+1)
                                     for j in range(16)), F(0))
        assert 0 < F(mass, 10**6) < probability


def private_covariance_controls():
    e, k2, k3 = ETA, F(109,25), F(561,50)
    base = F(1,12)-k2*e/2
    diagonal1 = base-(e/2)**2-(25*e)**2/F(47,50)
    diagonal2 = base-(k2*e/2)**2-(59*e)**2/F(47,50)
    a, c = F(4041,50000), F(1567,20000)
    assert diagonal1 > a and diagonal2 > c
    off = (1+k3)*e/4+k2*e**2/4+25*59*e**2/F(47,50)
    assert min(a-off**2/c,c-off**2/a) > F(7813,100000)
    assert F(84,11) == 24/F(22,7)
    assert F(12500,23439) == 1/(24*F(7813,100000))


def certificate(progress=None):
    signed_cutoff_controls(ETA, ROWS, DENSITY, SIGNED_NORMALIZER, PRIVATE)
    private_covariance_controls()
    shifted_private_covariance(ETA)
    prepared = prepare()
    strips = prepared[0]
    midpoint = midpoint_strip_table(ETA,strips)
    core_table = adaptive_core_table(ETA,strips)

    @lru_cache(maxsize=None)
    def private_choice(a,b):
        first,second = a*MARGIN_STEP,b*MARGIN_STEP
        return max(private_mass(first,second,code,strips) for code in range(40))

    records = []
    jobs = ((-1,(24,20),F(0),49426,F(225627942,10**14)),
            (-1,(24,20),F(49,2000),48625,F(197515123,10**14)),
            (1,(30,24),F(49,2000),9231,F(207450756,10**14)))
    for side, widths, threshold, count_expected, floor in jobs:
        total, count = policy_sum(side, widths, threshold, prepared, private_choice, midpoint, core_table)
        assert count == count_expected and total > floor
        row = {"side":side,"threshold_abs_max":str(threshold),"cells":count,
               "cell_widths":list(widths),"exact_lower_sum":str(total),"usable_floor":str(floor)}
        records.append(row)
        if progress is not None:
            progress(row)
    return {"scope":"exact signed envelopes on two threshold bands","eta":str(ETA),"records":records}


signed = SimpleNamespace(certificate=certificate)


# Fixed nine-key rational parameter factory; no driver-cell search.
from fractions import Fraction as Q
from math import factorial, prod, isqrt
ORDER_UNIT = 10 ** 9
RADIUS_ROWS = ((18000000, 763017, 640344), (25040000, 900043, 797437), (33650000, 972941, 913183), (28000000, 934161, 845818))
CALIBRATION = {'normal_prefactor_lower': '19947/25000', 'core_sd': '15221/1000000', 'terminal_sd': '393/20000'}

def root_ending_probability(root_cap, closing):
    c = min(300000000, closing)
    if c <= 0:
        return 0
    u = min(300000000, root_cap, 1000 * c // 1500)
    if u <= 0:
        return 0
    v = c - (1500 * u + 999) // 1000
    assert 0 <= v <= c
    moment = u * u * (c + 3 * v) + c * c * c + c * c * v + c * v * v + v * v * v
    denominator = 216 * 10 ** 18
    lower_area = u * (c + v) // 2 - (100 * u * moment + denominator - 1) // denominator
    return max(0, lower_area * 700 // 99)

def parameters(relaxed_endings=False):
    eta, threshold = (Q('9/20000'), Q('49/2000'))
    root2, root3 = (Q(283, 200), Q(1733, 1000))
    k = [Q(0, 1), Q(1, 1), Q(109, 25), Q(561, 50), Q(586, 25), Q(221, 5)]
    gamma = [2 * eta, root3 * Q(367, 100) * eta, root3 * Q(367, 100) * eta]
    low = [1 - 8 * k[2] * eta - gamma[0] ** 2, 1 - 6 * k[2] * eta - gamma[1] ** 2, 1 - 6 * k[2] * eta - gamma[2] ** 2]
    cross = root3 * Q(491, 20) * eta + gamma[0] * gamma[1]
    far = Q(353, 10) * eta + gamma[1] * gamma[2]
    lower = [[low[0], -cross, -cross], [-cross, low[1], -far], [-cross, -far, low[2]]]
    assert all((lower[i][i] > sum((abs(lower[i][j]) for j in range(3) if j != i)) for i in range(3)))
    unit = Q(63, 2500 * 400000 * ORDER_UNIT**5 * 1000000**6)
    assert ORDER_UNIT == S
    return {'strip': strip_table(), 'unit': unit, 'gamma': gamma, 'theta': eta,
            'threshold': threshold, 'paths': k, 'lower': lower}


nine_parameters = SimpleNamespace(parameters=parameters, RADIUS_ROWS=RADIUS_ROWS, CALIBRATION=CALIBRATION, root_ending_probability=root_ending_probability)


from math import factorial, isqrt
S=10**9
COUNT=30*60*820
def ceil(n,d=1):return (n+d-1)//d


def calibrate(p,data,radii,relaxed):
 assert p['theta']==F('9/20000') and p['threshold']==F('49/2000')
 assert len(radii)==4 and [row[0] for row in radii]==[18000000,25040000,33650000,28000000]
 pref=F(data['normal_prefactor_lower']);core_sd=F(data['core_sd']);terminal_sd=F(data['terminal_sd'])
 assert pref==F(19947,25000) and F(7071,5000)**2<2<F(283,200)**2 and F(173,50)**2<12
 atan5=sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(9)),F(0))
 atan239=F(1,239)-F(1,3*239**3)
 pi_upper=16*atan5-4*atan239
 assert pi_upper<F(355,113) and pref*pref*pi_upper<2
 variance=F(3)-2*F(7071,5000)
 lo2,hi2=F(141421,100000),F(70711,50000)
 assert lo2**2<2<hi2**2
 assert all((1+F(d,3))/lo2**d+d*d*hi2**(d-1)*(2*p['theta'])<1 for d in range(1,7))
 assert 9*F(943,1000)**2>8>9*F(9428,10000)**2
 assert 9*F(2828,6000)**2<2<9*F(2830,6000)**2
 assert 36*F(2828,12000)**2<2<36*F(2830,12000)**2
 assert 12*F(2885,10000)**2<1 and 6*F(4082,10000)**2<1
 assert 2*F(7072,10000)**2>1 and 2*(1-F(2929,10000))**2<1
 assert 9*(F(972,1000)-F(1,2))**2>2 and 2*(F(1208,1000)-F(1,2))**2>1
 assert 3*variance*p['theta']<core_sd**2 and 5*variance*p['theta']<terminal_sd**2
 for r,cm,tm in radii:
  assert all(isinstance(v,int) and 0<v<=10**6 for v in (cm,tm))
  for sd,mass in ((core_sd,cm),(terminal_sd,tm)):
   z=F(r,S)/sd
   polynomial=sum((F((-1)**j,2**j*factorial(j)*(2*j+1))*z**(2*j+1) for j in range(16)),F(0))
   assert pref*polynomial>F(mass,10**6)
 k2,k3,k4=p['paths'][2:5];assert (k2,k3,k4)==(F(109,25),F(561,50),F(586,25))
 theta=p['theta'];gamma=p['gamma'];r3=F(1733,1000)
 assert r3*r3>3
 eigen=min(p['lower'][i][i]-sum(abs(p['lower'][i][j]) for j in range(3) if i!=j) for i in range(3))
 assert eigen>F(93,100)
 first=k2+r3*(1+k3)+sum(gamma)/2
 p1,p2,p3,p4=(F(2),F(3,2)*k2,F(4,3)*k3,F(5,4)*k4) if relaxed else (F(1),k2,k3,k4)
 end=(p2+3*p4)/2+r3*((p1+p3)/2+(F(289,10) if relaxed else 26))+p3*sum(gamma)/2
 assert first<26 and end<113
 pair=(p2+p4)*theta/4+p3*theta*theta/4+26*113*theta*theta/F(93,100)
 assert pair<F(13,2000) and F(1,12)+F(13,2000)<F(9,100)


class Context:

    def __init__(self, relaxed):
        ns = vars(nine_parameters)
        p = ns['parameters'](relaxed)
        self.strip, self.unit = p['strip'], p['unit']
        self.theta, self.threshold = p['theta'], p['threshold']
        self.radii = ns['RADIUS_ROWS']
        self.pair = ns['root_ending_probability']
        calibrate(p, ns['CALIBRATION'], self.radii, relaxed)
        self.weights = [[self.radii[r][1] ** 4 * self.radii[e][2 if relaxed else 1] ** 2 for e in range(4)] for r in range(4)]
        self.mid = [ceil(972 * r[0], 1000) for r in self.radii]
        self.close = [[ceil(1208 * r[0] + 500 * e[0], 1000) for e in self.radii] for r in self.radii]

    def mass(self, d):
        return self.strip[min(30000, max(0, 346 * d // 10000000))]

    def candidate(self, g, r, e, root_kind, kind):
        first, branch, leaf, selected, el, er = g
        radius = self.radii[r][0]
        mid = self.mid[r]
        closing = self.close[r][e]
        if first <= ceil(radius, 3) or leaf <= mid:
            return 0
        L, R = (el - closing, er - closing)
        cap = min(branch - ceil(radius, 3), (selected - mid) * 1000 // 1415, min(L, R) * 2 // 3)
        if cap <= 0:
            return 0
        if root_kind == 0:
            root = cap // 4
        elif root_kind == 1:
            root = cap // 2
        elif root_kind == 2:
            root = 3 * cap // 4
        else:
            disc = L * L - L * R + R * R
            rad = isqrt(disc)
            rad += rad * rad < disc
            root = min(cap, max(0, 2 * (L + R - rad) // 9))
        l, rend = (L - ceil(3 * root, 2), R - ceil(3 * root, 2))
        common = self.mass(leaf - mid) * self.weights[r][e]
        if kind == 0:
            return common * self.mass(root) * self.mass(l) * self.mass(rend)
        if kind == 1:
            return common * self.pair(root, L) * self.mass(rend)
        return common * self.pair(root, R) * self.mass(l)


nine_kernel = SimpleNamespace(COUNT=COUNT, Context=Context)


NORMALIZER = F(247, 10000)

SQRT2 = (F(1414213562373, 1000000000000), F(707106781187, 500000000000))

SQRT3 = (F(108253175473, 62500000000), F(1732050807569, 1000000000000))

PATH = (0, 1, F(109, 25), F(561, 50), F(586, 25), F(221, 5))

EDGES = ((0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9), (6, 10), (6, 11), (8, 12), (8, 13))

DRIVERS = (((2, 1), (3, 1), (1, -2)), ((7, 1), (6, -1)), ((9, 1), (8, -1)))

GUARDS = frozenset({9, 2, 3, 7})

REFLECT = (0, 1, 3, 2, 4, 5, 8, 9, 6, 7, 12, 13, 10, 11)

GROUPS = ((1,), (2, 3), (4, 5), (6, 7, 8, 9), (10, 11), (12, 13))

ADDRESSES = tuple((z, h) for z in range(40) for h in range(z, 40))


def covariance_bounds(theta):
    lo2, hi2 = SQRT2
    lo3, hi3 = SQRT3
    assert lo2**2 < 2 < hi2**2 and lo3**2 < 3 < hi3**2
    gamma = (2*theta, hi3*F(367,100)*theta, hi3*F(367,100)*theta)
    rr = 1-8*PATH[2]*theta-gamma[0]**2
    qq = 1-6*PATH[2]*theta-gamma[1]**2
    cross = hi3*F(491,20)*theta+gamma[0]*gamma[1]
    far = F(353,10)*theta+gamma[1]**2
    assert min(rr-2*cross, qq-cross-far) > F(93,100)
    determinant = rr*qq**2-rr*far**2-2*qq*cross**2-2*cross**2*far
    atan5 = sum((F((-1)**j, (2*j+1)*5**(2*j+1)) for j in range(9)), F(0))
    pi_upper = F(355,113)
    assert 16*atan5-4*(F(1,239)-F(1,3*239**3)) < pi_upper
    assert NORMALIZER**2*(2*pi_upper)**4 < determinant
    density = tuple((value*S).__ceil__() for value in (1+2*PATH[2]*theta,cross,far))
    return gamma, density, determinant


def tree_paths():
    graph = [[] for _ in range(14)]
    for i, (a,b) in enumerate(EDGES):
        terminal = int(i >= 9)
        graph[a].append((b,terminal))
        graph[b].append((a,terminal))
    paths = []
    for root in range(14):
        entries, pending = {root:(0,0)}, [root]
        for a in pending:
            length, terminals = entries[a]
            for b, terminal in graph[a]:
                if b not in entries:
                    entries[b] = length+1, terminals+terminal
                    pending.append(b)
        assert len(entries) == 14
        paths.append(entries)
    return paths


def mean_rows(theta, extent, gamma, relaxed):
    paths = tree_paths()
    lo2, hi2 = SQRT2
    lo3, hi3 = SQRT3
    lo = tuple((1+F(d,3))/hi2**d for d in range(6))
    hi = tuple((1+F(d,3))/lo2**d for d in range(6))

    def signed_cov(v, w, sign):
        length, terminals = paths[v][w]
        assert terminals <= 1
        if sign < 0:
            return lo[length]
        error = PATH[length]*theta
        if relaxed and length:
            error *= 1+F(terminals,length)
        return hi[length]+error

    rows = []
    for v in range(1, 14):
        sign = -1 if v in GUARDS else 1
        depth = paths[v][0][0]
        root_cov = signed_cov(v,0,1)
        assert root_cov < 1
        values = [sign*signed_cov(v,0,sign)]
        for j, driver in enumerate(DRIVERS):
            upper = sum(sign*m*signed_cov(v,u,sign*m) for u,m in driver)
            if j:
                # The two endpoints are siblings in the other branch.
                if all(paths[v][u][0] == depth+2 for u,_ in driver):
                    costs = (F(785,100), F(74,5), F(289,10) if relaxed else F(26))
                    upper = min(upper, costs[depth-1]*theta)
                upper *= hi3 if upper >= 0 else lo3
            values.append(upper+root_cov*gamma[j])
        values.append(extent*(1-lo[depth]))
        rows.append(tuple((value*S).__ceil__() for value in values))
    for v, (x,y,z,h,bias) in enumerate(rows, 1):
        assert rows[REFLECT[v]-1] == (x,y,h,z,bias)
    return tuple(rows)


def grid():
    root = tuple((i*5000000, (i+1)*5000000) for i in range(20))
    root += tuple((i*10000000, (i+1)*10000000) for i in range(10,20))
    breadth = tuple((i*50000000, (i+1)*50000000) for i in range(60))
    ending = tuple((i*100000000, (i+1)*100000000) for i in range(40))
    return root,breadth,ending


def cap(coeff, interval):
    left, right = interval
    return coeff*(right if coeff >= 0 else left)


def margin_tables(rows, root, breadth, ending):
    xy = tuple(tuple(tuple(cap(a,x)+cap(b,y)+bias*S for a,b,_,_,bias in rows)
                     for y in breadth) for x in root)
    zh = tuple(tuple(cap(c,ending[z])+cap(d,ending[h]) for _,_,c,d,_ in rows)
               for z,h in ADDRESSES)
    return xy,zh


def shape(margins):
    return tuple(min(margins[v-1] for v in group) for group in GROUPS)


def density_tables(extent, coefficients, root, breadth, ending):
    extent_ticks = extent*S
    assert extent_ticks.denominator == 1
    shift = 2*extent_ticks.numerator
    yy, cross, far = coefficients
    xx = tuple(S*(3*(lo+hi+shift)**2+(hi-lo)**2) for lo,hi in root)
    y = tuple((yy*(3*(lo+hi)**2+(hi-lo)**2),6*cross*(lo+hi)) for lo,hi in breadth)
    zh = []
    for z,h in ADDRESSES:
        zl,zhigh = ending[z]
        hl,hhigh = ending[h]
        zm, hm = zl+zhigh, hl+hhigh
        weight = S*(3*(zm**2+hm**2)+(zhigh-zl)**2+(hhigh-hl)**2)
        zh.append((weight+6*far*zm*hm,zm+hm,1 if z==h else 2))
    return xx,y,zh


def lower_density(numerator):
    # The box-mean quadratic is numerator/(12*S^3). Apply
    # (1-q/4096)^4096 to half of that quadratic.
    denominator = 8192*12*S*S
    value = S-(numerator+denominator-1)//denominator
    assert value > 0
    for _ in range(12):
        value = value*value//S
    return value


def precision_nine_replay(strict, relaxed, codes, extent, rows, coefficients):
    root,breadth,ending = grid()
    xy,zh = margin_tables(rows[0],root,breadth,ending)
    xy_relaxed,zh_relaxed = margin_tables(rows[1],root,breadth,ending)
    xx,dense_y,dense_zh = density_tables(extent,coefficients,root,breadth,ending)
    sum_strict = sum_relaxed = count_strict = count_relaxed = 0
    for i, code in enumerate(codes):
        if code == 0:
            continue
        x_y,z_h = divmod(i,820)
        x,y = divmod(x_y,60)
        assert 1 <= code <= 192
        choice, kind = divmod(code-1,3)
        choice, root_kind = divmod(choice,4)
        r,e = divmod(choice,4)
        margins = [-((a+b+S-1)//S) for a,b in zip(xy[x][y],zh[z_h],strict=True)]
        # The four terminal leaves are the only scores whose paths to
        # a driver can contain a terminal edge.
        assert len(margins) == 13
        wider = margins[:9]+[-((a+b+S-1)//S) for a,b in
                              zip(xy_relaxed[x][y][9:],zh_relaxed[z_h][9:],strict=True)]
        sm = strict.candidate(shape(margins),r,e,root_kind,kind)
        rm = relaxed.candidate(shape(wider),r,e,root_kind,kind)
        if sm or rm:
            y_part,linear = dense_y[y]
            z_part,z_sum,multiple = dense_zh[z_h]
            density = lower_density(xx[x]+y_part+linear*z_sum+z_part)
            density *= multiple*(1 if x<20 else 2)
            sum_strict += density*sm
            sum_relaxed += density*rm
            count_strict += sm>0
            count_relaxed += rm>0
    return (sum_strict,count_strict),(sum_relaxed,count_relaxed)


def _nine():
    from fractions import Fraction as F
    import lzma
    COUNT = nine_kernel.COUNT
    Context = nine_kernel.Context
    LIMIT=65536

    def reorder_choices(codes, *, triangle_first):
        """Convert between x,y,z/h order and z/h,y,x order."""
        assert len(codes) == COUNT
        result = bytearray(COUNT)
        stride = 60*820
        for triangle in range(820):
            for y in range(60):
                start = y*820+triangle
                xy = slice(start, start+30*stride, stride)
                start = (triangle*60+y)*30
                triangular = slice(start, start+30)
                if triangle_first:
                    result[triangular] = codes[xy]
                else:
                    result[xy] = codes[triangular]
        return bytes(result)

    def decode(encoded):
        assert 0 < len(encoded) <= LIMIT and encoded.endswith(b"\n")
        data = encoded[:-1]
        assert len(data) % 11 == 0
        packed = bytearray()
        for lo in range(0, len(data), 11):
            value = 0
            for digit in data[lo:lo+11]:
                assert 33 <= digit <= 126
                value = value*94+digit-33
            assert value < 2**72
            packed.extend(value.to_bytes(9, "big"))
        decoder = lzma.LZMADecompressor(memlimit=64*1024**2)
        codes = decoder.decompress(packed, COUNT+1)
        assert decoder.eof and len(codes) == COUNT
        assert len(decoder.unused_data) < 9 and not any(decoder.unused_data)
        return reorder_choices(codes, triangle_first=False)

    def certificate(encoded, progress=None):
        codes = decode(encoded)
        strict, relaxed = Context(False), Context(True)
        assert strict.theta == relaxed.theta == price.ETA/2
        extent = F(49,2000)
        assert strict.threshold == relaxed.threshold == extent
        gamma, density, determinant = covariance_bounds(strict.theta)
        rows = tuple(mean_rows(strict.theta, extent, gamma, flag) for flag in (False,True))
        assert rows[0][:9] == rows[1][:9]
        totals = precision_nine_replay(strict,relaxed,codes,extent,rows,density)
        assert strict.unit == relaxed.unit == F(63,2500*400000*S**5*10**36)
        unit = strict.unit*NORMALIZER/F(63,2500)
        results = []
        for phase,(mass,count),floor in zip(("strict","relaxed"),totals,
                                               (F(1673,10**9),F(1673,10**9)),strict=True):
            exact = mass*unit
            assert exact > floor
            row = {"phase":phase, "eta":str(price.ETA), "threshold_abs_max":str(extent),
                   "exact_lower_sum":str(exact), "used_floor":str(floor),
                   "positive_codes":count,
                   "precision_determinant_lower":str(determinant)}
            results.append(row)
            if progress is not None:
                progress(row)
        return results

    return SimpleNamespace(certificate=certificate)


nine = _nine()


def combined_threshold_quota(zero, endpoint, extent, logarithm, radius, curvature):
    """Lower-bound the signed probability plus the capped quadratic deficit."""
    assert 0 < endpoint < zero and 0 < radius <= extent and curvature > 0
    assert 0 < logarithm < 1
    assert endpoint*sum(logarithm**j/factorial(j) for j in range(7)) > zero
    rate = logarithm/extent
    linear = zero*rate
    quadratic = curvature+zero*(F(1,2)-logarithm/6)*rate**2
    # Completing the square gives a lower bound throughout the inner range.
    inner = zero-linear**2/(4*quadratic)
    outer = endpoint+curvature*radius**2
    return {"zero_floor":str(zero), "endpoint_floor":str(endpoint),
            "log_ratio_upper":str(logarithm), "quadratic_radius":str(radius),
            "quadratic_coefficient":str(curvature), "inner_floor":str(inner),
            "outer_floor":str(outer), "combined_floor":str(min(inner, outer))}


def mixed_allocation(adjacent, baseline2, arms2, long_paths, size_four, d7, base, signed, nine, outer):
    A,k = price.allocation_inputs(
        adjacent, baseline2, arms2, long_paths, size_four, d7)
    assert F(signed["eta"]) == price.ETA
    assert len(nine) == 2 and {row["phase"] for row in nine} == {"strict", "relaxed"}
    assert all(F(row["eta"]) == price.ETA and F(row["threshold_abs_max"]) == F(49,2000) for row in nine)
    caps = {row["phase"]: F(row["used_floor"]) for row in nine}
    assert price.WEAK_WEIGHT == 0

    extent = F(49,2000)
    zero, = [F(row["usable_floor"]) for row in signed["records"] if F(row["threshold_abs_max"]) == 0]
    endpoints = {row["side"]:F(row["usable_floor"]) for row in signed["records"]
                 if F(row["threshold_abs_max"]) == extent}
    assert set(endpoints) == {-1,1} and len(signed["records"]) == 3
    diamond = outer["actual_base_deficit"]
    assert F(diamond["marking_cutoff"]) == price.ETA
    assert F(diamond["threshold_abs_max"]) == F(101,5000)
    quota = combined_threshold_quota(zero, min(endpoints.values()), extent,
                                    F(3327,25000), F(101,5000),
                                    F(diamond["deficit_coefficient"]))
    weights = (218978509, 167231461, 100, 100, 100, 1778578981, 394844746)
    assert F(weights[0],10**14) < F(quota["combined_floor"])
    # The existing q and d columns have the same constant and incidence debits.
    marking = price.marked_rows(A, weights, 335739986, 335740086, 205463189, k, caps, medium=(29115760,29115760), deficit=F(0))
    assert F(marking["credit"]) == F(2818956503,100000000000000)
    bands = [{"extent":str(extent), "signed_and_deficit_quota":quota,
              "marking":marking, "combined_credit":marking["credit"]}]
    assert F(outer["threshold_abs_max_of_central_band"]) == F(49,2000)
    assert F(base["path"]["gamma_upper"]) < price.PATH_UPPER
    saving = min(*(F(b["combined_credit"]) for b in bands), F(outer["outer_graph_credit"]))
    assert saving > price.SAVING and price.PATH_UPPER-saving < price.GAMMA
    return {
        'bands': bands,
        'price_allowance': str(A),
        'outer_credit': outer["outer_graph_credit"],
        'saving': str(saving),
        'circuit_coefficient': str(1+1/(2*price.GAMMA)),
    }


def main():
    import json
    from resource import RUSAGE_SELF, getrusage
    from time import monotonic, process_time

    started = monotonic()
    phases = []

    def progress(phase, value):
        row = {"phase": phase, "elapsed": monotonic()-started, "output": value}
        phases.append(row)
        print(json.dumps(row), file=sys.stderr, flush=True)

    def check(component, phase):
        result = component.certificate()
        progress(phase, result)
        return result

    progress("actual source intervals", source_intervals.certificate())
    base_result = base.base_components(progress)
    signed_result = signed.certificate(lambda row: progress("signed", row))
    choices = (here/"nine-choices.txt").read_bytes()
    assert sha256(choices).hexdigest() == "0b14a81944b97739be01f7f5d3a71d8e7def7dd2320c7fa8baabadfdf04f2bac"
    nine_result = nine.certificate(choices, lambda row: progress("nine " + row["phase"], row))
    d7_result = check(d7, "D7")
    outer_result = check(outer, "complete outer")
    adjacent = check(price.adjacent, "adjacent Price")
    baseline2 = check(price.baseline2, "D2 global baseline")
    arms2 = check(price.arms2, "D2 baseline and size-two arms")
    long_paths = check(price.long_paths, "D3-D6 Price")
    size_four = check(price.size_four, "D2-D3 size-four/six arms")
    final = mixed_allocation(adjacent, baseline2, arms2, long_paths, size_four,
                                              d7_result, base_result, signed_result, nine_result, outer_result)
    progress("exact final composition", final)
    print(json.dumps({"completed": "one-band numerical theorem", "elapsed": monotonic()-started,
                      "phases": phases,
                      "cpu_seconds": process_time(), "peak_rss_kib": getrusage(RUSAGE_SELF).ru_maxrss,
                      "circuit_coefficient": final["circuit_coefficient"]}), flush=True)


if __name__ == "__main__":
    main()
