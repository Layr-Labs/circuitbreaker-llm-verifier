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
                            adaptive_core_table, ROOT_TICKS, MID_STEP, MID_MAX_INDEX, root_policy, folded_private_bins)


"""Standalone directed quadrature for the 37-term crossing-path bound."""
from decimal import Context as DecimalContext, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction
from math import cos, factorial, pi
from typing import NamedTuple


Q_LOW = DecimalContext(prec=50, rounding=ROUND_FLOOR)
Q_HIGH = DecimalContext(prec=50, rounding=ROUND_CEILING)
Q_MID = DecimalContext(prec=50, rounding=ROUND_HALF_EVEN)
GRID = 1 << 128
NODE_COUNT = 62
WINDOW = Fraction(14,5)
PATH_UPPER = Fraction(".130647566423168")


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
  radius = Decimal('1e-43')
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
  for j in range(1,38):
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
    if j == 36:
      even_lower = total.low
    if j == 37:
      break

    next_state = []
    for middle in range(q):
            # B is symmetric, so B[last] indexes the old coordinate.
      active = [(i,(state[i][middle]*A[i][middle]) >> 128)
                for i in range(q) if state[i][middle]]
      active = [(i,a) for i,a in active if a]
      row = []
      for last, kernel in enumerate(B):
        contraction = sum(a*kernel[i] for i,a in active)
        row.append((C[middle][last]*contraction) >> 256)
      next_state.append(row)
    state = next_state
    debit += step_debit
    zero *= density_1

  assert Fraction(total.high) < PATH_UPPER
  return {"gamma_upper": str(total.high), "gamma_even_lower": str(even_lower),
          "initial_and_accumulated_debit_grid_units": debit,
          "per_step_debit_grid_units": step_debit}


from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial, isqrt, lcm, prod

UNIT = 1000000000000000000


def base_components(progress=None):
  path = path_series_enclosure()
  assert F(path["gamma_upper"]) < PATH_UPPER
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
ROWS = ((94395,-15833,1387,315),(83831,-21349,2646,917),(83831,47332,-27497,917),(71990,50917,-38606,1611),(-94280,-32338,724,315),(-83333,-44809,-27497,917))
ETA = F('57/50000')
DENSITY = (1003491987,1000000001,21835167)
SIGNED_NORMALIZERS = tuple(F(n,10**8) for n in (6136901,6150156,6163383))
SIGNED_NORMALIZER = SIGNED_NORMALIZERS[-1]
SIGNED_UNIT = SIGNED_NORMALIZER*F(1,4*10**85)
ROOT_EDGES = {-1: (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 95, 98, 101, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 149, 154, 159, 164, 169, 174, 180, 186, 192, 198, 204, 210, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270, 276, 282, 291, 300, 310, 320, 332, 344, 356, 368, 382, 396, 412, 428, 446, 464, 484, 504, 527, 550, 575, 600), 1: (0, 8, 16, 24, 32, 40, 48, 58, 68, 78, 88, 98, 108, 120, 132, 146, 162, 178, 194, 210, 230, 252, 278, 308, 344, 384, 430, 484, 544, 600)}
MARGIN_STEP = 1000000

def weight(q,scale=S):
 v=(scale*(1-q/4096)).__floor__()
 for _ in range(12):v=v*v//scale
 return v

@lru_cache(maxsize=1)
def strip_table():
 den=lcm(*[2**j*factorial(j)*(2*j+1)for j in range(16)])
 coeff=[(-1)**j*(den//(2**j*factorial(j)*(2*j+1)))*10000**(2*(15-j))for j in range(16)]
 divisor=100000*den*10000**31;best=0;result=[]
 for j in range(30001):
  p=0
  for c in reversed(coeff):p=p*j*j+c
  v=79788*S*j*p//divisor
  jj=min(j,16000);pol=1680*10000**6-280*jj*jj*10000**4+42*jj**4*10000**2-5*jj**6
  coarse=797*jj*pol*S//(10000*1680000*10000**6)
  best=max(best,v,coarse);assert 0<=best<S;result.append(best)
 return result

@lru_cache(maxsize=32768)
def clipped_diamond(u,c):
  v=c-(1415*u+999)//1000
  if not u or not c+v:
    return 0
  moment=u*u*(c+3*v)+(c+v)*(c*c+v*v)
  denominator=4629*(c+v)*10**18
  unit=10**12
  q=(5000*moment*unit+denominator-1)//denominator
  assert 0<=q<unit
  factor=unit-q+q*q*(3*unit-q)//(6*unit**2)
  return u*(c+v)*factor*42//(11*unit)

def private_policy(ordinary,shifted):
  @lru_cache(maxsize=None)
  def geometry(second):
    rows = []
    for half,tight,weak in ordinary:
      c = min(500000000,second-half)
      if c <= 0:
        break
      rows.append((c,1000*c//1415,tight,weak))
    c = min(500000000,second)
    return rows,c,1000*c//1415

  @lru_cache(maxsize=None)
  def private_choice(a,b):
    if a <= 0 or b <= 0:
      return 0
    first = a*MARGIN_STEP
    rows,c,limit = geometry(b*MARGIN_STEP)
    tight_sum = ordinary_sum = 0
    for width,cap,tight,weak in rows:
      area = clipped_diamond(first if first < cap else cap,width)
      tight_sum += tight*area
      ordinary_sum += weak*area
    one = min(tight_sum,ordinary_sum*50//47)*209//210

    tight_sum = ordinary_sum = 0
    for debit,tight,weak in shifted:
      width = first-debit
      if width <= 0:
        break
      area = clipped_diamond(width if width < limit else limit,c)
      tight_sum += tight*area
      ordinary_sum += weak*area
    two = min(tight_sum,ordinary_sum*50//47)*4169//4200
    return max(one,two)
  return private_choice


def prepare():
 strips=strip_table();joints={};a,b,c=DENSITY;den=4096*960000*10**9;U=10**18
 for dy,dz in ((8, 4), (30, 16)):
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


def policy_sum(side, widths, threshold, prepared, private_choice, midpoint, core_table, caps):
  assert side in (-1, 1)
  strips, joints = prepared
  dy, dz = widths
  upper_corr = tuple(ROWS[i][0] for i in (0,1,1,3,0,1))
  biases = [(threshold*h/F(11,200)).__ceil__() if (i >= 4) == (side > 0)
            else (-(100000-upper_corr[i])*threshold).__ceil__()
            for i,(*_,h) in enumerate(ROWS)]
  required = (1,1,MID_STEP,1,1,MID_STEP)
  edges = ROOT_EDGES[side]
  assert edges[0] == 0 and all(a < b for a, b in zip(edges, edges[1:]))
  shift = side*threshold
  bits = 144

  first_debits,first_shifts,second_shifts,root_payments = root_policy(
      private_choice,tuple(map(tuple,core_table)),MARGIN_STEP)
  first_cap = caps[0] + first_debits[-1]
  second_cap = caps[1] + (1415*first_debits[-1]+999)//1000
  mid_cap = MID_MAX_INDEX*MID_STEP
  last_radius = len(first_debits)-1
  full_guard = first_debits[-1]
  terminal_limit = len(midpoint)-1

  def face_mass(first, second, middle, terminal, root_guard, middle_guard):
    if middle >= mid_cap and middle_guard >= mid_cap:
      lower = upper = MID_MAX_INDEX
    else:
      lo,hi = (middle,middle_guard) if middle <= middle_guard else (middle_guard,middle)
      lower = lo//MID_STEP
      upper = MID_MAX_INDEX if hi >= mid_cap else hi//MID_STEP
    if first >= first_cap:
      first = first_cap
    else:
      first += first_shifts[first%MARGIN_STEP]
    if second >= second_cap:
      second = second_cap
    else:
      second += second_shifts[second%MARGIN_STEP]
    narrow = last_radius if root_guard >= full_guard else bisect_right(first_debits,root_guard)-1
    t = terminal//MARGIN_STEP_END
    terminal_pair = midpoint[t if t < terminal_limit else terminal_limit]
    if narrow < 0:
      if terminal_pair == 0:
        return 0
      best_root,_,_ = root_payments(first,second,lower,0,narrow)
      return best_root*terminal_pair >> bits
    s = 346*terminal//10**7
    sibling_only = 10**6*strips[s if s < 30000 else 30000]
    if terminal_pair == 0 and sibling_only == 0:
      return 0
    best_root,best_both,negative = root_payments(first,second,lower,upper,narrow)
    split_terminal = (best_both*sibling_only+negative*terminal_pair)//2
    return max(best_root*terminal_pair,split_terminal) >> bits

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
    left,right = left+shift,right+shift
    q = (left*left+left*right+right*right)/6
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
        if product:
          root_sum += joint_row[k]*isqrt(isqrt(isqrt(product)))
      count += max(0,end-start)
    total += (xhi-xlo)*xweight*root_sum
    lower = upper
  return SIGNED_UNIT*(total << bits)/2, count


def coupled_r(i, c_lo, c_hi, alpha, k2, k3):
  J,C3 = k2/2,k3/3
  return (
      max(0,2*J-2*c_lo)+max(0,2*J-alpha+c_hi)+max(0,c_hi-alpha),
      F(89,5)+2*c_hi,
      2-2*alpha**2+2*max(0,c_hi-alpha**2),
      4*J-4*alpha**3+2*max(0,c_hi-alpha**3),
      2*J+2*c_hi+2*max(0,J-c_lo),
      4*C3-2+2*c_hi+2*max(0,C3-c_lo),
  )[i]


def signed_cutoff_controls(eta, rows, density, normalizer):
  lo2, hi2 = F(1414213562373, 10**12), F(1414213562374, 10**12)
  lo3, hi3 = F(1732050807568, 10**12), F(1732050807569, 10**12)
  assert lo2**2 < 2 < hi2**2 and lo3**2 < 3 < hi3**2
  k2, k3 = F(109, 25), F(561, 50)
  paths = (F(0), F(1), k2, k3)
  lower = (F(0), 2*lo2/3, F(5, 6), lo2/2)
  upper = (F(0), 2*hi2/3, F(5, 6), hi2/2)
  r_abs, q_abs = 2*eta, hi3*F(367, 100)*eta
  alpha = lo2/2
  error_y = (2*k2-alpha,F(89,5),2-2*alpha**2,2*k2-4*alpha**3,2*k2,2*k3-2)
  from one_band_outer import paid_costs
  error_z = paid_costs()
  ideal_y = (-F(1, 6), -lo2/6, hi2/3, F(1, 2), -F(1, 3), -lo2/3)
  ideal_z = (F(0), F(0), -1/(2*hi3), -lo2/(2*hi3), F(0), -1/(2*hi3))
  for i, d in enumerate((1, 2, 2, 3, 1, 2)):
    covariance = upper[d]+paths[d]*eta
    assert covariance < 1
    bounds = (covariance if i < 4 else -lower[d],
              ideal_y[i]+min(error_y[i]+2*covariance,coupled_r(i,lower[d],covariance,alpha,k2,k3))*eta,
              ideal_z[i]+hi3*error_z[i]*eta+paths[d]*eta*q_abs,
              (1-lower[d])*F(11, 200))
    assert all(F(tick, 10**5) > bound for tick, bound in zip(rows[i], bounds, strict=True))

  off = hi3*F(1105, 100)*eta+r_abs*q_abs
  r_lo, q_lo = 1-(8*k2-2*lo2)*eta-r_abs**2, 1-6*k2*eta-q_abs**2
  assert min(r_lo, q_lo)-off > F(47,50)
  assert all(F(tick, 10**9) > bound for tick, bound in
             zip(density, (1+2*(k2-2*lo2)*eta, F(1), off), strict=True))
  pi_upper = F(355, 113)
  atan5 = sum((F((-1)**j, (2*j+1)*5**(2*j+1)) for j in range(9)), F(0))
  assert 16*atan5-4*(F(1, 239)-F(1, 3*239**3)) < pi_upper
  assert normalizer == SIGNED_NORMALIZERS[-1]
  for factor, norm in zip((3,F(5,2),2),SIGNED_NORMALIZERS,strict=True):
    grade_r = 1-((2*k2-lo2)*factor+(2*k2+lo2))*eta-r_abs**2
    assert grade_r > 0 and norm**2*(2*pi_upper)**3 < grade_r*q_lo-off**2


def certificate(progress=None):
  signed_cutoff_controls(ETA,ROWS,DENSITY,SIGNED_NORMALIZER)
  prepared = prepare()
  strips = prepared[0]
  extent = F('277/10000')
  h,b,c = F(3),F(5,2),F(2)
  variants = [tuple(folded_private_bins(ETA,cap,sd,shifted=shifted)
                    for cap,sd in ((h,24224),(b,22113),(c,19779))) for shifted in (False,True)]

  def ratio(weaker,stronger):
    return min(F(a,b) for a,b in zip(weaker,stronger,strict=True) if b)

  def half_ratio(rows):
    return ratio([mass for _,mass in rows[0]],[mass for _,mass in rows[2]])

  caps = (c,b,h)
  normalizers = tuple(reversed(SIGNED_NORMALIZERS))
  private = (F(1),F(47,50),min(map(half_ratio,variants)))

  @lru_cache(maxsize=None)
  def core(o,u):
    return tuple(map(tuple,adaptive_core_table(ETA,strips,(caps[o],caps[u]))))

  @lru_cache(maxsize=None)
  def terminal(o,u,w):
    return midpoint_strip_table(ETA,strips,(caps[o],caps[u],caps[w]))

  reference_core = tuple(value for row in core(0,1) for value in row)
  core_ratios = {(o,u):ratio(tuple(value for row in core(o,u) for value in row),reference_core)
                 for o,u in product(range(3),range(1,3))}
  terminal_ratios = {(o,u,w):ratio(terminal(o,u,w),terminal(0,1,1))
                     for o,u,w in product(range(3),range(1,3),range(1,3))}
  scales = {}
  for o,a,u,w in product(range(3),range(3),range(1,3),range(1,3)):
    scale = core_ratios[o,u]*terminal_ratios[o,u,w]*private[a]*normalizers[o]/normalizers[0]
    assert 0 < scale <= 1
    scales[o,a,u,w] = scale

  paired = []
  for shifted,rows in enumerate(variants):
    assert all([h for h,_ in row] == [h for h,_ in rows[2]] for row in rows)
    paired.append(tuple(((707107*h+999999)//1000000 if shifted else h,tight,weak)
                        for (h,tight),(_,weak) in zip(rows[2],rows[1],strict=True)))

  private_choice = private_policy(*paired)

  caps = tuple((value+MARGIN_STEP-1)//MARGIN_STEP*MARGIN_STEP for value in
               (1000*500000000//1415+max(row[0] for row in paired[1]),
                500000000+max(row[0] for row in paired[0])))

  jobs = ((-1,(8,4),extent,813458,323683082),(1,(30,16),extent,15113,350783352))
  records = []
  for side,widths,threshold,count_expected,ticks in jobs:
    total,count = policy_sum(side,widths,threshold,prepared,private_choice,terminal(0,1,1),core(0,1),caps)
    floor = F(ticks,10**14)
    assert count == count_expected and total > floor
    row = {"star_factor":"root-private", "side":side, "threshold_abs_max":str(threshold),
           "cells":count, "cell_widths":list(widths), "exact_lower_sum":str(total),
           "usable_floor":str(floor)}
    records.append(row)
    if progress is not None:
      progress(row)
  return {"eta":str(ETA), "records":records, "scales":scales}


signed = SimpleNamespace(certificate=certificate)


# Fixed nine-key rational parameter factory; no driver-cell search.
from fractions import Fraction as Q
from math import factorial, prod, isqrt
ORDER_UNIT = 10 ** 9
RADIUS_ROWS = ((18000000,706639,584309),(25040000,856187,742470),(33650000,950508,871882),(28000000,897852,794520))
CALIBRATION = {'normal_prefactor_lower': '19947/25000', 'core_sd': '17130/1000000', 'terminal_sd': '22115/1000000'}

def root_ending_probability(root_cap, closing):
  c = min(300000000, closing)
  if c <= 0:
    return 0
  u = min(root_cap, 2 * c // 3)
  if u <= 0:
    return 0
  v = c - (3 * u + 1) // 2
  assert 0 <= v <= c
  moment = u * u * (c + 3 * v) + (c + v) * (c * c + v * v)
  denominator = 216 * 10 ** 18
  lower_area = u * (c + v) // 2 - (100 * u * moment + denominator - 1) // denominator
  return max(0, lower_area * 700 // 99)

def parameters(relaxed_endings=False):
  eta, threshold = (Q('57/100000'), Q('277/10000'))
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
 assert p['theta']==F('57/100000') and p['threshold']==F('277/10000')
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
 return {'normal_prefactor_calibrated':True,'all_eight_defect_masses_calibrated':True,'pair_covariance_calibrated':True,'reconstruction_and_projection_guards_calibrated':True}

class Context:

  def __init__(self, relaxed):
    ns = vars(nine_parameters)
    p = ns['parameters'](relaxed)
    self.strip, self.unit = p['strip'], p['unit']
    self.theta, self.threshold = p['theta'], p['threshold']
    self.radii = ns['RADIUS_ROWS']
    self.pair = ns['root_ending_probability']
    self.guards = calibrate(p, ns['CALIBRATION'], self.radii, relaxed)
    self.weights = [[self.radii[r][1] ** 4 * self.radii[e][2 if relaxed else 1] ** 2 for e in range(4)] for r in range(4)]
    self.mid = [ceil(972 * r[0], 1000) for r in self.radii]
    self.root_charge = [ceil(r[0],3) for r in self.radii]
    self.close = [[ceil(1208 * r[0] + 500 * e[0], 1000) for e in self.radii] for r in self.radii]

  def mass(self, d):
        # The feasibility cap ensures each remaining margin is nonnegative.
    index = 173 * d // 5000000
    return self.strip[index if index < 30000 else 30000]

  def ending_mass(self, endings, closing, core_cap, root_kind, kind):
    left, right = endings
    L, R = left - closing, right - closing
    cap = min(core_cap, min(L, R) * 2 // 3)
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
      rad = 1 + isqrt(disc - 1)
      root = min(cap, max(0, 2 * (L + R - rad) // 9))
    debit = (3 * root + 1) // 2
    if kind == 0:
      return self.mass(root) * self.mass(L - debit) * self.mass(R - debit)
    if kind == 1:
      return self.pair(root, L) * self.mass(R - debit)
    return self.pair(root, R) * self.mass(L - debit)

  def candidates(self, g, strict_end, relaxed_end, recipe):
    first, branch, leaf, selected = g
    root_charge, mid, closing, root_kind, kind, strict_weight, relaxed_weight = recipe
    if first <= root_charge or leaf <= mid:
      return 0, 0
    core_cap = min(branch - root_charge, (selected - mid) * 1000 // 1415)
    if core_cap <= 0:
      return 0, 0
    common = self.mass(leaf - mid)
    strict = self.ending_mass(strict_end, closing, core_cap, root_kind, kind)
    relaxed = strict if strict_end == relaxed_end else self.ending_mass(
        relaxed_end, closing, core_cap, root_kind, kind)
    return common * strict * strict_weight, common * relaxed * relaxed_weight


nine_kernel = SimpleNamespace(COUNT=COUNT, Context=Context)


NORMALIZER = F(123513,5000000)

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
  rr = 1-(8*PATH[2]-2*lo2)*theta-gamma[0]**2
  qq = 1-6*PATH[2]*theta-gamma[1]**2
  cross = hi3*F(491,20)*theta+gamma[0]*gamma[1]
  far = F(353,10)*theta+gamma[1]**2
  assert min(rr-2*cross, qq-cross-far) > F(93,100)
  determinant = rr*qq**2-rr*far**2-2*qq*cross**2-2*cross**2*far
  atan5 = sum((F((-1)**j, (2*j+1)*5**(2*j+1)) for j in range(9)), F(0))
  pi_upper = F(355,113)
  assert 16*atan5-4*(F(1,239)-F(1,3*239**3)) < pi_upper
  assert NORMALIZER**2*(2*pi_upper)**4 < determinant
  density = tuple((value*S).__ceil__() for value in (1+2*(PATH[2]-2*lo2)*theta,cross,far))
  return gamma, density, determinant


def tree_paths():
  graph = [[] for _ in range(14)]
  for i, (a,b) in enumerate(EDGES):
    graph[a].append((b,i))
    graph[b].append((a,i))
  paths = []
  for root in range(14):
    entries,pending = {root:()},[root]
    for a in pending:
      for b,edge in graph[a]:
        if b not in entries:
          entries[b] = (*entries[a],edge)
          pending.append(b)
    assert len(entries) == 14
    paths.append(entries)
  return paths


from one_band_outer import revised_q_coefficient


def mean_rows(theta, extent, gamma, relaxed, side, start):
  paths = tree_paths()
  lo2, hi2 = SQRT2
  lo3, hi3 = SQRT3
  lo = tuple((1+F(d,3))/hi2**d for d in range(6))
  hi = tuple((1+F(d,3))/lo2**d for d in range(6))

  def signed_cov(v, w, sign):
    edges = paths[v][w]
    length,terminals = len(edges),sum(e>=9 for e in edges)
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
    depth = len(paths[v][0])
    root_cov = signed_cov(v,0,1)
    assert root_cov < 1
    values = [sign*signed_cov(v,0,sign)]
    for j, driver in enumerate(DRIVERS):
      upper = F(0)
      charges = [F(0)]*len(EDGES)
      for u,m in driver:
        signed = sign*m
        edges = paths[v][u]
        length = len(edges)
        upper += signed*(hi[length] if signed>0 else lo[length])
        if length:
          factor = PATH[length]/length if signed>0 else (lo2/2)**(length-1)
          for edge in edges:
            charges[edge] += signed*factor
      critical = upper
      upper += theta*sum(max(F(0),c)*(1+int(relaxed and i>=9))
                         for i,c in enumerate(charges))
      if j:
                # The two endpoints are siblings in the other branch.
        if all(len(paths[v][u]) == depth+2 for u,_ in driver):
          costs = (F(785,100), F(74,5), F(289,10) if relaxed else F(26))
          upper = min(upper, costs[depth-1]*theta)
        upper *= hi3 if upper >= 0 else lo3
      bound = upper+root_cov*gamma[j]
      if j:
        revised = revised_q_coefficient(v,j,depth,relaxed,critical,theta,PATH,lo3,hi3,REFLECT)
        if revised is not None:
          bound = min(bound,revised)
      values.append(bound)
    values.append(extent*(1-lo[depth]) if sign*side<0 else -start*(1-root_cov))
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


def margin_tables(rows, relaxed_rows, root, breadth, ending):
  def split(rows):
    xy = tuple((a,b,bias) for a,b,_,_,bias in rows)
    zh = tuple((c,d) for _,_,c,d,_ in rows)
    return xy, zh

  strict_xy, strict_zh = split(rows)
  relaxed_xy, relaxed_zh = split(relaxed_rows)
  assert strict_xy[:9] == relaxed_xy[:9] and strict_zh[:9] == relaxed_zh[:9]
  for xy, zh in ((strict_xy,strict_zh), (relaxed_xy,relaxed_zh)):
    assert xy[1] == xy[2] and xy[3] == xy[4] and zh[3] == zh[4]
    assert xy[5] == xy[7] and xy[6] == xy[8]
    assert zh[5] == zh[6] and zh[7] == zh[8]
    assert xy[9] == xy[10] == xy[11] == xy[12]
    assert zh[9] == zh[10] and zh[11] == zh[12]

  def root_breadth(x,y):
    def term(row):
      a,b,bias = row
      return cap(a,x)+cap(b,y)+bias*S

        # The four selected-core rows are a Cartesian product of
        # two root/breadth terms and two ending terms.
    return (term(strict_xy[0]), term(strict_xy[1]), term(strict_xy[3]),
            max(term(strict_xy[5]),term(strict_xy[6])),
            term(strict_xy[9]), term(relaxed_xy[9]))

  def endpoints(z,h):
    def term(row):
      c,d = row
      return cap(c,ending[z])+cap(d,ending[h])

    return (term(strict_zh[0]), max(term(strict_zh[1]),term(strict_zh[2])),
            term(strict_zh[3]), max(term(strict_zh[5]),term(strict_zh[7])),
            term(strict_zh[9]), term(strict_zh[11]),
            term(relaxed_zh[9]), term(relaxed_zh[11]))

  xy = tuple(tuple(root_breadth(x,y) for y in breadth) for x in root)
  zh = tuple(endpoints(z,h) for z,h in ADDRESSES)
  return xy,zh


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


THRESHOLD_BANDS = ((F(0),F(277,10000),-1),(F(0),F(277,10000),1))
NINE_TICKS = ((1859050,1424726),(2055903,1581244))


def precision_nine_replay(strict, relaxed, codes, bands, side_rows, coefficients):
  assert GROUPS == ((1,),(2,3),(4,5),(6,7,8,9),(10,11),(12,13))
  root,breadth,ending = grid()
  tables = []
  for (start,extent,side), rows in zip(bands,side_rows,strict=True):
    xy,zh = margin_tables(*rows,root,breadth,ending)
    shifts = tuple(2*side*S*t for t in (start,extent))
    assert all(t.denominator == 1 for t in shifts)
    lower,upper = (t.numerator for t in shifts)
    root_density = tuple(lower_density(S*(3*max((a+b+lower)**2,(a+b+upper)**2)+(b-a)**2))
                         for a,b in root)
    tables.append((xy,zh,root_density))
  _,dense_y,dense_zh = density_tables(F(0),coefficients,root,breadth,ending)
  drivers = tuple(tuple(lower_density(y_part+linear*z_sum+z_part)
                        for z_part,z_sum,_ in dense_zh) for y_part,linear in dense_y)
  sums = [[0,0] for _ in tables]
  counts = [[0,0] for _ in tables]
  assert strict.mid == relaxed.mid and strict.root_charge == relaxed.root_charge
  assert strict.close == relaxed.close and strict.strip is relaxed.strip
  assert strict.pair is relaxed.pair
    # Floor division by -S is the original negative ceiling.
  negative_scale = -S
  recipes = []
  for code in range(1,193):
    choice,kind = divmod(code-1,3)
    choice,root_kind = divmod(choice,4)
    r,e = divmod(choice,4)
    recipes.append((strict.root_charge[r],strict.mid[r],strict.close[r][e],root_kind,kind,
                    strict.weights[r][e],relaxed.weights[r][e]))
  for i,code in enumerate(codes):
    if code == 0:
      continue
    x_y,z_h = divmod(i,820)
    x,y = divmod(x_y,60)
    assert 1 <= code <= 192
    recipe = recipes[code-1]
    driver = drivers[y][z_h]
    multiple = dense_zh[z_h][2]*(1 if x<20 else 2)
    for slot,(xy,zh,weight) in enumerate(tables):
      a,b,c,d,ends,ends_relaxed = xy[x][y]
      u,v,w,k,left,right,left_relaxed,right_relaxed = zh[z_h]
      core = ((a+u)//negative_scale,(b+v)//negative_scale,
              (c+w)//negative_scale,(d+k)//negative_scale)
      strict_end = ((ends+left)//negative_scale,(ends+right)//negative_scale)
      relaxed_end = ((ends_relaxed+left_relaxed)//negative_scale,
                     (ends_relaxed+right_relaxed)//negative_scale)
      sm,rm = strict.candidates(core,strict_end,relaxed_end,recipe)
      density = (weight[x]*driver//S)*multiple
      sums[slot][0] += density*sm
      sums[slot][1] += density*rm
      counts[slot][0] += sm>0
      counts[slot][1] += rm>0
  return tuple(tuple(zip(s,c,strict=True)) for s,c in zip(sums,counts,strict=True))


def _nine():
  from fractions import Fraction as F
  import lzma
  COUNT = nine_kernel.COUNT
  Context = nine_kernel.Context
  LIMIT=65536

  def restore_order(codes):
    result = bytearray(COUNT)
    for triangle in range(820):
      for y in range(60):
        src = 30*(60*triangle+y)
        start = 820*y+triangle
        result[start::49200] = codes[src:src+30]
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
    return restore_order(codes)

  def certificate(encoded, progress=None):
    codes = decode(encoded)
    strict, relaxed = Context(False), Context(True)
    assert strict.theta == relaxed.theta == price.ETA/2
    assert strict.threshold == relaxed.threshold == F(277,10000)
    gamma,density,determinant = covariance_bounds(strict.theta)
    rows = tuple(tuple(mean_rows(strict.theta,hi,gamma,flag,side,lo) for flag in (False,True))
                 for lo,hi,side in THRESHOLD_BANDS)
    assert all(pair[0][:9] == pair[1][:9] for pair in rows)
    sums = precision_nine_replay(strict,relaxed,codes,THRESHOLD_BANDS,rows,density)
    assert strict.unit == relaxed.unit == F(63,2500*400000*S**5*10**36)
    unit = strict.unit*NORMALIZER/F(63,2500)
    results = []
    for i,phase in enumerate(("strict","relaxed")):
      band_rows = []
      for (lo,hi,side),pair,ticks in zip(THRESHOLD_BANDS,sums,NINE_TICKS,strict=True):
        mass,count = pair[i]
        exact,floor = mass*unit,F(ticks[i],10**12)
        assert exact > floor
        band_rows.append({"side":side,"band_abs_endpoints":[str(lo),str(hi)],
                          "exact_lower_sum":str(exact),"used_floor":str(floor),"positive_codes":count})
      row = {"phase":phase,"eta":str(price.ETA),"bands":band_rows,
             "calibration":(strict,relaxed)[i].guards,"precision_determinant_lower":str(determinant)}
      results.append(row)
      if progress is not None:
        progress(row)
    return results

  return SimpleNamespace(certificate=certificate)


nine = _nine()


EXPECTED_CREDIT = F('1804115787/50000000000000')
OUTER_CREDIT = F('3643068418759/100000000000000000')
GRAPH_COEFFICIENT = F('50000000/382814729')
BAND_WEIGHTS = ((323682897, 130665237, 54460609, 93, 185, 1976693348, 820661048), (801175715, 496676645, 496676460), (43176464, 19420985), (74917474, 59247915, 84281985, 9302527))


def mixed_allocation(adjacent, baseline2, arms2, long_paths, size_four, d7, base, signed, nine, outer):
  A,k,bilateral,star_cover,k0,groups,budget = price.allocation_inputs(
      adjacent, baseline2, arms2, long_paths, size_four, d7)
  assert F(signed["eta"]) == price.ETA
  assert price.WEAK_WEIGHT == 0
  assert len(nine) == 2 and {row["phase"] for row in nine} == {"strict","relaxed"}
  assert all(F(row["eta"]) == price.ETA and len(row["bands"]) == len(THRESHOLD_BANDS) for row in nine)
  global_share = k0*price.ETA/A
  assert 0 < global_share < 1 and price.ONE_LEAF_BONUS < global_share*(price.SIGMA-F(3,2))/2
  extent = F(277,10000)
  endpoints = {row["side"]:F(row["usable_floor"]) for row in signed["records"]
               if F(row["threshold_abs_max"]) == extent}
  assert set(endpoints) == {-1,1} and len(signed["records"]) == 2
  common = min(endpoints.values())
  expected = set(product(range(3),range(3),range(1,3),range(1,3)))
  assert common > 0 and signed["scales"].keys() == expected
  bands = []
  previous = {-1:F(0),1:F(0)}
  for slot,((lo,hi,side),(weights,source,medium_ticks,high_ticks)) in enumerate(zip(THRESHOLD_BANDS,(BAND_WEIGHTS,)*2,strict=True)):
    assert previous[side] == lo < hi <= extent
    previous[side] = hi
    ends = [str(lo),str(hi)]
    caps = {}
    for row in nine:
      bound = row["bands"][slot]
      assert bound["side"] == side and bound["band_abs_endpoints"] == ends
      caps[row["phase"]] = F(bound["used_floor"])
    strict_slack = []
    q = F(weights[0],10**14)
    medium_charges = tuple(F(t,10**14) for t in (*medium_ticks,0,0))
    high_charges = tuple(F(t,10**14) for t in high_ticks)
    for grades,scale in signed["scales"].items():
      paid = sum(d if g == 1 else h if g == 2 else F(0)
                 for g,d,h in zip(grades,medium_charges,high_charges,strict=True))
      strict_slack.append(common*scale+paid-q)
    assert q > 0 and min(strict_slack) > 0
    debit03,debit02 = 6*sum(high_charges),6*sum(medium_charges)
    marking = price.marked_rows(A,weights,*source,k,caps,medium=(100,100),
                                deficit=F(0),class_debit=debit03,class_debit02=debit02)
    bands.append({"side":side,"band_abs_endpoints":ends,"minimum_signed_grade_slack":str(min(strict_slack)),
                  "class_debits":[str(debit03),str(debit02)],"marking":marking,
                  "combined_credit":marking["credit"]})
  assert previous == {-1:extent,1:extent}
  assert F(outer["threshold_abs_max_of_central_band"]) == extent
  assert F(base["path"]["gamma_upper"]) < PATH_UPPER
  central = min(F(b["combined_credit"]) for b in bands)
  assert central == EXPECTED_CREDIT
  saving = min(central,F(outer["outer_graph_credit"]))
  assert F(outer["outer_graph_credit"]) == OUTER_CREDIT
  assert saving == central < OUTER_CREDIT
  assert PATH_UPPER - saving < GRAPH_COEFFICIENT
  return {"scope":"two signed threshold bands with a common ordering",
          "bilateral_medium":bilateral,"paired_star_cover":star_cover,
          "global_K0":str(k0),"raw_price_groups":{name:str(cap) for name,cap,*_ in groups},
          "price_guard_slacks":budget["coupling_guards"],"price_allowance":str(A),
          "bands":bands,"outer_credit":outer["outer_graph_credit"],
          "saving":str(saving),"circuit_coefficient":str(1+1/(2*GRAPH_COEFFICIENT))}


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

  from one_band_outer import sharp_scalar_dyad
  assert sharp_scalar_dyad() == (756, 12)
  progress("actual source intervals", source_intervals.certificate())
  base_result = base.base_components(progress)
  signed_result = signed.certificate(lambda row: progress("signed", row))
  choices = (here/"nine-choices.txt").read_bytes()
  assert sha256(choices).hexdigest() == "f12eb31c8bd49c76ba92749a83ad1601f8d3c8bb2f564ea57524e4924848a56a"
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
