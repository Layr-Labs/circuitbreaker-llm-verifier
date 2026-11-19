"""Verify all numerical premises of the one-band circuit lower bound."""
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import sys

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
import one_band_price as price
from one_band_outer import d7, outer


"""Standalone directed quadrature for the 47-term crossing-path bound."""
from decimal import Context as DecimalContext, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction
from math import cos, factorial, pi
from typing import NamedTuple


Q_LOW = DecimalContext(prec=100, rounding=ROUND_FLOOR)
Q_HIGH = DecimalContext(prec=100, rounding=ROUND_CEILING)
Q_MID = DecimalContext(prec=100, rounding=ROUND_HALF_EVEN)
GRID = 1 << 224
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
    radius = Decimal('1e-72')
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
            active = [(i,(state[i][middle]*A[i][middle]) >> 224)
                      for i in range(q) if state[i][middle]]
            active = [(i,a) for i,a in active if a]
            row = []
            for last, kernel in enumerate(B):
                contraction = sum(a*kernel[i] for i,a in active)
                row.append((C[middle][last]*contraction) >> 448)
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

def inverse(A):
    n=len(A)
    B=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        p=B[i][i]
        assert p>0
        B[i]=[x/p for x in B[i]]
        for j in range(n):
            if i!=j:
                p=B[j][i]
                B[j]=[x-p*y for x,y in zip(B[j],B[i])]
    return [row[n:] for row in B]

def path_distances():
    distances = [[0 if i == j else 100 for j in range(12)] for i in range(12)]
    for i, j in EDGES:
        distances[i][j] = distances[j][i] = 1
    for k in range(12):
        for i in range(12):
            for j in range(12):
                distances[i][j] = min(distances[i][j], distances[i][k]+distances[k][j])
    return distances

def signed_magnitude_bounds(terms):
    positive = negative = F(0)
    for term in terms:
        if term >= 0:
            positive += term
        else:
            negative -= term
    return positive, negative

def critical_path_multiplier(d, sqrt_two_upper):
    """Upper bound the exact weighted critical coefficient on d edges."""
    short = d // 2
    long = d - short
    factor = 4 - F(2)**(1-short) - F(2)**(1-long)
    power = 2**((d-1)//2) * (sqrt_two_upper if d % 2 == 0 else 1)
    return d * factor * power

THRESHOLD = F(9, 400)
THETA, FIRST_Z = F(19, 80000), F(11, 50)
SQRT2_LO, SQRT2_HI = F(141421, 100000), F(70711, 50000)
SQRT3_LO, SQRT3_HI = F(34641, 20000), F(86603, 50000)
LAMBDA_LO, LAMBDA_HI = 2*SQRT2_LO, 2*SQRT2_HI
DEFECT_STRIP, MAX_END = F(3, 40), F(9, 10)
X_STEP, G_STEP, FIRST_G = F(3, 100), F(3, 20), F(1, 20)


def dot(row, vector):
    return sum(x*y for x, y in zip(row, vector, strict=True))


def matrix_vector(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def regression_controls(theta):
    assert SQRT2_LO**2 < 2 < SQRT2_HI**2 and SQRT3_LO**2 < 3 < SQRT3_HI**2
    distances = path_distances()
    error = (F(0), *(theta*critical_path_multiplier(d, SQRT2_HI) for d in range(1, 7)))
    assert max(map(max, distances)) == 6
    assert all((1+F(d, 3))/SQRT2_LO**d+error[d] < 1 for d in range(1, 7))
    drivers = [{0: F(1)}, {0: 2*LAMBDA_HI, 1: F(-3), 2: F(-3)}]
    drivers.extend({i: SQRT3_HI, j: -SQRT3_HI}
                   for i, j in ((1, 2), (4, 5), (6, 7), (8, 9), (10, 11)))
    P, variances = [], []
    for i, row in enumerate(drivers):
        majorant = []
        for j, col in enumerate(drivers):
            positive, negative = signed_magnitude_bounds(
                a*b*error[distances[u][v]] for u, a in row.items() for v, b in col.items())
            majorant.append(max(positive, negative))
            if i == j:
                variances.append(1+positive)
        P.append(majorant)
    row_max = max(map(sum, P))
    assert row_max < 1
    W = inverse([[F(i == j)-P[i][j] for j in range(7)] for i in range(7)])
    assert all(a >= 0 for row in W for a in row)
    assert prod(variances)/636**2 < F(7, 44)**7
    defects = [tuple(max(signed_magnitude_bounds(
        a*3*theta*SQRT2_HI**distances[v][u] for u, a in row.items())) for row in drivers)
        for v in (0, 1, 2, 4, 6)]
    mean_rows = [tuple(dot(row, col) for col in zip(*W, strict=True)) for row in defects]
    return W, mean_rows, row_max


def slabs(theta):
    variance = 3*(3-LAMBDA_LO)*theta
    deviation = max(F(553, 50000), F(isqrt((variance*10**10).__floor__())+1, 10**5))
    assert variance < deviation**2 and F(797, 1000)**2 < F(7, 11)
    for index in range(44, 160):
        lo, hi = F(index, 200), F(index+1, 200)
        arg = min(2, lo*DEFECT_STRIP/deviation)
        strip = F(797, 1000)*sum(
            F((-1)**k, 2**k*factorial(k)*(2*k+1))*arg**(2*k+1) for k in range(12))
        assert 0 < strip < 1
        strip_floor = F((strip*10**9).__floor__(), 10**9)
        yield hi, (hi**7-lo**7)*strip_floor**5/7


def weak_controls(theta):
    assert theta == F(323, 800000)
    W, mean_rows, row_max = regression_controls(theta)
    integration = tuple(slabs(theta))

    denominator = lcm(*(moment.denominator for _, moment in integration))
    weights = [(hi**2, int(moment*denominator)) for hi, moment in integration]

    # A slab endpoint is i/200 and rate is in hundredths.
    squares = []
    for square, weight in weights:
        index_square = square*200**2
        assert index_square.denominator == 1
        squares.append((index_square.numerator, weight))
    exp_denominator = 64*100*200**2
    coefficients = tuple(c*exp_denominator**j for j,c in
                         enumerate((7,-42,210,-840,2520,-5040,5040),1))
    divisor = 5040*exp_denominator**7

    def exponential(n):
        assert 0 <= n <= exp_denominator
        value = -1
        for c in coefficients:
            value = n*value+c
        value = max(0,value*UNIT//divisor)
        for _ in range(6):
            value = value*value//UNIT
        return value

    @lru_cache(maxsize=None)
    def mass_integral(rate):
        numerator = sum(weight*exponential(rate*index_square) for index_square,weight in squares)
        return F(numerator,denominator*UNIT)

    w_denominator = lcm(*(entry.denominator for row in W for entry in row))
    w_integer = [[int(entry*w_denominator) for entry in row] for row in W]
    amplitude_denominator = 55000000
    rate_denominator = 2*amplitude_denominator**2*w_denominator

    shift = (3-LAMBDA_LO)*THRESHOLD/FIRST_Z
    assert 2*LAMBDA_LO-36*X_STEP-(6-2*LAMBDA_LO)*THRESHOLD/FIRST_Z > 0
    assert F(3, 4)*LAMBDA_HI < F(5, 2)
    total, kept, worst_rate = F(0), 0, F(0)
    for left, right, ga, gb in product(range(6), repeat=4):
        x_lo, x_hi = left*X_STEP, (left+1)*X_STEP
        y_lo, y_hi = right*X_STEP, (right+1)*X_STEP
        a_lo, a_hi = FIRST_G+ga*G_STEP, FIRST_G+(ga+1)*G_STEP
        b_lo, b_hi = FIRST_G+gb*G_STEP, FIRST_G+(gb+1)*G_STEP
        amplitudes = (1+THRESHOLD/FIRST_Z,
                      2*LAMBDA_HI-3*(x_lo+y_lo)+(6-2*LAMBDA_LO)*THRESHOLD/FIRST_Z,
                      SQRT3_HI*max(x_hi-y_lo, y_hi-x_lo), a_hi, b_hi, MAX_END, MAX_END)
        means = matrix_vector(mean_rows, amplitudes)
        if LAMBDA_LO-x_hi-y_hi-shift-means[0]-DEFECT_STRIP <= 0:
            continue
        radii = []
        for x, g_lo, g_hi, parent, child in ((x_hi, a_lo, a_hi, 1, 3), (y_hi, b_lo, b_hi, 2, 4)):
            midpoint = (LAMBDA_HI*x-1+shift+means[parent]+DEFECT_STRIP)/2
            guard = 1-F(3, 4)*(1+means[parent]+DEFECT_STRIP)+(F(3, 4)*LAMBDA_LO-F(5, 2))*x
            if midpoint+g_hi/(2*SQRT3_LO) >= 0 or guard <= F(3, 4)*shift:
                break
            selected = midpoint-g_lo/(2*SQRT3_HI)
            cap = (-LAMBDA_LO*selected-shift-means[child]-DEFECT_STRIP)*SQRT3_LO
            radius = min(MAX_END, F((cap*1000).__ceil__()-1, 1000))
            if radius <= 0:
                break
            radii.append(radius)
        if len(radii) != 2:
            continue
        actual_amplitudes = (*amplitudes[:5], *radii)
        scaled = [value*amplitude_denominator for value in actual_amplitudes]
        assert all(value.denominator == 1 for value in scaled)
        amplitudes_integer = [value.numerator for value in scaled]
        rate_numerator = 100*dot(amplitudes_integer, matrix_vector(w_integer, amplitudes_integer))
        rate = -(-rate_numerator//rate_denominator)
        volume = 24*SQRT3_LO*X_STEP**2*G_STEP**2*prod(radii)
        total += volume*mass_integral(rate)/636
        kept += 1
        worst_rate = max(rate, worst_rate)

    usable = F((total*10**12).__ceil__()-1, 10**12)
    assert 0 < usable < total and len(integration) == 116
    return {"exact_lower_sum": str(total), "usable_floor": str(usable),
            "boxes": kept, "integration_slabs": len(integration),
            "maximum_driver_row_sum": str(row_max)}


def base_components(progress=None):
    path = path_series_enclosure()
    assert F(path["gamma_upper"]) < F('1633094581/12500000000')
    if progress is not None:
        progress("path", path)
    weak = weak_controls(F('323/800000'))
    assert F(weak["usable_floor"]) >= F(12719, 10**12)
    if progress is not None:
        progress("weak", weak)
    return {"path": path, "weak": weak}


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
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, lcm
import json

S = 1000000000
ROWS = ((94363, -15810, 1583, 315), (83686, -21997, 2500, 917), (83686, 47438, -27827, 917), (71618, 50820, -38984, 1611), (-94280, -32476, 625, 315), (-83333, -45193, -27827, 917))
ETA = F('323/400000')
T = F('9/400')
CORE = (28000000, 31466000, 9333334, 769390973799)
PRIVATE = ((9000000, 622691), (10500000, 696994), (11750000, 750937), (12520000, 780613), (14000000, 830353), (15375000, 868501), (16825000, 901150), (18250000, 926583))
ENDING = ((48842000, 622691), (50342000, 696994), (51592000, 750937), (52362000, 780613), (53842000, 830353), (55217000, 868501), (56667000, 901150), (58092000, 926583))
DENSITY = (1007041401, 1000000001, 15463245)
WIDTHS = {-1: (3, 3), 1: (9, 4)}
SOURCE_SHA = 'b3e1b4d8df84ba9d80169553b41a2f1d706b7cb74ef1e79af336785ee486867a'
FLOOR = F('139358000/100000000000000')
COUNTS = {-1: 1239545, 1: 167778}
SIGNED_UNIT = F('1547/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000')
ROOT_EDGES = {-1: (0, 2, 4, 6, 8, 10, 12, 14, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 77, 82, 87, 93, 99, 105, 111, 117, 123, 129, 135, 141), 1: (0, 4, 8, 12, 16, 20, 24, 29, 34, 39, 44, 49, 54, 60, 66, 73, 81, 89, 97, 105)}
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
 for j in range(30001):
  p=0
  for c in reversed(coeff):p=p*j*j+c
  v=79788*S*j*p//divisor
  jj=min(j,16000);pol=1680*10000**6-280*jj*jj*10000**4+42*jj**4*10000**2-5*jj**6
  coarse=797*jj*pol*S//(10000*1680000*10000**6)
  best=max(best,v,coarse);assert 0<=best<S;result.append(best)
 return result

def diamond(first,second,third,half):
 c=min(300000000,second-half-(1415*third+999)//1000)
 if c<=0:return 0
 u=min(300000000,first-third,1000*c//1415)
 if u<=0:return 0
 v=c-(1415*u+999)//1000
 moment=u*u*(c+3*v)+c*c*c+c*c*v+c*v*v+v*v*v
 den=1182*10**18;pen=(625*u*moment+den-1)//den
 return max(0,(u*(c+v)//2-pen)*84//11)

def private_mass(first,second,code,strips):
 radius,cutoff,third,coremass=CORE
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
 for side,(dy,dz) in WIDTHS.items():
  joint=[]
  for j in range(0,1440,dy):
   row=[]
   for k in range(0,680,dz):
    N=a*(3*j*j+3*j*dy+dy*dy)+16*b*(3*k*k+3*k*dz+dz*dz)+6*c*(2*j+dy)*(2*k+dz)
    value=U-(N*U+den-1)//den
    for _ in range(12):value=value*value//U
    row.append(value*dy*dz)
   joint.append(row)
  joints[side]=joint
 return strips,joints


def policy_sum(side,prepared,private_choice,terminal_choice):
 assert side in(-1,1)
 strips,joints=prepared;joint=joints[side];dy,dz=WIDTHS[side];radius,cutoff,third,coremass=CORE
 biases=[((T/F(11,200))*h).__ceil__()if((i>=4)==(side>0))else 0 for i,(*_,h)in enumerate(ROWS)]
 total=0;pos=0;third_bias=(1415*third+999)//1000;area_den=1182*10**18
 reqs=(third+1,1,cutoff+1,ENDING[0][0]+1,radius//3+1,cutoff+1)
 edges=ROOT_EDGES[side]
 assert edges[0]==0 and all(a<b for a,b in zip(edges,edges[1:]))
 for xlo,xhi in zip(edges,edges[1:]):
  shift=side*T
  left,right=F(xlo,1000),F(xhi,1000)
  q=max((left*left+left*right+right*right)/3,((left+shift)**2+(left+shift)*(right+shift)+(right+shift)**2)/3)/2
  xweight=weight(q);root_sum=0
  for ji,j in enumerate(range(0,1440,dy)):
   aa=[];bb=[];low,high=0,680-dz
   for(a,b,c,_),bias,req in zip(ROWS,biases,reqs):
    at=-(10*a*(xhi if a>0 else xlo)+25*b*(j+dy*(b>0))+100*c*dz*(c>0)+10000*bias);sl=-100*c
    if sl>0:low=max(low,-((at-req)//sl))
    else:high=min(high,(at-req)//(-sl))
    aa.append(at);bb.append(sl)
   start=-(-low//dz)*dz
   first=aa[0]+bb[0]*start;second=aa[1]+bb[1]*start;terminal=aa[3]+bb[3]*start
   df,ds,dt=bb[0]*dz,bb[1]*dz,bb[3]*dz
   joint_row=joint[ji]
   for ki in range(start//dz,high//dz+1):
    privatecode=private_choice(first//MARGIN_STEP,second//MARGIN_STEP)
    endcode=terminal_choice(terminal//MARGIN_STEP);closing,endmass0=ENDING[endcode]
    if terminal>closing:
     index=346*(terminal-closing)//10**7
     if index>30000:index=30000
     endmass=endmass0*strips[index]
    else:endmass=0
    if privatecode%4:
     best=private_mass(first,second,privatecode,strips)
    else:
     half,mass=PRIVATE[privatecode//4]
     c=second-half-third_bias
     if c>300000000:c=300000000
     if c<=0:best=0
     else:
      u=first-third;limit=1000*c//1415
      if u>limit:u=limit
      if u<=0:best=0
      else:
       v=c-(1415*u+999)//1000
       moment=u*u*(c+3*v)+(c+v)*(c*c+v*v)
       area=u*(c+v)//2-(625*u*moment+area_den-1)//area_den
       best=mass*(area*84//11)if area>0 else 0
    pos+=1
    root_sum+=joint_row[ki]*best*endmass
    first+=df;second+=ds;terminal+=dt
  total+=(xhi-xlo)*xweight*root_sum
 assert pos==COUNTS[side]
 return SIGNED_UNIT*coremass*total,pos


def private_covariance_controls():
    e, k2, k3 = ETA, F(109,25), F(561,50)
    base = F(1,12)-k2*e/2
    diagonal1 = base-(e/2)**2-(25*e)**2/F(19,20)
    diagonal2 = base-(k2*e/2)**2-(59*e)**2/F(19,20)
    a, c = F(81,1000), F(79,1000)
    assert diagonal1 > a and diagonal2 > c
    off = (1+k3)*e/4+k2*e**2/4+25*59*e**2/F(19,20)
    assert min(a-off**2/c,c-off**2/a) > F(197,2500)
    assert F(84,11) == 24/F(22,7)
    assert F(625,1182) == 1/(24*F(197,2500))


def certificate(progress=None):
    private_covariance_controls()
    prepared = prepare()
    strips = prepared[0]

    @lru_cache(maxsize=None)
    def private_choice(a, b):
        first = a*MARGIN_STEP + MARGIN_STEP//2
        second = b*MARGIN_STEP + MARGIN_STEP//2
        return max(range(32), key=lambda code: private_mass(first, second, code, strips))

    @lru_cache(maxsize=None)
    def terminal_choice(z):
        margin = z*MARGIN_STEP + MARGIN_STEP//2
        return max(range(8), key=lambda i: ENDING[i][1]*strips[
            min(30000, max(0, 346*(margin-ENDING[i][0])//10**7))])

    records = []
    for side in (-1, 1):
        total, count = policy_sum(side, prepared, private_choice, terminal_choice)
        assert total > FLOOR
        row = {"side": side, "cells": count, "exact_lower_sum": str(total),
               "usable_floor": str(FLOOR)}
        records.append(row)
        if progress is not None:
            progress(row)

    return {"scope": "exact signed .0225 envelope with one retained core and computed choices",
            "eta": str(ETA), "threshold_abs_max": str(T), "source_sha256": SOURCE_SHA,
            "margin_bin_width": str(F(MARGIN_STEP, S)), "records": records,
            "private_choice_bins": private_choice.cache_info().currsize,
            "terminal_choice_bins": terminal_choice.cache_info().currsize}


signed = SimpleNamespace(certificate=certificate)


# Fixed nine-key rational parameter factory; no driver-cell search.
from fractions import Fraction as Q
from itertools import product
from math import factorial, prod, isqrt
ORDER_UNIT = 10 ** 9
RADIUS_ROWS = ((18000000, 788061, 666303), (25040000, 917515, 821302), (33650000, 980375, 929263), (28000000, 947826, 867351))
CALIBRATION = {'normal_prefactor_lower': '19947/25000', 'core_sd': '721/50000', 'terminal_sd': '931/50000'}

def rational_inverse(matrix):
    dimension = len(matrix)
    table = [list(row) + [Q(row_index == column) for column in range(dimension)] for row_index, row in enumerate(matrix)]
    for column in range(dimension):
        pivot = table[column][column]
        assert pivot > 0
        table[column] = [entry / pivot for entry in table[column]]
        for row_index in range(dimension):
            if row_index != column:
                multiplier = table[row_index][column]
                table[row_index] = [left - multiplier * right for left, right in zip(table[row_index], table[column])]
    return [row[dimension:] for row in table]

def tree_distances(order, edges):
    adjacent = [[] for unused in range(order)]
    for left, right in edges:
        adjacent[left].append(right)
        adjacent[right].append(left)
    result = []
    for root in range(order):
        distances = [-1] * order
        distances[root] = 0
        queue = [root]
        for vertex in queue:
            for neighbor in adjacent[vertex]:
                if distances[neighbor] < 0:
                    distances[neighbor] = distances[vertex] + 1
                    queue.append(neighbor)
        assert min(distances) >= 0
        result.append(distances)
    return result

def ceiling_integer(numerator, denominator=1):
    return -(-numerator // denominator)

def scaled_up(value):
    return ceiling_integer(ORDER_UNIT * value.numerator, value.denominator)

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
    eta, threshold = (Q(323, 800000), Q(9, 400))
    root2, root3 = (Q(283, 200), Q(1733, 1000))
    k = [Q(0, 1), Q(1, 1), Q(109, 25), Q(561, 50), Q(586, 25), Q(221, 5)]
    gamma = [2 * eta, root3 * Q(367, 100) * eta, root3 * Q(367, 100) * eta]
    low = [1 - 8 * k[2] * eta - gamma[0] ** 2, 1 - 6 * k[2] * eta - gamma[1] ** 2, 1 - 6 * k[2] * eta - gamma[2] ** 2]
    cross = root3 * Q(491, 20) * eta + gamma[0] * gamma[1]
    far = Q(353, 10) * eta + gamma[1] * gamma[2]
    lower = [[low[0], -cross, -cross], [-cross, low[1], -far], [-cross, -far, low[2]]]
    assert all((lower[i][i] > sum((abs(lower[i][j]) for j in range(3) if j != i)) for i in range(3)))
    inverse = rational_inverse(lower)
    assert all((entry >= 0 for row in inverse for entry in row))
    inverse_fixed = [[scaled_up(entry) for entry in row] for row in inverse]
    edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9), (6, 10), (6, 11), (8, 12), (8, 13)]
    distance = tree_distances(14, edges)

    def regression_bound(root_bound, driver_bounds):
        return [scaled_up(root_bound)] + [scaled_up(sum(((driver_bounds[j] + root_bound * gamma[j]) * inverse[j][i] for j in range(3)))) for i in range(3)]
    defect = []
    for vertex in (0, 1, 2, 3, 6, 8):
        ratios = [root2 ** d for d in distance[vertex]]
        drivers = [max(ratios[2] + ratios[3], 2 * ratios[1]), root3 * max(ratios[6], ratios[7]), root3 * max(ratios[8], ratios[9])]
        incident_budget = (5 if relaxed_endings and vertex in (6, 8) else 3) * eta
        defect.append(regression_bound(incident_budget * ratios[0], [incident_budget * x for x in drivers]))
    root_half = regression_bound(eta / 2, [eta * k[2], eta * root3 * (1 + k[3]) / 2, eta * root3 * (1 + k[3]) / 2])
    leaf_half = []
    for vertex in (4, 10, 12):
        ds = distance[vertex]
        path_error = [value * Q(d + 1, d) if relaxed_endings and vertex != 4 and d else value for d, value in enumerate(k)]
        leaf_half.append(regression_bound(eta * path_error[ds[0]] / 2, [eta * (path_error[ds[2]] + path_error[ds[3]] + 2 * path_error[ds[1]]) / 2, eta * root3 * (Q(74, 5) if vertex == 4 else (Q(289, 10) if relaxed_endings else Q(26, 1)) if vertex == 12 else (path_error[ds[6]] + path_error[ds[7]]) / 2), eta * root3 * (Q(74, 5) if vertex == 4 else (Q(289, 10) if relaxed_endings else Q(26, 1)) if vertex == 10 else (path_error[ds[8]] + path_error[ds[9]]) / 2)]))
    permutation = (0, 2, 1)
    assert gamma[1] == gamma[2]
    assert all((inverse_fixed[i][j] == inverse_fixed[permutation[i]][permutation[j]] for i in range(3) for j in range(3)))

    def reflected(row):
        return [row[0], row[1], row[3], row[2]]
    assert defect[0] == reflected(defect[0]) and defect[1] == reflected(defect[1])
    assert defect[2] == reflected(defect[3]) and defect[4] == reflected(defect[5])
    assert root_half == reflected(root_half)
    assert leaf_half[0] == reflected(leaf_half[0]) and leaf_half[1] == reflected(leaf_half[2])
    assert Q(63, 2500) ** 2 * (1 + 2 * k[2] * eta) < Q(7, 44) ** 4
    unit = ORDER_UNIT
    assert unit == S
    strip = strip_table()

    if relaxed_endings:
        assert 5 * (3 - Q(7071, 2500)) * eta < Q(931, 50000) ** 2
    mass_unit = Q(63, 2500 * 400000 * unit ** 5 * 1000000 ** 6)
    return {'strip': strip, 'unit': mass_unit, 'gamma': gamma, 'inverse': inverse_fixed, 'rows': defect + [root_half] + leaf_half, 'theta': eta, 'threshold': threshold, 'paths': k, 'lower': lower}


nine_parameters = SimpleNamespace(parameters=parameters, RADIUS_ROWS=RADIUS_ROWS, CALIBRATION=CALIBRATION, root_ending_probability=root_ending_probability)


from fractions import Fraction as F
from math import factorial, isqrt
S=10**9
COUNT=30*60*820
def ceil(n,d=1):return (n+d-1)//d


def calibrate(p,data,radii,relaxed):
 assert p['theta']==F(323,800000) and p['threshold']==F(9,400)
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
 assert eigen>F(19,20)
 first=k2+r3*(1+k3)+sum(gamma)/2
 p1,p2,p3,p4=(F(2),F(3,2)*k2,F(4,3)*k3,F(5,4)*k4) if relaxed else (F(1),k2,k3,k4)
 end=(p2+3*p4)/2+r3*((p1+p3)/2+(F(289,10) if relaxed else 26))+p3*sum(gamma)/2
 assert first<26 and end<113
 pair=(p2+p4)*theta/4+p3*theta*theta/4+26*113*theta*theta/F(19,20)
 assert pair<F(11,2500) and F(1,12)+F(11,2500)<F(9,100)
 return {'normal_prefactor_calibrated':True,'all_eight_defect_masses_calibrated':True,'pair_covariance_calibrated':True,'reconstruction_and_projection_guards_calibrated':True}

class Context:

    def __init__(self, relaxed):
        ns = vars(nine_parameters)
        p = ns['parameters'](relaxed)
        self.strip, self.unit, self.inv, self.rows = (p['strip'], p['unit'], p['inverse'], p['rows'])
        self.theta, self.threshold = (p['theta'], p['threshold'])
        gamma = p['gamma']
        self.radii = ns['RADIUS_ROWS']
        self.pair = ns['root_ending_probability']
        self.guards = calibrate(p, ns['CALIBRATION'], self.radii, relaxed)
        assert len(self.rows) == 10 and all((len(row) == 4 and min(row) >= 0 for row in self.rows))
        assert all((self.rows[i][2] == self.rows[i][3] for i in (0, 1, 6, 7)))
        assert all((self.rows[i][2] <= self.rows[i][3] for i in (2, 4, 8)))
        assert all((self.inv[i][j] == self.inv[j][i] for i in range(3) for j in range(3)))
        assert self.inv[0][1] == self.inv[0][2] and self.inv[1][1] == self.inv[2][2]
        self.weights = [[self.radii[r][1] ** 4 * self.radii[e][2 if relaxed else 1] ** 2 for e in range(4)] for r in range(4)]
        self.mid = [ceil(972 * r[0], 1000) for r in self.radii]
        self.close = [[ceil(1208 * r[0] + 500 * e[0], 1000) for e in self.radii] for r in self.radii]
        self.zterms = [tuple((row[2] * (j + 1) * 100000000 for row in self.rows)) for j in range(40)]
        self.hterms = [tuple((row[3] * (j + 1) * 100000000 for row in self.rows)) for j in range(40)]
        addresses = [(i, j) for i in range(40) for j in range(i, 40)]
        self.address = [[i*40-i*(i+1)//2+j for j in range(40)] for i in range(40)]
        self.defect_additions = [
            tuple(z+h+S-1 for z, h in zip(self.zterms[i], self.hterms[j]))
            for i, j in addresses
        ]
        self.density_parameters = (self.threshold, tuple(map(tuple, self.inv)), tuple(gamma))
        self.density_rows = []
        self.base = []
        ts = ceil(self.threshold.numerator * S, self.threshold.denominator)
        for ix in range(30):
            xl = ix * 5000000 if ix < 20 else (ix - 10) * 10000000
            xh = xl + (5000000 if ix < 20 else 10000000)
            root = xh + ts
            drift = tuple((ceil(g.numerator * root, g.denominator) for g in gamma))
            ends = [(2*i+1)*100000000 + 2*drift[1] for i in range(40)]
            diag, cross = self.inv[1][1], 2*self.inv[1][2]
            self.density_rows.append([
                (ends[i]+ends[j], diag*(3*(ends[i]**2+ends[j]**2)+2*100000000**2)+3*cross*ends[i]*ends[j],
                 (1 if ix < 20 else 2)*(1 if i == j else 2))
                for i, j in addresses
            ])
            xy = []
            for iy in range(60):
                yl, yh = (iy * 50000000, (iy + 1) * 50000000)
                D0 = yh + drift[0]
                bases = tuple((row[0] * root + row[1] * D0 + row[2] * drift[1] + row[3] * drift[2] for row in self.rows))
                constants = (yl // 3 - ceil(943 * xh, 1000) - ceil(572 * ts, 10000), 9428 * xl // 10000 + yl // 6 - ceil(572 * ts, 10000), 2828 * yl // 6000 - ceil(5 * xh, 6) - ceil(ts, 6), ceil(5 * xh, 6) + ceil(2830 * yh, 12000) + ceil(ts, 6), 5 * xl // 6 + 2828 * yl // 12000 - ceil(ts, 6), ceil(7072 * xh, 10000) + ceil(yh, 4) + ceil(2929 * ts, 10000))
                assert constants[3] + constants[4] >= 0
                # Squares and cross products are in units of 1/12.
                root_mid2 = xl+xh+2*ts
                y_mid2 = yl+yh+2*drift[0]
                root_square = 3*root_mid2**2+(xh-xl)**2
                y_square = 3*y_mid2**2+(yh-yl)**2
                dense = (root_square*S+self.inv[0][0]*y_square,
                         6*self.inv[0][1]*y_mid2, drift[1], 1 if ix < 20 else 2)
                xy.append((bases, constants, dense))
            self.base.append(xy)

    def geometry_pair(self, ix, iy, iz, ih, relaxed):
        bases, c, _ = self.base[ix][iy]
        addr = self.address[iz][ih]
        contributions = self.defect_additions[addr]
        d0, da, db, du, dw, dz, h0, ha, hw, hz = [
            (a+b)//S for a, b in zip(bases, contributions, strict=True)
        ]
        ob = ceil(d0, 3)
        bb = ceil(2830 * d0, 6000)
        f = c[0] - ob
        branch = c[1] - ob - h0
        leaf = c[2] - bb - ceil(da, 2) - ha
        common = c[3] + bb + ceil(1415 * h0, 1000)
        sl = 2885 * (iz * 100000000) // 10000 - common - ceil(db, 2)
        sr = 2885 * (ih * 100000000) // 10000 - common - ceil(du, 2)
        endcommon = c[5] + ceil(d0, 2) + ceil(3 * h0, 2)
        left = 4082 * (iz * 100000000) // 10000 - endcommon - ceil(7075 * db, 10000)
        right = 4082 * (ih * 100000000) // 10000 - endcommon - ceil(7075 * du, 10000)
        rb = relaxed.base[ix][iy][0]
        rc = relaxed.defect_additions[addr]
        ldw, ldz, lhw, lhz = ((rb[i]+rc[i])//S for i in (4, 5, 8, 9))
        shared = f, branch, leaf, min(sl, sr)
        return ((*shared, left-ceil(dw, 2)-hw, right-ceil(dz, 2)-hz),
                (*shared, left-ceil(ldw, 2)-lhw, right-ceil(ldz, 2)-lhz))

    def density(self, ix, iy, iz, ih):
        base, linear, _, _ = self.base[ix][iy][2]
        coordinate_sum, quadratic, multiplicity = self.density_rows[ix][self.address[iz][ih]]
        num = base + linear * coordinate_sum + quadratic
        exponent = ceil(num, 24 * S * S)
        v = S - ceil(exponent, 4096)
        assert v > 0
        for _ in range(12):
            v = v * v // S
        return v * multiplicity

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
        assert strict.density_parameters == relaxed.density_parameters
        assert all(strict.rows[i] == relaxed.rows[i] for i in (0, 1, 2, 3, 6, 7))
        addresses = [(z, h) for z in range(40) for h in range(z, 40)]
        strict_total = relaxed_total = count = 0
        for index, code in enumerate(codes):
            if code:
                xy, zh = divmod(index, 820)
                x, y = divmod(xy, 60)
                z, h = addresses[zh]
                if not 24 <= z + h <= 43:
                    continue
                density = strict.density(x, y, z, h)
                assert 1 <= code <= 192
                q, kind = divmod(code-1, 3)
                q, root_kind = divmod(q, 4)
                r, e = divmod(q, 4)
                sg, rg = strict.geometry_pair(x, y, z, h, relaxed)
                sm = strict.candidate(sg, r, e, root_kind, kind)
                rm = relaxed.candidate(rg, r, e, root_kind, kind)
                assert sm > 0 and rm > 0
                strict_total += density*sm
                relaxed_total += density*rm
                count += 1
        results = []
        for phase, context, total, floor in (
            ("strict", strict, strict_total, F(1241, 10**9)),
            ("relaxed", relaxed, relaxed_total, F(903, 10**9)),
        ):
            exact = context.unit*total
            assert exact > floor
            row = {"phase": phase, "eta": "323/400000", "threshold_abs_max": "9/400",
                   "exact_lower_sum": str(exact), "used_floor": str(floor),
                   "nonzero_codes": count, "calibration": context.guards}
            results.append(row)
            if progress is not None:
                progress(row)
        return results
    return SimpleNamespace(certificate=certificate)


nine = _nine()


def main():
    import json
    from resource import RUSAGE_SELF, getrusage
    from time import monotonic, process_time

    started = monotonic()

    def progress(phase, value):
        print(json.dumps({"phase": phase, "elapsed": monotonic()-started, "output": value}), flush=True)

    def check(component, phase):
        result = component.certificate()
        progress(phase, result)
        return result

    progress("actual source intervals", source_intervals.certificate())
    base_result = base.base_components(progress)
    signed_result = signed.certificate(lambda row: progress("signed", row))
    choices = (here/"nine-choices.txt").read_bytes()
    assert sha256(choices).hexdigest() == "c04c6170c36f7968dedc7b37ed0f8835bc2687a28d462950c6578f3725f2ca4b"
    nine_result = nine.certificate(choices, lambda row: progress("nine " + row["phase"], row))
    d7_result = check(d7, "D7")
    outer_result = check(outer, "complete outer")
    adjacent = check(price.adjacent, "adjacent Price")
    baseline2 = check(price.baseline2, "D2 global baseline")
    arms2 = check(price.arms2, "D2 baseline and size-two arms")
    long_paths = check(price.long_paths, "D3-D6 Price")
    size_four = check(price.size_four, "D2-D3 size-four arms")
    final = price.allocation.final_allocation(adjacent, baseline2, arms2, long_paths, size_four,
                                              d7_result, base_result, signed_result, nine_result, outer_result)
    progress("exact final composition", final)
    print(json.dumps({"completed": "one-band numerical theorem", "elapsed": monotonic()-started,
                      "cpu_seconds": process_time(), "peak_rss_kib": getrusage(RUSAGE_SELF).ru_maxrss,
                      "circuit_coefficient": final["circuit_coefficient"]}), flush=True)


if __name__ == "__main__":
    main()
