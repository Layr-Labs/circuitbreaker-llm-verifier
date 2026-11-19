#!/usr/bin/env python3
"""Standalone exact certificate; no network or repository imports."""


"""Standalone directed quadrature for the 47-term crossing-path bound."""
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction
from math import cos, factorial, pi
from typing import NamedTuple


Q_LOW = Context(prec=100, rounding=ROUND_FLOOR)
Q_HIGH = Context(prec=100, rounding=ROUND_CEILING)
Q_MID = Context(prec=100, rounding=ROUND_HALF_EVEN)
GRID = 1 << 224
NODE_COUNT = 93
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
    """Return (k! b^k P_k) at degrees 92 and 93, x=a/b."""
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
    for j in range(1,48):
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
        if j == 46:
            even_lower = total.low
        if j == 47:
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

    assert even_lower > Decimal('0.1306475664231676642')
    assert total.high < Decimal('0.130647566423167665')
    return {"gamma_upper": str(total.high), "gamma_even_lower": str(even_lower),
            "initial_and_accumulated_debit_grid_units": debit,
            "per_step_debit_grid_units": step_debit}

from fractions import Fraction

SCALE = 10**9
SCORE_SCALE = 10**5

# Coefficients of the six worst signed conditional means.
MEAN_ROWS = (
    (94329, -16039, 1921, 315),
    (83603, -21780, 4113, 917),
    (83603, 47316, -28012, 917),
    (71567, 50606, -39010, 1611),
    (-94280, -32705, 522, 315),
    (-83332, -45350, -28012, 917),
)


def gaussian_weight(exponent):
    """A downward fixed-point approximation to (1-exponent/4096)^4096."""
    start = (1 - exponent / 4096) * SCALE
    weight = start.numerator // start.denominator
    assert weight > 0
    for unused in range(12):
        weight = weight * weight // SCALE
    return weight


def strip_table():
    """Centered strip lower bounds at the standardized arguments j/10000."""
    result = []
    for j in range(16001):
        # v=j/10000; integrate 1-v^2/2+v^4/8-v^6/48 from 0 to v.
        den = 10000**6
        polynomial = 1680*den - 280*j*j*10000**4 + 42*j**4*10000**2 - 5*j**6
        result.append(797*j*polynomial*SCALE // (10000*1680000*den))
    return result


def rectangle_controls():
    return tuple((radius, 972*radius//1000, (2043*radius+999)//1000,
                  (radius+2)//3, (radius+1)//2, defect_mass**4)
                 for radius, defect_mass in ((18_000_000,748),
                                             (25_040_000,889), (33_650_000,967)))


def rectangle(margins, strips, controls):
    first, second, small_left, terminal, root, small_right = margins
    best = 0
    small = small_left if small_left<small_right else small_right
    for radius, cutoff, closing, third, half, residual_mass in controls:
        if small <= cutoff or 3*root <= radius:
            continue
        outer = terminal - closing
        if outer <= 0:
            continue
        capacity = first-third
        alternative = (second-half)*1000//1415-third
        if alternative<capacity:
            capacity=alternative
        if capacity <= 0:
            continue
        index=346*outer//10**7
        outer_mass = strips[index if index<16000 else 16000] * residual_mass
        for numerator in (5, 8):
            first_half_difference = capacity*numerator//10
            second_half_difference = second-half-(1415*(first_half_difference+third)+999)//1000
            index=346*first_half_difference//10**7
            first_mass = strips[index if index<16000 else 16000]
            index=346*second_half_difference//10**7
            second_mass = strips[0 if index<0 else index if index<16000 else 16000]
            candidate=first_mass*second_mass*outer_mass
            if candidate>best:
                best=candidate
    return best


def seven_key_floor(threshold):
    assert threshold in (Fraction(1, 50), Fraction(3, 100))
    bias_ratio = threshold / Fraction(11, 200)
    biases = tuple((bias_ratio*h).__ceil__() for _, _, _, h in MEAN_ROWS)
    x_step, yz_step = SCALE//250, SCALE//50
    assert x_step%SCORE_SCALE == yz_step%SCORE_SCALE == SCALE%SCORE_SCALE == 0
    x_weights = [gaussian_weight((threshold+Fraction(i+1,250))**2/2) for i in range(50)]
    y_weights = [gaussian_weight(Fraction(519*(j+1)**2, 1000*50**2)) for j in range(180)]
    z_weights = [gaussian_weight(Fraction(1033*(k+1)**2, 2000*50**2)) for k in range(340)]
    strips = strip_table()

    controls = rectangle_controls()
    mass_unit = Fraction(621, 10_000*250*50*50*SCALE**6*1000**4)
    target = Fraction(283,10**9) if threshold == Fraction(1,50) else Fraction(273249,10**12)
    stopping_integer = (target/mass_unit).__floor__()
    sum_cells = 0
    for i, x_weight in enumerate(x_weights):
        for j, y_weight in enumerate(y_weights):
            # Each rounded margin must be positive. Its numerator is affine
            # in the z index, so the retained cells form an integer interval.
            intercepts = []
            slopes = []
            low, high = 0, len(z_weights) - 1
            for (a,b,c,_), bias in zip(MEAN_ROWS, biases):
                at_zero = -(a*x_step*(i+(a>0)) + b*yz_step*(j+(b>0))
                            + c*yz_step*(c>0) + bias*SCALE)
                per_cell = -c*yz_step
                if per_cell > 0:
                    low = max(low, -((at_zero-SCORE_SCALE)//per_cell))
                else:
                    high = min(high, (at_zero-SCORE_SCALE)//(-per_cell))
                intercepts.append(at_zero//SCORE_SCALE)
                slopes.append(per_cell//SCORE_SCALE)
            a0,a1,a2,a3,a4,a5=intercepts
            b0,b1,b2,b3,b4,b5=slopes
            weighted_z_sum = 0
            for k in range(low, high+1):
                margins = (a0+b0*k,a1+b1*k,a2+b2*k,a3+b3*k,a4+b4*k,a5+b5*k)
                weighted_z_sum += z_weights[k]*rectangle(margins, strips, controls)
            sum_cells += x_weight*y_weight*weighted_z_sum
        if sum_cells > stopping_integer:
            break
    answer = mass_unit*sum_cells
    assert answer > target
    return answer

"""Original exact lower-sum implementation of the stated curvature integrals."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, isqrt

QUANTUM = 10**18
ROOT2_LOW, ROOT2_HIGH = Q(141421,100000), Q(70711,50000)
RESIDUAL_HIGH = Q(8579,50000)

def quantize(value):
    return max(0, value.numerator*QUANTUM//value.denominator)

from math import lcm
_CURV_SCALE = QUANTUM
_CURV_DEGREE = 23
_CURV_COMMON = lcm(*(factorial(j) * (2 * j + 1) for j in range(_CURV_DEGREE + 1)))
_CURV_COEFFICIENTS = [(-1) ** j * (_CURV_COMMON // (factorial(j) * (2 * j + 1))) * 2000000 ** (_CURV_DEGREE - j) for j in range(_CURV_DEGREE + 1)]
_CURV_DENOMINATOR = 5000000 * _CURV_COMMON * 2000000 ** _CURV_DEGREE

@lru_cache(None)
def interval_at_clamped_grid(t):
    value = _CURV_COEFFICIENTS[-1]
    for coefficient in reversed(_CURV_COEFFICIENTS[:-1]):
        value = value * t * t + coefficient
    return max(0, _CURV_SCALE * 3989 * t * value // _CURV_DENOMINATOR)

def gaussian_interval(thousandths):
    return interval_at_clamped_grid(min(3300, max(0, thousandths)))

@lru_cache(None)
def exponential_lower(value):
    assert 0 <= value <= 64
    argument = value / 64
    a, b = (argument.numerator, argument.denominator)
    result = -1
    power = b
    for j in range(6, -1, -1):
        result = result * a + (-1) ** j * (5040 // factorial(j)) * power
        power *= b
    result = max(0, _CURV_SCALE * result // (5040 * b ** 7))
    for unused in range(6):
        result = result * result // _CURV_SCALE
    return result

def interval_scaled(magnitude, interpolation, deviation):
    a, b = (magnitude.numerator, magnitude.denominator)
    c, d = (deviation.numerator, deviation.denominator)
    m, n = (interpolation.numerator, interpolation.denominator)
    numerator = 1000000 * a * a * d * d * (n - m)
    _CURV_DENOMINATOR = b * b * c * c * m
    return gaussian_interval(isqrt(numerator // _CURV_DENOMINATOR))



def round_argument_up(value):
    return -((-1000*value.numerator)//value.denominator)

def density_absolute(value, kind, rounded=False):
    if rounded:
        value=Q(round_argument_up(value),1000)
    rate,prefactor = {"path":(Q(6),27638),
                       "contrast":(Q(3,2),13819),
                       "middle":(Q(33,2),45834)}[kind]
    return prefactor*exponential_lower(rate*value*value)//10000



def rational_sqrt_lower(value):
    return Q(isqrt(value.numerator*QUANTUM**2//value.denominator),QUANTUM)

def scalar_normalizations():
    def arctangent(value,last):
        return sum(((-1)**k*value**(2*k+1)/(2*k+1)
                    for k in range(last+1)),Q(0))
    pi_upper=16*arctangent(Q(1,5),14)-4*arctangent(Q(1,239),15)
    pi_lower=16*arctangent(Q(1,5),15)-4*arctangent(Q(1,239),14)
    assert Q(6283,2000)<pi_lower<pi_upper<Q(3927,1250)
    assert ROOT2_LOW**2<2<ROOT2_HIGH**2
    for coefficient,squared_norm in [(Q(13819,5000),24),
                                     (Q(13819,10000),6),
                                     (Q(3989,5000),2),
                                     (Q(22917,5000),66)]:
        assert coefficient**2*pi_upper<squared_norm
    assert pi_upper*Q(477,500)<3
    assert pi_upper**2*11*Q(1439,5000)**2<9
    return pi_upper

def convolution_prefix(left, right, length):
    """Exact initial coefficients of a nonnegative integer convolution."""
    assert 0 <= length <= len(left)+len(right)-1
    assert min(left) >= 0 and min(right) >= 0
    maximum = min(len(left),len(right))*max(left)*max(right)
    width = max(1, (maximum.bit_length()+7)//8)
    def pack(values):
        return int.from_bytes(b"".join(value.to_bytes(width,"little") for value in values), "little")
    product = pack(left)*pack(right)
    raw = product.to_bytes(width*(len(left)+len(right)-1), "little")
    return [int.from_bytes(raw[i:i+width], "little") for i in range(0,width*length,width)]


def original_pair_integrals():
    """Fine size-(2,2) adjacent and no-attachment distance-two sums."""
    denominator=200
    roots=140
    children=180
    endings=roots+3*children
    first_density=[density_absolute(Q(i+1,denominator),"path") for i in range(roots)]
    terminal_density=[density_absolute(Q(i,denominator),"contrast") for i in range(endings)]
    cells=[]
    for root in range(roots):
        row=[]
        for child in range(children):
            argument=(2*Q(child+1,denominator)+Q(282843,100000)*Q(root+1,denominator))
            density=density_absolute(argument,"contrast",True)
            coefficient=2*(first_density[root]*density//QUANTUM)//denominator**2
            margin=int(Q(282842,100000)*child+root)
            row.append((margin,coefficient))
        cells.append(row)
    central_density=[density_absolute(Q(i+1,denominator),"middle") for i in range(160)]
    outside_density=[density_absolute(Q(i+1,denominator),"contrast") for i in range(160)]
    adjacent=Q(0);distance_two=Q(0);previous=Q(0)
    for index in range(1,80):
        interpolation=Q(index,80)**2
        first=[interval_scaled(Q(i,denominator),interpolation,Q(43,500)) for i in range(roots)]
        second=[interval_scaled(Q(i,denominator),interpolation,Q(26,125)) for i in range(children)]
        last=[interval_scaled(Q(i,denominator),interpolation,Q(84,125)) for i in range(endings)]
        convolution=[x//(QUANTUM*denominator) for x in
                     convolution_prefix([0]+terminal_density[1:],last,endings)]
        side=0
        for root,row in enumerate(cells):
            subtotal=sum(((coefficient*second[child]//QUANTUM)
                            *convolution[margin]//QUANTUM)
                         for child,(margin,coefficient) in enumerate(row))
            side+=first[root]*subtotal//QUANTUM
        adjacent+=(interpolation-previous)*Q(side,QUANTUM)**2*(1-interpolation)**3
        central=[interval_scaled(Q(i,denominator)*Q(282842,17158),interpolation,Q(1)) for i in range(160)]
        strip=[interval_scaled(Q(i,denominator)*Q(100000,17158),interpolation,Q(1)) for i in range(160)]
        middle_arms=[0]+convolution_prefix(outside_density,strip,159)
        middle_sum=0
        for middle,raw in enumerate(middle_arms):
            arm=raw//(QUANTUM*denominator)
            middle_sum+=((central_density[middle]*central[middle]//QUANTUM)
                          *(arm*arm//QUANTUM)//QUANTUM)//denominator
        ratio=interpolation/(1-interpolation)
        determinant=(1-interpolation)*rational_sqrt_lower(1/(1+Q(391,100)*ratio))
        distance_two+=(interpolation-previous)*Q(middle_sum,QUANTUM)*determinant
        previous=interpolation
    return Q(477,500)*adjacent,Q(1439,5000)*distance_two

@lru_cache(None)
def coarse_two_vertex_arms():
    """Shared positive arms for the distance-two and distance-three sums."""
    denominator=80
    def outward_density(q):
        return density_absolute(q,"contrast")
    roots=[[2*outward_density(Q(x+1+2*(a+1),denominator))//denominator
            for a in range(48)] for x in range(64)]
    transition=[[density_absolute((Q(b+1)+ROOT2_HIGH*(a+1))/denominator,"path",True)//denominator
                 for b in range(48)] for a in range(48)]
    minimum=[[int(2*ROOT2_LOW*b+a) for b in range(48)] for a in range(48)]
    maximum=max(63,max(map(max,minimum)))
    closing=[0]+[outward_density(Q(z,denominator))//denominator for z in range(1,maximum+1)]
    source_density=[density_absolute(Q(x+1,denominator),"middle")//denominator for x in range(64)]
    q_high=Q(70711,100000)
    deviations=[Q(0)]+[RESIDUAL_HIGH/2*sum(((j+1)*q_high**j for j in range(depth)),Q(0))
                       for depth in range(1,4)]
    data=[];credit=Q(0);previous=Q(0)
    for step in range(1,24):
        interpolation=Q(step,24)**2
        strips=[[interval_scaled(Q(x,denominator),interpolation,deviations[depth]) for x in range(48)]
                 for depth in (1,2)]
        terminal=[interval_scaled(Q(z,denominator),interpolation,2*deviations[3]) for z in range(maximum+1)]
        simple=[interval_scaled(Q(z,denominator),interpolation,RESIDUAL_HIGH) for z in range(maximum+1)]
        middle=[interval_scaled(Q(z,denominator),interpolation,RESIDUAL_HIGH/(2*ROOT2_LOW)) for z in range(64)]
        close_sum=[sum(closing[z]*terminal[m-z] for z in range(1,m+1))//QUANTUM for m in range(maximum+1)]
        baseline=[sum(closing[z]*simple[x-z] for z in range(1,x+1))//QUANTUM for x in range(64)]
        inner_by_a=[sum((transition[a][b]*strips[1][b]//QUANTUM)*close_sum[minimum[a][b]] for b in range(48))//QUANTUM for a in range(48)]
        arms=[]
        for x in range(64):
            total=0
            for a in range(48):
                total+=(roots[x][a]*strips[0][a]//QUANTUM)*inner_by_a[a]
            arms.append(total//QUANTUM)
        density_sum=0
        for x in range(64):
            # At least one attachment; this positive expansion excludes baseline^2.
            additional=2*(1-interpolation)*Q(baseline[x]*arms[x],QUANTUM**2)+(1-interpolation)**2*Q(arms[x]**2,QUANTUM**2)
            density_sum+=Q(source_density[x]*middle[x],QUANTUM**2)*additional
        ratio=interpolation/(1-interpolation)
        determinant=(1-interpolation)*rational_sqrt_lower(1/(1+Q(391,100)*ratio))
        credit+=Q(1439,5000)*(interpolation-previous)*determinant*density_sum
        data.append((interpolation,baseline,arms))
        previous=interpolation
    return credit,tuple(data)

def adjacent_larger_components():
    from operator import mul
    denominator=80
    dimensions=48
    assert denominator==80 and ROOT2_HIGH==Q(70711,50000)
    # Every original argument has denominator 50000*80. Rounding up
    # to a thousandth is exactly ceiling(numerator/4000).
    @lru_cache(None)
    def density_at_grid(index):
        return 27638*exponential_lower(Q(6*index*index,10**6))//10000//denominator
    def rounded_density(numerator):
        return density_at_grid((numerator+3999)//4000)
    root=[rounded_density(50000*(a+1)) for a in range(dimensions)]
    first=[[rounded_density(50000*(b+1)+70711*(a+1))
            for b in range(dimensions)] for a in range(dimensions)]
    kernel=[[[rounded_density(50000*(c+1)+70711*(b+1)+25000*(a+1))
              for a in range(dimensions)] for c in range(dimensions)] for b in range(dimensions)]
    minimum=[[int(2*ROOT2_LOW*b+a) for b in range(dimensions)] for a in range(dimensions)]
    maximum=max(map(max,minimum))
    closing=[0]+[density_absolute(Q(z,denominator),"contrast")//denominator for z in range(1,maximum+1)]
    q_high=Q(70711,100000)
    deviations=[Q(0)]+[RESIDUAL_HIGH/2*sum(((j+1)*q_high**j for j in range(depth)),Q(0))
                       for depth in range(1,10)]
    total=Q(0);previous=Q(0)
    for step in range(1,24):
        interpolation=Q(step,24)**2
        strips={depth:[interval_scaled(Q(a,denominator),interpolation,deviations[depth]) for a in range(dimensions)]
                for depth in range(1,9)}
        state=[[(root[a]*strips[1][a]//QUANTUM)*first[a][b]//QUANTUM*strips[2][b]//QUANTUM
                for b in range(dimensions)] for a in range(dimensions)]
        sizes={}
        for size in range(2,9):
            if size%2==0:
                terminal=[interval_scaled(Q(z,denominator),interpolation,2*deviations[size+1]) for z in range(maximum+1)]
                convolution=[sum(closing[z]*terminal[m-z] for z in range(1,m+1))//QUANTUM for m in range(maximum+1)]
                sizes[size]=Q(sum(state[a][b]*convolution[minimum[a][b]] for a in range(dimensions) for b in range(dimensions))//QUANTUM,QUANTUM)
            if size<8:
                state=[[(sum(map(mul, column, weights))//QUANTUM)*strip//QUANTUM
                        for weights,strip in zip(row,strips[size+1])]
                       for column,row in zip(zip(*state),kernel)]
        base=(1-interpolation)*sizes[2]
        extra=sum(((1-interpolation)**(size//2)*sizes[size] for size in (4,6,8)),Q(0))
        total+=Q(477,500)*(interpolation-previous)*(1-interpolation)*(2*base*extra+extra*extra)
        previous=interpolation
    return total

def distance_three_integral():
    extra_two,arms_by_time=coarse_two_vertex_arms()
    denominator=80
    density=[[81027*exponential_lower(18*(Q(a+1,denominator)**2+Q(b+1,denominator)**2
                           +ROOT2_HIGH*Q((a+1)*(b+1),denominator**2)))//(10000*denominator**2)
              for b in range(64)] for a in range(64)]
    margin=[[int(2*ROOT2_LOW*a+b) for b in range(64)] for a in range(64)]
    maximum=max(map(max,margin))
    total=Q(0);previous=Q(0)
    for interpolation,baseline,arms in arms_by_time:
        side=[interpolation.denominator*base+(interpolation.denominator-interpolation.numerator)*arm
              for base,arm in zip(baseline,arms)]
        strips=[interval_scaled(Q(m,denominator),interpolation,RESIDUAL_HIGH) for m in range(maximum+1)]
        summation=0
        for a in range(64):
            for b in range(64):
                weight=density[a][b]*strips[margin[a][b]]//QUANTUM
                weight=weight*strips[margin[b][a]]//QUANTUM
                summation+=weight*side[a]*side[b]
        probability=Q(summation,QUANTUM**3*interpolation.denominator**2)
        ratio=interpolation/(1-interpolation)
        determinant=(1-interpolation)*rational_sqrt_lower(1/((1+5*ratio)*(1+9*ratio)))
        total+=Q(9,40)*(interpolation-previous)*determinant*probability
        previous=interpolation
    return total

def curvature_replay():
    import sys,time
    started=time.monotonic()
    def progress(label,*values):
        print(label,round(time.monotonic()-started,3),*(str(value) for value in values),file=sys.stderr,flush=True)
    scalar_normalizations()
    adjacent_baseline,distance_two_baseline=original_pair_integrals()
    assert quantize(adjacent_baseline)==3607233186504135
    assert quantize(distance_two_baseline)==2203451254088009
    progress('curvature C6 baselines',adjacent_baseline,distance_two_baseline)
    adjacent_extra=adjacent_larger_components()
    assert adjacent_extra>Q(143,10**6)
    progress('curvature C7 adjacent extra',adjacent_extra)
    distance_two_extra,unused=coarse_two_vertex_arms()
    assert distance_two_extra>Q(7,5000)
    progress('curvature C7 distance-two extra',distance_two_extra)
    distance_three=distance_three_integral()
    assert distance_three>Q(27,125000)
    progress('curvature C7 distance-three',distance_three)
    assert quantize(adjacent_baseline)==3607233186504135
    assert quantize(distance_two_baseline)==2203451254088009
    assert adjacent_extra>Q(143,10**6)
    assert distance_two_extra>Q(7,5000)
    assert distance_three>Q(27,125000)
    return {"adjacent_baseline":str(adjacent_baseline),
            "distance_two_baseline":str(distance_two_baseline),
            "adjacent_extra":str(adjacent_extra),
            "distance_two_extra":str(distance_two_extra),
            "distance_three":str(distance_three)}

"""Original rational implementations of the nine-key and weak-pair boxes."""
from fractions import Fraction as Q
from itertools import product
from math import factorial, prod

ORDER_UNIT=10**9

def rational_inverse(matrix):
    dimension=len(matrix)
    table=[list(row)+[Q(row_index==column) for column in range(dimension)]
           for row_index,row in enumerate(matrix)]
    for column in range(dimension):
        pivot=table[column][column]
        assert pivot>0
        table[column]=[entry/pivot for entry in table[column]]
        for row_index in range(dimension):
            if row_index!=column:
                multiplier=table[row_index][column]
                table[row_index]=[left-multiplier*right
                                  for left,right in zip(table[row_index],table[column])]
    return [row[dimension:] for row in table]

def tree_distances(order,edges):
    adjacent=[[] for unused in range(order)]
    for left,right in edges:
        adjacent[left].append(right);adjacent[right].append(left)
    result=[]
    for root in range(order):
        distances=[-1]*order;distances[root]=0;queue=[root]
        for vertex in queue:
            for neighbor in adjacent[vertex]:
                if distances[neighbor]<0:
                    distances[neighbor]=distances[vertex]+1;queue.append(neighbor)
        assert min(distances)>=0
        result.append(distances)
    return result

def ceiling_integer(numerator,denominator=1):
    return -((-numerator)//denominator)

def scaled_up(value):
    return ceiling_integer(ORDER_UNIT*value.numerator,value.denominator)

def minimal_nine_key_floor(relaxed_endings=False):
    eta,threshold=Q(19,80000),Q(3,100)
    root2,root3=Q(283,200),Q(1733,1000)
    k=[Q(0),Q(1),Q(283,50),Q(18),Q(1132,25),Q(100)]
    gamma=[2*eta,root3*k[2]*eta,root3*k[2]*eta]
    low=[1-8*k[2]*eta-gamma[0]**2,
         1-6*k[2]*eta-gamma[1]**2,
         1-6*k[2]*eta-gamma[2]**2]
    cross=root3*(1+3*k[3])*eta+gamma[0]*gamma[1]
    far=6*k[4]*eta+gamma[1]*gamma[2]
    lower=[[low[0],-cross,-cross],[-cross,low[1],-far],[-cross,-far,low[2]]]
    assert all(lower[i][i]>sum(abs(lower[i][j]) for j in range(3) if j!=i) for i in range(3))
    inverse=rational_inverse(lower)
    assert all(entry>=0 for row in inverse for entry in row)
    inverse_fixed=[[scaled_up(entry) for entry in row] for row in inverse]
    edges=[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),
           (3,8),(3,9),(6,10),(6,11),(8,12),(8,13)]
    distance=tree_distances(14,edges)
    def regression_bound(root_bound,driver_bounds):
        return [scaled_up(root_bound)]+[
            scaled_up(sum((driver_bounds[j]+root_bound*gamma[j])*inverse[j][i]
                           for j in range(3))) for i in range(3)]
    defect=[]
    for vertex in (0,1,2,3,6,8):
        ratios=[root2**d for d in distance[vertex]]
        drivers=[max(ratios[2]+ratios[3],2*ratios[1]),
                 root3*max(ratios[6],ratios[7]),root3*max(ratios[8],ratios[9])]
        incident_budget=(5 if relaxed_endings and vertex in (6,8) else 3)*eta
        defect.append(regression_bound(incident_budget*ratios[0],[incident_budget*x for x in drivers]))
    root_half=regression_bound(eta/2,[eta*k[2],eta*root3*(1+k[3])/2,eta*root3*(1+k[3])/2])
    leaf_half=[]
    for vertex in (4,10,12):
        ds=distance[vertex]
        path_error=[value*Q(d+1,d) if relaxed_endings and vertex!=4 and d else value
                    for d,value in enumerate(k)]
        leaf_half.append(regression_bound(eta*path_error[ds[0]]/2,[
            eta*(path_error[ds[2]]+path_error[ds[3]]+2*path_error[ds[1]])/2,
            eta*root3*(path_error[ds[6]]+path_error[ds[7]])/2,
            eta*root3*(path_error[ds[8]]+path_error[ds[9]])/2]))
    permutation=(0,2,1)
    assert gamma[1]==gamma[2]
    assert all(inverse_fixed[i][j]==inverse_fixed[permutation[i]][permutation[j]]
               for i in range(3) for j in range(3))
    def reflected(row):
        return [row[0],row[1],row[3],row[2]]
    assert defect[0]==reflected(defect[0]) and defect[1]==reflected(defect[1])
    assert defect[2]==reflected(defect[3]) and defect[4]==reflected(defect[5])
    assert root_half==reflected(root_half)
    assert leaf_half[0]==reflected(leaf_half[0]) and leaf_half[1]==reflected(leaf_half[2])
    assert Q(63,2500)**2*(1+2*k[2]*eta)<Q(7,44)**4
    # Integer strip table encloses .797(v-v^3/6+v^5/40-v^7/336).
    unit=ORDER_UNIT
    strip=[]
    for index in range(16001):
        value=index*100000
        polynomial=1680*unit**6-280*value**2*unit**4+42*value**4*unit**2-5*value**6
        strip.append(797*value*polynomial//(1680000*unit**6))
    def mass(deviation):
        return strip[min(16000,max(0,(346*deviation//100)//100000))]
    def mean(coefficients,root,drivers):
        return ceiling_integer(coefficients[0]*root+sum(a*b for a,b in zip(coefficients[1:],drivers)),unit)
    if relaxed_endings:
        assert 5*(3-Q(7071,2500))*eta < Q(313,20000)**2
    total=0
    mass_unit=Q(63,2500*50000*unit**5*1000**6)
    target=Q(301,10**10) if relaxed_endings else Q(247011,10**12)
    early_integer_cutoff=(target/mass_unit).__floor__()
    threshold_scaled=scaled_up(threshold)
    residuals = tuple((residual, residual_mass, ceiling_integer(residual,3),
                       ceiling_integer(972*residual,1000), ceiling_integer(1708*residual,1000))
                      for residual,residual_mass in ((18000000,748),(25040000,889),(33650000,967)))
    for ix in range(10):
        xl=ix*20000000;xh=xl+20000000
        root=xh+threshold_scaled
        drift=[ceiling_integer(coefficient.numerator*root,coefficient.denominator) for coefficient in gamma]
        for iy in range(30):
            yl=iy*100000000;yh=yl+100000000
            first_before_defect=yl//3-ceiling_integer(943*xh,1000)-ceiling_integer(572*threshold_scaled,10000)
            branch_before_defect=9428*xl//10000+yl//6-ceiling_integer(572*threshold_scaled,10000)
            leaf_before_defect=2828*yl//6000-ceiling_integer(5*xh,6)-ceiling_integer(threshold_scaled,6)
            outer_bias=ceiling_integer(5*xh,6)+ceiling_integer(2830*yh,12000)+ceiling_integer(threshold_scaled,6)
            guard_base=5*xl//6+2828*yl//12000-ceiling_integer(threshold_scaled,6)
            closing_bias=ceiling_integer(7072*xh,10000)+ceiling_integer(yh,4)+ceiling_integer(2929*threshold_scaled,10000)
            drivers=[yh+drift[0],100000000+drift[1],100000000+drift[2]]
            minimum_d0=mean(defect[0],root,drivers)
            minimum_d1=mean(defect[1],root,drivers)
            minimum_ha=mean(leaf_half[0],root,drivers)
            if first_before_defect-ceiling_integer(minimum_d0,3)<=residuals[0][2]:
                continue
            if leaf_before_defect-ceiling_integer(2830*minimum_d0,6000)-ceiling_integer(minimum_d1,2)-minimum_ha<=residuals[0][3]:
                continue
            for iz,ih in ((z,h) for z in range(40) for h in range(z,40)):
                zl=iz*100000000;zh=zl+100000000
                hl=ih*100000000;hh=hl+100000000
                drivers=[yh+drift[0],zh+drift[1],hh+drift[2]]
                d0=mean(defect[0],root,drivers)
                origin_bias=ceiling_integer(d0,3)
                first_margin=first_before_defect-origin_bias
                if first_margin<=residuals[0][2]:
                    continue
                branch_bias=ceiling_integer(2830*d0,6000)
                da=mean(defect[1],root,drivers)
                ha=mean(leaf_half[0],root,drivers)
                leaf_a=leaf_before_defect-branch_bias-ceiling_integer(da,2)-ha
                if leaf_a<=residuals[0][3]:
                    continue
                db,du,dw,dz=[mean(row,root,drivers) for row in defect[2:]]
                h0=mean(root_half,root,drivers)
                hw,hz=[mean(row,root,drivers) for row in leaf_half[1:]]
                branch_margin=branch_before_defect-origin_bias-h0
                selected=[2885*value//10000-outer_bias-branch_bias-ceiling_integer(1415*h0,1000)-ceiling_integer(defect_at,2)
                          for value,defect_at in ((zl,db),(hl,du))]
                guards=[guard_base+2885*value//10000-branch_bias-ceiling_integer(1415*h0,1000)-ceiling_integer(defect_at,2)
                        for value,defect_at in ((zl,db),(hl,du))]
                endings=[4082*value//10000-closing_bias-ceiling_integer(d0,2)-ceiling_integer(3*h0,2)-ceiling_integer(7075*parent,10000)-ceiling_integer(last,2)-half
                         for value,parent,last,half in ((zl,db,dw,hw),(hl,du,dz,hz))]
                best=0
                for residual,residual_mass,third,midpoint,closing in residuals:
                    if first_margin<=third or leaf_a<=midpoint:
                        continue
                    cap=min(branch_margin-third,
                            (min(selected+guards)-midpoint)*1000//1415,
                            (min(endings)-closing)*2//3)
                    if cap<=0:
                        continue
                    for numerator in (1,2,3):
                        root_deviation=cap*numerator//4
                        last_deviations=[value-closing-ceiling_integer(3*root_deviation,2) for value in endings]
                        candidate=mass(root_deviation)*mass(leaf_a-midpoint)*prod(mass(value) for value in last_deviations)*residual_mass**6
                        best=max(best,candidate)
                if best:
                    exponent=ceiling_integer(root*root*unit+sum(drivers[i]*inverse_fixed[i][j]*drivers[j] for i in range(3) for j in range(3)),2*unit*unit)
                    density=unit-ceiling_integer(exponent,4096)
                    assert density>0
                    for unused in range(12):
                        density=density*density//unit
                    total+=density*best*(1 if iz==ih else 2)
                    if relaxed_endings and total>early_integer_cutoff:
                        break
            if relaxed_endings and total>early_integer_cutoff:
                break
        if relaxed_endings and total>early_integer_cutoff:
            break
    answer=total*mass_unit
    assert answer>target
    return answer

def weak_pair_floor():
    theta=Q(19,80000)
    lo2,hi2=Q(141421,10**5),Q(70711,50000)
    lo3,hi3=Q(34641,20000),Q(86603,50000)
    assert lo2**2<2<hi2**2 and lo3**2<3<hi3**2
    edges=[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(4,8),(4,9),(6,10),(6,11)]
    distances=tree_distances(12,edges)
    errors=[Q(0)]+[theta*d*d*hi2**(d-1) for d in range(1,7)]
    drivers=[{0:Q(1)},{0:4*hi2,1:Q(-3),2:Q(-3)}]
    drivers.extend({left:hi3,right:-hi3} for left,right in ((1,2),(4,5),(6,7),(8,9),(10,11)))
    def signed_bounds(values):
        plus=sum((value for value in values if value>=0),Q(0))
        minus=sum((-value for value in values if value<0),Q(0))
        return plus,minus
    majorant=[];variance=[]
    for i,left in enumerate(drivers):
        row=[]
        for j,right in enumerate(drivers):
            plus,minus=signed_bounds([a*b*errors[distances[u][v]] for u,a in left.items() for v,b in right.items()])
            row.append(max(plus,minus))
            if i==j:
                variance.append(1+plus)
        majorant.append(row)
    assert max(map(sum,majorant))<1
    inverse=rational_inverse([[Q(i==j)-majorant[i][j] for j in range(7)] for i in range(7)])
    assert all(value>=0 for row in inverse for value in row)
    assert prod(variance)/650**2<Q(7,44)**7
    amplitude=[Q(1131,1000),Q(5703,1000),Q(209,1000),Q(13,20),Q(13,20),Q(3,4),Q(3,4)]
    transformed=[sum(left*right for left,right in zip(row,amplitude)) for row in inverse]
    assert sum(left*right for left,right in zip(amplitude,transformed))/2<22
    caps=[Q(65,1000),Q(7,100),Q(7,100),Q(47,500),Q(47,500)]
    for vertex,cap in zip((0,1,2,4,6),caps):
        covariances=[max(signed_bounds([coefficient*3*theta*hi2**distances[vertex][point]
                                      for point,coefficient in row.items()])) for row in drivers]
        assert sum(left*right for left,right in zip(covariances,transformed))<cap
    shift=(3-2*lo2)*Q(3,100)/Q(23,100)
    chosen=(2*hi2*Q(3,25)-1+shift+caps[1]+Q(1,10))/2
    assert chosen+Q(13,20)/(2*lo3)<0
    assert 2*lo2*chosen+shift+caps[3]+Q(1,10)+Q(3,4)/lo3<0
    assert 2*lo2-2*Q(3,25)-shift-caps[0]-Q(1,10)>0
    assert 3*(3-2*lo2)*theta<(Q(23,100)/20)**2
    strip=sum((Q((-1)**k*2**(k+1),factorial(k)*(2*k+1)) for k in range(12)),Q(0))
    assert Q(797,1000)*strip>Q(19,20)
    integral=Q(0)
    for index in range(23,80):
        left,right=Q(index,100),Q(index+1,100)
        integral+=(right**7-left**7)/7*(1-22*right**2/256)**256
    volume=6*lo3*Q(3,25)**2*Q(13,20)**2*Q(3,2)**2
    answer=volume*integral*Q(19,20)**5/650
    assert answer>Q(4695,10**12)
    return answer

def local_source_controls():
    eta=Q(19,40000)
    a,b=Q(141421,100000),Q(70711,50000)
    root3=Q(173206,100000)
    ks=[Q(0),Q(1),4*b,Q(18),32*b]
    depths=[1,2,2,3,1,2]
    wlow=[Q(1),2*a/3,Q(5,6),a/2]
    whigh=[Q(1),2*b/3,Q(5,6),b/2]
    rows=[(94329,-16039,1921,315),(83603,-21780,4113,917),
          (83603,47316,-28012,917),(71567,50606,-39010,1611),
          (-94280,-32705,522,315),(-83332,-45350,-28012,917)]
    raw_r=[-Q(1,6),-a/6,b/3,Q(1,2),-Q(1,3),-a/3]
    raw_q=[Q(0),Q(0),-Q(1,2)/root3,-2*a/(4*root3),Q(0),-Q(1,2)/root3]
    r_errors=[2*ks[d] for d in (2,3,1,2,2,3)]
    q_errors=[root3*ks[d] for d in (3,4,2,3,1,2)]
    for index,((p,q,c,h),depth) in enumerate(zip(rows,depths)):
        root_upper=(whigh[depth]+ks[depth]*eta) if index<4 else -wlow[depth]
        regression=whigh[depth]+ks[depth]*eta
        assert Q(p,10**5)>root_upper
        assert Q(q,10**5)>raw_r[index]+eta*r_errors[index]+regression*2*eta
        assert Q(c,10**5)>raw_q[index]+eta*q_errors[index]+regression*eta*root3*ks[2]
        assert Q(h,10**5)>(1-wlow[depth])*Q(11,200)
    gamma0=2*eta;gamma1=root3*ks[2]*eta
    lower0=1-8*ks[2]*eta-gamma0**2
    lower1=1-6*ks[2]*eta-gamma1**2
    off=root3*(2+2*ks[3])*eta+gamma0*gamma1
    determinant=lower0*lower1-off**2
    assert Q(621,10000)**2*Q(44,7)**3<determinant
    assert 1+2*ks[2]*eta+off<Q(519,500)
    assert 1+off<Q(1033,1000)
    assert 3*(3-2*a)*eta<Q(313,20000)**2
    assert Q(173,50)**2<12
    for tolerance,mass,z in [(Q(18,1000),Q(748,1000),Q(23,20)),
                             (Q(2504,100000),Q(889,1000),Q(8,5)),
                             (Q(3365,100000),Q(967,1000),Q(43,20))]:
        assert tolerance/Q(313,20000)>=z
        lower=sum(((-Q(1,2))**degree*z**(2*degree+1)/(factorial(degree)*(2*degree+1))
                   for degree in range(12)),Q(0))
        assert Q(797,1000)*lower>mass
    return True

def outer_and_actual_threshold_controls():
    # A small independent finite rectangle sum for the adjacent private paths.
    @__import__('functools').lru_cache(None)
    def pivotal_strip(index):
        v=min(Q(index,100000),Q(5,4))
        return Q(11,8)*sum(((-Q(3,2))**degree*v**(2*degree+1)
                            /(factorial(degree)*(2*degree+1)) for degree in range(8)),Q(0))
    left_density=[Q(11,4)*(1-6*Q(i+1,100)**2/64)**64 for i in range(4,60)]
    density={};side=Q(0)
    for offset,i in enumerate(range(4,60)):
        for j in range(4,80):
            argument=2*Q(j+1,100)+Q(2829,1000)*Q(i+1,100)
            if argument not in density:
                density[argument]=Q(11,8)*(1-Q(3,2)*argument**2/64)**64
            margin_index=2828*j+1000*i-8000
            side+=Q(1,5000)*left_density[offset]*density[argument]*pivotal_strip(margin_index)
    pivotal=side*side
    assert pivotal>Q(18,1000)>Q(179,10000)
    # Distance-two moving-threshold cylinder, integrated on the fixed .005 grid.
    distance_two=Q(0)
    for cell in range(2,120):
        lower,upper=Q(cell,200),Q(cell+1,200)
        z=max(Q(0),lower-Q(1,125))
        strip=Q(69,50)*(z-z**3/2)
        density=Q(229,50)*(1-Q(33,2)*upper**2/128)**128
        distance_two+=density*strip**2/200
    assert distance_two>Q(45,1000)
    outer=(3*Q(39,100)*Q(56,1000)*Q(11,10)*Q(358,10000)
           +6*Q(39,100)*Q(165,1000)*Q(71,100)*Q(45,1000))*Q(3,100)**2/2
    assert outer==Q(671223267,10**14)>Q(67,10**7)
    # Actual-base conditional bands and the five-strip sign margins.
    rlow,rhigh=Q(8332,10000),Q(8361,10000)
    t1low,t1high=Q(9428,10000),Q(9433,10000)
    vlow=1-(2*t1high)**2/(2*(1+rlow))-(t1high-t1low)**2/(2*(1-rhigh))
    vhigh=1-(2*t1low)**2/(2*(1+rhigh))
    assert vlow>Q(29,1000) and vhigh<Q(32,1000)
    assert (Q(135,100000)+Q(52,100)*(Q(24,100000)+Q(428,100000)))/Q(29,1000)<Q(13,100)
    assert Q(3225,10000)*Q(1,5)-Q(1125,10000)*Q(3,100)-Q(6,100)>0
    assert Q(445,100)*(1-Q(45,29)/128)**128>Q(93,100)
    mass=Q(93,1000)*Q(1365,10000)**2*Q(78,100)**2*Q(999,1000)
    assert mass>Q(1,1000)
    coefficient=Q(39,100)*Q(163,1000)*Q(71,100)*Q(1,1000)/2
    assert coefficient==Q(451347,20000000000)
    return {"adjacent_zero_event":str(pivotal),"outer_distance_two_event":str(distance_two),
            "outer_saving":str(outer),"actual_threshold_unit":str(coefficient)}

"""Bind the finite inequalities to the single declared eventual theorem."""

def population_arithmetic():
    eta=Q(19,40000);K=Q(28238380378679347,10**18);A=eta*K
    p=Q(247011,10**12);u=Q(301,10**10);w=Q(4695,10**12);beta=Q(3,8)
    v11=((Q(31,20)+beta/2)*A-92*p-4*u-Q(232,3)*w)/26
    v9=p-u-v11
    z=13*v9+17*v11+12*u+Q(32,3)*w
    assert v9>0 and v11>0
    assert Q(7,32)*A<z<Q(11,32)*A
    c=A/40;g=Q(2,3)*(z-c-beta*A/4);B=c+g/2
    assert g>0
    coefficients={(0,0):Q(0),(0,1):A*(Q(1,2)+beta/4)+B,
                  (0,2):A/2-g,(0,3):3*A/4-3*B,(1,0):9*A/10+z,
                  (1,1):9*A/10-B,(1,2):17*A/20-2*B,
                  (2,0):7*A/10,(2,1):19*A/20-B,(3,0):21*A/20}
    overlap=[12*v9+16*v11+12*u-3*A/16,30*v9+40*v11+24*u-9*A/32,4*v9+6*v11-A/16]
    assert all(value>0 for value in overlap)
    ranges={}
    for name,q,d in [("inner",Q(283,10**9),Q(0)),
                      ("annulus",Q(273249,10**12),Q(451347,50000000000000))]:
        slacks={}
        for (a,b),payment in coefficients.items():
            required=q*[0,18,12,6][a]+v9*[0,35,22,9][a+b]+v11*[0,45,28,11][a+b]
            required+=w*[0,Q(88,3),Q(56,3),8][a+b]+u*([0,36,24,12][a+b]+8*(3-a-b)*(a>0))+14*d*a
            slacks[a,b]=payment-required
        assert all(value>=0 for value in slacks.values())
        assert {key for key,value in slacks.items() if value==0}=={(0,0),(0,1),(0,2)}
        saving=6*q+9*v9+11*v11+12*u+8*w+6*d
        assert saving>Q(4052,10**9)
        ranges[name]=saving
    assert ranges["inner"]==Q(165375882580092155427,41600000000000000000000000)+Q(35,13)*u
    assert ranges["annulus"]==Q(165195157204092155427,41600000000000000000000000)+Q(35,13)*u
    gamma=Q(26129513284633533,200000000000000000)-Q(4052,10**9)
    exact_coefficient=1+1/(2*gamma)
    declared=Q('4.82720873828033295823867219799062')
    assert Q((exact_coefficient*10**32).__floor__(),10**32)==declared<exact_coefficient
    assert Q(67,10**7)>ranges["inner"]
    return {"K":str(K),"inner_credit":str(ranges['inner']),"annular_credit":str(ranges['annulus']),
            "overlap_slacks":[str(value) for value in overlap],
            "minimum_slack":str(ranges['annulus']-Q(4052,10**9)),
            "graph_saving":"4052/1000000000","gamma":str(gamma),
            "coefficient":str(declared),"rounding_gap":str(exact_coefficient-declared)}

def certificate_main():
    import json,time,sys
    start=time.monotonic()
    def progress(label):
        print(label,round(time.monotonic()-start,3),file=sys.stderr,flush=True)
    path=path_series_enclosure();progress('path certificate')
    local_source_controls()
    q_inner=seven_key_floor(Fraction(1,50))
    q_annulus=seven_key_floor(Fraction(3,100))
    assert q_inner>Q(283,10**9) and q_annulus>Q(273249,10**12)
    nine=minimal_nine_key_floor();relaxed=minimal_nine_key_floor(True);weak=weak_pair_floor();progress('ordering probabilities')
    assert relaxed>Q(301,10**10) and weak>Q(4695,10**12)
    thresholds=outer_and_actual_threshold_controls();progress('threshold bounds')
    curvature=curvature_replay();progress('curvature integrals')
    population=population_arithmetic()
    print(json.dumps({"passed":True,"scope":"exact finite inequalities; general proofs in supplied manuscripts",
          "path":path,"q_inner":str(q_inner),"q_annulus":str(q_annulus),"minimal9":str(nine),
          "relaxed_minimal9":str(relaxed),"weak_pair_floor_1e18":(weak*10**18).__floor__(),
          "outer_pivotal_floor_1e12":(Q(thresholds['adjacent_zero_event'])*10**12).__floor__(),
          "outer_saving":thresholds['outer_saving'],"actual_threshold_unit":thresholds['actual_threshold_unit'],
          "curvature_floor_1e18":{name:(Q(value)*10**18).__floor__() for name,value in curvature.items()},
          "population":population,"wall_seconds":round(time.monotonic()-start,6)},sort_keys=True))

if __name__=='__main__':
    certificate_main()
