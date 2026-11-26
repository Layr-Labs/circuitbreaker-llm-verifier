"""Standalone directed quadrature for the 25-term crossing-path bound."""
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction
from functools import cache
from math import cos, factorial, pi, sumprod
from typing import NamedTuple


Q_LOW = Context(prec=100, rounding=ROUND_FLOOR)
Q_HIGH = Context(prec=100, rounding=ROUND_CEILING)
Q_MID = Context(prec=100, rounding=ROUND_HALF_EVEN)
GRID = 1 << 224
NODE_COUNT = 48
WINDOW = Fraction(5, 2)


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


@cache
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
            active = tuple((state[i][middle]*A[i][middle]) >> 224 for i in range(q))
            row = []
            for last, kernel in enumerate(B):
                contraction = sumprod(active, kernel)
                row.append((C[middle][last]*contraction) >> 448)
            next_state.append(row)
        state = next_state
        debit += step_debit
        zero *= density_1

    assert total.high < Decimal('0.13064756648')
    return {"gamma_upper": str(total.high), "gamma_even_lower": str(even_lower),
            "initial_and_accumulated_debit_grid_units": debit,
            "per_step_debit_grid_units": step_debit}


if __name__ == "__main__":
    import json
    print(json.dumps(path_series_enclosure(), sort_keys=True))
