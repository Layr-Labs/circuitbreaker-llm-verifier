"""Standalone outward-rounded quadrature for the crossing-path bound U and crossing density W1."""
from decimal import (
    Context as DecimalContext,
    Decimal,
    ROUND_CEILING,
    ROUND_FLOOR,
    ROUND_HALF_EVEN,
    localcontext,
)
from fractions import Fraction
from math import cos, factorial, pi
from typing import NamedTuple


Q_LOW = DecimalContext(prec=44, rounding=ROUND_FLOOR)
Q_HIGH = DecimalContext(prec=44, rounding=ROUND_CEILING)
Q_MID = DecimalContext(prec=44, rounding=ROUND_HALF_EVEN)
GRID_BITS = 124
GRID = 1 << GRID_BITS
NODE_COUNT = 58
J_ODD = 37
WINDOW = Fraction(13, 5)
PATH_UPPER = Fraction(1306475664232, 10**13)
CROSS_UPPER = Fraction("0.16226017196")


class Enclosure(NamedTuple):
    low: Decimal
    high: Decimal

    @staticmethod
    def exact(value):
        if isinstance(value, Enclosure):
            return value
        if isinstance(value, Fraction):
            n, d = Decimal(value.numerator), Decimal(value.denominator)
            return Enclosure(Q_LOW.divide(n, d), Q_HIGH.divide(n, d))
        point = Decimal(value)
        return Enclosure(point, point)

    def __add__(self, other):
        b = self.exact(other)
        return Enclosure(Q_LOW.add(self.low, b.low), Q_HIGH.add(self.high, b.high))

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
            return Enclosure(Q_LOW.multiply(self.low, b.low), Q_HIGH.multiply(self.high, b.high))
        endpoints = tuple((x, y) for x in self for y in b)
        return Enclosure(
            min(Q_LOW.multiply(x, y) for x, y in endpoints),
            max(Q_HIGH.multiply(x, y) for x, y in endpoints),
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.exact(other)
        assert b.high < 0 or b.low > 0
        return self * Enclosure(Q_LOW.divide(1, b.high), Q_HIGH.divide(1, b.low))

    def __rtruediv__(self, other):
        return self.exact(other) / self

    def square(self):
        magnitudes = (self.low.copy_abs(), self.high.copy_abs())
        minimum = Decimal(0) if self.low <= 0 <= self.high else min(magnitudes)
        maximum = max(magnitudes)
        return Enclosure(Q_LOW.multiply(minimum, minimum), Q_HIGH.multiply(maximum, maximum))

    def sqrt(self):
        assert self.low > 0
        lo = Q_LOW.next_minus(Q_MID.sqrt(self.low))
        hi = Q_HIGH.next_plus(Q_MID.sqrt(self.high))
        assert lo > 0 and Q_HIGH.multiply(lo, lo) <= self.low
        assert Q_LOW.multiply(hi, hi) >= self.high
        return Enclosure(lo, hi)

    def exp(self):
        if self.low == self.high:
            return _exp_point(self.low)
        return Enclosure(_exp_point(self.low).low, _exp_point(self.high).high)

    def fixed(self):
        nl, dl = self.low.as_integer_ratio()
        nh, dh = self.high.as_integer_ratio()
        return (nl * GRID // dl, -((-nh * GRID) // dh))


def _exp_point(x):
    """Elementary range-reduced Taylor enclosure of exp(x) using only Q_LOW/Q_HIGH +, -, *, /."""
    if x == 0:
        return Enclosure.exact(1)
    neg = x < 0
    t = x.copy_abs()
    k = 0
    scale = 1
    while t > scale:
        scale <<= 1
        k += 1
    k += 10
    scale <<= 10
    y = Enclosure(Q_LOW.divide(t, scale), Q_HIGH.divide(t, scale))  # 0 < y <= 2^-10
    term = Enclosure.exact(1)
    total = term
    for j in range(1, 14):
        term = (term * y) / j
        total = total + term
    # Remainder sum_{j>=12} y^j/j! <= 2*y^12/12! = (y^11/11!)*y/6 for 0 <= y <= 1/2.
    rem = (term * y) / 7
    res = Enclosure(total.low, Q_HIGH.add(total.high, rem.high))
    for _ in range(k):
        res = res.square()
    return (Enclosure.exact(1) / res) if neg else res


def pi_enclosure():
    def atan(inv, terms):
        # Alternating series with even term count: partial sum is lower, + next term is upper.
        lower = sum(Fraction((-1) ** k, (2 * k + 1) * inv ** (2 * k + 1)) for k in range(terms))
        upper = lower + Fraction(1, (2 * terms + 1) * inv ** (2 * terms + 1))
        return Enclosure(Enclosure.exact(lower).low, Enclosure.exact(upper).high)

    return 16 * atan(5, 38) - 4 * atan(239, 16)


def crossing_density_enclosure(pi_range):
    """Certify W_1 = (3/(2*pi)) * arcsin(1/3) via the positive Taylor series of arcsin(1/3)."""
    total = Fraction(0)
    term = Fraction(1, 3)
    for k in range(36):
        total += term
        term = term * (2 * k + 1) * (2 * k + 1) / ((2 * k + 2) * (2 * k + 3) * 9)
    # Consecutive ratio is < 1/9, so geometric tail is bounded by (9/8) * term.
    arcsin_enc = Enclosure(
        Enclosure.exact(total).low,
        Enclosure.exact(total + Fraction(9, 8) * term).high,
    )
    return (3 * arcsin_enc) / (2 * pi_range)


def legendre_value_at_fraction(a, b):
    """Return the last two scaled Legendre values k! b^k P_k(a/b) for k = q-1, q."""
    older, newer = 1, a
    bb = b * b
    for k in range(2, NODE_COUNT + 1):
        older, newer = newer, (2 * k - 1) * a * newer - ((k - 1) ** 2) * bb * older
    return older, newer


def approximate_root(index):
    with localcontext(Q_MID):
        q = NODE_COUNT
        x = Decimal(str(cos(pi * (index + 0.75) / (q + 0.5))))
        for _ in range(7):
            p, next_p = Decimal(1), x
            for k in range(2, q + 1):
                p, next_p = next_p, ((2 * k - 1) * x * next_p - (k - 1) * p) / k
            x -= next_p * (x * x - 1) / (q * (x * next_p - p))
        return x


def gauss_rule():
    rule = []
    last_left = Decimal(1)
    q = NODE_COUNT
    radius = Decimal("1e-36")
    for index in range(q):
        center = approximate_root(index)
        left = Q_LOW.subtract(center, radius)
        right = Q_HIGH.add(center, radius)
        assert -1 < left < right < last_left
        signs = [legendre_value_at_fraction(*bound.as_integer_ratio())[1] for bound in (left, right)]
        assert signs[0] * signs[1] < 0
        last_left = left
        a, b = center.as_integer_ratio()
        previous, current = legendre_value_at_fraction(a, b)
        derivative = Enclosure.exact(
            Fraction(q * b * (a * current - q * b * b * previous), factorial(q) * b**q * (a * a - b * b))
        )
        width = max(Q_HIGH.subtract(center, left), Q_HIGH.subtract(right, center))
        error = Q_HIGH.multiply(width, q**4)
        derivative += Enclosure(error.copy_negate(), error)
        root = Enclosure(left, right)
        weight = WINDOW / ((1 - root.square()) * derivative.square())
        rule.append(((root + 1) * WINDOW / 2, weight))
    return rule


def positive_arms(length):
    return 10 ** (length // 2) * 2 ** (length % 2)


def tensor_error(dimension, at_zero):
    q = NODE_COUNT
    if dimension == 2:
        precision_moment = 2 * 9**q
    elif dimension == 3:
        precision_moment = 2 * 12**q + 33**q
    else:
        precision_moment = 2 * 12**q + 2 * 36**q + (dimension - 4) * 39**q
    weight = Fraction(
        WINDOW ** (dimension + 2 * q) * factorial(q) ** 3 * precision_moment,
        (2 * q + 1) * 2**q * factorial(2 * q) ** 2 * (4 * WINDOW) ** max(0, dimension - 9),
    )
    return (at_zero * weight).high


def reference_enclosures():
    if not __debug__:
        raise RuntimeError("Run the certificate with assertions enabled.")
    pi_range = pi_enclosure()
    w1_range = crossing_density_enclosure(pi_range)
    assert Fraction(w1_range.high) < CROSS_UPPER

    root2 = Enclosure.exact(2).sqrt()
    density_1 = (6 / pi_range).sqrt()
    density_2 = 3 / (2 * pi_range)
    rule = gauss_rule()
    nodes = [x for x, _ in rule]
    weights = [w for _, w in rule]
    x2 = [x.square() for x in nodes]
    assert sum(w * (-s * Fraction(39, 2)).exp() for w, s in zip(weights, x2)).high < Decimal("0.25")
    z_unary = [w * density_1 * (-6 * s).exp() for w, s in zip(weights, x2)]
    assert sum(z_unary).high < 1

    q = NODE_COUNT
    plain = [[None] * q for _ in range(q)]
    rooted = [[None] * q for _ in range(q)]
    for i in range(q):
        for k in range(i + 1):
            xy = nodes[i] * nodes[k]
            plain[i][k] = plain[k][i] = (-6 * xy).exp()
            rooted[i][k] = rooted[k][i] = (-6 * root2 * xy).exp()
    unary_a = [(-s * Fraction(3, 2)).exp() for s in x2]
    unary_c = [(-12 * s).exp() for s in x2]
    pair = [w * (-s * Fraction(9, 2)).exp() for w, s in zip(weights, x2)]

    def tabulate(make):
        lower = []
        widest = 0
        for i in range(q):
            row = []
            for k in range(q):
                a, b = make(i, k).fixed()
                assert 0 <= a <= b
                row.append(a)
                if b - a > widest:
                    widest = b - a
            lower.append(row)
        return lower, widest

    A, debit_a = tabulate(lambda i, k: unary_a[i] * rooted[i][k])
    B, debit_b = tabulate(lambda i, k: plain[i][k])
    C, debit_c = tabulate(lambda i, k: unary_c[i] * z_unary[k] * rooted[i][k].square())
    state = []
    debit = 0
    for i in range(q):
        row = []
        for k in range(q):
            a, b = (density_2 * pair[i] * pair[k] * rooted[i][k]).fixed()
            assert 0 <= a <= b
            row.append(a)
            debit += b - a
        state.append(row)
    assert sum(map(sum, state)) + debit <= GRID
    step_debit = debit_a + debit_b + q * debit_c + 2 * q * q

    total = Enclosure.exact(0)
    pair_tail = density_2 * Enclosure.exact(Fraction(-9 * WINDOW**2, 2)).exp() / (54 * root2 * WINDOW**2)
    zero = density_2
    even_lower = None
    for j in range(1, J_ODD + 1):
        dim = j + 1
        mass = sum(map(sum, state))
        low = Q_LOW.divide(mass, GRID)
        high = Q_HIGH.divide(mass + debit, GRID)
        quad = tensor_error(dim, zero)
        arms = sum(
            Fraction(1, positive_arms(min(i, dim - 2)) * positive_arms(dim - 2 - min(i, dim - 2)))
            for i in range(dim)
        )
        error_high = Q_HIGH.add(quad, (pair_tail * arms).high)
        orthant = Enclosure(Q_LOW.subtract(low, quad), Q_HIGH.add(high, error_high))
        total += orthant * (3 * 2 ** (j - 1) * ((-1) ** (j + 1)))
        if j == J_ODD - 1:
            even_lower = total.low
        if j == J_ODD:
            break

        next_state = []
        for middle in range(q):
            active = [(i, (state[i][middle] * A[i][middle]) >> GRID_BITS) for i in range(q) if state[i][middle]]
            active = [(i, a) for i, a in active if a]
            row = []
            for last, kernel in enumerate(B):
                contraction = sum(a * kernel[i] for i, a in active)
                row.append((C[middle][last] * contraction) >> (2 * GRID_BITS))
            next_state.append(row)
        state = next_state
        debit += step_debit
        zero *= density_1

    assert Fraction(total.high) < PATH_UPPER
    return {
        "U": PATH_UPPER,
        "W1": CROSS_UPPER,
        "gamma_upper": str(total.high),
        "gamma_even_lower": str(even_lower),
        "w1_upper": str(w1_range.high),
        "w1_lower": str(w1_range.low),
        "accumulated_debit_grid_units": debit,
        "per_step_debit_grid_units": step_debit,
    }


def reference_upper_bounds():
    """Replay the certificate and return exact rational upper bounds (U, W1)."""
    res = reference_enclosures()
    return res["U"], res["W1"]


if __name__ == "__main__":
    import json
    import time

    t0 = time.time()
    bounds = reference_enclosures()
    dt = time.time() - t0
    print(
        json.dumps(
            {
                "U": str(bounds["U"]),
                "W1": str(bounds["W1"]),
                "gamma_upper": bounds["gamma_upper"],
                "gamma_even_lower": bounds["gamma_even_lower"],
                "w1_upper": bounds["w1_upper"],
                "elapsed_s": round(dt, 4),
            },
            indent=2,
        )
    )
