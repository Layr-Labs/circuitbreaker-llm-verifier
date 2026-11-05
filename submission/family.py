"""Uniform fixed polylog-entropy family for proof.md.

The source defines a total family with all small/failed shapes zero.
No hardware cutoff is part of the definition. The asymptotic polynomial degree
is fixed and enormous. Actual nontrivial instances are not practical tests.
"""
from functools import lru_cache
from fractions import Fraction
from math import comb
import json
import sys

LOG_ENTROPY_POWER = 24576
LONG_ORDER = 4096
C_LINEAR = 1 << 4110
RESILIENCE_DENOMINATOR = 4096
RESILIENCE_DEGREE = 4 * RESILIENCE_DENOMINATOR


def floor_root(n, degree):
    if n < 0 or degree < 1:
        raise ValueError('nonnegative integer and positive degree required')
    if n < 2 or degree == 1:
        return n
    if degree >= n.bit_length():
        return 1
    lo, hi = 0, 1 << ((n.bit_length() + degree - 1) // degree)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** degree <= n:
            lo = mid
        else:
            hi = mid
    return lo


def ceil_root(n, degree):
    r = floor_root(n, degree)
    return r if r ** degree == n else r + 1


def floor_root_power_two(n, exponent):
    """floor(n^(1/2^exponent)), without materializing an unnecessary degree."""
    if n < 2:
        return n
    if exponent >= n.bit_length().bit_length():
        return 1
    return floor_root(n, 1 << exponent)


def index_subsets(count, size):
    """All fixed-size subsets in increasing bit-mask order, using integers."""
    if size < 0 or size > count:
        return
    if size == 0:
        yield ()
        return
    mask, limit = (1 << size) - 1, 1 << count
    while mask < limit:
        yield tuple(j for j in range(count) if (mask >> j) & 1)
        lowest = mask & -mask
        following = mask + lowest
        mask = following | (((following ^ mask) >> 2) // lowest)


def poly_mod(x, modulus):
    while x and x.bit_length() >= modulus.bit_length():
        x ^= modulus << (x.bit_length() - modulus.bit_length())
    return x


def poly_gcd(x, y):
    while y:
        x, y = y, poly_mod(x, y)
    return x


def prime_divisors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


class PolynomialField:
    def __init__(self, modulus):
        if modulus < 2:
            raise ValueError('positive field degree required')
        self.modulus = modulus
        self.degree = modulus.bit_length() - 1

    def mul(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a >> self.degree:
                a ^= self.modulus
        return out

    def power(self, a, e):
        out = 1
        while e:
            if e & 1:
                out = self.mul(out, a)
            e >>= 1
            if e:
                a = self.mul(a, a)
        return out


def irreducible(modulus):
    if modulus < 2:
        return False
    d = modulus.bit_length() - 1
    field = PolynomialField(modulus)
    x = poly_mod(2, modulus)
    value = x
    divisors = {d // p for p in prime_divisors(d)}
    for j in range(1, d + 1):
        value = field.mul(value, value)
        if j in divisors and poly_gcd(value ^ x, modulus) != 1:
            return False
    return value == x


@lru_cache(None)
def small_field(degree):
    for low in range(1 << degree):
        modulus = (1 << degree) | low
        if irreducible(modulus):
            return PolynomialField(modulus)
    raise AssertionError('irreducible polynomial existence')


class TowerField:
    """Explicit degree-power-of-two tower; no irreducible search."""
    def __init__(self, degree, base=None):
        if degree < 1 or degree & (degree - 1):
            raise ValueError('tower degree must be a power of two')
        self.degree = degree
        if degree == 1:
            self.trace_one = 1
            self.base = None
        else:
            self.base = base if base is not None else tower_field(degree // 2)
            self.trace_one = self.base.trace_one << (degree // 2)

    def mul(self, a, b):
        # Explicit evaluation stack avoids a host-language recursion limit.
        pending, values = [('mul', self, a, b)], []
        while pending:
            task = pending.pop()
            if task[0] == 'mul':
                _, field, x, y = task
                if field.degree == 1:
                    values.append(x & y)
                    continue
                half = field.degree // 2
                mask = (1 << half) - 1
                x0, x1 = x & mask, x >> half
                y0, y1 = y & mask, y >> half
                pending.extend([('combine', field),
                                ('mul', field.base, x0 ^ x1, y0 ^ y1),
                                ('mul', field.base, x1, y1),
                                ('mul', field.base, x0, y0)])
            elif task[0] == 'combine':
                field = task[1]
                p2, p1, p0 = values.pop(), values.pop(), values.pop()
                pending.extend([('finish', p0, p2, field.degree // 2),
                                ('mul', field.base, field.base.trace_one, p1)])
            else:
                _, p0, p2, half = task
                values.append((p0 ^ values.pop()) | ((p2 ^ p0) << half))
        assert len(values) == 1
        return values[0]

    def power(self, a, e):
        out = 1
        while e:
            if e & 1:
                out = self.mul(out, a)
            e >>= 1
            if e:
                a = self.mul(a, a)
        return out

    def trace(self, a):
        out, value = 0, a
        for _ in range(self.degree):
            out ^= value
            value = self.mul(value, value)
        return out


@lru_cache(None)
def tower_field(degree):
    if degree < 1 or degree & (degree - 1):
        raise ValueError('tower degree must be a power of two')
    result, width = TowerField(1), 1
    while width < degree:
        width *= 2
        result = TowerField(width, result)
    return result


class LinearCode:
    def __init__(self, n, output, error_exp):
        self.n = n
        target = n.bit_length() + 2 * (8 * output).bit_length() + 2 * error_exp + 10
        self.w = 1 << (target - 1).bit_length()
        self.ell = 2 * self.w
        self.field = tower_field(self.w)

    def bit(self, x, index):
        mask = (1 << self.w) - 1
        a, b = index & mask, (index >> self.w) & mask
        value = 0
        for j in range((self.n + self.w - 1) // self.w - 1, -1, -1):
            value = self.field.mul(value, a) ^ ((x >> (j * self.w)) & mask)
        return (value & b).bit_count() & 1


def polynomial_design(ell, count, intersection):
    z = (ell - 1).bit_length()
    q = 1 << z
    if count > q ** (intersection + 1):
        raise ValueError('insufficient polynomial graphs')
    field = small_field(z)
    out = []
    for i in range(count):
        coefficients = []
        value = i
        for _ in range(intersection + 1):
            coefficients.append(value & (q - 1))
            value >>= z
        graph = []
        for x in range(ell):
            y = 0
            for c in reversed(coefficients):
                y = field.mul(y, x) ^ c
            graph.append(x * q + y)
        out.append(tuple(graph))
    return ell * q, out


def ceil_log2(n):
    if n < 1:
        raise ValueError('positive size required')
    return (n - 1).bit_length()


def capped_power(base, exponent, cap):
    """min(base**exponent, cap+1), with O(log cap) bits per temporary."""
    if base < 0 or exponent < 0 or cap < 0:
        raise ValueError('nonnegative arguments required')
    result = 1
    while exponent:
        if exponent & 1:
            result = min(cap + 1, result * base)
        exponent >>= 1
        if exponent:
            base = min(cap + 1, base * base)
    return min(result, cap + 1)


def primitive(field):
    order = (1 << field.degree) - 1
    factors = prime_divisors(order)
    candidate = 1
    while candidate <= order:
        if all(field.power(candidate, order // p) != 1 for p in factors):
            return candidate
        candidate += 1
    raise AssertionError('finite multiplicative group is cyclic')


def coefficient_list(value, width, count):
    mask = (1 << width) - 1
    return [(value >> (width * i)) & mask for i in range(count)]


def pack_coefficients(values, width):
    result = 0
    for i, value in enumerate(values):
        result |= value << (width * i)
    return result


def trim_coefficients(values):
    values = list(values)
    while values and values[-1] == 0:
        values.pop()
    return values


def field_poly_remainder(values, modulus, base):
    values = trim_coefficients(values)
    modulus = trim_coefficients(modulus)
    if not modulus:
        raise ValueError('zero modulus')
    inverse = base.power(modulus[-1], (1 << base.degree) - 2)
    while len(values) >= len(modulus):
        scale = base.mul(values[-1], inverse)
        shift = len(values) - len(modulus)
        for i, coefficient in enumerate(modulus):
            values[shift + i] ^= base.mul(scale, coefficient)
        values = trim_coefficients(values)
    return values


def field_poly_gcd(first, second, base):
    first, second = trim_coefficients(first), trim_coefficients(second)
    while second:
        first, second = second, field_poly_remainder(first, second, base)
    return first


class ExtensionField:
    """A monic polynomial extension of a binary field, in fixed power basis."""
    def __init__(self, base, lower_coefficients):
        self.base = base
        self.width = base.degree
        self.lower = tuple(lower_coefficients)
        self.dimension = len(self.lower)
        self.degree = self.width * self.dimension

    def mul(self, first, second):
        size, base = self.dimension, self.base
        x = coefficient_list(first, self.width, size)
        y = coefficient_list(second, self.width, size)
        product = [0] * (2 * size - 1)
        for i in range(size):
            for j in range(size):
                product[i + j] ^= base.mul(x[i], y[j])
        for i in range(2 * size - 2, size - 1, -1):
            top = product[i]
            if top:
                for j, coefficient in enumerate(self.lower):
                    product[i - size + j] ^= base.mul(top, coefficient)
        return pack_coefficients(product[:size], self.width)

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            exponent >>= 1
            if exponent:
                value = self.mul(value, value)
        return result


def extension_irreducible(base, lower):
    dimension = len(lower)
    if dimension < 2:
        return dimension == 1
    extension = ExtensionField(base, lower)
    x = 1 << base.degree
    value = x
    checkpoints = {dimension // p for p in prime_divisors(dimension)}
    modulus = list(lower) + [1]
    for j in range(1, dimension + 1):
        for _ in range(base.degree):
            value = extension.mul(value, value)
        if j in checkpoints:
            difference = coefficient_list(value ^ x, base.degree, dimension)
            if len(field_poly_gcd(modulus, difference, base)) != 1:
                return False
    return value == x


@lru_cache(None)
def extension_field(base_degree, dimension):
    base = small_field(base_degree)
    limit = 1 << (base_degree * dimension)
    mask = (1 << base_degree) - 1
    candidate = 1
    while candidate < limit:
        if candidate & mask:
            lower = coefficient_list(candidate, base_degree, dimension)
            if extension_irreducible(base, lower):
                return ExtensionField(base, lower)
        candidate += 1
    raise AssertionError('irreducible extension exists')


def compositions(total, width):
    """Lexicographic nonnegative compositions, with an explicit stack."""
    pending = [(0, total, ())]
    while pending:
        position, remaining, prefix = pending.pop()
        if position + 1 == width:
            yield prefix + (remaining,)
            continue
        for value in range(remaining, -1, -1):
            pending.append((position + 1, remaining - value, prefix + (value,)))


def monomial_prefix(dimension, degree_bound, count):
    result = []
    degree = 0
    while degree < degree_bound and len(result) < count:
        for powers in compositions(degree, dimension):
            result.append(powers)
            if len(result) == count:
                break
        degree += 1
    if len(result) != count:
        raise ValueError('insufficient multivariate coefficient capacity')
    return result


class ShortAffine:
    """Folded affine condenser followed by the supplied basic SU construction."""
    def __init__(self, n, entropy, output):
        self.n, self.entropy, self.output = n, entropy, output
        self.w1 = ceil_log2(n ** 16 * entropy)
        self.w2 = ceil_log2(n ** 128)
        self.input_small = entropy * self.w1
        self.su_degree = floor_root(entropy, 8)
        self.seed_length = self.w1 + 17 * self.w2
        self._ready = False

    def capacity_ok(self):
        if self.su_degree < 1:
            return False
        return comb(self.su_degree + 15, 16) * self.w2 >= self.input_small

    def prepare(self):
        if self._ready:
            return
        self.fold_field = small_field(self.w1)
        self.fold_generator = primitive(self.fold_field)
        self.su_field = small_field(self.w2)
        self.orbit = extension_field(self.w2, 16)
        self.orbit_generator = primitive(self.orbit)
        coefficient_count = (self.input_small + self.w2 - 1) // self.w2
        self.monomials = monomial_prefix(16, self.su_degree, coefficient_count)
        self._ready = True

    def fold(self, x, seed):
        field = self.fold_field
        coefficients = coefficient_list(x, self.w1, (self.n + self.w1 - 1) // self.w1)
        point, output = seed, 0
        for j in range(self.entropy):
            value = 0
            for coefficient in reversed(coefficients):
                value = field.mul(value, point) ^ coefficient
            output |= value << (j * self.w1)
            point = field.mul(point, self.fold_generator)
        return output

    def su_value(self, coefficients, vector):
        field = self.su_field
        coordinates = coefficient_list(vector, self.w2, 16)
        output = 0
        for coefficient, exponents in zip(coefficients, self.monomials):
            if coefficient:
                term = coefficient
                for value, exponent in zip(coordinates, exponents):
                    term = field.mul(term, field.power(value, exponent))
                output ^= term
        return output

    def __call__(self, x, seed):
        self.prepare()
        fold_seed = seed & ((1 << self.w1) - 1)
        vector = (seed >> self.w1) & ((1 << (16 * self.w2)) - 1)
        inner = seed >> (self.w1 + 16 * self.w2)
        condensed = self.fold(x, fold_seed)
        coefficients = coefficient_list(condensed, self.w2, len(self.monomials))
        output = 0
        for i in range(self.output):
            vector = self.orbit.mul(vector, self.orbit_generator)
            value = self.su_value(coefficients, vector)
            output |= ((value & inner).bit_count() & 1) << i
        return output


class LongLinear:
    def __init__(self, n, output, error_exponent):
        self.n, self.output, self.error_exponent = n, output, error_exponent
        self.code = LinearCode(n, output, error_exponent)
        self.seed_length = self.code.ell * (1 << ceil_log2(self.code.ell))
        self.entropy_needed = (2 * self.seed_length + (1 << (LONG_ORDER + 4)) * output
                               + 8 * (error_exponent + ceil_log2(output + 2) + 10))
        self._sets = None

    def capacity_ok(self):
        q = 1 << ceil_log2(self.code.ell)
        return LONG_ORDER < self.code.ell and self.output <= q ** (LONG_ORDER + 1)

    def __call__(self, x, seed):
        if self._sets is None:
            dimension, self._sets = polynomial_design(self.code.ell, self.output, LONG_ORDER)
            assert dimension == self.seed_length
        output = 0
        for i, positions in enumerate(self._sets):
            index = 0
            for j, position in enumerate(positions):
                index |= ((seed >> position) & 1) << j
            output |= self.code.bit(x, index) << i
        return output

class ResilientShape:
    @classmethod
    def build(cls, size):
        if size < 2:
            return None
        result = cls()
        result.size = size
        log_size = ceil_log2(size + 2)
        quotient = size // log_size
        if quotient < 2:
            return None
        result.z = quotient.bit_length() - 1
        result.q = 1 << result.z
        result.order = RESILIENCE_DEGREE
        quotient = result.q // (8 * result.order * (result.z + 1))
        if quotient < 1:
            return None
        result.w = quotient.bit_length() - 1
        if result.w < result.order + 1 or result.w > result.q:
            return None
        approximation = sum((Fraction(1, j << j) for j in range(1, 4 * result.z + 17)), Fraction(0))
        scaled = (1 << result.w) * (result.order * result.z + 1) * approximation
        result.blocks = scaled.numerator // scaled.denominator
        if not 1 <= result.blocks <= result.q or result.w * result.q > size:
            return None
        result.rows = 1 << (result.z * result.order)
        result.moment_order = (1 << 20) * (result.order * result.z + 1) + 4
        return result

    def evaluate(self, word):
        field = small_field(self.z)
        mask = self.q - 1
        row_index = 0
        while row_index < self.rows:
            coefficients = coefficient_list(row_index, self.z, self.order)
            row_value = 0
            offset = 0
            while offset < self.blocks:
                tribe = 1
                x = 0
                while x < self.w:
                    value = 0
                    for coefficient in reversed(coefficients):
                        value = field.mul(value, x) ^ coefficient
                    value = field.mul(value, x) ^ offset
                    if not ((word >> (x * self.q + value)) & 1):
                        tribe = 0
                        break
                    x += 1
                if tribe:
                    row_value = 1
                    break
                offset += 1
            if not row_value:
                return 0
            row_index += 1
        return 1


class FamilyShape:
    @classmethod
    def build(cls, n):
        if n < 2:
            return None
        log_n = ceil_log2(n + 2)
        entropy = capped_power(log_n, LOG_ENTROPY_POWER, n)
        if entropy > n:
            return None
        result = cls()
        result.n, result.log_n, result.k = n, log_n, entropy
        result.m = floor_root(entropy, 4)
        root_h = ceil_root(entropy, 4096)
        result.h = 1 << ceil_log2(root_h)
        result.ell = ceil_root(entropy, 1024)
        if result.h < 2 or result.ell < 1:
            return None
        result.l = result.h.bit_length() - 1
        result.qlen = (16 * result.h ** 2 + 16 * C_LINEAR) * result.ell
        result.d = (1 << 64) * (result.h ** 3 * result.ell + log_n) ** 2
        result.error_exponent = ceil_root(result.ell, 8)
        result.first_error_exponent = (1 << 20) * result.h ** 3 * result.ell + result.error_exponent
        result.short = ShortAffine(n, entropy, result.m)
        if not result.short.capacity_ok():
            return None
        result.s = result.short.seed_length
        result.b = (result.s + result.l - 1) // result.l
        result.row_count = 1 << result.s
        result.budget = ((1 << 20) * C_LINEAR * result.b * result.h
                         * (result.d + result.qlen + result.h * result.ell))
        if (result.m < 4 * result.budget
                or entropy < 4 * result.h * result.m + 4 * result.budget
                or result.h ** 2 < C_LINEAR
                or result.d < 4 * result.h * (result.qlen + (result.h + 1) * result.ell)
                or result.m < 4 * result.d):
            return None
        su_advice = 2 * result.m * (255 * (result.short.su_degree - 1) + 1) * result.short.w2
        if entropy < su_advice + ceil_log2(8 * result.m * n ** 8):
            return None
        if result.row_count > n ** 4096:
            return None
        result.x_refresh = LongLinear(n, result.ell, result.first_error_exponent)
        result.y_refresh = LongLinear(result.m, result.qlen, result.error_exponent)
        result.y_next = LongLinear(result.m, result.d, result.error_exponent)
        result.x_alt = LongLinear(n, result.ell, result.error_exponent)
        result.q_alt = LongLinear(result.qlen, result.ell, result.error_exponent)
        contracts = [(result.x_refresh, result.d, C_LINEAR * result.d),
                     (result.y_refresh, result.ell, C_LINEAR * result.qlen),
                     (result.y_next, result.ell, C_LINEAR * result.d),
                     (result.x_alt, result.ell, C_LINEAR * result.ell),
                     (result.q_alt, result.ell, C_LINEAR * result.ell)]
        for extractor, seed_length, entropy_budget in contracts:
            if (not extractor.capacity_ok()
                    or extractor.seed_length > seed_length
                    or extractor.entropy_needed > entropy_budget):
                return None
        result.resilient = ResilientShape.build(result.row_count)
        if result.resilient is None:
            return None
        resilience = result.resilient
        moment_order = resilience.moment_order
        if result.h < moment_order * resilience.w:
            return None
        charge = (20 + 2 * ceil_log2(result.b + 1) + 4 * ceil_log2(result.h + 1)
                  + ceil_log2(moment_order + 1)
                  + moment_order * (1 + ceil_log2(resilience.blocks + 1))
                  + 10 * resilience.order * resilience.z + 100)
        if result.error_exponent < charge:
            return None
        return result


class IndexedEngine:
    def __init__(self, shape, payload):
        self.shape, self.payload = shape, payload

    def lookahead(self, q_word):
        shape = self.shape
        mask = (1 << shape.ell) - 1
        seed = q_word & mask
        responses = []
        for j in range(shape.h):
            response = shape.x_alt(self.payload, seed)
            responses.append(response)
            if j + 1 < shape.h:
                seed = shape.q_alt(q_word, response)
        return responses

    def row(self, index):
        shape = self.shape
        original = shape.short(self.payload, index)
        current = original & ((1 << shape.d) - 1)
        index_word = index << (shape.b * shape.l - shape.s)
        selected = 0
        for j in range(shape.b):
            response = shape.x_refresh(self.payload, current)
            q_word = shape.y_refresh(original, response)
            responses = self.lookahead(q_word)
            digit = (index_word >> ((shape.b - j - 1) * shape.l)) & (shape.h - 1)
            selected = responses[digit]
            if j + 1 < shape.b:
                current = shape.y_next(original, selected)
        return selected

    def evaluate(self):
        word, index = 0, 0
        while index < self.shape.row_count:
            word |= (self.row(index) & 1) << index
            index += 1
        return self.shape.resilient.evaluate(word)


def evaluate(bits):
    shape = FamilyShape.build(len(bits))
    if shape is None:
        return 0
    payload = 0
    for i, bit in enumerate(bits):
        if bit not in (0, 1):
            raise ValueError('binary input required')
        payload |= bit << i
    return IndexedEngine(shape, payload).evaluate()


def binary_rank(words):
    pivots = {}
    for word in words:
        while word:
            leading = word.bit_length() - 1
            if leading in pivots:
                word ^= pivots[leading]
            else:
                pivots[leading] = word
                break
    return len(pivots)


def finite_audit():
    """Bounded evidence for primitives, not the universal disperser proof."""
    checks = {}
    count = 0
    for base in range(12):
        for exponent in range(14):
            for cap in (0, 1, 5, 97, 10000):
                assert capped_power(base, exponent, cap) == min(base ** exponent, cap + 1)
                count += 1
    checks['capped_power_cases'] = count
    count = 0
    for n in range(256):
        for degree in range(1, 13):
            value = floor_root(n, degree)
            assert value ** degree <= n < (value + 1) ** degree
            assert ceil_root(n, degree) == value + (value ** degree != n)
            count += 1
    checks['integer_root_cases'] = count

    base = small_field(2)
    count = 0
    for lower_word in range(16):
        lower = coefficient_list(lower_word, 2, 2)
        has_root = any((base.mul(x, x) ^ base.mul(lower[1], x) ^ lower[0]) == 0 for x in range(4))
        assert extension_irreducible(base, lower) == (not has_root)
        count += 1
    checks['extension_irreducibility_cases'] = count
    ext = extension_field(2, 2)
    for x in range(1, 16):
        assert ext.power(x, 15) == 1
        assert ext.mul(x, ext.power(x, 14)) == 1
    checks['extension_inverse_cases'] = 15

    orbit = extension_field(1, 16)
    generator = primitive(orbit)
    reference = PolynomialField((1 << 16) | pack_coefficients(orbit.lower, 1))
    count = 0
    state = 17
    for _ in range(128):
        state = (1103515245 * state + 12345) & 65535
        first = state
        state = (1103515245 * state + 12345) & 65535
        second = state
        assert orbit.mul(first, second) == reference.mul(first, second)
        count += 1
    checks['extension_multiplication_cases'] = count
    assert orbit.power(generator, 65535) == 1
    assert all(orbit.power(generator, 65535 // p) != 1 for p in prime_divisors(65535))
    checks['primitive_order_cases'] = 1

    toy = ShortAffine.__new__(ShortAffine)
    toy.w2, toy.su_field = 1, small_field(1)
    toy.monomials = monomial_prefix(16, 2, 8)
    count = 0
    for first in range(16):
        second = (first * 37 + 11) & 255
        vector = (first * 3191 + 101) & 65535
        a = toy.su_value(coefficient_list(first, 1, 8), vector)
        b = toy.su_value(coefficient_list(second, 1, 8), vector)
        c = toy.su_value(coefficient_list(first ^ second, 1, 8), vector)
        assert a ^ b == c
        direct = 0
        for i, exponents in enumerate(toy.monomials):
            monomial = 1
            for j, exponent in enumerate(exponents):
                if exponent:
                    monomial &= (vector >> j) & 1
            direct ^= ((first >> i) & 1) & monomial
        assert a == direct
        count += 1
    checks['orbit_polynomial_linearity_cases'] = count

    fold = small_field(5)
    gamma = primitive(fold)
    subspaces = set()
    for x in range(1, 16):
        for y in range(1, 16):
            if x != y:
                subspaces.add(tuple(sorted((x, y, x ^ y))))
    count = 0
    for space in subspaces:
        x, y = space[:2]
        bad = 0
        for seed in range(32):
            images = []
            for word in (x, y):
                image, point = 0, seed
                for j in range(2):
                    value = 0
                    for i in range(3, -1, -1):
                        value = fold.mul(value, point) ^ ((word >> i) & 1)
                    image |= value << (5 * j)
                    point = fold.mul(point, gamma)
                images.append(image)
            bad += binary_rank(images) < 2
            count += 1
        assert bad <= 2 * 4
    checks['folded_affine_rank_cases'] = count

    graphs, owners = [], []
    for row in range(16):
        coefficients = coefficient_list(row, 2, 2)
        for offset in range(4):
            points = set()
            for x in range(3):
                y = base.mul(base.mul(coefficients[1], x) ^ coefficients[0], x) ^ offset
                points.add(x * 4 + y)
            graphs.append(points)
            owners.append(row)
    count = 0
    for i in range(64):
        for j in range(i):
            if owners[i] != owners[j]:
                assert len(graphs[i] & graphs[j]) <= 2
                count += 1
            else:
                assert not graphs[i] & graphs[j]
    assert len({tuple(sorted(points)) for points in graphs}) == 64
    checks['polynomial_graph_intersections'] = count

    count = 0
    for size in range(2, 7):
        left_size = min(2, size // 2)
        right_size = min(2, size - left_size)
        def increasing(table, width):
            return all(not ((table >> x) & 1) or ((table >> (x | (1 << j))) & 1)
                       for x in range(1 << width) for j in range(width))
        left_tables = [t for t in range(1 << (1 << left_size)) if increasing(t, left_size)]
        right_tables = [t for t in range(1 << (1 << right_size)) if increasing(t, right_size)]
        for first in left_tables:
            for second in right_tables:
                sg = sh = sgh = 0
                for word in range((1 << size) - 1):
                    g = (first >> (word & ((1 << left_size) - 1))) & 1
                    h = (second >> ((word >> left_size) & ((1 << right_size) - 1))) & 1
                    sg += g
                    sh += h
                    sgh += g * h
                population = (1 << size) - 1
                assert population * sgh <= sg * sh
                count += 1
    checks['truncated_cube_covariance_cases'] = count
    count = 0
    for triggers in range(101):
        for order in range(2, 22, 2):
            upper = sum((-1) ** j * comb(triggers, j) for j in range(min(order, triggers) + 1))
            lower = sum((-1) ** j * comb(triggers, j) for j in range(min(order - 1, triggers) + 1))
            exact = int(triggers == 0)
            assert lower <= exact <= upper
            assert upper - lower == (comb(triggers, order) if triggers >= order else 0)
            count += 1
    checks['bonferroni_cases'] = count

    count = 0
    for n in (0, 1, 2, 8, 16, 32, 256, 10000):
        assert FamilyShape.build(n) is None
        assert evaluate([0] * n) == 0
        assert evaluate([1] * n) == 0
        count += 2
    checks['actual_zero_guard_cases'] = count
    return {'all_passed': True, 'checks': checks,
            'scope': 'Finite exact predicates only. Full-family tests exercise only the explicit zero guard. Universal affine extraction, transcript conditioning, resilience, uniformity and the circuit theorem rely on the supplied mathematical proofs.'}


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared finite audit takes {}')
    print(json.dumps(finite_audit(), sort_keys=True))
