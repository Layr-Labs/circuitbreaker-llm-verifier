"""Fixed total certificate-prefix evaluator for proof.md.

All fixed integer parameters are part of one uniform algorithm. Its polynomial
runtime has an enormous exponent. At physically practical lengths the explicit
shape guards return zero; the theorem concerns every sufficiently large integer
length and proves that those guards are eventually false.
"""
from functools import lru_cache

STRENGTH = 1024
LOG_STRENGTH = 10
SHORT_ALPHABET_BITS = 1024 * STRENGTH


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


def greedy_design(ell, count, intersection, alphabet_bits=SHORT_ALPHABET_BITS):
    q = 1 << alphabet_bits
    words, graphs = [], []
    for value in range(q ** ell):
        word = tuple((value >> (alphabet_bits * j)) & (q - 1) for j in range(ell))
        if all(sum(x == y for x, y in zip(word, previous)) <= intersection for previous in words):
            words.append(word)
            graphs.append(tuple(j * q + word[j] for j in range(ell)))
            if len(words) == count:
                return ell * q, graphs
    raise ValueError('greedy design exhausted')


def extract_code_bits(code, sets, x, seed):
    out = 0
    for i, positions in enumerate(sets):
        index = sum(((seed >> position) & 1) << j for j, position in enumerate(positions))
        out |= code.bit(x, index) << i
    return out


class SeededExtractor:
    """Linear code-bit evaluator; parameters are checked by the family shape."""
    def __init__(self, n, entropy, output, error_exp, short=False):
        self.n, self.entropy, self.output = n, entropy, output
        self.error_exp, self.short = error_exp, short
        self.code = LinearCode(n, output, error_exp)
        if short:
            self.intersection = (entropy.bit_length() - 1) // 4
            self.seed_length = self.code.ell << SHORT_ALPHABET_BITS
        else:
            self.intersection = 4096
            self.seed_length = self.code.ell * (1 << (self.code.ell - 1).bit_length())
        self._sets = None

    def capacity_ok(self):
        if self.short:
            return self.output * (1 << self.code.ell) < 1 << (SHORT_ALPHABET_BITS * (self.intersection + 1))
        q = 1 << (self.code.ell - 1).bit_length()
        return self.output <= q ** (self.intersection + 1)

    def sets(self):
        if self._sets is None:
            if self.short:
                dimension, self._sets = greedy_design(self.code.ell, self.output, self.intersection)
            else:
                dimension, self._sets = polynomial_design(self.code.ell, self.output, self.intersection)
            assert dimension == self.seed_length
        return self._sets

    def __call__(self, x, seed):
        return extract_code_bits(self.code, self.sets(), x, seed)


class FamilyShape:
    """Exact fixed parameter shapes; failed guards specify output zero."""
    B_EXP = SHORT_ALPHABET_BITS + LOG_STRENGTH + 40
    TREE_EXP = SHORT_ALPHABET_BITS + 3 * LOG_STRENGTH + 90
    B = 1 << B_EXP
    H = 1 << (TREE_EXP + 2)
    Q = 1 << (SHORT_ALPHABET_BITS + 2 * LOG_STRENGTH + 80)
    J = 1 << (14 + LOG_STRENGTH)

    @classmethod
    def build(cls, n):
        t = floor_root(n, 1 << cls.TREE_EXP)
        if t < cls.B + 2:
            return None
        result = cls()
        result.n, result.t = n, t
        result.cutoff = floor_root(n, 4 * STRENGTH)
        if result.cutoff < 2 * t:
            return None
        result.k = ceil_root(n, STRENGTH)
        result.K = floor_root(n, 2 * STRENGTH)
        result.R = floor_root(result.K, 8)
        result.pmax = floor_root(n, STRENGTH * (1 << 24))
        result.ell = cls.H * t ** cls.B + 1
        if result.ell ** (6 * cls.J + 10) > result.pmax:
            return None
        result.short = SeededExtractor(n, result.K, result.R, 2, short=True)
        if result.ell * result.pmax > result.R:
            return None
        # Explicit greedy-design existence guard, eventually satisfied.
        if not result.short.capacity_ok():
            return None
        result.group_log = (result.R.bit_length() - 1) // (2 * (1 << 20))
        if result.group_log < 1:
            return None
        result.G = 1 << result.group_log
        if (result.short.seed_length + result.group_log - 1) // result.group_log >= cls.B:
            return None
        result.block_long = SeededExtractor(n, result.K, result.R, ceil_root(result.K, 1024))
        if not result.block_long.capacity_ok():
            return None
        result.rao_output = floor_root(result.R, 2)
        if result.block_long.seed_length > result.rao_output:
            return None
        result.stages = []
        r = result.R
        for i in range(result.group_log):
            next_r = r >> 4106
            if next_r < 1 or ((result.G >> i) * r) > n:
                return None
            w = floor_root(r, 10)
            degree = 2 * w - 1
            if degree < 3:
                return None
            ext = SeededExtractor(n, r // 2, next_r, ceil_root(r // 2, 1024))
            if not ext.capacity_ok():
                return None
            if ext.seed_length > floor_root(degree, 16):
                return None
            result.stages.append((r, next_r, w, degree, ext))
            r = next_r
        if r < result.rao_output:
            return None
        # One modulus per listed stage, even if degrees happen to repeat.
        result.prefix = sum(stage[3] for stage in result.stages)
        if result.prefix >= n:
            return None
        return result


class ExtractionEngine:
    def __init__(self, shape, payload, moduli):
        self.shape, self.payload = shape, payload
        self.fields = [PolynomialField(modulus) for modulus in moduli]
        self._challenges, self._selected = {}, {}

    def rao(self, rows):
        s = self.shape
        rows = list(rows) + [0] * (s.G - len(rows))
        for i, (r, next_r, w, degree, ext) in enumerate(s.stages):
            matrix = sum(row << (j * r) for j, row in enumerate(rows))
            next_rows = []
            mask = (1 << w) - 1
            for j in range(0, len(rows), 2):
                pair = ((rows[j] & mask) | ((rows[j + 1] & mask) << w)) & ((1 << degree) - 1)
                seed = self.fields[i].power(pair, 7) & ((1 << ext.seed_length) - 1)
                next_rows.append(ext(matrix, seed))
            rows = next_rows
        assert len(rows) == 1
        return rows[0] & ((1 << s.rao_output) - 1)

    def block(self, blocks):
        s = self.shape
        rows = [s.short(blocks[0], seed) for seed in range(1 << s.short.seed_length)]
        j = 1
        while len(rows) > 1:
            next_rows = []
            for start in range(0, len(rows), s.G):
                seed = self.rao(rows[start:start + s.G]) & ((1 << s.block_long.seed_length) - 1)
                next_rows.append(s.block_long(blocks[j], seed))
            rows = next_rows
            j += 1
        return rows[0]

    @lru_cache(None)
    def children(self, node):
        start, end, depth = node
        if end - start < self.shape.cutoff:
            return ()
        length, extra = divmod(end - start, self.shape.t)
        children = []
        left = start
        for j in range(self.shape.t):
            right = left + length + int(j < extra)
            children.append((left, right, depth + 1))
            left = right
        return tuple(children)

    def vector(self, node):
        start, end, depth = node
        return self.payload & (((1 << (end - start)) - 1) << start)

    @lru_cache(None)
    def slice_mask(self, confidence):
        s = self.shape
        return sum(((1 << confidence) - 1) << (j * s.pmax) for j in range(s.ell))

    def cr(self, node, descendant, confidence):
        s = self.shape
        challenge = sum(row << (j * s.pmax) for j, row in enumerate(self.challenge(descendant)))
        mask = self.slice_mask(confidence)
        wanted = challenge & mask
        x = self.vector(node)
        for seed in range(1 << s.short.seed_length):
            if s.short(x, seed) & mask == wanted:
                return False
        return True

    def selected(self, node):
        self.ensure(node)
        return self._selected[node]

    def challenge(self, node):
        self.ensure(node)
        return self._challenges[node]

    def ensure(self, node):
        if node in self._challenges:
            return
        s = self.shape
        pending = [(node, False)]
        while pending:
            current, ready = pending.pop()
            if current in self._challenges:
                continue
            children = self.children(current)
            if not children:
                self._challenges[current] = ()
                self._selected[current] = None
            elif not ready:
                pending.append((current, True))
                pending.extend((child, False) for child in reversed(children))
            else:
                # All child challenges are now cached. CR cannot recurse down.
                chosen_child = next((child for child in children
                                     if self.cr(current, child, s.pmax)), children[0])
                mask = (1 << s.pmax) - 1
                # Own rows FIRST give fixed indices despite mixed leaf depths.
                own = [self.block([self.vector(children[j]) for j in chosen]) & mask
                       for chosen in index_subsets(len(children), s.B)]
                self._challenges[current] = tuple(own) + self._challenges[chosen_child]
                self._selected[current] = chosen_child

    def output(self):
        s = self.shape
        path = [(0, s.n, 0)]
        while self.children(path[-1]):
            path.append(self.selected(path[-1]))
        for j in range(0, len(path) - s.Q, s.Q):
            u = j // s.Q
            p1 = s.ell ** (6 * (s.J - u) + 3)
            if not self.cr(path[j], path[j + s.Q], p1):
                p2 = s.ell ** (6 * (s.J - u) + 6)
                return int(self.cr(path[j], path[j + s.Q], p2))
        return 0


def evaluate(bits):
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError('the input alphabet is {0,1}')
    shape = FamilyShape.build(len(bits))
    if shape is None:
        return 0
    offset, moduli = 0, []
    for stage in shape.stages:
        degree = stage[3]
        modulus = (1 << degree) | sum(bits[offset + j] << j for j in range(degree))
        if not irreducible(modulus):
            return 0
        moduli.append(modulus)
        offset += degree
    payload = sum(bit << j for j, bit in enumerate(bits[offset:]))
    return ExtractionEngine(shape, payload, moduli).output()
