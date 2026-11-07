from math import gcd as integer_gcd
"""Total tree disperser evaluator. Proofs in the supplied evidence."""
from math import isqrt, comb
from functools import lru_cache
from itertools import combinations

def is_prime(n):
    return n >= 2 and (n == 2 or (n % 2 and all((n % d for d in range(3, isqrt(n) + 1, 2)))))

@lru_cache(None)
def prime_at_least(n):
    n = max(3, n)
    if n % 2 == 0:
        n += 1
    while not is_prime(n):
        n += 2
    return n

def ceil_log2(n):
    return (n - 1).bit_length() if n > 1 else 0

def root(n, d):
    if n < 0 or d < 1:
        raise ValueError('nonnegative radicand, positive root')
    if n < 2:
        return n
    lo = 0
    hi = 1 << (n.bit_length() + d - 1) // d
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** d <= n:
            lo = mid
        else:
            hi = mid
    return lo

def ceil_root(n, d):
    a = root(n, d)
    return a + (a ** d < n)

def pack(bits):
    return sum((b << j for j, b in enumerate(bits)))

def power(a, e, f):
    v = 1
    while e:
        if e & 1:
            v = mul(v, a, f)
        a = mul(a, a, f)
        e >>= 1
    return v

def trace(a, f):
    v = 0
    for _ in range(f.bit_length() - 1):
        v ^= a
        a = mul(a, a, f)
    assert v in (0, 1)
    return v

def mod(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << a.bit_length() - b.bit_length()
    return a

def gcd(a, b):
    while b:
        a, b = (b, mod(a, b))
    return a

def quotient(a, b):
    q = 0
    while a and a.bit_length() >= b.bit_length():
        s = a.bit_length() - b.bit_length()
        q ^= 1 << s
        a ^= b << s
    assert not a
    return q

def kernel(rows, n):
    rows = list(rows)
    pivots = []
    rank = 0
    for col in range(n):
        pivot = next((i for i in range(rank, n) if rows[i] >> col & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = (rows[pivot], rows[rank])
        for i in range(n):
            if i != rank and rows[i] >> col & 1:
                rows[i] ^= rows[rank]
        pivots.append(col)
        rank += 1
    basis = []
    for free in range(n):
        if free in pivots:
            continue
        v = 1 << free
        for row, pivot in enumerate(pivots):
            if rows[row] >> free & 1:
                v ^= 1 << pivot
        basis.append(v)
    return basis

def factor_one(f):
    while True:
        d = f.bit_length() - 1
        cols = [mod(1 << 2 * j, f) ^ 1 << j for j in range(d)]
        rows = [sum(((cols[j] >> i & 1) << j for j in range(d))) for i in range(d)]
        basis = kernel(rows, d)
        if len(basis) == 1:
            return f
        b = next((b for b in basis if b not in (0, 1)))
        g = gcd(f, b)
        h = quotient(f, g)
        assert g not in (1, f) and h not in (1, f)
        f = min(g, h)

def mul(a, b, f):
    out = 0
    degree = f.bit_length() - 1
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> degree & 1:
            a ^= f
    return out

@lru_cache(None)
def irreducible_prime(k):
    assert k >= 3 and all((k % d for d in range(2, isqrt(k) + 1)))
    residue = 2 % k
    m = 1
    while residue != 1:
        residue = 2 * residue % k
        m += 1
    order = (1 << m) - 1
    valuation = 0
    while order % k == 0:
        order //= k
        valuation += 1
    f = factor_one((1 << k) - 1)
    assert f.bit_length() - 1 == m
    for _ in range(1, valuation):
        lifted = sum(((f >> j & 1) << j * k for j in range(m + 1)))
        f = factor_one(lifted)
        assert f.bit_length() - 1 == m
    a = 2

    def kmul(x, y):
        return mul(x, y, f)
    zero = (0,) * k
    one = (1,) + (0,) * (k - 1)
    beta = (0, 1) + (0,) * (k - 2)

    def add(x, y):
        return tuple((a ^ b for a, b in zip(x, y)))

    def square(x):
        out = [0] * k
        for j, c in enumerate(x):
            if not c:
                continue
            c = kmul(c, c)
            if 2 * j >= k:
                c = kmul(c, a)
            out[2 * j % k] ^= c
        return tuple(out)

    def emul(x, y):
        tmp = [0] * (2 * k - 1)
        for i, c in enumerate(x):
            if c:
                for j, d in enumerate(y):
                    if d:
                        tmp[i + j] ^= kmul(c, d)
        for j in range(2 * k - 2, k - 1, -1):
            if tmp[j]:
                tmp[j - k] ^= kmul(tmp[j], a)
        return tuple(tmp[:k])
    current = beta
    gamma = zero
    for _ in range(m):
        gamma = add(gamma, current)
        for _ in range(k):
            current = square(current)
    assert current == beta and any(gamma[1:])
    poly = [one]
    root = gamma
    for _ in range(k):
        out = [zero] * (len(poly) + 1)
        for j, c in enumerate(poly):
            out[j + 1] = add(out[j + 1], c)
            out[j] = add(out[j], emul(c, root))
        poly = out
        root = square(root)
    assert root == gamma
    assert all((not any(c[1:]) and c[0] in (0, 1) for c in poly))
    return sum((c[0] << j for j, c in enumerate(poly)))

@lru_cache(None)
def code_parameters(n, m, z):
    if n < 1 or m < 1 or z < 0:
        return None
    bound = 1024 * (n + 1) * (2 * m) ** 6 * (1 << 6 * z)
    s = prime_at_least(bound.bit_length())
    return (s, 2 * s)

def coefficients(bits, s):
    return tuple((pack(bits[j:j + s]) for j in range(0, len(bits), s)))

def code_bit(coeff, s, index):
    f = irreducible_prime(s)
    mask = (1 << s) - 1
    a = index & mask
    b = index >> s & mask
    value = 0
    for c in reversed(coeff):
        value = mul(value, a, f) ^ c
    return trace(mul(b, value, f), f)

def prime_divisors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return tuple(out)

def digits(a, s, length=None):
    q = 1 << s
    out = []
    while a:
        out.append(a & q - 1)
        a >>= s
    if length is not None:
        if len(out) > length:
            raise ValueError('coefficient range')
        out.extend([0] * (length - len(out)))
    return out

def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a

def poly_remainder(a, b, f):
    a = trim(a)
    b = trim(b)
    if not b:
        raise ZeroDivisionError
    inverse = power(b[-1], (1 << f.bit_length() - 1) - 2, f)
    while len(a) >= len(b):
        c = mul(a[-1], inverse, f)
        j = len(a) - len(b)
        for i, v in enumerate(b):
            a[j + i] ^= mul(c, v, f)
        a = trim(a)
    return a

def poly_gcd(a, b, f):
    a = trim(a)
    b = trim(b)
    while b:
        a, b = (b, poly_remainder(a, b, f))
    if not a:
        return []
    inverse = power(a[-1], (1 << f.bit_length() - 1) - 2, f)
    return [mul(v, inverse, f) for v in a]

def extension_mul(a, b, s, modulus, f):
    d = len(modulus) - 1
    aa = digits(a, s, d)
    bb = digits(b, s, d)
    tmp = [0] * (2 * d - 1)
    for i, x in enumerate(aa):
        if x:
            for j, y in enumerate(bb):
                if y:
                    tmp[i + j] ^= mul(x, y, f)
    for j in range(2 * d - 2, d - 1, -1):
        c = tmp[j]
        if c:
            for i in range(d):
                tmp[j - d + i] ^= mul(c, modulus[i], f)
    return sum((v << s * i for i, v in enumerate(tmp[:d])))

def extension_power(a, e, s, modulus, f):
    out = 1
    while e:
        if e & 1:
            out = extension_mul(out, a, s, modulus, f)
        a = extension_mul(a, a, s, modulus, f)
        e >>= 1
    return out

def extension_irreducible(modulus, s, f):
    d = len(modulus) - 1
    if d < 1 or modulus[-1] != 1:
        return False
    q = 1 << s
    z = 0 if d == 1 else q
    if d == 1:
        z = modulus[0]
    current = z
    check = {d // p for p in prime_divisors(d)}
    for j in range(1, d + 1):
        current = extension_power(current, q, s, modulus, f)
        if j in check:
            diff = digits(current ^ z, s)
            if len(poly_gcd(modulus, diff, f)) != 1:
                return False
    return current == z

@lru_cache(None)
def singer(s, d):
    if s < 3 or not is_prime(s) or d < 1:
        raise ValueError('field parameters')
    q = 1 << s
    Q = q ** d
    f = irreducible_prime(s)
    if d == 1:
        modulus = (0, 1)
    else:
        modulus = None
        for a in range(Q):
            candidate = tuple(digits(a, s, d) + [1])
            if extension_irreducible(candidate, s, f):
                modulus = candidate
                break
        assert modulus is not None
    divisors = prime_divisors(Q - 1)
    for g in range(1, Q):
        if all((extension_power(g, (Q - 1) // p, s, modulus, f) != 1 for p in divisors)):
            return (f, modulus, g)
    raise AssertionError('primitive element existence')
def rational(numerator,denominator):
    divisor=integer_gcd(numerator,denominator)
    return numerator//divisor,denominator//divisor

def valid_mode(n,mode):
    if type(n)is not int or n<1 or type(mode)is not tuple or len(mode)!=2:return False
    p,q=mode
    return (type(p)is int and type(q)is int and p>0 and q>0
            and (3*q<=20*p and 25*p<=13*q) and max(p,q)<=100*ceil_log2(n+1))

@lru_cache(None)
def regime(mode):
    if type(mode)is not tuple or len(mode)!=2:raise ValueError('rational rate')
    p,q=mode
    if type(p)is not int or type(q)is not int or p<1 or q<1 or not (3*q<=20*p and 25*p<=13*q):
        raise ValueError('rate outside3/20..13/25')
    an=100000*q+q-p;ad=100000*(q-p)
    sn=100001*ad*q+400000*an*p;sd=100000*ad*q
    cn=9999*(q-p)*sd;cd=10000*q*sn
    return dict(beta=rational(p,q),alpha=rational(an,ad),
                zeta=rational(1000000*p-4*q,4000000*q),
                group=rational(q-p,q),slice=rational(1000000*p-q,1000000*q),
                seed_guard=rational(sn,sd),arity=rational(cn,cd))

def field_cap(mode):
    p,q=mode
    if 100*p<=19*q:return 2
    if 10*p<=3*q:return 3
    if 5*p<=2*q:return 4
    if 25*p<=12*q:return 5
    return 6

@lru_cache(None)
def optimal_parameters(n, K, mode=(1,2), take=None):
    if not valid_mode(n,mode):return None
    if not 2 <= K <= n:
        return None
    L = ceil_log2(n + 1)
    w = K.bit_length() - 1
    bn, bd = regime(mode)['beta']
    an, ad = regime(mode)['alpha']
    m = root(K ** bn, bd)
    if take is not None:
        if type(take)is not int or not 1<=take<=m:return None
        m=take
    d = (an * L + ad * w - 1) // (ad * w)
    h = d * ceil_root(n, d) + 1
    bound = (1 << 30) * h * d * (32 * m) ** 4
    s = prime_at_least(bound.bit_length())
    while s <= (L+1)**2 and (1<<s) <= bound*s:
        s = prime_at_least(s+1)
    q = 1 << s
    if s > (L + 1) ** 2:
        return None
    Q = q ** d
    if Q > (n + 1) ** field_cap(mode) or h ** d > (n + 1) ** field_cap(mode):
        return None
    if q <= (1 << 30) * h * d * s * (32 * m) ** 4 or h >= q - 1 or m >= Q:
        return None
    if K < 400 * m * h * d * s * s + ceil_log2(32 * m) + 3:
        return None
    return (m, s, d, h, (d + 1) * s)

def exact_degree_exponents(d, total):
    for bars in combinations(range(total + d - 1), d - 1):
        ends = (-1,) + bars + (total + d - 1,)
        yield tuple((ends[j + 1] - ends[j] - 1 for j in range(d)))

def monomial_exponents(d, h):
    for total in range(h):
        yield from exact_degree_exponents(d, total)

def evaluate_message(bits, vector, s, d, h, f):
    coords = digits(vector, s, d)
    powers = []
    for c in coords:
        row = [1]
        for _ in range(1, h):
            row.append(mul(row[-1], c, f))
        powers.append(row)
    value = 0
    for bit, exponents in zip(bits, monomial_exponents(d, h)):
        if bit:
            term = 1
            for j, e in enumerate(exponents):
                term = mul(term, powers[j][e], f)
            value ^= term
    return value

def optimal_extract_with_parameters(bits, seed, s, d, h, take):
    f, modulus, g = singer(s, d)
    Q = 1 << s * d
    v = seed & Q - 1
    b = seed >> s * d & (1 << s) - 1
    out = []
    for _ in range(take):
        v = extension_mul(v, g, s, modulus, f)
        value = evaluate_message(bits, v, s, d, h, f)
        out.append(trace(mul(b, value, f), f))
    return tuple(out)

def optimal_extract(bits, seed, K, take, mode=(1,2)):
    p = optimal_parameters(len(bits), K, mode=mode)
    if p is None or not 0 <= take <= p[0]:
        raise ValueError('invalid optimal extractor parameters')
    m, s, d, h, seedlen = p
    if not 0 <= seed < 1 << seedlen:
        raise ValueError('seed range')
    return optimal_extract_with_parameters(bits, seed, s, d, h, take)
constant_parameters = optimal_parameters
constant_extract = optimal_extract

def response_parameters(n, K, take, mode=(1,2)):
    return optimal_parameters(n,K,mode=mode,take=take)

def response_extract(bits, seed, K, take, mode=(1,2)):
    p=response_parameters(len(bits),K,take,mode=mode)
    if p is None:raise ValueError('invalid response parameters')
    m,s,d,h,seedlen=p
    if not 0<=seed<(1<<seedlen):raise ValueError('seed range')
    return optimal_extract_with_parameters(bits,seed,s,d,h,take)

def restrict_seed(seed, indices):
    return sum(((seed >> j & 1) << i for i, j in enumerate(indices)))

@lru_cache(None)
def small_parameters(n, r, mode=(1,2)):
    if not valid_mode(n,mode):return None
    if not 1 <= r <= n or ceil_log2(n + 1) > root(r, 2000):
        return None
    m = root(r ** 999999, 1000000)
    z = root(r ** regime(mode)['zeta'][0], regime(mode)['zeta'][1])
    if z < 2:
        return None
    s, ell = code_parameters(n, m, z)
    p = prime_at_least(max(3, ceil_log2(ell)))
    Q = 1 << p
    count = 1
    A = 0
    while count < m:
        count *= Q
        A += 1
    A = max(1, A)
    d = Q * Q
    need = d + m * (1 << A - 1) + s + z + ceil_log2(4 * m)
    if ell > Q or A > ell or 2 * need > r:
        return None
    return (m, z, s, ell, p, Q, A, d)

@lru_cache(None)
def small_set(n, r, i, mode=(1,2)):
    pp = small_parameters(n, r, mode=mode)
    if pp is None:
        raise ValueError('invalid small-error design')
    m, z, s, ell, p, Q, A, d = pp
    if not 0 <= i < m:
        raise ValueError('design row')
    f = irreducible_prime(p)
    coeff = []
    value = i
    for _ in range(A):
        coeff.append(value % Q)
        value //= Q
    assert value == 0
    result = []
    for alpha in range(ell):
        beta = 0
        for c in reversed(coeff):
            beta = mul(beta, alpha, f) ^ c
        result.append(alpha * Q + beta)
    return tuple(result)

def small_extract(bits, seed, r, take, mode=(1,2)):
    p = small_parameters(len(bits), r, mode=mode)
    if p is None or not 0 <= take <= p[0]:
        raise ValueError('invalid small-error extractor')
    coeff = coefficients(bits, p[2])
    return tuple((code_bit(coeff, p[2], restrict_seed(seed, small_set(len(bits), r, i, mode=mode))) for i in range(take)))

@lru_cache(None)
def near_parameters(n, r, mode=(1,2)):
    p = small_parameters(n, r, mode=mode)
    M = r // 4
    if p is None or M < 1:
        return None
    J = (M + p[0] - 1) // p[0]
    return (M, J, J * p[-1])

def near_extract(bits, seed, r, mode=(1,2)):
    pp = near_parameters(len(bits), r, mode=mode)
    if pp is None:
        raise ValueError('invalid iterated extractor')
    M, J, dtotal = pp
    p = small_parameters(len(bits), r, mode=mode)
    d = p[-1]
    mask = (1 << d) - 1
    out = []
    for j in range(J):
        take = min(p[0], M - len(out))
        out.extend(small_extract(bits, seed >> j * d & mask, r, take, mode=mode))
    return tuple(out)

@lru_cache(None)
def cubic_parameters(w):
    if w < 1:
        return None
    k = prime_at_least(2 * w)
    out = isqrt(w) // 100
    d0 = (2 * k + 4) // 5 + 10
    if out < 1 or 5 * k > 11 * w or 10 * (w - d0) < w:
        return None
    return (k, out)

def cubic_extract(bits, take):
    if len(bits) % 2:
        raise ValueError('two equal slices')
    w = len(bits) // 2
    p = cubic_parameters(w)
    if p is None or not 0 <= take <= p[1]:
        raise ValueError('invalid cubic seed')
    f = irreducible_prime(p[0])
    value = power(pack(bits), 7, f)
    return tuple((value >> i & 1 for i in range(take)))

@lru_cache(None)
def reservoir_parameters(n, K, R, t, take, mode=(1,2)):
    if not valid_mode(n,mode):return None
    if not 1 <= K <= n or R < 1 or t < 1 or (t > root(K ** regime(mode)['group'][0], regime(mode)['group'][1])):
        return None
    pp = near_parameters(n, K, mode=mode)
    if pp is None:
        return None
    M, J, seedlen = pp
    w = root(K ** regime(mode)['slice'][0], regime(mode)['slice'][1])
    if w < 1 or w > min(R, M) or (not 0 <= take <= min(R, M)):
        return None
    if (2 * t + ceil_log2(t)) * w > K // 4:
        return None
    cp = cubic_parameters(w)
    if cp is None or cp[1] < seedlen:
        return None
    counts = []
    while t > 1:
        counts.append(t)
        t = (t + 1) // 2
    return (w, M, seedlen, tuple(counts))

def reservoir_extract(bits, rows, K, take, mode=(1,2)):
    R = len(rows[0]) if rows else 0
    if any((len(row) != R for row in rows)):
        raise ValueError('ragged reservoir rows')
    pp = reservoir_parameters(len(bits), K, R, len(rows), take, mode=mode)
    if pp is None:
        raise ValueError('invalid entropy-reservoir parameters')
    w, M, seedlen, counts = pp
    rows = tuple((tuple(row) for row in rows))
    for count in counts:
        assert len(rows) == count
        nextrows = []
        for i in range(0, count, 2):
            pair = rows[i][:w] + rows[min(i + 1, count - 1)][:w]
            seed = pack(cubic_extract(pair, seedlen))
            nextrows.append(near_extract(bits, seed, K, mode=mode))
        rows = tuple(nextrows)
    return rows[0][:take]

@lru_cache(None)
def block_parameters(n, K, mode=(1,2)):
    if not valid_mode(n,mode):return None
    cp = constant_parameters(n, K, mode=mode)
    sp = small_parameters(n, K, mode=mode)
    R = root(K ** regime(mode)['beta'][0], regime(mode)['beta'][1])
    if cp is None or sp is None or R < 1 or (R > min(cp[0], sp[0])):
        return None
    g = root(K ** regime(mode)['group'][0], regime(mode)['group'][1])
    if g < 2:
        return None
    rp = reservoir_parameters(n, K, R, g, sp[-1], mode=mode)
    if rp is None:
        return None
    v = 1 << cp[-1]
    groups = 1
    b = 1
    while groups < v:
        groups *= g
        b += 1
    return (R, g, b, sp[-1])

def block_extract(blocks, K, mode=(1,2)):
    n = len(blocks[0]) if blocks else 0
    p = block_parameters(n, K, mode=mode)
    if p is None or len(blocks) != p[2] or any((len(x) != n for x in blocks)):
        raise ValueError('invalid block-source parameters')
    R, g, b, seedlen = p
    cp = constant_parameters(n, K, mode=mode)
    rows = tuple((constant_extract(blocks[0], seed, K, R, mode=mode) for seed in range(1 << cp[-1])))
    previous = blocks[0]
    for x in blocks[1:]:
        new = []
        for i in range(0, len(rows), g):
            group = rows[i:i + g]
            group += ((0,) * R,) * (g - len(group))
            seed = pack(reservoir_extract(previous, group, K, seedlen, mode=mode))
            new.append(small_extract(x, seed, K, R, mode=mode))
        rows = tuple(new)
        previous = x
    assert len(rows) == 1
    return rows[0]

def interval_bounds(n, t, level, index):
    total = t ** level
    return (index * n // total, (index + 1) * n // total)

def level_interval_bounds(n, log_products, level, index):
    total = 1 << log_products[level]
    return (index * n // total, (index + 1) * n // total)

def active_level_indices(log_products, gap):
    active = [0]
    crossing = gap
    for j, v in enumerate(log_products[1:], 1):
        if v >= crossing:
            active.append(j)
            crossing = (v // gap + 1) * gap
    return tuple(active)


from fractions import Fraction as F
from itertools import product
RANGES = ((F(3, 50), F(2, 25)), (F(1, 2), F(3, 5)), (F(1), F(13, 10)), (F(3), F(18, 5)), (F(1), F(6, 5)), (F(3, 25), F(3, 20)))
BASE = (F(7, 100), F(27, 50), F(527, 500), F(33743, 10000), F(25927, 25000), F(269991, 2000000))


def rationals(lower, upper, denominator_cap):
    """Unique increasing rationals in a closed interval with bounded denominator."""
    if denominator_cap < 1:
        return ()
    values = []
    for denominator in range(1, denominator_cap + 1):
        start = -((-lower.numerator * denominator) // lower.denominator)
        stop = (upper.numerator * denominator) // upper.denominator
        for numerator in range(start, stop + 1):
            if integer_gcd(numerator, denominator) == 1:
                values.append(F(numerator, denominator))
    return tuple(sorted(values))


def analytic_prefilter(config):
    if len(config)!=6:return False
    h, a, wm, w0, gap, confidence = config
    return (all(lo <= value <= hi for value, (lo, hi) in zip(config, RANGES))
            and a / wm < F(13, 25) and a / w0 > F(3, 20)
            and wm > gap
            and F(2000001, 1000000) * h / gap < confidence
            and confidence < a / 4 - F(102, 100000000) * w0
            )


def configurations(denominator_cap):
    """BASE first, followed by filtered grid tuples in fixed lexicographic order."""
    yield BASE
    axes = tuple(rationals(lo, hi, denominator_cap) for lo, hi in RANGES)
    for config in product(*axes):
        if config != BASE and analytic_prefilter(config):
            yield config


def primitive_error(z, near_rounds, slice_length, group_size, block_count):
    """Exact conservative error arithmetic, CONDITIONAL on primitive hypotheses.

    This function does not establish those hypotheses. The eventual integration
    must first check the actual small/near/cubic/reservoir/block constructors.
    floor(sqrt(w)/50) weakens the cubic bound, preserving its upper-bound direction.
    """
    from math import isqrt
    if min(z, near_rounds, slice_length, group_size, block_count) < 1:
        return None
    epsilon = F(1, 1 << z)
    near = near_rounds * epsilon
    cubic = 2 * F(1, 1 << (isqrt(slice_length) // 50))
    rounds = (group_size - 1).bit_length()
    group = rounds * (2 * near + cubic)
    block = (block_count - 1) * (2 * epsilon + group)
    return dict(small=epsilon, near=near, cubic=cubic, group=group, block=block)

def response_rank_budget(params):
    """At most t_j-1 earlier siblings are forced at each visited depth."""
    return sum((1<<a)-1 for a in params['arities'][:params['H']])*params['ell']*params['P']

def geometry(p, dimension):
    H, ell, width = p['H'], p['ell'], p['P']
    active, V, arities = p['active'], p['products'], p['arities']
    if H < 1 or len(active) < 2 or active[0] != 0 or active[-1] != H:
        return None
    if any(a >= b for a, b in zip(active, active[1:])):
        return None
    if not 0 < dimension <= p['n']:
        return None
    largest_gap = max(V[b] - V[a] for a, b in zip(active, active[1:]))
    if largest_gap < 1:
        return None
    # A conservative lower bound for BOTH useful active-node entropies.
    coefficient = 1 + sum((1 << (V[b]-V[a]))-1 for a,b in zip(active,active[1:]))
    initial = (dimension+coefficient-1)//coefficient
    rank = response_rank_budget(p)
    resiliency = (H + 1) * ell * width
    if initial - rank < max(p['thresholds'][:H]) + 1:
        return None
    factor = 1
    hidden = []
    for j in range(active[1],H):
        next_factor = factor * 2 * p['blocks'][j]
        if initial // next_factor - rank < p['thresholds'][j] + 1:
            return None
        margin = (initial // factor - rank) // (2 << arities[j])
        if margin < p['thresholds'][j] + resiliency + 1:
            return None
        hidden.append(margin)
        factor = next_factor
    largest_leaf = (p['n'] + (1 << V[H]) - 1) // (1 << V[H])
    if initial // factor - rank <= largest_leaf:
        return None
    return dict(initial=initial, rank=rank, resiliency=resiliency,
                hidden=tuple(hidden), final=initial // factor - rank,
                largest_leaf=largest_leaf)

def output_events(H, responses, ell, width, confidences, block_error):
    if min(H, responses, ell, width) < 1 or not 0 <= block_error <= 1:
        return None
    if any(not 1 <= first <= second <= width for first, second in confidences):
        return None
    eta = (2 * responses + 2) ** H * (block_error + F(1, 1 << width))
    path = H * responses * (F(1, 1 << width) + eta)
    ancestors = F(0)
    margins = []
    for first, second in confidences:
        zero = F(1, 1 << (ell * second))
        one = F(1, 1 << (ell * first)) - responses * (F(1, 1 << second) + eta)
        margin = min(zero, one) - path - ancestors
        if margin <= 0:
            return None
        margins.append(margin)
        ancestors += responses * (F(1, 1 << first) + eta)
    return dict(eta=eta, path=path, margins=tuple(margins))


def ceiling(value):
    return -((-value.numerator) // value.denominator)

def config_mode(W, Z, rate):
    value = max(F(3, 20), min(F(13, 25), rate * Z / W))
    return (value.numerator, value.denominator)

def config_entropy_log(Z, B, j, minimum, maximum):
    """Historical fixed-decline arithmetic witness; not the current planner."""
    return max(ceiling(minimum * Z), ceiling(maximum * Z) - (25001*j*B+99999)//100000)

def next_entropy_log(current, minimum, block_count):
    """Pay the actual integer log of the public ordinary entropy cost."""
    return max(minimum, current-ceil_log2(2*block_count))

@lru_cache(None)
def parameters_for_config(N, config):
    if len(config) != 6 or not all(lo <= value <= hi for value, (lo, hi) in zip(config, RANGES)):
        return None
    hs, rate, wm, w0, spacing, confidence = config
    if N < 2:
        return None
    L = ceil_log2(N)
    B = ceil_log2(L + 1)
    Z = ceil_root(L ** 3, 4) * ceil_root(B, 4)
    base = ceil_root(L, 4) * ceil_root(B ** 3, 4)
    Wmin = ceiling(wm * Z)
    Wmax = ceiling(w0 * Z)
    entropy_log = 13 * Z
    Lpad = ceil_root(Wmax ** 20, 19)
    floorlog = N.bit_length() - 1
    if entropy_log > floorlog or Lpad > floorlog or Wmin ** 3 < (L + 1) ** 2:
        return None
    k = 1 << entropy_log
    ell_log = ceiling(hs * ceil_root(L, 2) * ceil_root(B, 2))
    P_log = ceiling(confidence * Z)
    if min(ell_log, P_log) < 1 or ell_log + P_log > floorlog:
        return None
    ell = 1 << ell_log
    P = 1 << P_log
    gap = ceiling(spacing * Z)
    products = [0]
    arities = []
    ambients = []
    blocks = []
    seeds = []
    thresholds = []
    modes = []
    W = Wmax
    while True:
        v = products[-1]
        count = 1 << v
        j = len(products) - 1
        K = 1 << W
        mode = config_mode(W, Z, rate)
        cn, cd = regime(mode)['arity']
        sn, sd = regime(mode)['seed_guard']
        ambient = max((N + count - 1) // count, 1 << Lpad)
        Lj = ceil_log2(ambient)
        A = max(1, hs.numerator * cn * W * base * L // (hs.denominator * cd * Z * Lj))
        arities.append(A)
        if v + A > floorlog:
            break
        cp = response_parameters(ambient, K, ell*P, mode=mode)
        if cp is None or sd * cp[-1] > sn * Lj or ell * P > cp[0]:
            return None
        bp = block_parameters(ambient, K, mode=mode)
        if bp is None or bp[0] < P:
            return None
        ambients.append(ambient)
        blocks.append(bp[2])
        seeds.append(cp[-1])
        thresholds.append(K)
        modes.append(mode)
        products.append(v + A)
        W = next_entropy_log(W, Wmin, bp[2])
    Amax = max(arities)
    if Amax >= gap:
        return None
    all_active = active_level_indices(products, gap)
    accepted = all_active
    if len(accepted) < 2:
        return None
    H = accepted[-1]
    active = tuple((j for j in all_active if j <= H))
    hp = len(active) - 1
    total_rows = 0
    for j in range(H):
        t = 1 << arities[j]
        b = blocks[j]
        if t < b:
            return None
        total_rows += comb(t, b)
        if total_rows > ell:
            return None
    confidence_unit = (ell_log + 1999999) // 2000000
    confidences = tuple(((1 << 4000002 * (hp - u) * confidence_unit, 1 << 4000002 * (hp - u) * confidence_unit + 2000001 * confidence_unit) for u in range(hp)))
    if any((p2 > P for p1, p2 in confidences)):
        return None
    return dict(config=config, n=N, k=k, thresholds=tuple(thresholds), modes=tuple(modes), ell=ell, P=P, H=H, products=tuple(products), arities=tuple(arities), ambients=tuple(ambients), blocks=tuple(blocks), seeds=tuple(seeds), active=active, confidences=confidences)


def finite_certificate(params):
    if geometry(params, params['dimension']) is None:
        return False
    worst_error = F(0)
    for j in range(params['H']):
        ambient, K, mode = params['ambients'][j], params['thresholds'][j], params['modes'][j]
        bp = block_parameters(ambient, K, mode=mode)
        sp = small_parameters(ambient, K, mode=mode)
        np = near_parameters(ambient, K, mode=mode)
        if bp is None or sp is None or np is None:
            return False
        rp = reservoir_parameters(ambient, K, bp[0], bp[1], sp[-1], mode=mode)
        if rp is None:
            return False
        errors = primitive_error(sp[1], np[1], rp[0], bp[1], bp[2])
        worst_error = max(worst_error, errors['block'])
    responses = 1 << max(params['seeds'][:params['H']])
    return output_events(params['H'], responses, params['ell'], params['P'],
                         params['confidences'], worst_error) is not None

def minimum_geometry_dimension(params):
    H,ell,P=params['H'],params['ell'],params['P']
    active,V,arity=params['active'],params['products'],params['arities']
    if H<1 or len(active)<2:return None
    coefficient=1+sum((1<<(V[b]-V[a]))-1 for a,b in zip(active,active[1:]))
    rank=response_rank_budget(params)
    resiliency=(H+1)*ell*P
    needed=rank+max(params['thresholds'][:H])+1
    factor=1
    for j in range(params['active'][1],H):
        nxt=factor*2*params['blocks'][j]
        needed=max(needed,nxt*(rank+params['thresholds'][j]+1),factor*(rank+2*(1<<arity[j])*(params['thresholds'][j]+resiliency+1)))
        factor=nxt
    leaf=(params['n']+(1<<V[H])-1)//(1<<V[H])
    needed=max(needed,factor*(rank+leaf+1))
    return coefficient*(needed-1)+1

@lru_cache(None)
def tree_parameters(N):
    if N < 2:
        return None
    L = ceil_log2(N)
    B = ceil_log2(L + 1)
    Z = ceil_root(L ** 3, 4) * ceil_root(B, 4)
    if 13 * Z > N.bit_length() - 1:
        return None
    best = None
    for config in configurations(B):
        plan = parameters_for_config(N, config)
        if plan is None:
            continue
        dimension = minimum_geometry_dimension(plan)
        if dimension is None or dimension > N:
            continue
        if best is not None and dimension >= best['dimension']:
            continue
        candidate = plan | dict(dimension=dimension)
        if finite_certificate(candidate):
            best = candidate
    return best

def evaluate(bits):
    if any((type(b) is not int or b not in (0, 1) for b in bits)):
        raise ValueError('binary input')
    p = tree_parameters(len(bits))
    if p is None:
        return 0
    n = p['n']
    thresholds = p['thresholds']
    ell = p['ell']
    P = p['P']
    H = p['H']
    products = p['products']
    arities = p['arities']
    ambients = p['ambients']
    blocks = p['blocks']
    data = tuple(bits)
    challenges = {}
    choices = {}

    def ambient(level, index):
        start, end = level_interval_bounds(n, products, level, index)
        return data[start:end] + (0,) * (ambients[level] - end + start)

    def test(level, index, challenge, confidence):
        x = ambient(level, index)
        K = thresholds[level]
        target = challenge + ((0,) * P,) * (ell - len(challenge))
        for seed in range(1 << p['seeds'][level]):
            response = response_extract(x, seed, K, ell * P, mode=p['modes'][level])
            if all((response[j * P:j * P + confidence] == target[j][:confidence] for j in range(ell))):
                return 0
        return 1
    for index in range(1 << products[H]):
        challenges[H, index] = ()
    for level in range(H - 1, -1, -1):
        t = 1 << arities[level]
        b = blocks[level]
        K = thresholds[level]
        for index in range(1 << products[level]):
            selected = 0
            for j in range(t):
                if test(level, index, challenges[level + 1, index * t + j], P):
                    selected = j
                    break
            choices[level, index] = selected
            rows = list(challenges[level + 1, index * t + selected])
            for indices in combinations(range(t), b):
                source_blocks = []
                for j in indices:
                    child = index * t + j
                    start, end = level_interval_bounds(n, products, level + 1, child)
                    source_blocks.append(data[start:end] + (0,) * (ambients[level] - end + start))
                rows.append(block_extract(tuple(source_blocks), K, mode=p['modes'][level])[:P])
            challenges[level, index] = tuple(rows)
    path = [0]
    for level in range(H):
        path.append((path[-1] << arities[level]) + choices[level, path[-1]])
    for u, (p1, p2) in enumerate(p['confidences']):
        level = p['active'][u]
        nextlevel = p['active'][u + 1]
        challenge = challenges[nextlevel, path[nextlevel]]
        if not test(level, path[level], challenge, p1):
            return test(level, path[level], challenge, p2)
    return 0

# BASE compatibility interfaces used solely by the scoped component audit.
def level_mode(W,Z):
    return config_mode(W,Z,BASE[1])

def level_entropy_log(Z,B,j):
    return config_entropy_log(Z,B,j,BASE[2],BASE[3])
