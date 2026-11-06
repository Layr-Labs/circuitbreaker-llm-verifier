"""Exact bounded supplementary predicates; universal claims rest on written proofs."""
import importlib.util
import itertools
from bisect import bisect_left
from hashlib import sha256
import json
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt, gcd
from pathlib import Path
from random import Random


def load_family():
    path = Path('/snapshot/files/family.py')
    if not path.exists():
        path = Path(__file__).with_name('family.py')
    spec = importlib.util.spec_from_file_location('candidate_family', path)
    f = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f)
    return f


def literal_dimensions(k, r, t):
    H = sum(range(2, r+2))
    T = sum(max(k-j, 0) for j in range(2, t+1))
    U = sum(max(k-j, 0) for j in range(1, r+1))
    old = H+max(T, U//(t+1))+1
    if k < 2*r+1 or 2*k < (t+1)**2:
        return old, old
    h = next(h for h in range(2, r+1) if h*(t+1) >= 2*r)
    F = Fraction(r*(2*k-2*r-1+h), 2*(t+1))
    new = H+max(T, F.numerator//F.denominator)+1
    assert Fraction(U, t+1)-F == Fraction(r*(r-h), 2*(t+1))
    assert new <= old
    return old, new


def tail_sum(k, A, t, count, x, literal):
    if literal:
        return sum(min(x, (k-j-A)//2) for j in range(t, count+1))
    # Sum consecutive floor-halves by pairing adjacent integers. This uses
    # a product prefix independently of family.py's squared-prefix formula.
    end = min(count, k-A-2*x)
    first, last = k-A-end-1, k-A-count-1
    first_sum = (first//2)*((first+1)//2) if first >= 0 else 0
    last_sum = (last//2)*((last+1)//2) if last >= 0 else 0
    return (end-t+1)*x+first_sum-last_sum


@lru_cache(None)
def literal_count_maximum(k, t, count):
    # Keep the independent count loop; share its identical exact A scan
    # across the positive/zero oracles and their different r prefixes.
    best = -1
    for A in range(t-1, k-count-1):
        upper = min(A//(t-1), (k-t-A)//2,
                    (2*(t-1)*k-t*(t-1)-2*t*A)//(2*(t-1)))
        if upper >= 1:
            best = max(best, A+tail_sum(k, A, t, count, upper, False))
    return best


@lru_cache(None)
def literal_integer_bound(k, r, t, witness=False):
    best = sum(k-j for j in range(2, t+1))
    best_dims = [k if i < t-1 else i+2 for i in range(r)]
    if not witness:
        return max([best]+[literal_count_maximum(k, t, count) for count in range(t, r+1)]), best_dims
    for count in range(t, r+1):
        for A in range(t-1, k-count-1):
            upper = min(A//(t-1), (k-t-A)//2,
                        (2*(t-1)*k-t*(t-1)-2*t*A)//(2*(t-1)))
            if upper < 1:
                continue
            # Small witness cases enumerate EVERY pivot, not only the optimizer.
            pivots = range(1, upper+1) if witness else (upper,)
            for x in pivots:
                value = A+tail_sum(k, A, t, count, x, witness)
                if value <= best:
                    continue
                best = value
                if witness:
                    top = [x]*(t-1)
                    left = A-sum(top)
                    for j in range(1, t):
                        add = min(left, k-A-x-j-top[j-1])
                        assert add >= 0
                        top[j-1] += add
                        left -= add
                    assert left == 0
                    deficits = top+[min(x, (k-j-A)//2) for j in range(t, count+1)]
                    best_dims = [d+i+2 for i, d in enumerate(deficits)]+[i+2 for i in range(count, r)]
    return best, best_dims



@lru_cache(None)
def literal_zero_bound(k, r, t, witness=False):
    H = sum(range(2, r+2))
    ds = [k if i < t-1 else i+1 for i in range(r)]
    best = sum(ds)-H
    if best < 0:
        best, ds = 0, [i+2 for i in range(r)]
    for count in range(1, t):
        A = min(k-count-1, (2*count*k-count*(count+1))//(2*(count+1)))
        if A < count or A <= best:
            continue
        best = A
        if witness:
            deficits = [1]*count
            left = A-count
            for j in range(1, count+1):
                add = min(left, k-A-j-1)
                deficits[j-1] += add
                left -= add
            assert left == 0
            ds = [x+i+2 for i, x in enumerate(deficits)]+[i+2 for i in range(count, r)]
    if not witness:
        return max([best]+[literal_count_maximum(k, t, count) for count in range(t, r+1)]), ds
    for count in range(t, r+1):
        for A in range(t-1, k-count-1):
            upper = min(A//(t-1), (k-t-A)//2,
                        (2*(t-1)*k-t*(t-1)-2*t*A)//(2*(t-1)))
            for x in (range(1, upper+1) if witness else (upper,)):
                if x < 1:
                    continue
                value = A+tail_sum(k, A, t, count, x, witness)
                if value <= best:
                    continue
                best = value
                if witness:
                    top = [x]*(t-1)
                    left = A-sum(top)
                    for j in range(1, t):
                        add = min(left, k-A-x-j-top[j-1])
                        assert add >= 0
                        top[j-1] += add
                        left -= add
                    assert left == 0
                    deficits = top+[min(x, (k-j-A)//2) for j in range(t, count+1)]
                    ds = [x+i+2 for i, x in enumerate(deficits)]+[i+2 for i in range(count, r)]
    return best, ds


def zero_good(ds, t, k):
    eligible = [i for i, d in enumerate(ds) if d >= i+2]
    for chosen in itertools.combinations(eligible, t):
        positive = [i for i in chosen if ds[i] > i+2]
        if positive and (sum(ds[i]-i-2 for i in chosen)
                         +max(ds[i] for i in positive) > k+1):
            return True
    return False

def oracle(n, mode):
    best = None
    data = None
    for k in range(2, n):
        if any(k%d == 0 for d in range(2, isqrt(k)+1)):
            continue
        m = n-k
        r = (m+k-1)//k
        if r < 3 or k <= r+1:
            continue
        H = sum(range(2, r+2))
        for t in range(2, r):
            old, previous = literal_dimensions(k, r, t)
            if mode == 'old':
                d = old
            elif mode == 'weighted':
                d = previous
            elif mode == 'exact':
                capacity = sum(k-j for j in range(2, t+1))
                if H+capacity+1 >= m:
                    continue
                lower = H+capacity+1
                B = sum(range(2, t+2))
                if best is not None and (4*k+8*lower-4*B-3, k, t) >= best:
                    continue
                d = previous if previous == H+capacity+1 else H+literal_integer_bound(k, r, t)[0]+1
            else:
                capacity = sum(k-j for j in range(2, t+1))
                lower = H+max(0, capacity-r+t-1)+1
                if lower >= m:
                    continue
                B = sum(range(2, t+2))
                if best is not None and (4*k+8*lower-4*B-3, k, t) >= best:
                    continue
                d = H+literal_zero_bound(k, r, t)[0]+1
            B = sum(range(2, t+2))
            if not (B <= d < m):
                continue
            key = (4*k+8*d-4*B-3, k, t)
            if best is None or key < best:
                best = key
                data = (k, m, r, t)
    return best, data



def constructor_oracle(n):
    best = (4*n, 0, 0, 0)
    result = None
    for k in range(2, n+1):
        if any(k % divisor == 0 for divisor in range(2, isqrt(k)+1)):
            continue
        top_count = (n+k-1)//k
        for count in (top_count, top_count-1):
            if count < 3 or k <= count+1:
                continue
            m = min(n, count*k)
            ignored = n-m
            H = sum(range(2, count+2))
            for t in range(2, count):
                B = sum(range(2, t+2))
                capacity = sum(k-j for j in range(2, t+1))
                lower = H+max(0, capacity-count+t-1)+1
                if lower >= m or (4*ignored+8*lower-4*B-3, k, t, ignored) >= best:
                    continue
                d = H+literal_zero_bound(k, count, t)[0]+1
                if not B <= d < m:
                    continue
                key = (4*ignored+8*d-4*B-3, k, t, ignored)
                if key < best:
                    best = key
                    result = k, m, count, t, ignored
    return best, result

LOWER_STEP = [(1, 1), (15, -1), (17, 1), (19, 1), (21, -1), (26, -1), (29, 1), (31, 1), (33, -1), (43, 1), (44, -1), (61, 1), (63, -1), (65, -1), (71, 1), (73, 1), (75, -1), (103, 1), (229, 1), (242, -1), (271, 1), (285, -1), (323, 1), (385, -2), (421, 1), (439, 1), (440, -1), (493, 1), (494, -1), (571, 1), (770, -1), (841, 1), (1155, -1), (1273, 1), (1891, 1), (1924, -1), (3001, 1), (3003, -1), (3080, -1), (3361, 1), (3781, 1), (4004, -1), (4861, 1), (5005, -1), (6511, 1), (6864, -1), (7981, 1), (7995, -1), (8821, 1), (8853, -1), (9931, 1), (10725, -1), (11933, 1), (11934, -1), (13441, 1), (13442, -1), (14923, 1), (17654, -1), (19951, 1)]
UPPER_STEP = [(1, 1), (15, -1), (17, 1), (19, 1), (21, -1), (22, -1), (23, 1), (26, -1), (29, 1), (35, -1), (41, 1), (45, -1), (47, 1), (52, -1), (59, 1), (65, -1), (66, -1), (67, 1), (78, -1), (79, 1), (107, 1), (135, -1), (210, 1), (275, -1), (289, 1), (385, -1), (418, -1), (419, 1), (521, 1), (585, -1), (629, 1), (795, -1), (839, 1), (936, -1), (1049, 1), (1144, -1), (1717, 1), (1925, -1), (2068, -1), (2081, 1), (3135, -1), (3149, 1), (5005, -1), (5039, 1), (7007, -1), (7349, 1), (8086, -1), (8087, 1), (9009, -1), (9011, 1), (10036, -1), (10079, 1), (12376, -1), (15107, 1), (16588, -1), (16589, 1), (18096, -1)]


def check_floor_scheme():
    positive = (1, 6, 10, 14, 105)
    negative = (2, 3, 5, 7, 11, 13, 385, 1001)
    period = 30030
    assert all(period % v == 0 for v in positive+negative)
    assert sum((Fraction(1, v) for v in positive), Fraction()) == sum((Fraction(1, v) for v in negative), Fraction())
    assert len(positive)+len(negative) == 13
    assert LOWER_STEP[0] == UPPER_STEP[0] == (1, 1)
    for rows in (LOWER_STEP, UPPER_STEP):
        assert len(dict(rows)) == len(rows)
        assert all(type(j) is int and type(v) is int and 1 <= j < period and v != 0 for j, v in rows)
        assert all(j >= 15 for j, v in rows if j != 1)
    lower, upper = dict(LOWER_STEP), dict(UPPER_STEP)
    p = q = 0
    minimum, maximum = 100, -100
    for j in range(period):
        p += lower.get(j, 0)
        q += upper.get(j, 0)
        M = sum(j//v for v in positive)-sum(j//v for v in negative)
        assert q <= M <= p
        minimum, maximum = min(minimum, M), max(maximum, M)
    assert (minimum, maximum) == (-1, 4)
    assert sum(lower.values()) == 4 and sum(upper.values()) == -1
    c, d = Fraction(96863, 100000), Fraction(105557, 100000)
    def signed_sums(rows):
        plus = sum((Fraction(v, j) for j, v in rows if j >= 2 and v > 0), Fraction())
        minus = sum((Fraction(v, j) for j, v in rows if j >= 2 and v < 0), Fraction())
        return plus, minus
    pP, nP = signed_sums(LOWER_STEP)
    pQ, nQ = signed_sums(UPPER_STEP)
    assert (1+nP)*c+pP*d <= 1
    assert pQ*c+(1+nQ)*d >= 1
    def root10(n):
        lo, hi = 0, 1
        while hi**10 <= n:
            hi *= 2
        while hi-lo > 1:
            mid = (lo+hi)//2
            if mid**10 <= n:
                lo = mid
            else:
                hi = mid
        assert lo >= 1 and lo**10 <= n < (lo+1)**10
        return lo
    for rows in (LOWER_STEP, UPPER_STEP):
        weight = sum((Fraction(abs(v), root10(j**9)) for j, v in rows if j >= 2), Fraction())
        assert weight < Fraction(189, 200)
    assert Fraction(39, 2)+400*Fraction(189, 200) < 400
    factors = {}
    for value, sign in [(v, -1) for v in positive]+[(v, 1) for v in negative]:
        n, prime = value, 2
        while prime*prime <= n:
            while n % prime == 0:
                factors[prime] = factors.get(prime, Fraction())+Fraction(sign, value)
                n //= prime
            prime += 1
        if n > 1:
            factors[n] = factors.get(n, Fraction())+Fraction(sign, value)
    expected = {2: Fraction(17, 105), 3: Fraction(11, 70), 5: Fraction(43, 462),
                7: Fraction(281, 4290), 11: Fraction(43, 455), 13: Fraction(6, 77)}
    assert factors == expected and all(v > 0 for v in factors.values())
    def log_lower(v):
        z = Fraction(v-1, v+1)
        return 2*sum(z**(2*i+1)/(2*i+1) for i in range(20))
    assert sum((factors[v]*log_lower(v) for v in factors), Fraction()) > Fraction(49, 50)
    assert sum((factors[v]*(v-1) for v in factors), Fraction()) < 5
    assert c*Fraction(49, 50) > Fraction(9, 10)
    ratio = Fraction(109, 100)
    assert ratio-d/c == Fraction(2367, 9686300)
    E0 = 400*(1+ratio)*2**198+4*2**111+884*2**74
    assert E0 < Fraction(21303, 96863000)*2**220
    a, b = Fraction(25289, 20000), Fraction(2756503, 2000000)
    sa, sb = Fraction(5622388283, 5000000000), Fraction(1467486957, 1250000000)
    assert b == ratio*a+Fraction(1, 1000000) and ratio/2**220 < Fraction(1, 1000000)
    assert sa*sa > a and sb*sb > b
    assert 1/(2*a*a)+sa < Fraction(143721, 100000)
    assert 1/(2*(ratio*a)**2)+sb < Fraction(143721, 100000)
    assert Fraction(25, 12)+a*sa/2 < 3 and Fraction(25, 13)+b*sb/2 < 3
    assert Fraction(143721, 100000)+Fraction(3, 2**146)+Fraction(3, 2**292) < Fraction(114977, 80000)
    a0, b0 = Fraction(129, 100), Fraction(135451, 100000)
    s0, s1 = Fraction(113578167, 100000000), Fraction(116383419, 100000000)
    assert s0*s0 > a0 and s1*s1 > b0
    assert 1/(2*a0*a0)+s0 < Fraction(3591, 2500)
    assert 1/(2*(Fraction(21, 20)*a0)**2)+s1 < Fraction(3591, 2500)
    assert Fraction(25, 12)+a0*s0/2 < 3 and Fraction(25, 13)+b0*s1/2 < 3
    assert Fraction(21, 20)/2**18 < Fraction(1, 100000)
    assert Fraction(387, 200)*2**220+1 < 2**221
    assert Fraction(3591, 2500)+Fraction(3, 4096)+Fraction(3, 2**24) < Fraction(114977, 80000)
    return {'period_cells': period, 'kernel_min': minimum, 'kernel_max': maximum,
            'lower_terms': len(LOWER_STEP), 'upper_terms': len(UPPER_STEP),
            'certificate_sha256': sha256(json.dumps([positive, negative, LOWER_STEP, UPPER_STEP], separators=(',', ':')).encode()).hexdigest()}

RATIONAL_SCHEMES = [{'kind': 'lower', 'denominator': 123, 'coefficients': [[1, 123], [2, -123], [3, -123], [5, -123], [6, 123], [7, -123], [10, 123], [11, -123], [13, -123], [14, 123], [15, 123], [21, 123], [22, 123], [26, 123], [28, -91], [30, -246], [33, 123], [35, 123], [39, 123], [42, -246], [44, -123], [52, -123], [55, 123], [65, 123], [66, -246], [70, -246], [77, 123], [78, -246], [84, 91], [91, 123], [105, -246], [110, -155], [130, -246], [132, 123], [140, -32], [143, 91], [154, -123], [156, 123], [165, -278], [182, -246], [195, -155], [210, 278], [231, -246], [260, 123], [273, -246], [286, -337], [308, -32], [330, 369], [364, 91], [385, -251], [390, 278], [429, -401], [455, -369], [462, 246], [546, 278], [572, 155], [660, -305], [715, -91], [770, 160], [780, -310], [858, 460], [910, 333], [924, -214], [1001, -337], [1092, -91], [1155, 342], [1365, 256], [1430, 91], [1540, -77], [1716, -460], [1820, -247], [2002, 273], [2145, 96], [2310, -475], [2730, -470], [2860, -359], [3003, 251], [4004, -210], [4290, -492], [4620, 215], [5005, 96], [5460, 251], [6006, -474], [8580, 492], [10010, -224], [12012, 288], [15015, -276], [20020, 101], [30030, 492], [60060, 295]], 'steps': [[1, 123], [17, 123], [19, 123], [23, 123], [29, 32], [37, 123], [41, 123], [44, -123], [47, 123], [56, -91], [59, 123], [73, 123], [83, 123], [84, -123], [90, -123], [97, 123], [101, 123], [103, 123], [120, -123], [139, 123], [140, -246], [149, 91], [163, 123], [165, -32], [193, 123], [208, -123], [210, -214], [229, 123], [251, 123], [294, -123], [330, -64], [347, 123], [349, 123], [383, 123], [385, -5], [409, 32], [420, -64], [450, -123], [468, -91], [503, 91], [532, -91], [607, 123], [613, 123], [643, 32], [647, 32], [659, 123], [660, -123], [769, 123], [780, -246], [823, 123], [829, 91], [858, -155], [907, 123], [910, -159], [924, -214], [977, 91], [980, -64], [1031, 123], [1033, 123], [1037, 32], [1039, 123], [1073, 64], [1087, 123], [1091, 123], [1140, -123], [1144, -59], [1155, -64], [1217, 123], [1223, 123], [1229, 123], [1247, 123], [1320, -123], [1364, -123], [1365, -113], [1517, 123], [1537, 123], [1609, 73], [1637, 32], [1699, 123], [1709, 32], [1760, -32], [1792, -91], [1794, -123], [1804, -10], [1949, 59], [1974, -123], [1980, -214], [2087, 123], [2089, 123], [2099, 91], [2100, -91], [2145, -246], [2273, 32], [2309, 64], [2413, 113], [2420, -5], [2693, 123], [2713, 118], [2730, -186], [2929, 10], [2939, 123], [2964, -32], [2990, -91], [3003, -118], [3080, -355], [3120, -32], [3343, 64], [3349, 123], [3401, 123], [3407, 123], [3413, 91], [3416, -91], [3432, -91], [3617, 123], [3623, 123], [3629, 123], [3637, 59], [3640, -406], [3822, -59], [3861, -64], [3990, -91], [4004, -151], [4223, 69], [4241, 27], [4247, 123], [4253, 123], [4259, 32], [4283, 123], [4523, 91], [4607, 60], [4619, 123], [4620, -391], [4680, -37], [5148, -214], [5429, 32], [5443, 32], [5445, -32], [5717, 91], [5720, -22], [6323, 32], [6329, 32], [6359, 32], [6367, 123], [6370, -36], [6380, -155], [6707, 96], [6719, 91], [6765, -32], [6776, -123], [7254, -106], [7619, 32], [7639, 123], [7643, 123], [7644, -59], [7700, -141], [8257, 27], [8293, 32], [8399, 91], [8447, 91], [8448, -59], [8580, -423], [9433, 37], [9437, 32], [9701, 123], [9719, 32], [9799, 123], [9827, 123], [9851, 32], [9853, 123], [10289, 77], [10709, 32], [10710, -72], [11000, -32], [11010, -123], [11032, -91], [11044, -32], [11088, -5], [12011, 50], [12012, -145], [12679, 24], [12829, 54], [12839, 123], [12863, 123], [12869, 123], [12870, -337], [13104, -5], [13640, -32], [13650, -63], [13860, -300], [14000, -91], [14891, 113], [14893, 123], [14897, 32], [14909, 123], [14947, 32], [14957, 123], [15652, -128], [15834, -70], [16484, -91], [16500, -50], [17639, 6], [17933, 74], [18259, 35], [18443, 64], [18467, 91], [18479, 123], [18480, -177], [18564, -91], [18590, -133], [19793, 92], [19937, 32], [19939, 123], [19949, 123], [19963, 32], [19997, 123], [20017, 123], [20020, -533], [21839, 15], [21840, -192], [23087, 33], [23593, 10], [23999, 27], [24007, 123], [24023, 123], [24853, 47], [24860, -32], [24882, -170], [26723, 84], [26737, 59], [26869, 64], [26879, 123], [28541, 96], [28543, 123], [28600, -126], [29480, -32], [29484, -5], [30030, -491], [30264, -59], [30316, -59], [30360, -123], [30415, -5], [31973, 79], [31979, 32], [32759, 4], [33427, 101], [33431, 123], [33653, 54], [33659, 32], [33823, 28], [33833, 123], [33877, 32], [34007, 105], [34297, 64], [34313, 123], [34320, -364], [36740, -16], [39373, 32], [39449, 59], [39461, 32], [39463, 123], [39467, 123], [39468, -396], [40040, -99], [42677, 114], [42689, 91], [42877, 23], [42881, 123], [42893, 123], [42899, 32], [42900, -501], [45427, 22], [45439, 91], [45443, 123], [45864, -155], [45890, -49], [48048, -415], [50483, 10], [50587, 37], [50593, 123], [50773, 27], [50789, 32], [50801, 32], [50803, 123], [50819, 123], [51479, 64], [51480, -131], [52612, -91], [52624, -14], [54285, -22], [54600, -240], [56056, -1], [59119, 24], [59123, 123], [60007, 18], [60023, 123], [60029, 123], [60037, 32], [60041, 123], [60043, 123], [60059, 123]]}, {'kind': 'upper', 'denominator': 12, 'coefficients': [[1, 12], [2, -12], [3, -12], [5, -12], [6, 12], [7, -12], [10, 12], [11, -12], [13, -12], [14, 12], [15, 12], [21, 12], [22, 12], [26, 12], [30, -24], [33, 12], [35, 12], [39, 12], [42, -24], [44, -12], [55, 12], [65, 12], [66, -24], [70, -24], [77, 12], [78, -24], [91, 12], [105, -36], [110, -36], [130, -24], [132, 12], [140, -12], [143, 12], [154, -12], [165, -24], [182, -24], [195, -36], [210, 36], [220, 12], [231, -36], [273, -24], [286, -36], [330, 48], [364, -12], [385, -36], [390, 48], [429, -24], [455, -24], [462, 36], [546, 24], [572, 12], [660, -36], [715, -24], [770, 35], [780, -24], [858, 24], [924, -24], [1001, -36], [1092, -24], [1155, 37], [1365, 35], [1430, 34], [1540, -23], [1716, -24], [1820, 15], [2002, 47], [2145, 12], [2310, -48], [2730, -34], [2860, -22], [3003, 13], [4004, -24], [4290, -48], [4620, 36], [5005, 24], [5460, 9], [6006, -48], [8580, 36], [10010, -20], [12012, 48], [15015, -22], [20020, 4], [30030, 48], [60060, -11]], 'steps': [[1, 12], [17, 12], [19, 12], [23, 12], [31, 12], [37, 12], [47, 12], [53, 12], [61, 12], [66, -12], [71, 12], [73, 12], [78, -12], [79, 12], [88, -12], [97, 12], [113, 12], [120, -12], [126, -12], [157, 12], [199, 12], [210, -36], [220, -12], [223, 12], [241, 12], [280, -12], [315, -12], [337, 12], [367, 12], [390, -12], [397, 12], [420, -12], [431, 12], [433, 12], [491, 12], [546, -24], [550, -12], [585, -12], [589, 12], [601, 12], [631, 12], [661, 12], [726, -12], [728, -12], [770, -1], [780, -24], [799, 12], [924, -24], [936, -12], [1003, 12], [1021, 12], [1031, 12], [1144, -12], [1155, -23], [1181, 12], [1213, 12], [1219, 12], [1277, 12], [1279, 12], [1365, -13], [1430, -14], [1471, 12], [1577, 12], [1597, 12], [1601, 12], [1717, 12], [1820, -45], [1849, 12], [1873, 12], [2002, -1], [2184, -36], [2190, -12], [2200, -24], [2207, 12], [2227, 12], [2251, 12], [2263, 12], [2341, 12], [2347, 12], [2381, 12], [2383, 12], [2530, -12], [2548, -12], [2574, -24], [2663, 1], [3003, -11], [3120, -24], [3150, -12], [3276, -12], [3301, 12], [3319, 1], [3510, -12], [3529, 12], [3533, 12], [3541, 12], [3718, -9], [3721, 11], [4004, -1], [4095, -1], [4097, 12], [4201, 10], [4550, -12], [4620, -23], [4621, 12], [4633, 12], [4950, -12], [5170, -12], [5171, 12], [5460, -60], [5461, 12], [5471, 12], [5479, 12], [5503, 1], [5777, 12], [5839, 12], [6053, 12], [6435, -12], [6490, -12], [6552, -24], [6600, -12], [6637, 11], [6943, 1], [7722, -12], [8190, -11], [8580, -13], [8581, 12], [8587, 12], [9130, -8], [9156, -12], [9240, -11], [9241, 12], [9253, 12], [9257, 12], [9473, 11], [10013, 12], [10033, 12], [10081, 12], [10091, 1], [10920, -23], [11235, -12], [11237, 12], [12012, -60], [12013, 12], [12101, 12], [12107, 11], [12391, 12], [12497, 4], [13170, -3], [13200, -12], [13860, -22], [14196, -1], [14213, 12], [14303, 2], [14633, 12], [14701, 12], [14711, 3], [15481, 9], [16695, -11], [17160, -60], [17161, 12], [17167, 1], [17383, 12], [17473, 12], [17557, 10], [18041, 12], [18097, 2], [19030, -12], [19320, -24], [19740, -23], [19741, 9], [21340, -11], [21341, 10], [21583, 12], [22932, -12], [23323, 12], [23333, 2], [25200, -24], [25213, 12], [25217, 12], [25307, 3], [26600, -1], [26620, -12], [26623, 12], [26629, 12], [26633, 11], [27301, 12], [28393, 12], [29820, -22], [29835, -12], [30030, -7], [30576, -12], [30580, -12], [31200, -8], [31243, 12], [31249, 10], [32341, 11], [33181, 11], [35178, -9], [35200, -12], [35201, 12], [35207, 12], [36037, 12], [36053, 12], [36131, 1], [38220, -11], [39160, -1], [39186, -12], [39795, -1], [40040, -26], [41581, 12], [41593, 12], [41597, 2], [43680, -47], [45103, 12], [49504, -8], [49507, 1], [49771, 12], [49781, 12], [49783, 12], [49787, 12], [52416, -20], [53690, -24], [57552, -9]]}]

def check_rational_schemes():
    period = 60060
    def logs(v):
        z = Fraction(v-1, v+1)
        a = 2*sum((z**(2*i+1)/(2*i+1) for i in range(60)), Fraction())
        return a, a+2*z**121/(121*(1-z*z))
    lb = {p: logs(p) for p in (2,3,5,7,11,13)}
    def root10(n):
        lo, hi = 0, 1
        while hi**10 <= n:hi *= 2
        while hi-lo > 1:
            m = (lo+hi)//2
            if m**10 <= n:lo=m
            else:hi=m
        assert lo >= 1 and lo**10 <= n < (lo+1)**10
        return lo
    report=[]
    L,U=Fraction(977463,1000000),Fraction(1022562,1000000)
    for scheme in RATIONAL_SCHEMES:
        kind, den = scheme['kind'], scheme['denominator']
        coeff, steps = scheme['coefficients'], scheme['steps']
        assert kind in ('lower','upper') and type(den) is int and den>0
        assert coeff[0] == [1,den] and steps[0] == [1,den]
        assert len(dict(coeff))==len(coeff) and len(dict(steps))==len(steps)
        assert all(type(d) is int and type(v) is int and d>0 and period%d==0 and v!=0 for d,v in coeff)
        assert all(type(j) is int and type(v) is int and 1<=j<=period and v!=0 for j,v in steps)
        assert all(steps[i][0]<steps[i+1][0] for i in range(len(steps)-1))
        assert sum((Fraction(v,d) for d,v in coeff), Fraction()) == 0
        jumps=dict(steps);value=0;minimum=10**9;maximum=-10**9
        for j in range(period):
            value += jumps.get(j,0)
            M=sum(v*(j//d) for d,v in coeff)
            assert value>=M if kind=='lower' else value<=M
            minimum,maximum=min(minimum,M),max(maximum,M)
        value += jumps.get(period,0)
        assert value>=maximum if kind=='lower' else value<=minimum
        factors={p:Fraction() for p in lb}
        for d,v in coeff:
            n=d
            for p in lb:
                exponent=0
                while n%p==0:n//=p;exponent+=1
                factors[p]-=Fraction(v*exponent,den*d)
            assert n==1
        A0=sum((w*lb[p][int(w<0)] for p,w in factors.items()),Fraction())
        A1=sum((w*lb[p][int(w>=0)] for p,w in factors.items()),Fraction())
        plus=sum((Fraction(v,den*j) for j,v in steps if j>1 and v>0),Fraction())
        minus=sum((Fraction(v,den*j) for j,v in steps if j>1 and v<0),Fraction())
        error=Fraction(3,2)*sum(Fraction(abs(v),den) for d,v in coeff)
        error+=2000*sum((Fraction(abs(v),den*root10(j**9)) for j,v in steps if j>1),Fraction())
        assert error<2000
        if kind=='lower':assert (1+minus)*L+plus*U < A0
        else:assert plus*L+(1+minus)*U > A1
        report.append({'kind':kind,'kernel_terms':len(coeff),'step_terms':len(steps),'period_cells':period})
    L,U=Fraction(977463,1000000),Fraction(1022562,1000000)
    q=Fraction(5231,5000)
    assert L*q-U>0
    E=2000*(1+q)*2**243+4*2**136+1084*2**91
    assert E<(L*q-U)*2**270
    a=Fraction(6451,5000);b=q*a+Fraction(1,1000000)
    def sqrt_upper(x):
        scale=10**10;n=x.numerator*scale*scale//x.denominator
        return Fraction(isqrt(n)+1,scale)
    sa,sb=sqrt_upper(a),sqrt_upper(b)
    assert sa*sa>a and sb*sb>b
    assert 1/(2*a*a)+sa < Fraction(35906,25000)
    assert 1/(2*(q*a)**2)+sb < Fraction(35906,25000)
    assert Fraction(25,12)+a*sa/2<3 and Fraction(25,13)+b*sb/2<3
    assert q/2**270<Fraction(1,1000000)
    assert Fraction(35906,25000)+Fraction(3,2**180)+Fraction(3,2**360)<Fraction(1149,800)
    a0=Fraction(6501,5000);ratio=Fraction(103,100);b0=ratio*a0+Fraction(1,100000)
    s0,s1=sqrt_upper(a0),sqrt_upper(b0)
    assert s0*s0>a0 and s1*s1>b0
    assert 1/(2*a0*a0)+s0<Fraction(35901,25000)
    assert 1/(2*(ratio*a0)**2)+s1<Fraction(35901,25000)
    assert Fraction(25,12)+a0*s0/2<3 and Fraction(25,13)+b0*s1/2<3
    assert ratio/2**18 < Fraction(1,100000)
    assert 2**14 == 16384
    assert Fraction(35901,25000)+Fraction(3,16384)+Fraction(3,2**28)<Fraction(1149,800)
    assert a0*2**270+1<2**271
    return {'schemes':report,'certificate_sha256':sha256(json.dumps(RATIONAL_SCHEMES,separators=(',',':')).encode()).hexdigest()}

def generate_interval_cover():
    def fifth_ceil(n):
        if n <= 1:
            return n
        x = 1 << ((n.bit_length()+4)//5)
        while True:
            y = (4*x+n//(x**4))//5
            if y >= x:
                return x+(x**5 < n)
            x = y
    sieve=bytearray(b'\x01')*3000001;sieve[:2]=b'\x00\x00'
    for p in range(2,isqrt(3000000)+1):
        if sieve[p]:sieve[p*p::p]=b'\x00'*(((3000000-p*p)//p)+1)
    primes=[p for p in range(5,3000001) if sieve[p]]
    def dims(k,r,t):
        T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
        return r*(r+3)//2+max(T,U//(t+1))+1,t*(t+3)//2
    
    def valid(n,cost):return (10000*cost)**5<=114900**5*n**4
    
    def choose(n):
        root=fifth_ceil(n**3);middle=bisect_left(primes,(132*root)//100)
        best=None;result=None
        for k in primes[max(0,middle-2):middle+2]:
            r0=(n+k-1)//k
            for r in (r0,r0-1):
                if r<3 or k<=r+1:continue
                m=min(n,r*k);ignored=n-m;q=isqrt(r)
                for t in sorted({2,q,q+1}):
                    if not 2<=t<r:continue
                    d,B=dims(k,r,t)
                    if not B<=d<m:continue
                    loss=4*ignored+8*d-4*B-3
                    key=(loss,k,t,ignored)
                    if best is None or key<best:best=key;result=k,r,t,d,B
        assert result is not None and valid(n,best[0]),(n,best,result)
        return result

    rows = []
    current = 114900**5//40000**5+1
    limit = 2**35-1
    while current <= limit:
        k, r, t, dimension, degree = choose(current)
        base = 8*dimension-4*degree-3
        def cost(n):
            return base+4*max(0, n-r*k)
        maximum = min(limit, r*k+k)
        assert (r-1)*k < current <= maximum
        assert dimension < min(current, r*k)
        if valid(maximum, cost(maximum)):
            end = maximum
        else:
            lo, hi = current, maximum
            while hi-lo > 1:
                middle = (lo+hi)//2
                if valid(middle, cost(middle)):
                    lo = middle
                else:
                    hi = middle
            end = lo
        assert valid(current, cost(current)) and valid(end, cost(end))
        rows.append((current, end, k, r, t))
        current = end+1
    return rows


def check_interval_cover(f, rows):
    assert 40000**5*195 <= 114900**5 < 40000**5*196
    expected_left = 196
    endpoints = 0
    # Independent trial division uses a fresh tiny sieve, not the generator's sieve.
    maximum = isqrt(max(row[2] for row in rows))
    small = bytearray(b'\x01')*(maximum+1)
    small[:2] = b'\x00\x00'
    for divisor in range(2, isqrt(maximum)+1):
        if small[divisor]:
            small[divisor*divisor::divisor] = b'\x00'*((maximum-divisor*divisor)//divisor+1)
    divisors = [q for q in range(2, maximum+1) if small[q]]
    def checked_prime(n):
        if n < 2:
            return False
        for q in divisors:
            if q*q > n:
                break
            if n % q == 0:
                return False
        return True
    for left, right, k, r, t in rows:
        assert all(type(value) is int for value in (left, right, k, r, t))
        assert left == expected_left and left <= right
        assert checked_prime(k) and k >= 5
        assert r >= 3 and k > r+1 and 2 <= t < r
        H = (r+1)*(r+2)//2-1
        T = sum(k-j for j in range(2, t+1))
        U = r*(2*k-r-1)//2
        d_ref = H+max(T, U//(t+1))+1
        degree = sum(range(2, t+2))
        assert left > (r-1)*k and right <= r*k+k
        assert degree <= d_ref < min(left, r*k)
        base = 8*d_ref-4*degree-3
        assert base >= 0
        for endpoint in (left, right):
            cost = base+4*max(0, endpoint-r*k)
            assert (10000*cost)**5 <= 114900**5*endpoint**4
            endpoints += 1
        expected_left = right+1
    assert rows and expected_left == 2**35
    return {'intervals': len(rows), 'endpoint_checks': endpoints,
            'first': rows[0], 'last': rows[-1],
            'row_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()}

def chain_trial(n):
 return n>=2 and (n==2 or n%2 and all(n%d for d in range(3,isqrt(n)+1,2)))
def chain_jacobi(a,n):
 result=1;a%=n
 while a:
  while a%2==0:
   a//=2
   if n%8 in (3,5):result=-result
  a,n=n,a
  if a%4==n%4==3:result=-result
  a%=n
 return result if n==1 else 0
def chain_next_trial(n):
 n=max(3,n|1)
 while not chain_trial(n):n+=2
 return n

def chain_proth(start,limit):
 m=(start.bit_length()+1)//2
 while True:
  a=(start-1+(1<<m)-1)//(1<<m)
  if a%2==0:a+=1
  if a < 1<<m:break
  m+=1
 while (p:=a*(1<<m)+1)<=limit and a<1<<m:
  if gcd(p, 6541380665835015) != 1:
   a+=2
   continue
  for b in range(2,65):
   symbol=chain_jacobi(b,p)
   if symbol==-1:
    if pow(b,(p-1)//2,p)==p-1:return [p,m,a,b]
    break
  a+=2
 return None

def generate_prime_chain():
 low=2**18;high=2**271;p=chain_next_trial(low);rows=[[p,0,0,0]]
 while p<high:
  start=(102*p+99)//100;limit=103*p//100
  if p>=2**30:start=(1029*p+999)//1000
  row=None
  if p>=2**30:
   row=chain_proth(start,limit)
   if row is None:row=chain_proth(p+1,limit)
  if row is None:
   assert p<2**32,('no Proth prime',p,start,limit)
   q=chain_next_trial(start);assert q<=limit,(p,q,limit);row=[q,0,0,0]
  assert p<row[0] and 100*row[0]<=103*p
  rows.append(row);p=row[0]
 return rows

def check_prime_chain(rows):
    assert rows and 2**18 <= rows[0][0] and 100*rows[0][0] <= 103*2**18
    small_rows = [row[0] for row in rows if row[1] == 0]
    bound = isqrt(max(small_rows))
    sieve = bytearray(b'\x01')*(bound+1)
    sieve[:2] = b'\x00\x00'
    for v in range(2, isqrt(bound)+1):
        if sieve[v]:
            sieve[v*v::v] = b'\x00'*((bound-v*v)//v+1)
    divisors = [v for v in range(2, bound+1) if sieve[v]]
    def independent_trial(n):
        if n < 2:
            return False
        for v in divisors:
            if v*v > n:
                break
            if n % v == 0:
                return False
        return True
    previous = None
    proths = trials = 0
    for row in rows:
        assert len(row) == 4 and all(type(v) is int for v in row)
        p, m, a, b = row
        if m:
            assert m >= 1 and 0 < a < 2**m and a % 2 == 1
            assert p == a*2**m+1 and 1 < b < p
            assert pow(b, (p-1)//2, p) == p-1
            proths += 1
        else:
            assert a == b == 0 and p < 2**32 and independent_trial(p)
            trials += 1
        if previous is not None:
            assert previous < p and 100*p <= 103*previous
        previous = p
    assert previous >= 2**271
    return {'primes': len(rows), 'proth_certificates': proths, 'trial_primes': trials,
            'first': rows[0], 'last': rows[-1],
            'row_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()}

def run():
    f = load_family()
    floor_scheme = check_floor_scheme()
    rational_schemes = check_rational_schemes()
    interval_cover = check_interval_cover(f, generate_interval_cover())
    prime_chain = check_prime_chain(generate_prime_chain())
    prefixes = 0
    for k in (2, 3, 5, 7):
        for low in range(1 << k):
            q = (1 << k) | low
            reducible = any(f.poly_mod(q, (1 << d) | a) == 0
                            for d in range(1, k//2+1)
                            for a in range(1 << d))
            assert f.irreducible_prime_degree(q, k) == (not reducible)
            prefixes += 1

    factor_cases = 0
    for degree in range(1, 9):
        for low in range(1 << degree):
            polynomial = (1 << degree) | low
            derivative = sum(((polynomial >> j) & 1) << (j-1) for j in range(1, degree+1, 2))
            if f.poly_gcd(polynomial, derivative) != 1:
                continue
            factor = f.factor_one(polynomial)
            assert f.poly_mod(polynomial, factor) == 0
            fd = factor.bit_length()-1
            assert fd >= 1
            assert not any(f.poly_mod(factor, (1 << j) | q) == 0
                           for j in range(1, fd//2+1) for q in range(1 << j))
            factor_cases += 1

    constructed_fields = []
    for k in (3, 5, 7, 11, 13, 17, 19, 23, 31, 43):
        polynomial = f.irreducible_prime(k)
        assert polynomial.bit_length() == k+1
        assert f.irreducible_prime_degree(polynomial, k)
        constructed_fields.append({'k': k, 'polynomial': polynomial})

    compressed_polynomials = 0
    compressed_irreducibles = []
    for k in (5, 7, 11, 13, 17):
        seen = set()
        valid = 0
        for encoded in range(1 << (k-5)):
            bits = [(encoded >> j) & 1 for j in range(k-5)]
            q = f.decode_prefix(k, bits)
            assert q not in seen
            seen.add(q)
            # Independent literal binary polynomial product with Z^2(Z+1)^2.
            R = (1 << (k-4)) | encoded
            product = 0
            for i in range(k-3):
                if (R >> i) & 1:
                    for j in (2, 4):
                        product ^= 1 << (i+j)
            assert q == product ^ 1
            assert q.bit_length() == k+1 and q & 3 == 1
            assert (q >> (k-1)) & 1 == 0
            assert q.bit_count() % 2 == 1
            assert sum((q >> j) & 1 for j in range(1, k+1, 2)) % 2 == 0
            valid += f.irreducible_prime_degree(q, k)
            compressed_polynomials += 1
        assert valid > 0
        compressed_irreducibles.append({'k': k, 'valid': valid})

    zero_minor_cases = 0
    rng_minor = Random(20261107)
    for k, q in ((3, 11), (5, 37), (7, 131)):
        assert f.irreducible_prime_degree(q, k)
        for h in range(2, k):
            for _ in range(8):
                basis = [1 << j for j in range(k)]
                for _ in range(4*k):
                    i, j = rng_minor.sample(range(k), 2)
                    basis[i] ^= basis[j]
                rng_minor.shuffle(basis)
                matrix = [[f.field_pow(v, 1 << u, q) for v in basis[:h]] for u in range(h)]
                for col in range(h):
                    pivot = next(i for i in range(col, h) if matrix[i][col])
                    matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
                    inverse = f.field_pow(matrix[col][col], (1 << k)-2, q)
                    matrix[col] = [f.field_mul(x, inverse, q) for x in matrix[col]]
                    for i in range(col+1, h):
                        factor = matrix[i][col]
                        matrix[i] = [x ^ f.field_mul(factor, y, q) for x, y in zip(matrix[i], matrix[col])]
                zero_minor_cases += 1

    distributions = 0
    for m in range(3, 9):
        for t in range(2, m):
            h = next(h for h in range(2, m+1) if h*(t+1) >= 2*m)
            high = tuple(range(m-h, m))
            low = tuple(range(m-h))
            lam = Fraction(h*(t+1), m)-2
            assert 0 <= lam <= 1 and len(low) >= t-1
            coefficients = [Fraction(0) for _ in range(m)]
            marked_weight = Fraction(0)
            mass = Fraction(0)
            for extra_high, case_mass in ((0, 1-lam), (1, lam)):
                if not case_mass:
                    continue
                for j in high:
                    for hs in itertools.combinations([i for i in high if i != j], extra_high):
                        for ls in itertools.combinations(low, t-1-extra_high):
                            weight = case_mass/h/comb(h-1, extra_high)/comb(len(low), t-1-extra_high)
                            chosen = (j,)+hs+ls
                            assert len(chosen) == len(set(chosen)) == t
                            mass += weight
                            marked_weight += weight*(j+2)
                            for i in chosen:
                                coefficients[i] += weight
                            coefficients[j] += weight
            assert mass == 1
            assert coefficients == [Fraction(t+1, m)]*m
            assert marked_weight == Fraction(sum(i+2 for i in high), h)
            distributions += 1

    dimension_cases = 0
    attaining_vectors = 0
    for r in range(3, 21):
        for k in (r+2, 2*r, 2*r+1, 2*r+3):
            for t in range(2, r):
                old, previous = literal_dimensions(k, r, t)
                d, B = f.certificate_data(k, r, t)
                assert B == sum(range(2, t+2)) and d <= previous <= old
                if r <= 10:
                    assert d <= sum(range(2, r+2))+f.isolation_deficit_bound(k, r, t)+1
                if r <= 10:
                    E, ds = literal_zero_bound(k, r, t, witness=True)
                    assert d == sum(range(2, r+2))+E+1
                    assert sum(ds) == d-1 and all(0 <= x <= k for x in ds)
                    assert not zero_good(ds, t, k)
                    attaining_vectors += 1
                dimension_cases += 1

    vectors = 0
    maxima_candidates = 0
    exact_maxima = 0
    for r in range(3, 6):
        k = 2*r+1
        for t in range(2, r):
            actual = -1
            for ds in itertools.product(range(k+1), repeat=r):
                vectors += 1
                total = sum(ds)
                if total <= actual:
                    continue
                maxima_candidates += 1
                good = zero_good(ds, t, k)
                if not good:
                    actual = total
            assert f.certificate_data(k, r, t)[0] == actual+1
            exact_maxima += 1

    comparisons = []
    for n, expected in ((500, (1461, 1437, 1429, 1413)), (1000, (2533, 2493, 2485, 2485)),
                        (2000, (4513, 4449, 4433, 4417)), (10000, (16809, 16601, 16577, 16521))):
        old, _ = oracle(n, 'old')
        previous, _ = oracle(n, 'weighted')
        exact, _ = oracle(n, 'exact')
        new, params = oracle(n, 'zero')
        assert (old[0], previous[0], exact[0], new[0]) == expected
        assert f.reference_parameters(n) == params
        k, _, r, t = params
        d, B = f.certificate_data(k, r, t)
        assert new == (4*k+8*d-4*B-3, k, t)
        comparisons.append({'N': n, 'old_loss': old[0], 'weighted_loss': previous[0], 'positive_exact_loss': exact[0], 'new_loss': new[0]})

    compressed_comparisons = []
    for n, expected in ((500, 1393), (1000, 2465), (2000, 4397), (10000, 16501)):
        best = (4*n, 0, 0, 0)
        params = None
        for saved in range(6):
            refkey, refdata = oracle(n+saved, 'zero')
            if refdata is None:
                continue
            k, m, r, t = refdata
            key = (refkey[0]-4*saved, k, t, 5-saved)
            if key < best:
                best = key
                params = k, m, r, t, 5-saved
        assert best[0] == expected and f.compressed_parameters(n) == params
        compressed_comparisons.append({'N': n, 'loss': best[0], 'parameters': params})

    constructor_comparisons = []
    for n, expected in ((500, 1257), (1000, 2265), (2000, 4061), (10000, 15365)):
        key, params = constructor_oracle(n)
        assert key[0] == expected and f.parameters(n) == params
        constructor_comparisons.append({'N': n, 'loss': key[0], 'parameters': params})

    prefix_reduction_cases = 0
    for n in range(10, 81):
        actual = (4*n, 0, 0, 0)
        for k in range(2, n+1):
            if not f.is_prime(k):
                continue
            for ignored in range(min(k, n)+1):
                m = n-ignored
                r = (m+k-1)//k
                if r < 3 or k <= r+1:
                    continue
                H = sum(range(2, r+2))
                for t in range(2, r):
                    dimension = H+literal_zero_bound(k, r, t)[0]+1
                    B = sum(range(2, t+2))
                    if B <= dimension < m:
                        actual = min(actual, (4*ignored+8*dimension-4*B-3, k, t, ignored))
                    prefix_reduction_cases += 1
        assert constructor_oracle(n)[0] == actual

    recurrence_cases = 0
    for values in itertools.product(range(8), repeat=3):
        q = 0b1011
        e = [1, 0, 0]
        for value in values:
            for j in (2, 1):
                e[j] ^= f.field_mul(value, e[j-1], q)
        literal = 0
        for a, b in itertools.combinations(values, 2):
            literal ^= f.field_mul(a, b, q)
        assert e[2] == literal
        recurrence_cases += 1

    assert f.parameters(94) == (31, 93, 3, 2, 1)
    k, m, r, t, ignored = f.parameters(94)
    q = f.irreducible_prime(k)
    assert f.irreducible_prime_degree(q, k)
    prefix = [0]*ignored
    rng = Random(20261106)
    for _ in range(16):
        payload = [rng.randrange(2) for _ in range(m)]
        padded = payload+[0]*(r*k-m)
        ys = []
        for i in range(r):
            value = sum(bit << j for j, bit in enumerate(padded[i*k:(i+1)*k]))
            ys.append(f.field_pow(value, (1 << (r+1-i))-1, q))
        expected = 0
        for subset in itertools.combinations(ys, t):
            product = 1
            for value in subset:
                product = f.field_mul(product, value, q)
            expected ^= product
        assert f.evaluate(prefix+payload) == (expected & 1)
        assert f.evaluate([1]*ignored+payload) == (expected & 1)
    for n in range(25):
        assert f.evaluate([0]*n) == 0
    return {'all_passed': True, 'floor_scheme': floor_scheme, 'rational_schemes': rational_schemes, 'interval_cover': interval_cover, 'prime_chain': prime_chain, 'constructed_fields': constructed_fields, 'factor_cases': factor_cases,
            'constructor_comparisons': constructor_comparisons, 'prefix_reduction_cases': prefix_reduction_cases, 'compressed_polynomials': compressed_polynomials,
            'compressed_irreducibles': compressed_irreducibles, 'compressed_comparisons': compressed_comparisons, 'prime_degree_prefixes': prefixes,
            'zero_deficit_minor_cases': zero_minor_cases, 'marked_distributions': distributions, 'dimension_cases': dimension_cases,
            'exhaustive_dimension_vectors': vectors, 'maxima_candidates': maxima_candidates,
            'exact_maximum_cases': exact_maxima, 'attaining_vectors': attaining_vectors,
            'global_comparisons': comparisons, 'symmetric_recurrences': recurrence_cases,
            'valid_prefix_evaluator_cases': 16, 'zero_inputs': 25}


if __name__ == '__main__':
    json.load(__import__('sys').stdin)
    print(json.dumps(run(), sort_keys=True))
