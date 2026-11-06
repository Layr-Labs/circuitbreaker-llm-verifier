"""Exact bounded supplementary predicates; universal claims rest on written proofs."""
import importlib.util
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt
from pathlib import Path
from random import Random


FINITE_LOSS_COVER = [(299, 365, 73, 4, 2), (366, 423, 47, 8, 3), (424, 477, 53, 8, 3), (478, 531, 59, 8, 3), (532, 603, 67, 8, 3), (604, 670, 67, 9, 3), (671, 747, 83, 8, 3), (748, 830, 83, 9, 3), (831, 927, 103, 8, 3), (928, 1030, 103, 9, 3), (1031, 1143, 127, 8, 3), (1144, 1270, 127, 9, 3), (1271, 1413, 157, 8, 3), (1414, 1570, 157, 9, 3), (1571, 1773, 197, 8, 3), (1774, 1921, 113, 16, 4), (1922, 2032, 127, 15, 4), (2033, 2159, 127, 16, 4), (2160, 2329, 137, 16, 4), (2330, 2512, 157, 15, 4), (2513, 2669, 157, 16, 4), (2670, 2839, 167, 16, 4), (2840, 3043, 179, 16, 4), (3044, 3247, 191, 16, 4), (3248, 3438, 191, 17, 4), (3439, 3664, 229, 15, 4), (3665, 3893, 229, 16, 4), (3894, 4097, 241, 16, 4), (4098, 4369, 257, 16, 4), (4370, 4688, 293, 15, 4), (4689, 4981, 293, 16, 4), (4982, 5287, 311, 16, 4), (5288, 5486, 211, 25, 5), (5487, 5725, 229, 24, 5), (5726, 5954, 229, 25, 5), (5955, 6214, 239, 25, 5), (6215, 6526, 251, 25, 5), (6527, 6838, 263, 25, 5), (6839, 7202, 277, 25, 5), (7203, 7618, 293, 25, 5), (7619, 7982, 307, 25, 5), (7983, 8242, 317, 25, 5), (8243, 8606, 331, 25, 5), (8607, 9022, 347, 25, 5), (9023, 9334, 359, 25, 5), (9335, 9698, 373, 25, 5), (9699, 10114, 389, 25, 5), (10115, 10525, 421, 24, 5), (10526, 10946, 421, 25, 5), (10947, 11414, 439, 25, 5), (11415, 11882, 457, 25, 5), (11883, 12454, 479, 25, 5), (12455, 12974, 499, 25, 5), (12975, 13546, 521, 25, 5), (13547, 14066, 541, 25, 5), (14067, 14638, 563, 25, 5), (14639, 15133, 409, 36, 6), (15134, 15577, 421, 36, 6), (15578, 16021, 433, 36, 6), (16022, 16613, 449, 36, 6), (16614, 17057, 461, 36, 6), (17058, 17518, 461, 37, 6), (17519, 18019, 487, 36, 6), (18020, 18611, 503, 36, 6), (18612, 19277, 521, 36, 6), (19278, 20017, 541, 36, 6), (20018, 20609, 557, 36, 6), (20610, 21127, 571, 36, 6), (21128, 21719, 587, 36, 6), (21720, 22459, 607, 36, 6), (22460, 23066, 607, 37, 6), (23067, 23717, 641, 36, 6), (23718, 24383, 659, 36, 6), (24384, 25049, 677, 36, 6), (25050, 25937, 701, 36, 6), (25938, 26603, 719, 36, 6), (26604, 27343, 739, 36, 6), (27344, 28157, 761, 36, 6), (28158, 29119, 787, 36, 6), (29120, 29933, 809, 36, 6), (29934, 30549, 599, 50, 7), (30550, 31263, 613, 50, 7), (31264, 32050, 641, 49, 7), (32051, 32691, 641, 50, 7), (32692, 33303, 653, 50, 7), (33304, 34150, 683, 49, 7), (34151, 34833, 683, 50, 7), (34834, 35751, 701, 50, 7), (35752, 36650, 733, 49, 7), (36651, 37383, 733, 50, 7), (37384, 38301, 751, 50, 7), (38302, 39350, 787, 49, 7), (39351, 40137, 787, 50, 7), (40138, 41050, 821, 49, 7), (41051, 41950, 839, 49, 7), (41951, 42850, 857, 49, 7), (42851, 43850, 877, 49, 7), (43851, 44727, 877, 50, 7), (44728, 45550, 911, 49, 7), (45551, 46450, 929, 49, 7), (46451, 47350, 947, 49, 7), (47351, 48350, 967, 49, 7), (48351, 49550, 991, 49, 7), (49551, 50650, 1013, 49, 7), (50651, 51650, 1033, 49, 7), (51651, 53050, 1061, 49, 7), (53051, 54350, 1087, 49, 7), (54351, 55450, 1109, 49, 7), (55451, 56559, 1109, 50, 7), (56560, 57650, 1153, 49, 7), (57651, 59050, 1181, 49, 7), (59051, 60231, 1181, 50, 7), (60232, 61450, 1229, 49, 7), (61451, 62502, 947, 65, 8), (62503, 63850, 1277, 49, 7), (63851, 64878, 983, 65, 8), (64879, 65802, 997, 65, 8), (65803, 66858, 1013, 65, 8), (66859, 68046, 1031, 65, 8), (68047, 69234, 1049, 65, 8), (69235, 70554, 1069, 65, 8), (70555, 71742, 1087, 65, 8), (71743, 72798, 1103, 65, 8), (72799, 74118, 1123, 65, 8), (74119, 75595, 1163, 64, 8), (75596, 76758, 1163, 65, 8), (76759, 77946, 1181, 65, 8), (77947, 79266, 1201, 65, 8), (79267, 80718, 1223, 65, 8), (80719, 82434, 1249, 65, 8), (82435, 83785, 1289, 64, 8), (83786, 85074, 1289, 65, 8), (85075, 86262, 1307, 65, 8), (86263, 87582, 1327, 65, 8), (87583, 89245, 1373, 64, 8), (89246, 90618, 1373, 65, 8), (90619, 92495, 1423, 64, 8), (92496, 94055, 1447, 64, 8), (94056, 95615, 1471, 64, 8), (95616, 97045, 1493, 64, 8), (97046, 98538, 1493, 65, 8), (98539, 100295, 1543, 64, 8), (100296, 101855, 1567, 64, 8), (101856, 103805, 1597, 64, 8), (103806, 105365, 1621, 64, 8), (105366, 106986, 1621, 65, 8), (106987, 108485, 1669, 64, 8), (108486, 110305, 1697, 64, 8), (110306, 111995, 1723, 64, 8), (111996, 113461, 1367, 82, 9), (113462, 115505, 1777, 64, 8), (115506, 116947, 1409, 82, 9), (116948, 118441, 1427, 82, 9), (118442, 120101, 1447, 82, 9), (120102, 122093, 1471, 82, 9), (122094, 123587, 1489, 82, 9), (123588, 125413, 1511, 82, 9), (125414, 127073, 1531, 82, 9), (127074, 128567, 1549, 82, 9), (128568, 130393, 1571, 82, 9), (130394, 132551, 1597, 82, 9), (132552, 134377, 1619, 82, 9), (134378, 135871, 1637, 82, 9), (135872, 137531, 1657, 82, 9), (137532, 139318, 1699, 81, 9), (139319, 141017, 1699, 82, 9), (141018, 142843, 1721, 82, 9), (142844, 144503, 1741, 82, 9), (144504, 146534, 1787, 81, 9), (146535, 148321, 1787, 82, 9), (148322, 150313, 1811, 82, 9), (150314, 151973, 1831, 82, 9), (151974, 153914, 1877, 81, 9), (153915, 155791, 1877, 82, 9), (155792, 157783, 1901, 82, 9), (157784, 159818, 1949, 81, 9), (159819, 161767, 1949, 82, 9), (161768, 163759, 1973, 82, 9), (163760, 165751, 1997, 82, 9), (165752, 168241, 2027, 82, 9), (168242, 170399, 2053, 82, 9), (170400, 172723, 2081, 82, 9), (172724, 175213, 2111, 82, 9), (175214, 177371, 2137, 82, 9), (177372, 179363, 2161, 82, 9), (179364, 181466, 2213, 81, 9), (181467, 183679, 2213, 82, 9), (183680, 186058, 2269, 81, 9), (186059, 188354, 2297, 81, 9), (188355, 190651, 2297, 82, 9), (190652, 192678, 1889, 101, 10), (192679, 194567, 1889, 102, 10), (194568, 196962, 1931, 101, 10), (196963, 199002, 1951, 101, 10), (199003, 201246, 1973, 101, 10), (201247, 203286, 1993, 101, 10), (203287, 205734, 2017, 101, 10), (205735, 207978, 2039, 101, 10), (207979, 210426, 2063, 101, 10), (210427, 212466, 2083, 101, 10), (212467, 214549, 2083, 102, 10), (214550, 217158, 2129, 101, 10), (217159, 219606, 2153, 101, 10), (219607, 222258, 2179, 101, 10), (222259, 224706, 2203, 101, 10), (224707, 226909, 2203, 102, 10), (226910, 229602, 2251, 101, 10), (229603, 231846, 2273, 101, 10), (231847, 234294, 2297, 101, 10), (234295, 236591, 2297, 102, 10), (236592, 239394, 2347, 101, 10), (239395, 241842, 2371, 101, 10), (241843, 244698, 2399, 101, 10), (244699, 247146, 2423, 101, 10), (247147, 249594, 2447, 101, 10), (249595, 252246, 2473, 101, 10), (252247, 255306, 2503, 101, 10), (255307, 258162, 2531, 101, 10), (258163, 260814, 2557, 101, 10), (260815, 262143, 2591, 101, 10)]

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
    def paired_prefix(n):
        return (n//2)*((n+1)//2) if n >= 0 else 0
    return ((end-t+1)*x+paired_prefix(k-A-end-1)
            -paired_prefix(k-A-count-1))


@lru_cache(None)
def literal_integer_bound(k, r, t, witness=False):
    best = sum(k-j for j in range(2, t+1))
    best_dims = [k if i < t-1 else i+2 for i in range(r)]
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

def run():
    f = load_family()
    cover_left = 299
    cover_endpoints = 0
    for left, right, k, r, t in FINITE_LOSS_COVER:
        assert left == cover_left and left <= right
        assert all(k % q for q in range(2, isqrt(k)+1)) and k >= 5
        assert r >= 3 and k > r+1 and 2 <= t < r
        H = sum(range(2, r+2))
        T = sum(k-j for j in range(2, t+1))
        U = sum(k-j for j in range(1, r+1))
        d_ref = H+max(T, U//(t+1))+1
        B = sum(range(2, t+2))
        assert left > (r-1)*k and right <= r*k+k
        assert B <= d_ref < min(left, r*k)
        A = 8*d_ref-4*B-3
        assert A >= 0
        for endpoint in (left, right):
            cost = A+4*max(0, endpoint-r*k)
            assert (2*cost)**5 <= 25**5*endpoint**4
            cover_endpoints += 1
        cover_left = right+1
    assert len(FINITE_LOSS_COVER) == 220 and cover_left == 2**18
    assert 8**5*298 <= 25**5

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
    return {'all_passed': True, 'finite_cover_intervals': len(FINITE_LOSS_COVER),
            'finite_cover_endpoint_checks': cover_endpoints, 'finite_cover_last': cover_left-1, 'constructed_fields': constructed_fields, 'factor_cases': factor_cases,
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
