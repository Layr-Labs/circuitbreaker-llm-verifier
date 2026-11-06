"""Exact bounded supplementary predicates; universal claims rest on written proofs."""
import importlib.util
import itertools
from bisect import bisect_left
from hashlib import sha256
import json
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt
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
    sieve=bytearray(b'\x01')*1000001;sieve[:2]=b'\x00\x00'
    for p in range(2,isqrt(1000000)+1):
        if sieve[p]:sieve[p*p::p]=b'\x00'*(((1000000-p*p)//p)+1)
    primes=[p for p in range(5,1000001) if sieve[p]]
    def dims(k,r,t):
        T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
        return r*(r+3)//2+max(T,U//(t+1))+1,t*(t+3)//2
    
    def valid(n,cost):return (10000*cost)**5<=115213**5*n**4
    
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
    current = 115213**5//40000**5+1
    limit = 2**30-1
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
    assert 40000**5*198 <= 115213**5 < 40000**5*199
    expected_left = 199
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
            assert (10000*cost)**5 <= 115213**5*endpoint**4
            endpoints += 1
        expected_left = right+1
    assert rows and expected_left == 2**30
    M = lambda v: v-v//2-v//3-v//5+v//30
    assert all(M(v) in (0, 1) for v in range(30))
    assert all(M(v) == 1 for v in range(1, 6))
    assert all(M(v) == 1 for v in range(7, 10))
    assert all(M(v) == 0 for v in range(24, 29))
    ell, upper = Fraction(51072, 50999), Fraction(59595, 50999)
    assert ell == 1+ell/24-upper/29
    assert upper == 1+upper/6-ell/7+upper/10
    assert upper/ell == Fraction(19865, 17024)
    for d, cap in [(24, 21), (29, 19), (6, 41), (7, 38), (10, 32)]:
        assert Fraction(cap, 100)**2*d > 1
    assert Fraction(15, 2)+32*Fraction(21, 100)+80*Fraction(19, 100) < 32
    assert Fraction(15, 2)+80*Fraction(41, 100)+32*Fraction(38, 100)+80*Fraction(32, 100) < 80
    def log_lower(q, count):
        z = Fraction(q-1, q+1)
        return 2*sum(z**(2*j+1)/(2*j+1) for j in range(count))
    assert (14*log_lower(2, 3)+9*log_lower(3, 3)+5*log_lower(5, 5))/30 > Fraction(9, 10)
    q = Fraction(11669, 10000)
    assert q-upper/ell == Fraction(191, 10640000)
    assert 11863284**2 > 2**47
    assert 12815092**2 > q*2**47 and 54763**3 > q*2**47
    assert 1420920528 < Fraction(1719, 106400000)*2**47
    a, b = Fraction(30583, 25000), Fraction(1427492109, 1000000000)
    sa, sb = Fraction(11060379741, 10000000000), Fraction(11947770123, 10000000000)
    assert b == q*a+Fraction(1, 10**9) and q/2**47 < Fraction(1, 10**9)
    assert sa*sa > a and sb*sb > b
    assert 1/(2*a*a)+sa < Fraction(28803, 20000)
    assert 1/(2*(q*a)**2)+sb < Fraction(28803, 20000)
    assert Fraction(25, 12)+a*sa/2 < 3 and Fraction(25, 14)+b*sb/2 < 3
    assert Fraction(28803, 20000)+Fraction(3, 2**31)+Fraction(3, 2**62) < Fraction(115213, 80000)
    a0, b0 = Fraction(33, 25), Fraction(133321, 100000)
    s0, s1 = Fraction(1148913, 1000000), Fraction(115465, 100000)
    assert s0*s0 > a0 and s1*s1 > b0
    assert 1/(2*a0*a0)+s0 < Fraction(359, 250)
    assert 1/(2*Fraction(13332, 10000)**2)+s1 < Fraction(359, 250)
    assert Fraction(25, 12)+a0*s0/2 < 3 and Fraction(25, 13)+b0*s1/2 < 3
    assert Fraction(101, 100*2**18) < Fraction(1, 100000)
    assert Fraction(99, 50)*2**47+1 < 2**48
    assert Fraction(359, 250)+Fraction(3, 4096)+Fraction(3, 2**24) < Fraction(115213, 80000)
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
  for b in range(2,65):
   symbol=chain_jacobi(b,p)
   if symbol==-1:
    if pow(b,(p-1)//2,p)==p-1:return [p,m,a,b]
    break
  a+=2
 return None

def generate_prime_chain():
 low=2**18;high=2**48;p=chain_next_trial(low);rows=[[p,0,0,0]]
 while p<high:
  start=(1008*p+999)//1000;limit=101*p//100
  row=None
  if p>=2**30:
   row=chain_proth(start,limit)
   if row is None:row=chain_proth(p+1,limit)
  if row is None:
   assert p<2**32,('no Proth prime',p,start,limit)
   q=chain_next_trial(start);assert q<=limit,(p,q,limit);row=[q,0,0,0]
  assert p<row[0] and 100*row[0]<=101*p
  rows.append(row);p=row[0]
 return rows



def check_prime_chain(rows):
    assert rows and 2**18 <= rows[0][0] and 100*rows[0][0] <= 101*2**18
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
            assert previous < p and 100*p <= 101*previous
        previous = p
    assert previous >= 2**48
    return {'primes': len(rows), 'proth_certificates': proths, 'trial_primes': trials,
            'first': rows[0], 'last': rows[-1],
            'row_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()}

def run():
    f = load_family()
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
    return {'all_passed': True, 'interval_cover': interval_cover, 'prime_chain': prime_chain, 'constructed_fields': constructed_fields, 'factor_cases': factor_cases,
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
