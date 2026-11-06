"""Exact bounded supplementary predicates; universal claims rest on written proofs."""
import importlib.util
import sys
sys.dont_write_bytecode = True
import itertools
from bisect import bisect_left
from hashlib import sha256
import json
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt, gcd, prod
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

def load_certificates():
    path=Path('/snapshot/files/certificates.py')
    if not path.exists():path=Path(__file__).with_name('certificates.py')
    spec=importlib.util.spec_from_file_location('candidate_certificates',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.RATIONAL_SCHEMES

RATIONAL_SCHEMES=load_certificates()

def check_rational_schemes():
    def logs(v):
        z = Fraction(v-1, v+1)
        a = 2*sum((z**(2*i+1)/(2*i+1) for i in range(60)), Fraction())
        return a, a+2*z**121/(121*(1-z*z))
    lb = {p: logs(p) for p in (2,3,5,7,11,13,17)}
    def root8(n):
        lo, hi = 0, 1
        while hi**8 <= n:hi *= 2
        while hi-lo > 1:
            m = (lo+hi)//2
            if m**8 <= n:lo=m
            else:hi=m
        assert lo >= 1 and lo**8 <= n < (lo+1)**8
        return lo
    report=[]
    L,U=Fraction(986117,1000000),Fraction(1013811,1000000)
    assert [o['kind'] for o in RATIONAL_SCHEMES]==['lower','upper']
    for scheme in RATIONAL_SCHEMES:
        period=scheme['period']
        assert type(period) is int and period in (3063060,)
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
        M=0;boundary_checks=0
        for left in range(0,period,65536):
            right=min(period,left+65536)
            if left:
                assert M==sum(v*((left-1)//d) for d,v in coeff)
                boundary_checks+=1
            differences=[0]*(right-left)
            for d,v in coeff:
                first=max(d,((left+d-1)//d)*d)
                for j in range(first,right,d):differences[j-left]+=v
            for offset,delta in enumerate(differences):
                j=left+offset
                value+=jumps.get(j,0);M+=delta
                assert value>=M if kind=='lower' else value<=M
                minimum,maximum=min(minimum,M),max(maximum,M)
            assert M==sum(v*((right-1)//d) for d,v in coeff)
            boundary_checks+=1
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
        error+=18000*sum((Fraction(abs(v),den*root8(j**7)) for j,v in steps if j>1),Fraction())
        assert error<18000
        if kind=='lower':assert (1+minus)*L+plus*U < A0
        else:assert plus*L+(1+minus)*U > A1
        report.append({'kind':kind,'kernel_terms':len(coeff),'step_terms':len(steps),'period_cells':period,'direct_boundary_checks':boundary_checks})
    L,U=Fraction(986117,1000000),Fraction(1013811,1000000)
    q=Fraction(20563,20000)
    assert L*q-U>0
    E=18000*(1+q)*2**210+4*2**121+964*2**81
    assert E<(L*q-U)*2**240
    a=Fraction(2602763,2000000);b=q*a+Fraction(1,100000000)
    def sqrt_upper(x):
        scale=10**20;n=x.numerator*scale*scale//x.denominator
        return Fraction(isqrt(n)+1,scale)
    sa,sb=sqrt_upper(a),sqrt_upper(b)
    assert sa*sa>a and sb*sb>b
    assert 1/(2*a*a)+sa < Fraction(14360113,10000000)
    assert 1/(2*(q*a)**2)+sb < Fraction(14360113,10000000)
    assert Fraction(25,12)+a*sa/2<3 and Fraction(25,13)+b*sb/2<3
    assert q/2**240<Fraction(1,100000000)
    assert 8*Fraction(14360113,10000000)<Fraction(11488091,1000000)
    a0=Fraction(130157859,100000000);ratio=Fraction(51,50);b0=ratio*a0+Fraction(1,4000)
    s0,s1=sqrt_upper(a0),sqrt_upper(b0)
    assert s0*s0>a0 and s1*s1>b0
    assert 1/(2*a0*a0)+s0<Fraction(14360083,10000000)
    assert 1/(2*(ratio*a0)**2)+s1<Fraction(14360083,10000000)
    assert Fraction(25,12)+a0*s0/2<3 and Fraction(25,13)+b0*s1/2<3
    assert ratio/2**17 < Fraction(1,100000)
    assert 2**14 == 16384
    assert 8*Fraction(14360083,10000000)<Fraction(11488091,1000000)
    assert Fraction(4,3)**5>4
    assert a0*Fraction(4,3)*2**240+1<2**241
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
    sieve=bytearray(b'\x01')*1000001;sieve[:2]=b'\x00\x00'
    for p in range(2,isqrt(1000000)+1):
        if sieve[p]:sieve[p*p::p]=b'\x00'*(((1000000-p*p)//p)+1)
    primes=[p for p in range(5,1000001) if sieve[p]]
    def dims(k,r,t):
        T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
        return r*(r+3)//2+max(T,U//(t+1))+1,t*(t+3)//2
    
    def valid(n,cost):return (1000000*cost)**5<=11488091**5*n**4
    
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
    current = 11488091**5//4000000**5+1
    limit = 2**28-1
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
    assert 4000000**5*195 <= 11488091**5 < 4000000**5*196
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
            assert (1000000*cost)**5 <= 11488091**5*endpoint**4
            endpoints += 1
        expected_left = right+1
    assert rows and expected_left == 2**28
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

CHAIN_WHEEL=prod(p for p in range(3,998,2) if chain_trial(p))

def chain_proth(start,limit):
 m=(start.bit_length()+1)//2
 while True:
  a=(start-1+(1<<m)-1)//(1<<m)
  if a%2==0:a+=1
  if a < 1<<m:break
  m+=1
 while (p:=a*(1<<m)+1)<=limit and a<1<<m:
  if gcd(p, CHAIN_WHEEL) != 1:
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
 low=2**12;high=2**241;p=chain_next_trial(low);rows=[[p,0,0,0]]
 while p<high:
  start=(1019*p+999)//1000;limit=51*p//50
  if p>=2**30:start=(1019*p+999)//1000
  if p>=2**64:start=(10199*p+9999)//10000
  row=None
  if p>=2**30:
   row=chain_proth(start,limit)
   if row is None:row=chain_proth(p+1,limit)
  if row is None:
   assert p<2**32,('no Proth prime',p,start,limit)
   q=chain_next_trial(start)
   if q>limit:q=chain_next_trial(p+1)
   assert q<=limit,(p,q,limit);row=[q,0,0,0]
  assert p<row[0] and 50*row[0]<=51*p
  rows.append(row);p=row[0]
 return rows

def check_prime_chain(rows):
    assert rows and 2**12 <= rows[0][0] and 50*rows[0][0] <= 51*2**12
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
            assert previous < p and 50*p <= 51*previous
        previous = p
    assert previous >= 2**241
    return {'primes': len(rows), 'proth_certificates': proths, 'trial_primes': trials,
            'first': rows[0], 'last': rows[-1],
            'row_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()}


def check_nearest_rounding():
    cases=0
    for r in range(3,513):
        q=isqrt(r);t=q+(4*r>(2*q+1)**2)
        assert 2<=t<r
        for k in range(r+2,r+17):
            T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
            E=max(T,U//(t+1))
            assert 3*E+k>=0 and (3*E+k)**2<=9*k*k*r
            cases+=1
    assert Fraction(43,100)-Fraction(3,64)-Fraction(3,2**18)>Fraction(3,8)
    sharper_cases=0
    for r in range(2049,2305):
        q=isqrt(r);t=q+(4*r>(2*q+1)**2)
        for k in range(r+2,r+17):
            T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
            E=max(T,U//(t+1));z=91*E+45*k
            assert z>=0 and z*z<=91*91*k*k*r
            sharper_cases+=1
    assert sharper_cases==3840
    assert Fraction(45,91)*Fraction(129,100)-Fraction(3,64)-Fraction(3,2**18)>Fraction(59,100)
    assert 8*Fraction(14360113,10000000)<Fraction(11488091,1000000)
    assert 8*Fraction(14360083,10000000)<Fraction(11488091,1000000)
    return {'original_cases':cases,'sharper_regime_cases':sharper_cases}

def check_biased_rounding():
    cases=0
    for r in range(2049,2305):
        j=isqrt(r);k0=isqrt(r**3);k0+=k0*k0<r**3
        for k in range(k0,k0+15):
            es=[]
            for t in (j,j+1):
                assert 2<=t<r and k>r+1
                T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
                es.append(max(T,U//(t+1)))
            E=min(es);left=4*r*E+2*r*k-r*r;right=4*r*k-r*r+2*k
            assert left>=0 and right>0 and left*left<=right*right*r
            cases+=1
    assert (1+Fraction(1,4096))**3<Fraction(129,100)**2
    assert Fraction(7,6)**2>Fraction(27,20)
    assert Fraction(8,7)**2>Fraction(129,100)
    assert Fraction(275,129)+Fraction(63,40)<4
    assert Fraction(129,200)+Fraction(175,1032)-Fraction(4,64)-Fraction(13,4*64**3)>Fraction(3,4)
    assert Fraction(6639300809,5000000000)<Fraction(27,20)
    assert Fraction(20563,20000)*Fraction(2602763,2000000)+Fraction(1,100000000)<Fraction(27,20)
    return cases

def check_refined_thresholds():
    v0=Fraction(649,500);su=Fraction(2278596059,2000000000)
    assert su*su>v0
    K0=v0/2+1/(4*v0*su);G=Fraction(275,129)+Fraction(63,40)
    assert 55**5<2**29
    assert Fraction(10,9)+Fraction(2,55**2)<Fraction(9,10)*55
    assert (1+Fraction(1,55**2))**3<Fraction(129,100)**2
    assert K0-G/55-Fraction(13,4*55**3)>Fraction(3,4)
    assert K0-Fraction(4,8192)-Fraction(13,4*8192**3)>Fraction(327,400)
    assert Fraction(130157859,100000000)>v0 and Fraction(2602763,2000000)>v0
    cases=0
    for r in range(1513,1769):
        j=isqrt(r);k0=isqrt(r**3);k0+=k0*k0<r**3
        for k in range(k0,k0+15):
            es=[]
            for t in (j,j+1):
                assert 2<=t<r and k>r+1
                T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
                es.append(max(T,U//(t+1)))
            E=min(es);left=4*r*E+2*r*k-r*r;right=4*r*k-r*r+2*k
            assert left>=0 and right>0 and left*left<=right*right*r
            cases+=1
    return cases

def check_biased_threshold():
    lo=Fraction(13,10);hi=Fraction(67,50);u=Fraction(57009,50000)
    assert u*u>lo and Fraction(579,500)**2>hi
    K=lo/2+1/(4*lo*u);G=Fraction(55,26)+Fraction(116379,100000)
    assert 48**5<2**28
    assert Fraction(10,9)+Fraction(2,48**2)<Fraction(9,10)*48
    assert (1+Fraction(1,48**2))**3<Fraction(129,100)**2
    assert K-G/48-Fraction(13,4*48**3)>Fraction(3,4)
    assert K-Fraction(4,8192)-Fraction(13,4*8192**3)>Fraction(1309,1600)
    assert 8*K>Fraction(6549,1000)
    assert 8*G-Fraction(100,67)<25
    assert Fraction(6549,1000)-Fraction(25,48)-Fraction(23,48**3)>6
    assert Fraction(6549,1000)-Fraction(25,8192)-Fraction(23,2**39)>Fraction(1309,200)
    a=Fraction(2602763,2000000);q=Fraction(20563,20000)
    a0=Fraction(130157859,100000000);ratio=Fraction(51,50)
    assert a>lo and a0>lo
    assert q*a+Fraction(1,100000000)<hi
    assert ratio*a0+Fraction(1,100000)<hi
    assert ratio/48**3<Fraction(1,100000)
    assert a0*48**3>2**17
    cases=0
    for r in range(1153,1409):
        j=isqrt(r);k0=isqrt(r**3);k0+=k0*k0<r**3
        for k in range(k0,k0+15):
            es=[]
            for t in (j,j+1):
                assert 2<=t<r and k>r+1 and t*(t+3)>=r
                T=(t-1)*k-t*(t+1)//2+1;U=r*k-r*(r+1)//2
                es.append(max(T,U//(t+1)))
            E=min(es);left=4*r*E+2*r*k-r*r;right=4*r*k-r*r+k
            assert left>=0 and right>0 and left*left<=right*right*r
            cases+=1
    return cases

def check_shifted_choice():
    lo=Fraction(13,10);hi=Fraction(67,50);c=Fraction(13,50)
    assert lo**5/(1+hi/48**2)**3>Fraction(19,10)**2
    assert hi**5<Fraction(21,10)**2
    assert Fraction(1,2)-1/(4*Fraction(21,10))+c/33==Fraction(1497,3850)<Fraction(39,100)
    assert Fraction(19,10)*(Fraction(39,100)**2-2*c)+Fraction(39,200)+Fraction(1,2)==-Fraction(401,100000)
    G3=Fraction(20,13)+Fraction(19,25)*hi*Fraction(579,500)
    assert G3==Fraction(22081871,8125000)
    assert 8*G3-Fraction(100,67)<Fraction(81,4)
    cases=0
    for r in range(1153,1409):
        j=isqrt(r)
        for multiplier in range(192,209):
            k=isqrt(r**3*multiplier**2)//100
            assert 100*k*k>361*r**3 and 100*k*k<441*r**3
            E=min(max((t-1)*k-t*(t+1)//2+1,(r*k-r*(r+1)//2)//(t+1)) for t in (j,j+1))
            left=100*E+50*k+50*r;right=100*k*r-25*r*r+26*k
            assert left>0 and right>0 and left*left*r<=right*right
            cases+=1
    return cases

def check_earlier_continuous():
    lo=Fraction(13,10);hi=Fraction(67,50);ratio=Fraction(51,50)
    assert 36**5<2**26 and 42**5<2**27
    assert 36**2/hi>30**2
    assert lo**5/(1+hi/36**2)**3>Fraction(19,10)**2
    assert hi**5<Fraction(21,10)**2
    assert Fraction(1,2)-1/(4*Fraction(21,10))+Fraction(13,50)/30<Fraction(39,100)
    assert Fraction(1,2)+1/(4*Fraction(19,10))<Fraction(3,4)
    assert 1/lo+Fraction(2,36**2)<lo*36
    assert Fraction(130157859,100000000)*36**3>2**15
    assert ratio/36**3<Fraction(1,40000)
    assert Fraction(6549,1000)-Fraction(81,4*42)-Fraction(21,42**3)>6
    assert Fraction(14360113,10000000)+Fraction(22081871,8125000*36**2)+Fraction(3,36**4)<2
    return True

def check_twenty_continuous():
    lo=Fraction(13,10);hi=Fraction(67,50);c=Fraction(7,25)
    assert 16**5==2**20 and 16**2/hi>13**2
    assert lo**5/(1+hi/16**2)**3>Fraction(19,10)**2
    assert hi**5<Fraction(21,10)**2
    assert Fraction(1,2)-1/(4*Fraction(21,10))+c/13<Fraction(41,100)
    assert Fraction(19,10)*(Fraction(41,100)**2-2*c)+Fraction(41,200)+Fraction(1,2)==-Fraction(3961,100000)
    assert 2/lo+Fraction(39,50)*lo*Fraction(57009,50000)<Fraction(2703,1000)
    assert 2/hi+Fraction(39,50)*hi*Fraction(579,500)<Fraction(2703,1000)
    assert 8*Fraction(2703,1000)-Fraction(100,67)<Fraction(81,4)
    assert Fraction(51,50)/16**3<Fraction(1,4000)
    assert Fraction(130157859,100000000)*16**3>2**12
    assert 1/lo+Fraction(2,16**2)<lo*16
    assert Fraction(14360113,10000000)+Fraction(2703,1000*16**2)+Fraction(3,16**4)<2
    cases=0
    for r in range(192,448):
        j=isqrt(r)
        for multiplier in range(192,209):
            k=isqrt(r**3*multiplier**2)//100
            assert 100*k*k>361*r**3 and 100*k*k<441*r**3
            E=min(max((t-1)*k-t*(t+1)//2+1,(r*k-r*(r+1)//2)//(t+1)) for t in (j,j+1))
            left=100*E+50*k+50*r;right=100*k*r-25*r*r+28*k
            assert left>0 and right>0 and left*left*r<=right*right
            cases+=1
    return cases

def check_sharp_constants():
    a=Fraction(2602763,2000000);b=Fraction(20563,20000)*a+Fraction(1,100000000)
    assert b==Fraction(53520615969,40000000000)
    assert Fraction(130157859,100000000)>a
    assert Fraction(6639300809,5000000000)<b
    def upper(v):
        scale=10**10
        return Fraction(isqrt(v.numerator*scale**2//v.denominator)+1,scale)
    ua,ub=upper(a),upper(b)
    assert ua==Fraction(11407810921,10**10) and ua**2>a and ub**2>b
    Kstar=a/2+1/(4*a*ua);Gstar=Fraction(2678755,1000000)
    assert 8*Kstar>Fraction(655269,100000)
    assert 2/a+Fraction(153,200)*a*ua<Gstar
    assert 2/b+Fraction(153,200)*b*ub<Gstar
    assert 8*Gstar-2/b<Fraction(997,50)
    assert Fraction(655269,100000)-Fraction(997,50*8192)-Fraction(21,2**39)>Fraction(26201,4000)
    c=Fraction(53,200);bshift=Fraction(2007,5000)
    assert Fraction(1,2)-1/(4*Fraction(21,10))+c/13<bshift
    assert Fraction(19,10)*(bshift**2-2*c)+bshift/2+Fraction(1,2)==-Fraction(42069,250000000)
    assert c*Fraction(19,10)-Fraction(1,2)==Fraction(7,2000)
    cases=0
    for r in range(192,448):
        j=isqrt(r)
        for multiplier in range(192,209):
            k=isqrt(r**3*multiplier**2)//100
            E=min(max((t-1)*k-t*(t+1)//2+1,(r*k-r*(r+1)//2)//(t+1)) for t in (j,j+1))
            left=200*E+100*k+100*r;right=200*k*r-50*r*r+53*k
            assert left>0 and right>0 and left*left*r<=right*right
            cases+=1
    return cases

def check_enclosing_constants():
    a=Fraction(2602763,2000000);b=Fraction(53520615969,40000000000)
    wl=Fraction(19169,10000);wh=Fraction(20709,10000)
    c=Fraction(21,80);cap=Fraction(9987,25000)
    assert wl**2<a**5/(1+b/256)**3 and b**5<wh**2
    assert Fraction(1,2)-1/(4*wh)+c/13<cap
    assert wl*(cap**2-2*c)+cap/2+Fraction(1,2)<0 and c*wl>Fraction(1,2)
    def root(z):
        u=Fraction(isqrt(z.numerator*10**20//z.denominator)+1,10**10)
        assert (u-Fraction(1,10**10))**2<=z<u**2
        return u
    ua,ub=root(a),root(b);G=Fraction(1337443,500000)
    assert 2/a+Fraction(61,80)*a*ua<G
    assert 2/b+Fraction(61,80)*b*ub<G
    K=Fraction(3276349,500000);D=Fraction(995217,50000)
    assert 8*(a/2+1/(4*a*ua))>K
    assert 8*G-2/b<D
    assert K-D/8192-Fraction(21,2**39)>Fraction(655026,100000)
    return {'all_passed':True,'K':str(K),'G':str(G),'D':str(D),'shift':str(c)}


def check_division_remainder():
    wl=Fraction(19169,10000);lo=Fraction(36958,100000)
    hi=Fraction(9987,25000);h=Fraction(133,200)
    assert lo<Fraction(1,2)-1/(4*wl)
    assert (1+lo**2)/2+wl*lo**3>h and h*hi<1
    b=Fraction(53520615969,40000000000)
    ub=Fraction(isqrt(b.numerator*10**20//b.denominator)+1,10**10)
    assert ub**2>b and 8*h/ub>Fraction(45991,10000)
    assert Fraction(7,10)*13**2>h and Fraction(7,4)*13>h
    cases=0
    for r in range(170,426):
        for mult in range(19170,20709,96):
            k=isqrt(r**3*mult**2)//10000
            if not wl**2*r**3<k*k<Fraction(20709,10000)**2*r**3:continue
            j=isqrt(r)
            E=min(max((t-1)*k-t*(t+1)//2+1,(r*k-r*(r+1)//2)//(t+1)) for t in (j,j+1))
            L=400*E+200*k+200*r-400
            R=400*k*r-100*r*r+105*k-266*r
            assert L>0 and R>0 and L*L*r<=R*R
            cases+=1
    return {'all_passed':True,'cases':cases,'h':str(h)}


def check_late_shift():
    a=Fraction(2602763,2000000);b=Fraction(53520615969,40000000000)
    wl=Fraction(19169,10000);c=Fraction(247,1000)
    assert 8192**2>7000**2*b
    A=Fraction(1,2)+c/7000
    assert A*A-2*c+1/(16*wl*wl)<0
    assert wl*(A*A-2*c)+Fraction(1,2)-1/(16*wl)<0
    assert Fraction(1,2)-1/(4*Fraction(20709,10000))+c/7000<Fraction(2,5)
    G=Fraction(26508961,10000000)
    for v in (a,b):
        u=Fraction(isqrt(v.numerator*10**20//v.denominator)+1,10**10)
        assert u*u>v and 2/v+Fraction(747,1000)*v*u<G
    D=Fraction(985621,50000)
    assert 8*G-2/b<D
    assert Fraction(3276349,500000)-D/8192-Fraction(29,2**39)>Fraction(655029,100000)
    return {'all_passed':True,'shift':str(c),'D':str(D)}


def check_adaptive_shift():
    wl=Fraction(19169,10000);wh=Fraction(20709,10000);h=Fraction(191,1000)
    def c(w):return Fraction(1,8)+1/(4*w)-1/(32*w*w)
    def H(w):return w/8+Fraction(7,16)-1/(32*w)+1/(64*w*w)
    assert 0<c(wh)<c(wl)<Fraction(1,4)
    assert Fraction(3,16)+1/(64*wh*wh)>h
    assert H(Fraction(21,10))<Fraction(7,10)
    assert Fraction(7,10)*Fraction(2,5)+Fraction(21,160)<1
    assert Fraction(1,2)-1/(4*wh)+Fraction(1,4*13)<Fraction(2,5)
    a=Fraction(2602763,2000000);b=Fraction(53520615969,40000000000)
    G=Fraction(264437261,100000000)
    for v in (a,b):
        u=Fraction(isqrt(v.numerator*10**20//v.denominator)+1,10**10)
        assert u*u>v
        assert 9/(4*v)+Fraction(5,8)*v*u-1/(32*v**3*u)<G
    ub=u
    assert 8*G-2/b<Fraction(1966023,100000)
    assert 8*h/ub>Fraction(4128,3125)
    assert Fraction(24411,100000)*16>Fraction(45991,10000)-Fraction(4128,3125)
    assert Fraction(5219,100000)*8192>Fraction(45991,10000)-Fraction(4128,3125)
    return True


def check_shifted_cut():
    wl=Fraction(19169,10000);wh=Fraction(20709,10000);z=Fraction(339,1000)
    assert 1/(4*wh)-Fraction(1,52)-z/(wl*13**2)>0
    assert Fraction(2,5)+z/(wl*13**2)<Fraction(1,2)
    b=Fraction(53520615969,40000000000)
    ub=Fraction(isqrt(b.numerator*10**20//b.denominator)+1,10**10)
    assert Fraction(106,25)/ub>Fraction(366551,100000)
    K=Fraction(3276349,500000);D=Fraction(1966023,100000);P=Fraction(366551,100000)
    assert K*16>D and P*16>37
    assert K-D/40-Fraction(37,40**3)>6
    assert K-D/8192-Fraction(37,2**39)>Fraction(655029,100000)
    cases=0
    for r in range(170,426):
        for mult in range(19170,20709,96):
            k=isqrt(r**3*mult**2)//10000
            if not wl**2*r**3<k*k<wh**2*r**3:continue
            j=isqrt(r)
            E=min(max((t-1)*k-t*(t+1)//2+1,(r*k-r*(r+1)//2)//(t+1)) for t in (j,j+1))
            L=800*k*E+400*k*k+200*k*r-1600*k
            R=800*k*k*r-200*k*r*r+100*k*k-25*r**3-424*k*r
            assert L>0 and R>0 and L*L*r<=R*R
            cases+=1
    return cases


def run():
    shifted_cut = check_shifted_cut()
    adaptive_shift = check_adaptive_shift()
    late_shift = check_late_shift()
    division_remainder = check_division_remainder()
    enclosing_constants = check_enclosing_constants()
    sharp_constants = check_sharp_constants()
    twenty_continuous = check_twenty_continuous()
    earlier_continuous = check_earlier_continuous()
    shifted_choice = check_shifted_choice()
    f = load_family()
    nearest_rounding = check_nearest_rounding()
    biased_rounding = check_biased_rounding()
    refined_rounding = check_refined_thresholds()
    biased_threshold = check_biased_threshold()
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
    return {'all_passed': True, 'shifted_cut': shifted_cut, 'adaptive_shift': adaptive_shift, 'late_shift': late_shift, 'division_remainder': division_remainder, 'enclosing_constants': enclosing_constants, 'sharp_constants': sharp_constants, 'twenty_continuous_cases': twenty_continuous, 'earlier_continuous_guards': earlier_continuous, 'shifted_choice_cases': shifted_choice, 'biased_threshold_cases': biased_threshold, 'refined_threshold_cases': refined_rounding, 'biased_rounding_cases': biased_rounding, 'nearest_rounding_cases': nearest_rounding, 'floor_scheme': floor_scheme, 'rational_schemes': rational_schemes, 'interval_cover': interval_cover, 'prime_chain': prime_chain, 'constructed_fields': constructed_fields, 'factor_cases': factor_cases,
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
