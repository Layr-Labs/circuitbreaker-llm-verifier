"""Uniform evaluator for the deterministic-field zero-deficit family in proof.md.

This module constructs the selected prime-degree field deterministically.
Integers encode binary polynomials with bit j the coefficient of Z^j.
"""
from math import isqrt
from functools import lru_cache


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def old_certificate_data(k, r, t):
    """Exact original positive-deficit dimension and degree."""
    capacity = (t-1)*k-t*(t+1)//2+1
    total = r*k-r*(r+1)//2
    dimension = r*(r+3)//2+max(capacity, total//(t+1))+1
    return dimension, t*(t+3)//2


def weighted_certificate_data(k, r, t):
    """Refine the dimension by weighted marked-subset averaging when guarded."""
    old, degree = old_certificate_data(k, r, t)
    if k < 2*r+1 or 2*k < (t+1)**2:
        return old, degree
    marked = (2*r+t)//(t+1)
    averaged = r*(2*k-2*r-1+marked)//(2*(t+1))
    capacity = (t-1)*k-t*(t+1)//2+1
    dimension = r*(r+3)//2+max(capacity, averaged)+1
    return dimension, degree


def isolation_deficit_bound(k, r, t):
    """Exact largest positive-deficit sum failing the original isolation test."""
    best = (t-1)*k-t*(t+1)//2+1
    for count in range(t, r+1):
        for top_sum in range(t-1, k-count-1):
            pivot = min(
                top_sum//(t-1),
                (k-t-top_sum)//2,
                ((t-1)*k-t*(t-1)//2-t*top_sum)//(t-1))
            if pivot < 1:
                continue
            residual = k-top_sum
            last_full = min(count, residual-2*pivot)
            value = (top_sum+(last_full-t+1)*pivot
                     +(residual-last_full-1)**2//4
                     -(residual-count-1)**2//4)
            best = max(best, value)
    return best


def zero_isolation_deficit_bound(k, r, t):
    """Exact excess above H for failed isolation allowing zero-deficit blocks."""
    capacity = (t-1)*k-t*(t+1)//2+1
    best = max(0, capacity-r+t-1)
    for count in range(1, t):
        value = min(k-count-1,
                    (2*count*k-count*(count+1))//(2*(count+1)))
        if value >= count:
            best = max(best, value)
    # Section13: at fixed top sum the largest feasible count dominates.
    for top_sum in range(t-1, k-t-1):
        count = min(r, k-top_sum-2)
        pivot = min(top_sum//(t-1), (k-t-top_sum)//2,
                    ((t-1)*k-t*(t-1)//2-t*top_sum)//(t-1))
        if pivot < 1:
            continue
        residual = k-top_sum
        last_full = min(count, residual-2*pivot)
        value = (top_sum+(last_full-t+1)*pivot
                 +(residual-last_full-1)**2//4
                 -(residual-count-1)**2//4)
        best = max(best, value)
    return best


@lru_cache(None)
def certificate_data(k, r, t):
    """Section12 threshold: selected blocks may have zero deficit."""
    H = r*(r+3)//2
    return H+zero_isolation_deficit_bound(k, r, t)+1, t*(t+3)//2


def certificate_score(k, r, t):
    dimension, degree = certificate_data(k, r, t)
    return 8*dimension - 4*degree - 3


@lru_cache(None)
def reference_parameters(n):
    best = None
    result = None
    for k in range(2, n):
        if not is_prime(k):
            continue
        m = n-k
        r = (m+k-1)//k
        if r < 3 or k <= r+1:
            continue
        for t in range(2, r):
            capacity = (t-1)*k-t*(t+1)//2+1
            if r*(r+3)//2+max(0, capacity-r+t-1)+1 >= m:
                continue  # Section12 lower bound proves no admissible pair is omitted.
            lower = r*(r+3)//2+max(0, capacity-r+t-1)+1
            B = t*(t+3)//2
            if best is not None and (4*k+8*lower-4*B-3, k, t) >= best:
                continue  # A lower loss key cannot improve the incumbent.
            dimension, degree = certificate_data(k, r, t)
            if not (degree <= dimension < m):
                continue
            key = (4*k+8*dimension-4*degree-3, k, t)
            if best is None or key < best:
                best = key
                result = k, m, r, t
    return result



def decode_prefix(k, bits):
    """Q=1+(Z^4+Z^2)R, monic R degree k-4 and next coefficient zero."""
    if len(bits) != k-5:
        raise ValueError("wrong encoded prefix length")
    R = (1 << (k-4)) | sum(bit << j for j, bit in enumerate(bits))
    return 1 ^ (R << 4) ^ (R << 2)


@lru_cache(None)
def compressed_parameters(n):
    """Optimize six actual-prefix lengths, including the constant-zero option."""
    best = (4*n, 0, 0, 0)
    result = None
    for saved in range(6):
        ref = reference_parameters(n+saved)
        if ref is None:
            continue
        k, m, r, t = ref
        loss = 4*(k-saved)+certificate_score(k, r, t)
        ignored = 5-saved
        key = (loss, k, t, ignored)
        if key < best:
            best = key
            result = k, m, r, t, ignored
    return result


@lru_cache(None)
def parameters(n):
    best=(4*n,0,0,0);result=None
    for k in range(2,n+1):
        if not is_prime(k):continue
        r0=(n+k-1)//k
        for r in (r0,r0-1):
            if r<3 or k<=r+1:continue
            m=min(n,r*k);ignored=n-m
            assert (m+k-1)//k==r and 0<=ignored<=k
            H=r*(r+3)//2
            for t in range(2,r):
                T=(t-1)*k-t*(t+1)//2+1;B=t*(t+3)//2
                lower=H+max(0,T-r+t-1)+1
                if lower>=m or (4*ignored+8*lower-4*B-3,k,t,ignored)>=best:continue
                d,B=certificate_data(k,r,t)
                if not B<=d<m:continue
                key=(4*ignored+8*d-4*B-3,k,t,ignored)
                if key<best:best=key;result=(k,m,r,t,ignored)
    return result

def poly_mod(a, b):
    if b == 0:
        raise ZeroDivisionError("zero polynomial divisor")
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_mod(a, b)
    return a


def field_mul(a, b, modulus):
    degree = modulus.bit_length() - 1
    high_bit = 1 << degree
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & high_bit:
            a ^= modulus
    return out


def field_pow(a, exponent, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = field_mul(result, a, modulus)
        exponent >>= 1
        if exponent:
            a = field_mul(a, a, modulus)
    return result


def irreducible_prime_degree(modulus, k):
    """Exact for prime k, as established in proof.md section 2."""
    if modulus.bit_length() != k + 1 or not is_prime(k):
        return False
    x = 2
    for _ in range(k):
        x = field_mul(x, x, modulus)
    return x == 2 and poly_gcd(modulus, 0b110) == 1


def evaluate(bits):
    """Return F_N on the sequence of N bits, in written input order."""
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("the input alphabet is {0,1}")
    params = parameters(len(bits))
    if params is None:
        return 0
    k, m, r, t, ignored = params
    modulus = irreducible_prime(k)
    payload = list(bits[ignored:]) + [0] * (r * k - m)
    symmetric = [1] + [0] * t
    for index in range(r):
        block = payload[index * k:(index + 1) * k]
        value = sum(bit << j for j, bit in enumerate(block))
        # i=index+1, e_i=r+1-i=r-index, w_i=e_i+1.
        power = (1 << (r - index + 1)) - 1
        value = field_pow(value, power, modulus)
        for j in range(min(index + 1, t), 0, -1):
            symmetric[j] ^= field_mul(value, symmetric[j - 1], modulus)
    return symmetric[t] & 1


# Deterministic prime-degree field construction proved in proof.md Appendix A.
def mod(a,b):
    while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a

def gcd(a,b):
    while b:a,b=b,mod(a,b)
    return a

def quotient(a,b):
    q=0
    while a and a.bit_length()>=b.bit_length():
        s=a.bit_length()-b.bit_length();q^=1<<s;a^=b<<s
    assert not a
    return q

def kernel(rows,n):
    rows=list(rows);pivots=[];rank=0
    for col in range(n):
        pivot=next((i for i in range(rank,n) if (rows[i]>>col)&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(n):
            if i!=rank and (rows[i]>>col)&1:rows[i]^=rows[rank]
        pivots.append(col);rank+=1
    basis=[]
    for free in range(n):
        if free in pivots:continue
        v=1<<free
        for row,pivot in enumerate(pivots):
            if (rows[row]>>free)&1:v^=1<<pivot
        basis.append(v)
    return basis

def factor_one(f):
    while True:
        d=f.bit_length()-1
        cols=[mod(1<<(2*j),f)^(1<<j) for j in range(d)]
        rows=[sum(((cols[j]>>i)&1)<<j for j in range(d)) for i in range(d)]
        basis=kernel(rows,d)
        if len(basis)==1:return f
        b=next(b for b in basis if b not in (0,1))
        g=gcd(f,b);h=quotient(f,g)
        assert g not in (1,f) and h not in (1,f)
        f=min(g,h)

def mul(a,b,f):
    out=0;degree=f.bit_length()-1
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if (a>>degree)&1:a^=f
    return out

@lru_cache(None)
def irreducible_prime(k):
    assert k>=3 and all(k%d for d in range(2,isqrt(k)+1))
    residue=2%k;m=1
    while residue!=1:residue=(2*residue)%k;m+=1
    order=(1<<m)-1;valuation=0
    while order%k==0:order//=k;valuation+=1
    f=factor_one((1<<k)-1)
    assert f.bit_length()-1==m
    for _ in range(1,valuation):
        lifted=sum(((f>>j)&1)<<(j*k) for j in range(m+1))
        f=factor_one(lifted)
        assert f.bit_length()-1==m
    a=2
    def kmul(x,y):return mul(x,y,f)
    zero=(0,)*k;one=(1,)+(0,)*(k-1)
    beta=(0,1)+(0,)*(k-2)
    def add(x,y):return tuple(a^b for a,b in zip(x,y))
    def square(x):
        out=[0]*k
        for j,c in enumerate(x):
            if not c:continue
            c=kmul(c,c)
            if 2*j>=k:c=kmul(c,a)
            out[(2*j)%k]^=c
        return tuple(out)
    def emul(x,y):
        tmp=[0]*(2*k-1)
        for i,c in enumerate(x):
            if c:
                for j,d in enumerate(y):
                    if d:tmp[i+j]^=kmul(c,d)
        for j in range(2*k-2,k-1,-1):
            if tmp[j]:tmp[j-k]^=kmul(tmp[j],a)
        return tuple(tmp[:k])
    current=beta;gamma=zero
    for _ in range(m):
        gamma=add(gamma,current)
        for _ in range(k):current=square(current)
    assert current==beta and any(gamma[1:])
    poly=[one];root=gamma
    for _ in range(k):
        out=[zero]*(len(poly)+1)
        for j,c in enumerate(poly):
            out[j+1]=add(out[j+1],c)
            out[j]=add(out[j],emul(c,root))
        poly=out;root=square(root)
    assert root==gamma
    assert all(not any(c[1:]) and c[0] in (0,1) for c in poly)
    return sum(c[0]<<j for j,c in enumerate(poly))

