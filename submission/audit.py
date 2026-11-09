"""Unchanged family evaluator and bounded parameter consistency checks."""
import json,sys
from math import isqrt

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def certificate_data(k,m,r,t):
    ell=m-(r-1)*k
    h=r*(r+3)//2-(ell==1)
    capacity=(t-1)*k-t*(t+1)//2+1-min(k-ell,t-1)
    total=r*k-r*(r+1)//2-(k-1)*(ell<=2)
    degree=t*(t+3)//2+t*(ell<=2)
    return h+max(capacity,total//(t+1))+1,degree


def parameters(n):
    best=result=None
    for k in range(2,n):
        if not is_prime(k):continue
        m=n-k;r=(m+k-1)//k
        if r<3 or k<=r+1:continue
        for t in range(2,r):
            d,b=certificate_data(k,m,r,t)
            if not b<=d<m:continue
            key=(3*k+8*d-4*b-3,k,t)
            if best is None or key<best:best=key;result=k,m,r,t
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
    k, m, r, t = params
    modulus = (1 << k) | sum(bit << j for j, bit in enumerate(bits[:k]))
    if not irreducible_prime_degree(modulus, k):
        return 0
    payload = list(bits[k:]) + [0] * (r * k - m)
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


def old_parameters(n):
    best=result=None
    for k in range(2,n):
        if not is_prime(k):continue
        m=n-k;r=(m+k-1)//k
        if r<3 or k<=r+1:continue
        for t in range(2,r):
            h=r*(r+3)//2;capacity=(t-1)*k-t*(t+1)//2+1
            total=r*k-r*(r+1)//2;b=t*(t+3)//2
            d=h+max(capacity,total//(t+1))+1
            if not b<=d<m:continue
            key=(4*k+8*d-4*b-3,k,t)
            if best is None or key<best:best=key;result=key
    return result


def secondary_reference(k):
    assert is_prime(k)
    lo,hi=0,1
    while 2*hi**3<=k:hi*=2
    while lo+1<hi:
        mid=(lo+hi)//2
        if 2*mid**3<=k:lo=mid
        else:hi=mid
    t=lo;r=(4*t*t+t+2)//4;n=k*(r+1);m=n-k
    d,b=certificate_data(k,m,r,t)
    assert 2*t**3<=k<2*(t+1)**3
    assert r>=3 and k>r+1 and 2<=t<r and b<=d<m
    assert m==r*k
    return {'N':n,'k':k,'t':t,'r':r,'D':d,'B':b,'loss':3*k+8*d-4*b-3}


def aligned_reference(n):
    lo,hi=0,1
    while 2*hi**5<=n:hi*=2
    while lo+1<hi:
        mid=(lo+hi)//2
        if 2*mid**5<=n:lo=mid
        else:hi=mid
    t=lo;denominator=(8*t+1)**2;numerator=64*n
    k=(numerator+denominator-1)//denominator
    while not is_prime(k):k+=1
    m=n-k;r=(m+k-1)//k;d,b=certificate_data(k,m,r,t)
    assert 2*t**5<=n<2*(t+1)**5
    assert r>=3 and k>r+1 and 2<=t<r and b<=d<m
    return {'N':n,'k':k,'t':t,'r':r,'D':d,'B':b,
            'target_numerator':numerator,'target_denominator':denominator,
            'gap_numerator':k*denominator-numerator,'loss':3*k+8*d-4*b-3}


def run():
    isolated=0
    for k in (2,3,5,7,11):
        for low in range(1<<k):
            q=(1<<k)|low
            if not irreducible_prime_degree(q,k):continue
            assert q.bit_count()%2==1
            for j in range(k):assert not irreducible_prime_degree(q^(1<<j),k)
            isolated+=1
    for n in (0,1,2,8,32,128,256):assert evaluate([0]*n)==0
    examples=[]
    for n in (128,192,256,512,1024,2048,4096,8192):
        selected=parameters(n);previous=old_parameters(n)
        assert selected and previous
        k,m,r,t=selected;d,b=certificate_data(k,m,r,t);ell=m-(r-1)*k
        weights=list(range(r+1,1,-1));widths=[k]*(r-1)+[ell]
        eligible=sorted(w for w,l in zip(weights,widths) if l>w)
        h=sum(min(w,l) for w,l in zip(weights,widths))
        capacity=sum(sorted((max(l-w,0) for w,l in zip(weights,widths)),reverse=True)[:t-1])
        total=sum(k+1-w for w in eligible)
        assert (d,b)==(h+max(capacity,total//(t+1))+1,sum(eligible[:t]))
        loss=3*k+8*d-4*b-3
        assert loss<=previous[0]-previous[1]
        examples.append({'N':n,'k':k,'r':r,'t':t,'ell':ell,'D':d,'B':b,
                         'loss':loss,'old_loss':previous[0],'old_k':previous[1],
                         'old_t':previous[2],'saving':previous[0]-loss})
    references=[secondary_reference(k) for k in (101,1009,10007,100003,1000003,10000019,100000007,1000000007)]
    aligned=[aligned_reference(n) for n in (128,1024,16384,262144,4194304,67108864,1000000000,1000000000000)]
    return {'all_passed':True,'aligned_references':aligned,'secondary_references':references,'isolated_valid_prefixes':isolated,'zero_inputs':7,
            'parameter_examples':examples,
            'scope':'Bounded exact consistency only; all universal claims are supplied mathematical proofs.'}

if __name__=='__main__':
    if json.load(sys.stdin)!={}:raise ValueError('expected empty input')
    print(json.dumps(run(),sort_keys=True))
