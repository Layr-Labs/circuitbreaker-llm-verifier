"""Unchanged family evaluator and bounded parameter consistency checks."""
import json,sys
from math import isqrt

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def parameters(n):
    if n == 0:
        return None
    target = 32 * n**3
    lo, hi = 0, 1
    while hi**5 <= target:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**5 > target:
            hi = mid
        else:
            lo = mid
    k = hi
    while not is_prime(k):
        k += 1
    if k >= n:
        return None
    m = n - k
    r = (m + k - 1) // k
    t = isqrt(r)
    if t * t != r:
        t += 1
    if t >= r or k * t < r * (r + 1):
        return None
    return k, m, r, t


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


def certificate_parameters(n):
    selected=parameters(n)
    if selected is None:return None
    k,m,r,t=selected;ell=m-(r-1)*k
    weights=list(range(r+1,1,-1));widths=[k]*(r-1)+[ell]
    eligible=sorted(w for w,l in zip(weights,widths) if l>w)
    if len(eligible)<t:return None
    h=sum(min(w,l) for w,l in zip(weights,widths))
    capacity=sum(sorted((max(l-w,0) for w,l in zip(weights,widths)),reverse=True)[:t-1])
    total=sum(max(k+1-w,0) for w in eligible)
    dimension=h+max(capacity,total//(t+1))+1;degree=sum(eligible[:t])
    full=r*(r+3)//2+max(sum(max(k-j,0) for j in range(2,t+1)),sum(max(k-j,0) for j in range(1,r+1))//(t+1))+1
    old_degree=t*(t+3)//2
    return {'N':n,'k':k,'r':r,'t':t,'ell':ell,'D':dimension,'B':degree,
            'full_D':full,'full_B':old_degree,'loss':3*k+8*dimension-4*degree-3,
            'previous_loss':3*k+8*full-4*old_degree-3}


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
    for n in (128,192,256,1024,16384,262144,4194304,67108864):
        row=certificate_parameters(n)
        assert row and 4<=row['B']<=row['D']<=row['full_D']
        assert row['B']>=row['full_B'] and row['loss']<=row['previous_loss']
        examples.append(row)
    return {'all_passed':True,'isolated_valid_prefixes':isolated,
            'zero_inputs':7,'parameter_examples':examples,
            'scope':'Bounded implementation consistency only; all universal claims are supplied mathematical proofs.'}

if __name__=='__main__':
    if json.load(sys.stdin)!={}:raise ValueError('expected empty input')
    print(json.dumps(run(),sort_keys=True))
