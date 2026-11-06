"""Uniform evaluator for the zero-deficit-selection certificate-prefix family in proof.md.

This module does not choose a field polynomial: those bits are part of the
input. Integers encode polynomials over F_2 with bit j the coefficient of Z^j.
"""
from math import isqrt


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
    for count in range(t, r+1):
        for top_sum in range(t-1, k-count-1):
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


def certificate_data(k, r, t):
    """Section12 threshold: selected blocks may have zero deficit."""
    H = r*(r+3)//2
    return H+zero_isolation_deficit_bound(k, r, t)+1, t*(t+3)//2


def certificate_score(k, r, t):
    dimension, degree = certificate_data(k, r, t)
    return 8*dimension - 4*degree - 3


def parameters(n):
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
