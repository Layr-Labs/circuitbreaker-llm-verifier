"""Exact check of the length guards and remainder bound of proof.md Sections 4-5."""
import math, random, sys
def isprime(k):
    if k < 2: return False
    i = 2
    while i * i <= k:
        if k % i == 0: return False
        i += 1
    return True
def params(N):
    B = 1
    lo, hi = 1, 4 * N + 4
    while lo < hi:  # least B with B^5 > 32 N^3
        mid = (lo + hi) // 2
        if mid ** 5 > 32 * N ** 3: hi = mid
        else: lo = mid + 1
    k = lo
    while not isprime(k): k += 1
    m = N - k; r = -(-m // k); t = math.isqrt(r); t += t * t < r
    return k, m, r, t
def check(N):
    k, m, r, t = params(N)
    assert k < N and 1 <= t < r and k * t >= r * (r + 1), N
    R = 4 * k + 9 * ((k * r) // t + k * t + r * (r + 1) + 2)
    assert R ** 5 <= 182 ** 5 * N ** 4, (N, R)  # R <= 182 N^(4/5), exact
    return R / N
rng = random.Random(0)
Ns = list(range(2**20, 2**20 + 20000)) + [rng.randrange(2**20, 2**60) for _ in range(20000)] + [2**e + d for e in range(20, 61) for d in (-1, 0, 1) if 2**e + d >= 2**20]
worst = max(check(N) for N in Ns)
print("guards and R(N)<=182 N^(4/5) hold on", len(Ns), "lengths; max R/N =", round(worst, 4))
