"""Complete exact enumeration for proof.md Section 5a.

For every integer N with 112 <= N < 10**7 this recomputes the evaluator's parameters
exactly (B least with B^5 > 32 N^3, k least prime >= B, m=N-k, r=ceil(m/k),
t=ceil(sqrt(r))), checks the guards k<N, t<r, k*t>=r(r+1), and checks the integer
inequality (2R)^5 <= 76^5 * N^4 for R = (7/2)k + 7 d_N, d_N = floor(kr/t)+kt+r(r+1)+1.
That inequality is equivalent to R <= 38 * N^(4/5).
"""
import math, sys

def isprime(q):
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    i = 3
    while i * i <= q:
        if q % i == 0:
            return False
        i += 2
    return True

def least_B(N):
    lo, hi = 1, 1
    while hi**5 <= 32 * N**3:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if mid**5 > 32 * N**3:
            hi = mid
        else:
            lo = mid + 1
    return lo

def check(lo=112, hi=10**7):
    """Cover every integer N in [lo,hi) exactly, segment by segment.

    B (hence k) and r=ceil((N-k)/k) are nondecreasing step functions of N, and t, the
    guards t<r, k*t>=r(r+1) and d_N depend only on (k,r). On a maximal run of N with
    fixed (k,r), k<N and R<=38*N^(4/5) are monotone in N, so they hold on the whole
    run iff they hold at its first element. We therefore test the first N of every
    run exactly, and assert that the runs tile [lo,hi) with no gap.
    """
    failures, first, runs = 0, [], 0
    max_ratio, argmax = 0.0, None
    N = lo
    while N < hi:
        B = least_B(N)
        k = B
        while not isprime(k):
            k += 1
        # B stays the least valid value until 32 N^3 >= B^5; find that N by binary search.
        a, b = N, hi
        while a < b:
            mid = (a + b) // 2
            if 32 * mid**3 >= B**5:
                b = mid
            else:
                a = mid + 1
        endB = a                                  # first N with a larger B (or hi)
        while N < endB:
            ok = k < N
            if ok:
                m = N - k
                r = -(-m // k)
                t = math.isqrt(r)
                t += (t * t < r)
                ok = t < r and k * t >= r * (r + 1)
                endR = min(endB, k * r + k + 1)   # next N gives r+1
            else:
                endR = N + 1
            if ok:
                d = (k * r) // t + k * t + r * (r + 1) + 1
                twoR = 7 * k + 14 * d
                ok = twoR**5 <= (76**5) * N**4
                ratio = twoR / 2 / N**0.8
                if ratio > max_ratio:
                    max_ratio, argmax = ratio, N
            if not ok:
                failures += 1
                if len(first) < 5:
                    first.append(N)
            runs += 1
            N = endR
    return {"range": [lo, hi], "runs_checked": runs, "failures": failures,
            "first_failures": first, "max_R_over_N45": round(max_ratio, 4), "argmax_N": argmax}

if __name__ == "__main__":
    print(check())
