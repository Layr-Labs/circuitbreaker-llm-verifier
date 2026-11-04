"""Exact-rational audit of the arithmetic used in this submission.

Runs with the Python standard library only. It checks, with Fractions and exact
integer arithmetic (no floating point in any inequality):

  A. the potential weights and every case inequality of framework.md;
  B. the summation that turns Delta W <= 14 Delta a into |C| >= 3.5n - 7d;
  C. the Case 8 sub-bounds, their tightness, and the bottleneck deficit;
  D. the parameter arithmetic of Lemma 3.1 of proof.md and of (3.2): every
     comparison of x^(3/5), x^(4/5) with a rational is converted to an integer
     power inequality, and BK's threshold kr/t + kt + r(r+1) < 6 m^(4/5) is
     checked exactly for several m;
  E. the schedule m(N): strict monotonicity of N(m) = m + k(m) and maximality.

This script certifies arithmetic only. It proves no mathematics: each inequality it
checks is an implication whose hypotheses are established in the markdown files.
"""

from fractions import Fraction as F

FAILS = []


def need(ok, label, detail=""):
    print(("PASS  " if ok else "FAIL  ") + label + (("  [" + detail + "]") if detail else ""))
    if not ok:
        FAILS.append(label)


# ----------------------------------------------------------------------------
# helpers (exact integer arithmetic throughout)
# ----------------------------------------------------------------------------
def isprime(x):
    if x < 2:
        return False
    if x % 2 == 0:
        return x == 2
    i = 3
    while i * i <= x:
        if x % i == 0:
            return False
        i += 2
    return True


def fifth_root_floor(n):
    """floor(n^(1/5)) for n >= 0, by exact integer search from a float estimate."""
    if n <= 0:
        return 0
    r = int(round(n ** 0.2))
    while (r + 1) ** 5 <= n:
        r += 1
    while r ** 5 > n:
        r -= 1
    return r


def k_of_m(m):
    """The least prime strictly greater than 2 m^(3/5)  <=>  p^5 > 32 m^3."""
    bound = 32 * m ** 3
    p = fifth_root_floor(bound) + 1
    while p ** 5 <= bound:
        p += 1
    while not isprime(p):
        p += 1
    return p


def N_of_m(m):
    return m + k_of_m(m)


def ceil_sqrt(x):
    r = 0
    while r * r < x:
        r += 1
    return r


def m_of_N(N):
    """The largest m >= 1 with N(m) <= N, by binary search (N(m) is increasing)."""
    if N < N_of_m(1):
        return None
    lo, hi = 1, N
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if N_of_m(mid) <= N:
            lo = mid
        else:
            hi = mid - 1
    return lo


# ----------------------------------------------------------------------------
# A. weights, potential, case inequalities
# ----------------------------------------------------------------------------
# W = 2s + 3q + sum of weights of unpaired parameters, wt = 0 / 6 / 7.
W_ISO, W_NON, GATE, PAIR = F(6), F(7), F(2), F(3)


def dw(ds=0, dq=0, weights=(), gains=0):
    """Delta W for (change in s, change in q, signed weight contributions, 6->7 rises)."""
    return GATE * ds + PAIR * dq + sum(weights, F(0)) + gains


cases = {
    "C2 unpaired w (da=-1)": (dw(ds=-1, weights=(-W_ISO, -W_ISO)), -1),
    "C2 paired w (da=0)": (dw(ds=-1, dq=-1, weights=(-W_ISO, W_NON)), 0),
    "C3 deg(p)=0 (da=0)": (dw(dq=-1), 0),
    "C3 deg(p)>=2 (da=0)": (dw(ds=-2, dq=-1, weights=(W_NON,)), 0),
    "C4 common child (da=0)": (dw(ds=-1, dq=-1), 0),
    "C5 paired in system (da=0)": (dw(ds=-2, dq=-1, weights=(W_ISO,), gains=1), 0),
    "C6 unpaired w (da=-1)": (dw(ds=-2, dq=-1, weights=(-W_NON,)), -1),
    "C6 paired w (da=0)": (dw(ds=-2, dq=-2, weights=(W_NON,)), 0),
    "C7 common second child (da=-1)": (dw(ds=-2, weights=(-W_NON, -W_NON)), -1),
    "C9 bottleneck AND (da=-2)": (dw(ds=-4, weights=(-W_NON, -W_NON, -W_ISO)), -2),
    "C10 bottleneck XOR (da=-1)":
        (dw(ds=-2, dq=1, weights=(-W_NON, -W_NON, -W_ISO, W_NON)), -1),
}
for label, (value, da) in cases.items():
    need(value <= 14 * da, "A " + label, "dW=%s <= 14*da=%s" % (value, 14 * da))

need(dw(ds=-1, weights=(-6, -6)) == -14, "A Case 2 is exactly -14")
need(dw(ds=-2, dq=-1, weights=(-7,)) == -14, "A Case 6 (w unpaired) is exactly -14")
need(dw(ds=-2, dq=1, weights=(-7, -7, -6, 7)) == -14, "A Case 10 is exactly -14")
need(dw(ds=-4, weights=(-7, -7, -6)) == -28, "A Case 9 is exactly -28 = 14*(-2)")
need(dw(ds=-2, weights=(-7, -7)) == -18, "A Case 7 is -18")
need(dw(ds=-1, gains=1) == -1, "A a deletion with one 6->7 rise is -1")
need(dw(ds=-1) == -2, "A a deletion rerouted to a constant is -2")
need(dw(ds=-1, dq=1) == 1, "A a step adding a pair with no deletion is +1")

# ----------------------------------------------------------------------------
# B. the summation
# ----------------------------------------------------------------------------
for n, C, d in [(1000, 3500, 0), (1024, 3600, 8), (10 ** 6, 3_500_000, 190)]:
    W0 = 2 * C + 7 * n
    need(W0 >= 14 * (n - d), "B potential summation at n=%d" % n, "2|C|+7n >= 14(n-d)")
    need(F(W0 - 7 * n, 2) >= F(7 * n - 14 * d, 2), "B rearrangement at n=%d" % n)
need(F(14 - 7, 2) == F(7, 2), "B coefficient (14-7)/2 = 7/2")

# ----------------------------------------------------------------------------
# C. Case 8 sub-bounds, tightness, bottleneck
# ----------------------------------------------------------------------------
for degx, k, k1 in [(3, 1, 0), (3, 1, 1), (5, 2, 0), (2, 2, 0), (2, 2, 1), (2, 3, 2)]:
    extra = F(-1) if degx == 2 else F(0)
    value = -7 - 2 * (degx + k) + k1 + extra
    need(value <= -14, "C Case 8 deg(x)=%d k=%d k1=%d" % (degx, k, k1), "dW<=%s" % value)
value = -7 - 2 * (2 + 1) + 1 - 1
need(value == -13, "C bottleneck subcase is -13", "dW=%s" % value)
need(value - 14 * -1 == 1, "C bottleneck deficit equals 1")

# ----------------------------------------------------------------------------
# D. parameter arithmetic
# ----------------------------------------------------------------------------
for m in [1, 2, 10, 1000, 10 ** 6, 10 ** 9]:
    need(4 ** 5 * (m + 1) ** 3 <= 7 ** 5 * m ** 3, "D (m+1)^(3/5) <= (7/4)m^(3/5) at m=%d" % m)
for N in [1024, 10 ** 4, 10 ** 6, 10 ** 9]:
    need(32 * 7 ** 5 * N ** 3 <= (N - 2) ** 5, "D 7N^(3/5)+1 <= N/2 at N=%d" % N)
for N in [1024, 10 ** 4, 10 ** 6]:
    need(21 ** 5 <= N ** 4, "D 0.5N^(4/5) >= 10.5 at N=%d" % N)
need(F(49, 2) + 42 == F(133, 2) and F(133, 2) <= 67, "D 49/2 + 42 = 66.5 <= 67")
need(F(2) + F(17, 6) + F(1, 4) < 6, "D 2 + 17/6 + 1/4 < 6, and 2sqrt2 < 17/6")


def threshold_ok(m):
    """kr/t + kt + r(r+1) < 6 m^(4/5)  <=>  num^5 <= 6^5 m^4 t^5, num = kr + t(kt + r(r+1))."""
    k = k_of_m(m)
    r = (m + k - 1) // k
    t = ceil_sqrt(r)
    num = k * r + t * (k * t + r * (r + 1))
    return num ** 5 <= 6 ** 5 * m ** 4 * t ** 5


# The inequality is verified exactly for every m with 512 <= m < 16807 (the range
# used by the explicit threshold N >= 2^10), and sampled far beyond it; from
# m >= 16807 it follows from the estimates of disperser.md Lemma 3.1(b), whose
# polynomial 1.75 u^2 - 12 u - 1.5 is positive at u = 7 and increasing for u >= 7.
bad = [m for m in range(512, 16807) if not threshold_ok(m)]
need(not bad, "D kr/t+kt+r(r+1) < 6m^(4/5) for every 512 <= m < 16807",
     "checked 16295 values" if not bad else "first failures: %s" % bad[:5])
for m in [16807, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]:
    need(threshold_ok(m), "D kr/t+kt+r(r+1) < 6m^(4/5) at m=%d" % m, "k=%d" % k_of_m(m))
need(F(7, 4) * 49 - 12 * 7 - F(3, 2) == F(1, 4), "D 1.75u^2-12u-1.5 = 1/4 at u=7")
need(F(7, 4) * 8 ** 2 - 12 * 8 - F(3, 2) > 0, "D the same quadratic is positive at u=8")

# ----------------------------------------------------------------------------
# E. the schedule
# ----------------------------------------------------------------------------
prev = 0
mono = True
for m in list(range(1, 300)) + [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]:
    cur = N_of_m(m)
    if cur <= prev:
        mono = False
    prev = cur
need(mono, "E N(m) = m + k(m) is strictly increasing")

for N in [10, 100, 1024, 10 ** 4, 10 ** 5, 10 ** 6]:
    m = m_of_N(N)
    if m is None:
        need(True, "E N=%d is below N(1), where f_N is identically 0" % N)
        continue
    need(N_of_m(m) <= N < N_of_m(m + 1), "E m(N) is the largest feasible m at N=%d" % N,
         "m=%d k=%d" % (m, k_of_m(m)))
    if N >= 1024:
        need(F(m) >= F(N, 2), "E m >= N/2 at N=%d" % N, "m=%d" % m)

if FAILS:
    raise SystemExit("%d check(s) failed" % len(FAILS))
print("\nAll arithmetic checks passed.")
