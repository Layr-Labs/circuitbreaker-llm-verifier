"""Exact finite check of the explicit remainder in proof.md section 5.

For each tested length N it recomputes the family parameters with family.py,
forms the exact loss R(N)=(7/2)k+7*d_N with d_N=floor(kr/t)+kt+r(r+1)+1,
and compares it with the declared explicit remainder

    E(N) = (10927/300)N^(4/5) + 42N^(3/5) + (7/2)N^(2/5) + (119/24)N^(1/5) + 35/2.

Each comparison R(N) < E(N) is decided exactly: the irrational powers of
x=N^(1/5) are replaced by rational lower bounds x_lo<=x, obtained by
integer bisection on x_lo^5<=N, so E(x_lo)<=E(N) and R(N)<E(x_lo) suffices.
It also runs the unchanged supplementary audit in audit.py, so that one
declared experiment covers both. A finite run does not prove the universal
statement; proof.md does.
"""
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys

THRESHOLD = 16807  # 7^5, the analytic threshold of proof.md section 5
COEFFICIENTS = (Fraction(10927, 300), Fraction(42), Fraction(7, 2),
                Fraction(119, 24), Fraction(35, 2))


def load_module(name):
    path = Path('/snapshot/files') / name
    if not path.exists():
        path = Path(__file__).with_name(name)
    spec = importlib.util.spec_from_file_location('candidate_' + path.stem, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fifth_root_lower(n, scale=10**6):
    """Largest q/scale with (q/scale)^5 <= n."""
    lo, hi = 0, scale * (n + 1)
    target = n * scale**5
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**5 <= target:
            lo = mid
        else:
            hi = mid
    return Fraction(lo, scale)


def explicit_remainder(x):
    c4, c3, c2, c1, c0 = COEFFICIENTS
    return c4 * x**4 + c3 * x**3 + c2 * x**2 + c1 * x + c0


def exact_loss(params):
    k, r, t = params[0], params[2], params[3]
    d = (k * r) // t + k * t + r * (r + 1) + 1
    return Fraction(7, 2) * k + 7 * d


def run(dense_limit=100000):
    f = load_module('family.py')
    guard_failures = [n for n in range(1, THRESHOLD) if f.parameters(n) is None]
    lengths = list(range(THRESHOLD, dense_limit + 1))
    n = dense_limit
    while n < 10**15:
        n = n * 3 // 2 + 7
        lengths.append(n)
    worst = Fraction(0)
    worst_n = None
    for n in lengths:
        params = f.parameters(n)
        assert params is not None, n
        k, _, r, t = params
        assert k < n and 1 <= t < r and k * t >= r * (r + 1), n
        loss = exact_loss(params)
        bound = explicit_remainder(fifth_root_lower(n))
        assert loss < bound, n
        ratio = loss / bound
        if ratio > worst:
            worst, worst_n = ratio, n
    audit = load_module('audit.py').run()
    assert audit['all_passed']
    return {
        'all_passed': True,
        'audit': audit,
        'largest_guard_failure_below_threshold': max(guard_failures),
        'lengths_checked': len(lengths),
        'dense_range': [THRESHOLD, dense_limit],
        'largest_sparse_length': lengths[-1],
        'max_loss_over_declared_remainder': float(worst),
        'max_ratio_length': worst_n,
        'scope': 'Finite exact check only; the all-length statement is proved in proof.md section 5.',
    }


if __name__ == '__main__':
    request = json.load(sys.stdin)
    if request != {}:
        raise ValueError('the declared check takes an empty input object')
    print(json.dumps(run(), sort_keys=True))
