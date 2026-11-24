"""Replay the radial upper bound, the four savings, and the strict transfer."""
from fractions import Fraction
import json
from pathlib import Path
import resource
import sys
from time import process_time

ENTRYPOINT = Path(__file__).resolve()
EVIDENCE = ENTRYPOINT.parent
if ENTRYPOINT == Path("/snapshot/entrypoint.py"):
    EVIDENCE = EVIDENCE / "files"
sys.path.insert(0, str(EVIDENCE))

import cavity_arithmetic as branch
from radial_path import path_series_enclosure, pi_enclosure


def score(rational):
    return Fraction(int(rational["numerator"]), int(rational["denominator"]))


def main():
    if not __debug__:
        raise RuntimeError("The certificate requires assertions.")
    start = process_time()
    phases = []

    def report(name, **values):
        row = {"phase": name, "cpu_seconds": round(process_time() - start, 4), **values}
        phases.append(row)
        print(json.dumps(row, sort_keys=True), file=sys.stderr, flush=True)

    request = sys.stdin.read()
    assert not request.strip() or json.loads(request) == {}
    scale = branch.SCALE
    pi_hi = Fraction(pi_enclosure().high)
    assert branch.UNIT == scale
    assert branch.ROOT_LO**2 < 2*scale**2 < branch.ROOT_HI**2
    assert 2*branch.INV_R2_LO**2 < scale**2
    assert branch.CG_NUM**2 * pi_hi < 24*scale**2
    assert 2*branch.C2_NUM*pi_hi < 3*scale
    assert 8*branch.C3_NUM**2 * pi_hi**3 < 108*scale**2
    gamma = Fraction(path_series_enclosure()["gamma_upper"])
    upper = Fraction(13064756648, 10**11)
    assert gamma < upper
    report("radial crossing cover", upper=str(gamma))

    first, den1, second, den2, pair_totals, pair_denominator = branch.cavity(96, 18, 18, progress=report)
    assert den1 == 1061683200000000 and den2 == 45298483200000000
    roots = tuple(Fraction(count, den1) for count in first) + tuple(Fraction(count, den2) for count in second)
    assert roots == (
        Fraction(164715102419, 88473600000000),
        Fraction(431822033417, 353894400000000),
        Fraction(2898007771877, 3774873600000000),
        Fraction(1382763477463, 2516582400000000),
    )
    saving = sum(roots)
    assert saving == Fraction(99636594286381, 22649241600000000)
    report("four disjoint root savings", values=[str(value) for value in roots],
           total=str(saving))

    penalty = Fraction(4757, 10000)
    split_penalty = Fraction(577, 400000)
    argument = Fraction(353553391, 10**9)
    assert 0 < argument < 1 and 8 * argument**2 > 1
    arctan_upper = sum((-1)**k * argument**(2*k+1) / (2*k+1) for k in range(9))
    edge_density_upper = Fraction(8113009, 50000000)
    assert 3 * arctan_upper / (2 * Fraction(pi_enclosure().low)) < edge_density_upper

    assert Fraction(3, 7) <= penalty < Fraction(1, 2)
    assert 0 < split_penalty <= penalty and 1 - penalty - 2 * split_penalty > 0
    quad_lead = Fraction(3, 4) - penalty / 2
    excess = quad_lead - penalty
    correction = penalty * split_penalty / (1 + split_penalty)
    discriminant = (excess + 3*correction/2)**2 - 4*(quad_lead + correction/2)*correction
    assert discriminant == Fraction(-57251147, 802309664645000)
    assert discriminant < 0

    degree_one = roots[0] + roots[1]
    degree_two = roots[2] + roots[3]
    assert pair_denominator == 4529848320000000000000000
    pair_floors = tuple(Fraction(value, pair_denominator) for value in pair_totals)
    assert pair_floors == (
        Fraction(7645194187106450777, 62914560000000000000000),
        Fraction(14286193147524308791, 377487360000000000000000),
        Fraction(150566082781961777, 20971520000000000000000),
    )
    pair_weights = tuple(1 - 2*penalty + (2*k-3)*split_penalty for k in range(3))
    assert all(weight > 0 for weight in pair_weights)
    rebate = sum(weight * floor for weight, floor in zip(pair_weights, pair_floors, strict=True))
    assert rebate == Fraction(1136062787521397377607753, 150994944000000000000000000000)
    report("isolated-center improvement", floors=[str(floor) for floor in pair_floors],
           rebate=str(rebate))
    width = ((1 - split_penalty) * upper + split_penalty * edge_density_upper
             - (1 - penalty - 2 * split_penalty) * degree_one
             - (1 - penalty) * degree_two - rebate)
    assert width == Fraction(2153994546625377104919643583, 16777216000000000000000000000)
    report("unit-factor interpolation", discriminant=str(discriminant),
           edge_density_upper=str(edge_density_upper))
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(2447221, 500000)
    assert coefficient == score(bound["guarantee"]["leading"]["coefficient"])
    assert score(bound["projection"]["exponent"]) == 1
    assert score(bound["guarantee"]["leading"]["exponent"]) == 1
    assert bound["guarantee"]["remainder"] == {"kind": "explicit", "terms": []}
    threshold = 1 / (2 * (coefficient - 1))
    assert 0 < width < threshold
    limit = 1 + 1 / (2 * width)
    report("strict cutwidth transfer", width=str(width), required_maximum=str(threshold),
           limiting_coefficient=str(limit))
    print(json.dumps({"completed": "cubic crossing-cover certificate",
                      "circuit_coefficient": str(coefficient), "phases": phases,
                      "cpu_seconds": process_time() - start,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
