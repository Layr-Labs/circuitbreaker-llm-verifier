"""Replay the radial upper bound, the five savings, and the strict transfer."""
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

    request = json.loads(sys.stdin.read())
    assert request in ({"part": "positive"}, {"part": "remaining"})
    part = request["part"]
    scale = branch.SCALE
    pi_hi = Fraction(pi_enclosure().high)
    assert branch.UNIT == scale
    assert branch.ROOT_LO**2 < 2*scale**2 < branch.ROOT_HI**2
    assert 2*branch.INV_R2_LO**2 < scale**2
    assert branch.CG_NUM**2 * pi_hi < 24*scale**2
    assert 2*branch.C2_NUM*pi_hi < 3*scale
    assert 8*branch.C3_NUM**2 * pi_hi**3 < 108*scale**2
    upper = Fraction(1633094581, 12500000000)
    positive1, negative1, positive2, negative2 = (
        Fraction(3066820985769, 1638400000000000),
        Fraction(1352332129989, 1107558400000000),
        Fraction(323900959788639, 419430400000000000),
        Fraction(31681288333347, 57593036800000000),
    )
    expected_bonus = (
        Fraction(789595681841030773071, 2879651840000000000000000),
        Fraction(101073994950203183679, 719912960000000000000000),
        Fraction(5683607378432641323, 143982592000000000000000),
        Fraction(85128686549922520527, 719912960000000000000000),
        Fraction(3144098781649314927, 359956480000000000000000),
    )
    if part == "positive":
        gamma = Fraction(path_series_enclosure()["gamma_upper"])
        assert gamma < upper
        report("radial crossing cover", upper=str(gamma))
        first, den1, second, den2, _, _ = branch.cavity(
            320, 40, 24, progress=report, positive_only=True)
        assert den1 == 26214400000000000 and den2 == 1677721600000000000
        assert first[1] == 0 and second[1] == 0
        result = (Fraction(first[0], den1), Fraction(second[0], den2))
        assert result == (positive1, positive2)
        report("positive root pair", values=[str(item) for item in result])
    else:
        first, den1, second, den2, pair_totals, pair_denominator = branch.cavity(
            104, 16, 18, progress=report, negative_only=True)
        assert den1 == 1107558400000000 and den2 == 57593036800000000
        assert first[0] == 0 and second[0] == 0
        result = (Fraction(first[1], den1), Fraction(second[1], den2))
        assert result == (negative1, negative2)
        report("negative root pair", values=[str(item) for item in result])
        assert pair_denominator == 5759303680000000000000000
        pair_floors = tuple(Fraction(val, pair_denominator) for val in pair_totals)
        assert pair_floors == expected_bonus
        report("isolated-center improvement", floors=[str(item) for item in pair_floors])

    penalty = Fraction(236349, 500000)
    split_penalty = Fraction(18403, 10000000)
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
    assert discriminant == Fraction(-666946946853, 50184199335204500000)
    assert discriminant < 0

    degree_one = positive1 + negative1
    degree_two = positive2 + negative2
    pair_weights = tuple(1 - 2*penalty + (2*k-3)*split_penalty for k in range(3))
    assert all(weight > 0 for weight in pair_weights)
    rebate = sum(weight * floor for weight, floor in zip(pair_weights, expected_bonus[:3], strict=True))
    additional_target = Fraction(8963129, 20000000)
    assert pair_weights[0] + additional_target <= 1 - penalty
    assert pair_weights[0] + 2 * additional_target <= 2 * penalty
    rebate += pair_weights[0] * expected_bonus[3] + additional_target * expected_bonus[4]
    assert rebate == Fraction(944900049264572164110577953, 28796518400000000000000000000000)
    width = ((1 - split_penalty) * upper + split_penalty * edge_density_upper
             - (1 - penalty - 2 * split_penalty) * degree_one
             - (1 - penalty) * degree_two - rebate)
    assert width == Fraction(3696211249941850775474886159547, 28796518400000000000000000000000)
    report("unit-factor interpolation", discriminant=str(discriminant),
           edge_density_upper=str(edge_density_upper))
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(489541, 100000)
    assert coefficient == score(bound["guarantee"]["leading"]["coefficient"])
    assert score(bound["projection"]["exponent"]) == 1
    assert score(bound["guarantee"]["leading"]["exponent"]) == 1
    assert bound["guarantee"]["remainder"] == {"kind": "explicit", "terms": []}
    threshold = 1 / (2 * (coefficient - 1))
    assert 0 < width < threshold
    limit = 1 + 1 / (2 * width)
    assert coefficient < Fraction("4.8954101") < limit
    report("combined consistency check", width=str(width), required_maximum=str(threshold),
           limiting_coefficient=str(limit))
    print(json.dumps({"completed": "one part of the two-part cubic certificate",
                      "part": part, "requires_other_part": True,
                      "combined_width": str(width), "verified_floors": list(map(str, result if part == "positive" else result + pair_floors)),
                      "circuit_coefficient": str(coefficient), "phases": phases,
                      "cpu_seconds": process_time() - start,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
