"""Replay the radial upper bound, the seven savings, and the strict transfer."""
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
from root_integration import positive_ray_floors


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
    assert request in tuple({"part": name} for name in ("positive", "negative"))
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
        Fraction(6132802254699, 3276800000000000),
        Fraction(8044239057, 6553600000000),
        Fraction(647713169760441, 838860800000000000),
        Fraction(39692631200919, 71932313600000000),
    )
    expected_bonus = (
        Fraction(2559403719770627073, 8991539200000000000000),
        Fraction(493995166845699416907, 1798307840000000000000000),
        Fraction(33078437798832419751, 359661568000000000000000),
        Fraction(1120461605982560727, 7193231360000000000000),
        Fraction(2897097909552259911, 224788480000000000000000),
    )
    expected_rays = (
        Fraction(3225870319753189, 1474560000000000000000),
        Fraction(3932004641225701, 368640000000000000000),
    )
    if part == "positive":
        gamma = Fraction(path_series_enclosure()["gamma_upper"])
        assert gamma < upper
        report("radial crossing cover", upper=str(gamma))
        first, den1, second, den2, _, _ = branch.cavity(
            320, 40, 24, progress=report, positive_only=True)
        assert first[1] == second[1] == 0
        base = (Fraction(first[0], den1), Fraction(second[0], den2))
        assert base == (positive1, positive2)
        report("positive root pair", floors=[str(v) for v in base])
        ray_values, denominator = positive_ray_floors(branch)
        adjusted = tuple(Fraction(value, denominator) for value in ray_values)
        assert adjusted == expected_rays
        report("positive ray count floors", floors=[str(v) for v in adjusted])
    else:
        first, den1, second, den2, totals, denominator = branch.cavity(
            112, 20, 28, progress=report, negative_only=True)
        assert first[0] == second[0] == 0
        base = (Fraction(first[1], den1), Fraction(second[1], den2))
        assert base == (negative1, negative2)
        report("negative root pair", floors=[str(v) for v in base])
        adjusted = tuple(Fraction(value, denominator) for value in totals)
        assert adjusted == expected_bonus
        report("isolated-center improvement", floors=[str(v) for v in adjusted])
    verified = base + adjusted

    penalty = Fraction(47091, 100000)
    split_penalty = Fraction(2102599, 1000000000)
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
    assert discriminant == Fraction(-105331260465423, 100420961892255480100000)
    assert discriminant < 0

    degree_one = positive1 + negative1
    degree_two = positive2 + negative2
    pair_weights = tuple(1 - 2*penalty + (2*k-3)*split_penalty for k in range(3))
    assert all(weight > 0 for weight in pair_weights)
    rebate = sum(weight * floor for weight, floor in zip(pair_weights, expected_bonus[:3], strict=True))
    additional_target = Fraction(889947797, 2000000000)
    assert pair_weights[0] + additional_target <= 1 - penalty
    assert pair_weights[0] + 2 * additional_target <= 2 * penalty
    rebate += pair_weights[0] * expected_bonus[3] + additional_target * expected_bonus[4]
    rebate += sum(weight * floor for weight, floor in zip(pair_weights[:2], expected_rays, strict=True))
    assert rebate == Fraction(406564145190519541121747715463, 8092385280000000000000000000000000)
    width = ((1 - split_penalty) * upper + split_penalty * edge_density_upper
             - (1 - penalty - 2 * split_penalty) * degree_one
             - (1 - penalty) * degree_two - rebate)
    assert width == Fraction(148364259508455484820287437447791, 1156055040000000000000000000000000)
    report("unit-factor interpolation", discriminant=str(discriminant),
           edge_density_upper=str(edge_density_upper))
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(612, 125)
    assert coefficient == score(bound["guarantee"]["leading"]["coefficient"])
    assert score(bound["projection"]["exponent"]) == 1
    assert score(bound["guarantee"]["leading"]["exponent"]) == 1
    assert bound["guarantee"]["remainder"] == {"kind": "explicit", "terms": []}
    threshold = 1 / (2 * (coefficient - 1))
    assert 0 < width < threshold
    limit = 1 + 1 / (2 * width)
    assert limit == Fraction(726391779508455484820287437447791, 148364259508455484820287437447791)
    assert coefficient < limit
    report("combined consistency check", width=str(width), required_maximum=str(threshold),
           limiting_coefficient=str(limit))
    print(json.dumps({"completed": "one part of the two-part cubic certificate",
                      "part": part,
                      "required_other_part": {"positive": "negative", "negative": "positive"}[part],
                      "verified_floors": [str(v) for v in verified],
                      "combined_width": str(width),
                      "circuit_coefficient": str(coefficient), "phases": phases,
                      "cpu_seconds": process_time() - start,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
