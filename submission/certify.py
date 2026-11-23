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

    first, den1, second, den2 = branch.cavity(76, 24, 24, progress=report)
    assert den1 == 887193600000000 and den2 == 22475571200000000
    roots = tuple(Fraction(count, den1) for count in first) + tuple(Fraction(count, den2) for count in second)
    assert roots == (
        Fraction(550213338463, 295731200000000),
        Fraction(45197628179, 36966400000000),
        Fraction(906255874791, 1182924800000000),
        Fraction(3081964908519, 5618892800000000),
    )
    saving = sum(roots)
    assert saving == Fraction(790744743289, 179804569600000)
    report("four disjoint root savings", values=[str(value) for value in roots],
           total=str(saving))

    penalty = Fraction(49249, 100000)
    multiplier = 1 + 9*(1-2*penalty)**2 / (16*penalty*(3-2*penalty))
    width = multiplier * (upper - (1-penalty)*saving)
    assert multiplier == Fraction(39700164001, 39695087992)
    assert width == Fraction(458332807490279849354381269309, 3568679105817844121600000000000)
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(1223277, 250000)
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
