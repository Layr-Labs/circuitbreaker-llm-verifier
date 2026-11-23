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


ROOTS = {
    'first': (Fraction(26199709529, 14745600000000), Fraction(1920189511, 1658880000000), Fraction(0, 1), Fraction(0, 1)),
    'second': (Fraction(45582353, 552960000000), Fraction(2091107119, 33177600000000), Fraction(152374370417, 199065600000000), Fraction(174410504053, 318504960000000)),
    'all': (Fraction(82245716827, 44236800000000), Fraction(13498299113, 11059200000000), Fraction(152374370417, 199065600000000), Fraction(174410504053, 318504960000000)),
}


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

    request = json.loads(sys.stdin.read() or "{}")
    part = request["part"]
    assert set(request) == {"part"} and part in ROOTS
    if part != 'second':
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

    first, den1, second, den2 = branch.cavity(72, 24, 24, progress=report, part=part)
    assert den1 == 796262400000000 and den2 == 19110297600000000
    roots = tuple(Fraction(count, den1) for count in first) + tuple(Fraction(count, den2) for count in second)
    assert roots == ROOTS[part]
    assert tuple(a + b for a, b in zip(ROOTS["first"], ROOTS["second"], strict=True)) == ROOTS["all"]
    saving = sum(ROOTS["all"])
    assert saving == Fraction(155458852481, 35389440000000)
    penalty = Fraction(49249, 100000)
    multiplier = 1 + 9*(1-2*penalty)**2 / (16*penalty*(3-2*penalty))
    width = multiplier * (Fraction(13064756648, 10**11) - (1-penalty)*saving)
    assert width == Fraction(2004702698372600349660503321, 15608743719862272000000000000)
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == score(bound["guarantee"]["leading"]["coefficient"]) == Fraction(611629, 125000)
    assert score(bound["projection"]["exponent"]) == score(bound["guarantee"]["leading"]["exponent"]) == 1
    assert bound["guarantee"]["remainder"] == {"kind": "explicit", "terms": []}
    assert 0 < width < 1 / (2 * (coefficient - 1))
    report("replayed root slice", part=part, values=list(map(str, roots)))
    print(json.dumps({"part": part, "verified_roots": list(map(str, roots)),
                      "combined_width": str(width), "circuit_coefficient": str(coefficient),
                      "requires_other_slice": part != "all", "phases": phases,
                      "cpu_seconds": process_time() - start,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
