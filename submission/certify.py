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

from cavity_arithmetic import branch_tables
from radial_path import path_series_enclosure
from root_integration import integrate_roots


def score(rational):
    return Fraction(int(rational["numerator"]), int(rational["denominator"]))


def main():
    start = process_time()
    phases = []

    def report(name, **values):
        row = {"phase": name, **values}
        phases.append(row)
        print(json.dumps(row, sort_keys=True), file=sys.stderr, flush=True)

    request = sys.stdin.read()
    assert not request.strip() or json.loads(request) == {}
    gamma = Fraction(path_series_enclosure()["gamma_upper"])
    upper = Fraction(13064756648, 10**11)
    assert gamma < upper
    report("radial crossing cover", upper=str(gamma))

    tables = branch_tables(report)
    roots = integrate_roots(*tables)
    assert roots == (
        Fraction(77243471358243, 41943040000000000),
        Fraction(48843291382437, 41943040000000000),
        Fraction(254019548270397, 335544320000000000),
        Fraction(714166460848791, 1342177280000000000),
    )
    saving = sum(roots)
    assert saving == Fraction(5765021061632139, 1342177280000000000)
    report("four disjoint root savings", values=[str(value) for value in roots],
           total=str(saving))

    width = Fraction(80089, 80080) * (upper - Fraction(36, 71) * saving)
    assert width == Fraction(94277575837016466171469, 733768318976000000000000)
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(489153, 100000)
    assert coefficient == score(bound["guarantee"]["leading"]["coefficient"])
    assert score(bound["projection"]["exponent"]) == 1
    assert score(bound["guarantee"]["leading"]["exponent"]) == 1
    assert bound["guarantee"]["remainder"] == {"kind": "explicit", "terms": []}
    threshold = 1 / (2 * (coefficient - 1))
    assert 0 < width < threshold
    limit = 1 + 1 / (2 * width)
    report("strict cutwidth transfer", width=str(width), required_maximum=str(threshold),
           limiting_coefficient=str(limit))
    print(json.dumps({"completed": "one-band numerical theorem",
                      "circuit_coefficient": str(coefficient), "phases": phases,
                      "cpu_seconds": process_time() - start,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


if __name__ == "__main__":
    main()
