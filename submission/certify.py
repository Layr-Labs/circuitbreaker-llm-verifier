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
        Fraction(7774910480781, 4194304000000000),
        Fraction(24881981336331, 20971520000000000),
        Fraction(128009116150311, 167772160000000000),
        Fraction(144620874738381, 268435456000000000),
    )
    saving = sum(roots)
    assert saving == Fraction(5827595462269497, 1342177280000000000)
    report("four disjoint root savings", values=[str(value) for value in roots],
           total=str(saving))

    penalty = Fraction(49257, 100000)
    multiplier = 1 + 9*(1-2*penalty)**2 / (16*penalty*(3-2*penalty))
    width = multiplier * (upper - (1-penalty)*saving)
    assert multiplier == Fraction(4411483561, 4410931512)
    assert width == Fraction(76051814369729864129395244548969, 592025205904244736000000000000000)
    bound = json.loads((EVIDENCE / "bound.json").read_text(encoding="utf-8"))
    coefficient = score(bound["projection"]["coefficient"])
    assert coefficient == Fraction(489223, 100000)
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
